# Error-adaptive tessellation: the glue between FEData and the isubd core
# (src/isubd.jl), and the compute-graph wiring for `solutionplot(...;
# adaptive=true)` and `meshplot(...; adaptive=true)` (#161).
#
# Refinement is driven by two interpolation-error estimators, both instances
# of DeviationLoD along a triangle's edges, refining when *either* asks:
#   - geometry error: the exact geometry (dofhandler interpolation + warps)
#     vs the flat triangle, relative to the grid's bounding-box diagonal;
#   - solution error: the exact field polynomial vs the linear vertex-color
#     interpolation, relative to the field's value span.
# The camera is deliberately not an input: keys recompute only when the
# solution, a warp or a tolerance changes. Refinement is always conforming
# (watertight) — the vertices are duplicated across cell boundaries anyway,
# which is what keeps DG fields plottable, so there is nothing to gain from
# the unconstrained variant a plot could expose.
#
# Everything derived from the dataset alone — the base domain, its adjacency,
# the continuous geometry mapping, the warp and field evaluators with their
# coefficient buffers — lives in one `IsubdSubstrate`, built lazily and cached
# on the FEData (`_substrate`), shared by every adaptive plot of it. The
# per-plot part is the key buffer and the decode buffers.
#
# Graph topology per plot (all nodes prefixed subd_ to avoid clashes):
#
#   ds_u, geometry_tol, solution_tol, max_depth ──> subd_keys ─┐
#   warp scales and fields ──┘                                 ├─> subd_positions, subd_ξ,
#   ds_u, warp scales ─────────────────────────────────────────┘   subd_faces, subd_color
#
# The persistent LEB key buffer lives in plot-local state captured by the
# subd_keys computation; the node returns `nothing` when an update does not
# change the key set, so the decode is skipped. Geometry, faces and colors
# are emitted by ONE edge, so they can never be resolved against different
# key sets (mismatched buffers rendered as dropped triangles).

# One dof field, evaluated at arbitrary (cell, ξ) — the continuous-geometry
# and per-vertex-colour workhorse, and the hot path of everything adaptive
# (the error estimators query it far more often than the rendering does).
#
# The value is summed straight from the reference shape functions, the way
# `geometric_map` does, rather than through `PointValues`: for H1 (Lagrange)
# fields the value mapping is the identity, so `reinit!`'s machinery — which
# evaluates every shape function into a buffer for `function_value` to read
# back — is pure overhead, worth a measured 3x. Fields whose values are
# mapped (Piola, i.e. H(curl)/H(div)) would need the real thing; they are
# unsupported here as in the rest of the package (#151).
#
# `Tv` is the sample type the coefficient buffers hold, Float32 by default:
# it is what the renderer draws (GLMakie uploads Float32), so sampling in it
# means the estimators measure the field the viewer actually sees — and the
# matching eps(Tv)-scaled tolerance floors keep refinement from chasing
# deviations below that resolution. Pass Float64 where evaluation itself is
# the point (tests, diagnostics).
struct FieldEvaluator{Tv,IP,PF}
    ips::Vector{IP}                   # one interpolation per subdofhandler with the field
    sdh_of_cell::Vector{Int}          # cell -> index into ips, 0 when the field is absent
    # cell -> global dofs of the field, an empty view when absent. CSR-style
    # (one flat data vector + offsets) so the hot loops read contiguously and
    # the map uploads to a device as two plain buffers.
    celldofs_field::Ferrite.ArrayOfVectorViews{Int,1}
    ncomps::Int
    # Monomial coefficients per cell, when the interpolation admits them (see
    # polyeval.jl). `prepare!` fills them for the cells that will be sampled;
    # `nothing` means every evaluation sums shape functions instead.
    poly::PF
end

function FieldEvaluator(dh::Ferrite.DofHandler, field::Symbol,
                        ::Type{Tv} = Float32) where {Tv<:AbstractFloat}
    sdhs = getsubdofhandlers(dh, field)
    isempty(sdhs) && error("field :$field not found in the DofHandler")
    ip_field = Ferrite.getfieldinterpolation(first(sdhs), field)
    ξ0 = Ferrite.Vec(ntuple(d -> 0.0, Ferrite.getrefdim(ip_field)))
    ncomps = length(Ferrite.reference_shape_value(ip_field, ξ0, 1))
    ips = [Ferrite.getfieldinterpolation(sdh, field) for sdh in sdhs]
    ncells = Ferrite.getncells(Ferrite.get_grid(dh))
    sdh_of_cell = zeros(Int, ncells)
    sizehint = length(Ferrite.dof_range(first(sdhs), field))
    celldofs_field = Ferrite.ArrayOfVectorViews(Int[], (ncells,); sizehint) do buf
        for (si, sdh) in enumerate(sdhs)
            rng = Ferrite.dof_range(sdh, field)
            for cell_idx in sdh.cellset
                sdh_of_cell[cell_idx] = si
                for d in view(Ferrite.celldofs(dh, cell_idx), rng)
                    Ferrite.push_at_index!(buf, d, cell_idx)
                end
            end
        end
    end
    # only a single interpolation can share one coefficient layout
    basis = length(ips) == 1 ? PolyBasis(only(ips)) : nothing
    T = ncomps == 1 ? Tv : Tensors.Vec{ncomps,Tv}
    poly = basis === nothing ? nothing : PolyField(basis, ncells, T)
    return FieldEvaluator{Tv,eltype(ips),typeof(poly)}(ips, sdh_of_cell, celldofs_field,
                                                       ncomps, poly)
end

_sample_type(::FieldEvaluator{Tv}) where {Tv} = Tv

# Fill the coefficients of the cells that will be evaluated. Cheap (one small
# matvec per cell) next to the sampling that follows, and skipped entirely
# when the field has no polynomial form.
prepare!(::FieldEvaluator{<:Any,<:Any,Nothing}, cells, u::AbstractVector) = nothing
function prepare!(ev::FieldEvaluator, cells, u::AbstractVector)
    pf = ev.poly
    nodal = _nodal_buffer(pf)
    for cell in cells
        ev.sdh_of_cell[cell] == 0 && continue
        dofs = ev.celldofs_field[cell]
        _gather_nodal!(nodal, dofs, u, ev.ncomps)
        refresh_cell!(pf, cell, nodal)
    end
    return nothing
end

_nodal_buffer(pf::PolyField{refdim,N,P,T}) where {refdim,N,P,T} = Vector{T}(undef, N)

# Epoch-guarded refresh: at most one `prepare!` per evaluator per solution
# epoch (see `IsubdSubstrate.epoch`), however many plots and graph nodes
# sample the evaluator. All callers pass the substrate's `used_cells`, so a
# skipped refresh never means missing cells.
function _refresh!(ev::FieldEvaluator, cells, u::AbstractVector, epoch::Int)
    pf = ev.poly
    pf === nothing && return nothing
    pf.epoch[] == epoch && return nothing
    prepare!(ev, cells, u)
    pf.epoch[] = epoch
    return nothing
end

function _gather_nodal!(nodal::Vector{T}, dofs, u, ::Int) where {T<:AbstractFloat}
    @inbounds for k in eachindex(nodal)
        nodal[k] = u[dofs[k]]
    end
    return nodal
end
function _gather_nodal!(nodal::Vector{Tensors.Vec{vdim,T}}, dofs, u,
                        ::Int) where {vdim,T<:AbstractFloat}
    @inbounds for k in eachindex(nodal)
        o = (k - 1) * vdim
        nodal[k] = Tensors.Vec{vdim,T}(ntuple(c -> u[dofs[o + c]], vdim))
    end
    return nodal
end

# Evaluate at one reference point of one cell; `nothing` outside the field's
# subdomain.
@inline function evaluate_at(ev::FieldEvaluator, cell_idx::Int, ξ, u::AbstractVector)
    si = ev.sdh_of_cell[cell_idx]
    si == 0 && return nothing
    pf = ev.poly
    if pf !== nothing && pf.filled[cell_idx]
        return evaluate(pf, cell_idx, ξ)
    end
    return _sum_shape_values(ev.ips[si], ev.celldofs_field[cell_idx], ξ, u)
end

# function barrier: the interpolation is only abstractly typed in the vector
# above, and this loop must specialize on it
function _sum_shape_values(ip, dofs::AbstractVector{Int}, ξ, u::AbstractVector)
    val = Ferrite.reference_shape_value(ip, ξ, 1) * u[dofs[1]]
    @inbounds for i in 2:length(dofs)
        val += Ferrite.reference_shape_value(ip, ξ, i) * u[dofs[i]]
    end
    return val
end

# A vectorized interpolation's basis functions are the scalar ones repeated
# once per component, and asking for each of them re-evaluates the scalar
# underneath — so a displacement field costs `vdim` times what it needs to.
# Evaluate the scalar basis once instead and gather the components (Ferrite
# numbers them consecutively per scalar base function).
function _sum_shape_values(ip::Ferrite.VectorizedInterpolation{vdim}, dofs::AbstractVector{Int},
                           ξ, u::AbstractVector) where {vdim}
    sip = ip.ip
    val = zero(Tensors.Vec{vdim,Float64})
    @inbounds for i in 1:Ferrite.getnbasefunctions(sip)
        N = Ferrite.reference_shape_value(sip, ξ, i)
        o = (i - 1) * vdim
        val += N * Tensors.Vec{vdim}(ntuple(c -> u[dofs[o + c]], vdim))
    end
    return val
end

# Global identity of every tessellation vertex of a cell, so that shared base
# edges can be matched *exactly* — by integer id, never by comparing floating
# point coordinates. A vertex coinciding with one of the cell's geometric nodes
# takes that node's global id (which is how two cells recognize their shared
# edge, and how two facets of one cell recognize theirs); anything else is
# cell-interior (a fan centre) and gets a fresh negative id that nobody else
# can collide with.
function _vertex_gids(cell, coords, counter::Base.RefValue{Int})
    refcoords = Ferrite.reference_coordinates(Ferrite.geometric_interpolation(typeof(cell)))
    nodes = cell.nodes
    return map(coords) do ξ
        j = findfirst(rc -> isapprox(rc, ξ; atol=1e-12), refcoords)
        j === nothing ? (counter[] -= 1) : Int(nodes[j])
    end
end

# The isubd base domain of a dataset: one LEB-ordered corner triple per base
# triangle, the triangle -> cell map, and a global vertex id per corner (for
# the adjacency table). Static per plot, like the visibility mask.
#
# In 2D the base is built as a fan from each cell's centre over its *element
# edges*, so a base triangle's split edge is always an element edge. That is
# what makes the base compatible in the sense conforming refinement needs:
# element edges are shared by exactly two base triangles which both treat them
# as their split edge (a "diamond"), while the fan's interior edges are legs on
# both sides. For quadrilaterals this reproduces the existing centre fan; for
# triangles it replaces the single base triangle, whose split edge (the longest
# reference edge) would generally meet a neighbour's leg.
#
# In 3D the same construction applies per *surface facet*: a facet is on the
# surface when the cell across it is missing or not part of the body (see
# `FEData.solid`), and it is fanned from its centre, so a base triangle's split
# edge is again an element edge — now shared by exactly the two surface facets
# meeting there. Note this draws strictly less than the static path, which
# tessellates every facet of every visible cell including the ones buried
# inside the body.
# Is this facet of the cell part of the drawn surface? It is when nothing sits
# across it, or what sits there is not part of the body (an interior cell of a
# solid mesh is, a cell a CrinkleClip removed is not). Without topology every
# facet is taken, which is the 2D case and the pre-topology fallback.
function _is_surface_facet(ds::FEData, cell_id::Int, facet::Int)
    ds.topology === nothing && return true
    neighbours = ds.topology.face_face_neighbor[cell_id, facet]
    isempty(neighbours) && return true
    return !ds.solid[first(neighbours)[1]]
end

# One centre fan of a convex polygon given by its rim vertices in order:
# (v[i+1], centre, v[i]) keeps the rim's winding while making the rim edge the
# triangle's split edge, which is what the conforming refinement needs.
function _push_fan!(corners, cornergids, cellmap, cell_id, rim, rimgids, counter)
    centre = sum(rim) / length(rim)
    centre_gid = (counter[] -= 1)
    for i in eachindex(rim)
        j = mod1(i + 1, length(rim))
        push!(corners, (rim[j], centre, rim[i]))
        push!(cornergids, (rimgids[j], centre_gid, rimgids[i]))
        push!(cellmap, cell_id)
    end
    return nothing
end

# The base domain: cell fans normally, Voronoi-region fans on a
# quadrature-point partition. Returns corners, corner gids, and the per-base-
# triangle maps (cell, sharing group, quadrature point, element-edge mask).
function _isubd_base_triangles(ds::FEData)
    return ds.qp_partition === nothing ? _isubd_base_cells(ds) : _isubd_base_qp(ds, ds.qp_partition)
end

function _isubd_base_cells(ds::FEData)
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    refdim = Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(first(cells))))
    corners = NTuple{3,Ferrite.Vec{refdim,Float64}}[]
    cornergids = NTuple{3,Int}[]
    cellmap = Int[]
    counter = Ref(0)
    for (cell_id, cell) in enumerate(cells)
        ds.visible[cell_id] || continue
        Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(cell))) == refdim ||
            error("adaptive tessellation requires a single reference dimension across the grid")
        refshape = getrefshape(cell)
        tess = reference_tessellation(refshape)
        isempty(tess.triangles) && continue     # e.g. line cells carry no surface
        if refdim == 2 && !isempty(tess.edges)
            gids = _vertex_gids(cell, tess.coords, counter)
            # the element edges, in whatever order the tessellation lists them
            rim = unique(Iterators.flatten(tess.edges))
            centre = sum(tess.coords[i] for i in rim) / length(rim)
            centre_gid = (counter[] -= 1)
            for (a, b) in tess.edges
                # (b, centre, a) keeps the fan's winding while making the
                # element edge (b, a) the triangle's split edge
                push!(corners, (tess.coords[b], centre, tess.coords[a]))
                push!(cornergids, (gids[b], centre_gid, gids[a]))
                push!(cellmap, cell_id)
            end
        elseif refdim == 3
            vcoords = Ferrite.reference_coordinates(Ferrite.Lagrange{refshape,1}())
            vgids = _vertex_gids(cell, vcoords, counter)
            for (f, face) in enumerate(Ferrite.reference_faces(refshape))
                _is_surface_facet(ds, cell_id, f) || continue
                _push_fan!(corners, cornergids, cellmap, cell_id,
                           [vcoords[v] for v in face], [vgids[v] for v in face], counter)
            end
        else
            gids = _vertex_gids(cell, tess.coords, counter)
            for tri in tess.triangles
                c = leb_order((tess.coords[tri[1]], tess.coords[tri[2]], tess.coords[tri[3]]))
                # recover the permutation leb_order applied, to keep the ids aligned
                perm = map(x -> findfirst(i -> tess.coords[i] === x, tri), c)
                push!(corners, c)
                push!(cornergids, (gids[tri[perm[1]]], gids[tri[perm[2]]], gids[tri[perm[3]]]))
                push!(cellmap, cell_id)
            end
        end
    end
    isempty(corners) && error("adaptive tessellation found no visible cells to tessellate")
    # sharing per cell, no quadrature points, every split edge an element edge
    return corners, cornergids, cellmap, cellmap, zeros(Int, length(cellmap)),
           fill(true, length(cellmap))
end

# The base of a quadrature-point partition: every Voronoi region of every
# drawn face, fanned from its centroid so the region's rim — the piecewise-
# constant jumps live there — consists of split edges, which longest-edge
# bisection subdivides but never crosses. Region-region rims meet as split
# edges on both sides (a diamond), so the conforming machinery applies
# unchanged; each region is its own vertex-sharing group, keeping the jumps
# crisp under refinement.
#
# Corner identity: a rim vertex is a geometric node (its global id), a cut
# point on an *element edge* (an id keyed by the edge's node pair and the
# parametric position along it, so the two cells sharing the edge — or the
# two surface facets meeting there in 3D — agree on it whenever their rules
# cut at the same spots; differing rules degrade to an unpaired boundary,
# which costs conformity along that edge and nothing else), or interior (a
# per-cell id from the clipped coordinates, which both adjacent regions of
# one face compute to well below the rounding).
function _isubd_base_qp(ds::FEData, qp::QPPartition)
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    refdim = Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(first(cells))))
    corners = NTuple{3,Ferrite.Vec{refdim,Float64}}[]
    cornergids = NTuple{3,Int}[]
    cellmap = Int[]
    groupmap = Int[]
    qpmap = Int[]
    edgemask = Bool[]
    counter = Ref(0)
    edge_cuts = Dict{Tuple{Int,Int,Float64},Int}()
    regions_cache = Dict{Type,Any}()
    ngroups = 0
    for (cell_id, cell) in enumerate(cells)
        ds.visible[cell_id] || continue
        Ferrite.getrefdim(Ferrite.geometric_interpolation(typeof(cell))) == refdim ||
            error("adaptive tessellation requires a single reference dimension across the grid")
        refshape = getrefshape(cell)
        regions = get!(() -> _qp_face_regions(refshape, _qr_for(qp.qr, refshape)),
                       regions_cache, typeof(cell))
        isempty(regions) && continue                 # e.g. line cells carry no surface
        refcorners = Ferrite.reference_coordinates(Ferrite.Lagrange{refshape,1}())
        refedges = Ferrite.reference_edges(refshape)
        local_pool = Dict{NTuple{refdim,Float64},Int}()
        for (fi, qpi, poly) in regions
            refdim == 3 && !_is_surface_facet(ds, cell_id, fi) && continue
            rimgids = [_qp_vertex_gid(ξ, refcorners, refedges, cell.nodes, local_pool,
                                      edge_cuts, counter) for ξ in poly]
            ngroups += 1
            _push_fan!(corners, cornergids, cellmap, cell_id, poly, rimgids, counter)
            append!(groupmap, fill(ngroups, length(poly)))
            append!(qpmap, fill(qpi, length(poly)))
            for i in eachindex(poly)
                j = mod1(i + 1, length(poly))
                push!(edgemask, _on_element_edge((poly[i] + poly[j]) / 2, refcorners, refedges))
            end
        end
    end
    isempty(corners) && error("adaptive tessellation found no visible cells to tessellate")
    return corners, cornergids, cellmap, groupmap, qpmap, edgemask
end

# Nearest classification of a reference point: node id, element-edge cut id,
# or cell-local interior id (see `_isubd_base_qp`).
function _qp_vertex_gid(ξ, refcorners, refedges, nodes, local_pool, edge_cuts, counter)
    j = findfirst(rc -> isapprox(rc, ξ; atol=1e-9), refcorners)
    j === nothing || return Int(nodes[j])
    for (a, b) in refedges
        A, B = refcorners[a], refcorners[b]
        AB = B - A
        t = ((ξ - A) ⋅ AB) / sum(abs2, AB)
        1e-9 < t < 1 - 1e-9 || continue
        sum(abs2, ξ - (A + t * AB)) <= 1e-18 || continue
        ga, gb = Int(nodes[a]), Int(nodes[b])
        key = ga < gb ? (ga, gb, round(t; digits=9) + 0.0) : (gb, ga, round(1 - t; digits=9) + 0.0)
        return get!(() -> (counter[] -= 1), edge_cuts, key)
    end
    key = ntuple(d -> round(ξ[d]; digits=9) + 0.0, length(ξ))
    return get!(() -> (counter[] -= 1), local_pool, key)
end

function _on_element_edge(ξ, refcorners, refedges)
    for (a, b) in refedges
        A, B = refcorners[a], refcorners[b]
        AB = B - A
        t = clamp(((ξ - A) ⋅ AB) / sum(abs2, AB), 0.0, 1.0)
        sum(abs2, ξ - (A + t * AB)) <= 1e-18 && return true
    end
    return false
end

# Pair up base triangles along shared edges. Only exact 2-triangle matches
# become neighbours, and only when both sides agree on the edge's role: a split
# edge may pair with a split edge and a leg with a leg, never across. Anything
# else stays a boundary, which costs conformity along that edge and nothing
# else (see `key_neighbour`).
function _base_adjacency(cornergids::Vector{NTuple{3,Int}})
    edges = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int,Bool}}}()
    for (t, g) in enumerate(cornergids)
        for (e, (i, j)) in enumerate(((1, 3), (1, 2), (2, 3)))   # EDGE_S, EDGE_L, EDGE_R
            a, b = g[i], g[j]
            forward = a < b
            push!(get!(Vector{Tuple{Int,Int,Bool}}, edges, forward ? (a, b) : (b, a)),
                  (t, e, forward))
        end
    end
    adjacency = [(NO_NEIGHBOR, NO_NEIGHBOR, NO_NEIGHBOR) for _ in cornergids]
    incompatible = 0
    for entries in values(edges)
        length(entries) == 2 || continue                     # boundary, or a non-manifold edge
        (t1, e1, f1), (t2, e2, f2) = entries
        if (e1 == EDGE_S) != (e2 == EDGE_S)
            incompatible += 1
            continue
        end
        reversed = f1 != f2
        adjacency[t1] = Base.setindex(adjacency[t1], (t2, e2, reversed), e1)
        adjacency[t2] = Base.setindex(adjacency[t2], (t1, e1, reversed), e2)
    end
    return adjacency, incompatible
end

# One upstream WarpByVector, ready to evaluate: the `Deformation` record with
# the FieldEvaluator built for its field — against the dof handler and
# solution *captured at warp time*, so a later Gradient rebinding the
# dataset's handler does not orphan the warp. The evaluator sits in a Ref
# because switching the warp's field observable rebuilds it; the Ref is typed
# so the mapping's hot loop stays concretely dispatched, which is also why a
# switch to a differently interpolated field cannot be followed.
struct WarpEval{EV<:FieldEvaluator,UO<:Makie.Observable,SO<:Makie.Observable}
    dh::Ferrite.AbstractDofHandler   # cold: only touched on a rebuild
    u::UO
    name::Makie.Observable{Symbol}
    scale::SO
    ev::Base.RefValue{EV}
    built_for::Base.RefValue{Symbol}
end

function WarpEval(d::Deformation, ::Type{Tv} = Float32) where {Tv<:AbstractFloat}
    ev = FieldEvaluator(d.dh, d.field[], Tv)
    return WarpEval{typeof(ev),typeof(d.u),typeof(d.scale)}(
        d.dh, d.u, d.field, d.scale, Ref(ev), Ref(d.field[]))
end

function _refresh_warp!(w::WarpEval{EV}, cells, epoch::Int) where {EV}
    if w.name[] !== w.built_for[]
        newev = FieldEvaluator(w.dh, w.name[], _sample_type(w.ev[]))
        newev isa EV ||
            error("switching the warp field from :$(w.built_for[]) to :$(w.name[]) changes the " *
                  "field's interpolation, which an adaptive plot cannot follow — recreate the plot")
        w.ev[] = newev
        w.built_for[] = w.name[]
    end
    _refresh!(w.ev[], cells, w.u[], epoch)
    return nothing
end

# Continuous geometry x(ξ) [+ Σ scaleᵢ · fieldᵢ(ξ) for upstream warps] of one
# cell, as the isubd mapping closure. Reads the current solution and warp
# scales non-reactively — reactivity is the graph nodes' concern (they list
# the warp observables as inputs, see `_warp_inputs!`). Evaluated in the
# substrate's sample type `T` (Float32 by default — what the renderer draws),
# so the geometry estimator measures the surface the viewer actually sees;
# the eps(T)-scaled tolerance floors keep it from chasing deviations below
# that resolution.
function _isubd_mapping(ds::FEData{dim}, cellmap::Vector{Int}, cellcoords,
                        warps, ::Type{T}) where {dim,T}
    grid = Ferrite.get_grid(ds.dh)
    cells = Ferrite.getcells(grid)
    gips = [Ferrite.geometric_interpolation(typeof(c)) for c in cells]
    # The geometry's own coefficients never change — the node coordinates are
    # fixed — so they are computed once here rather than per update.
    geo = _geometry_poly(gips, cellcoords, cellmap, Val(dim), T)
    function mapping(base_id::Int, ξ)
        cell_id = cellmap[base_id]
        x = (geo !== nothing && geo.filled[cell_id]) ? evaluate(geo, cell_id, ξ) :
            geometric_map(gips[cell_id], cellcoords[cell_id], ξ)
        for w in warps
            d = evaluate_at(w.ev[], cell_id, ξ, w.u[])
            d === nothing && continue
            x += T(w.scale[]) * d
        end
        # the warp fallback path may accumulate in Float64; pin the closure's
        # return type to the sample type either way
        return Tensors.Vec{dim,T}(Tuple(x))
    end
    return mapping
end

function _geometry_poly(gips, cellcoords, cellmap, ::Val{dim}, ::Type{T}) where {dim,T}
    basis = PolyBasis(first(gips))
    basis === nothing && return nothing
    poly = PolyField(basis, length(gips), Tensors.Vec{dim,T})
    for cell in unique(cellmap)
        # a differently interpolated cell keeps the shape-function path
        PolyBasis(gips[cell]) === nothing && continue
        length(cellcoords[cell]) == size(poly.coeffs, 1) || continue
        refresh_cell!(poly, cell, cellcoords[cell])
    end
    return poly
end

# Positions leave the graph as Float32 points for rendering (a plain copy
# when the substrate already samples in Float32, the default). Filled into a
# buffer that is handed out as is, see `_adaptive_solutionplot!`.
function _render_positions!(out::Vector{GeometryBasics.Point{dim,Float32}}, ps) where {dim}
    resize!(out, length(ps))
    @inbounds for i in eachindex(ps)
        out[i] = GeometryBasics.Point{dim,Float32}(ps[i]...)
    end
    return out
end

# Everything the adaptive path derives from the dataset alone — independent
# of any plot's tolerances: the base domain with its adjacency and continuous
# mapping, the warp evaluators, and one FieldEvaluator per sampled dof field
# with its coefficient buffers. Built once per FEData (lazily, by
# `_substrate`) and shared by every adaptive plot of it; the per-plot state is
# the key buffer and the decode buffers.
struct IsubdSubstrate{B<:IsubdBase,CC,W<:Vector}
    base::B
    cellmap::Vector{Int}                # base triangle -> cell
    # base triangle -> vertex-sharing group. Vertices are shared within a
    # group and duplicated across, so the drawn field may jump at group
    # boundaries: per cell normally (inter-element jumps of DG fields), per
    # Voronoi region on a quadrature-point partition (piecewise-constant
    # regions stay crisp under refinement).
    groupmap::Vector{Int}
    qpmap::Vector{Int}                  # base triangle -> quadrature point (0 without a partition)
    # base triangle -> does its split edge lie on an *element* edge? Always
    # true for cell fans (their split edges are element edges by
    # construction); on a quadrature-point partition the region rims also cut
    # through cell interiors, and the meshplot wireframe draws only the
    # element-edge part.
    edgemask::Vector{Bool}
    cellcoords::CC
    warps::W                            # WarpEvals, in application order
    used_cells::Vector{Int}             # cells any base triangle samples
    diag::Float64                       # grid bounding-box diagonal
    evaluators::Dict{Symbol,FieldEvaluator}  # per sampled dof field, lazily
    # Solution epoch: bumped whenever the dataset's (or a warp stage's) dof
    # vector, or a warp scale or field, changes. `_refresh!` uses it to run
    # each evaluator's `prepare!` at most once per update — not once per graph
    # node per plot — and `_dev_cache` to invalidate the deviation memos.
    epoch::Base.RefValue{Int}
    # Per-term deviation memos (`:geometry`, or a field name), shared by every
    # plot of this dataset and handed to DeviationLoD. A deviation is a pure
    # function of (key, term, epoch) — it contains no tolerance and no
    # refinement state — so within one epoch the second plot's refinement, a
    # re-tolerated first plot, or the wireframe next to the surface all decide
    # from lookups (~2ns) instead of re-sampling the fields (~500ns/key).
    dev_caches::Dict{Symbol,Dict{UInt64,Float64}}
    dev_epoch::Base.RefValue{Int}       # epoch the memos are valid for
end

# The id bookkeeping of the base build works on exact Float64 reference
# coordinates; the corners the pipeline then samples at are rounded once into
# the sample type here.
function _corners_as(corners::Vector{NTuple{3,Ferrite.Vec{refdim,Float64}}},
                     ::Type{T}) where {refdim,T}
    T === Float64 && return corners
    return [ntuple(i -> Tensors.Vec{refdim,T}(Tuple(tri[i])), 3) for tri in corners]
end

function _build_substrate(ds::FEData{dim}, ::Type{T}) where {dim,T}
    corners, cornergids, cellmap, groupmap, qpmap, edgemask = _isubd_base_triangles(ds)
    adjacency, _ = _base_adjacency(cornergids)
    grid = Ferrite.get_grid(ds.dh)
    cellcoords = [[Tensors.Vec{dim,T}(Tuple(x)) for x in Ferrite.getcoordinates(grid, i)]
                  for i in 1:Ferrite.getncells(grid)]
    warps = [WarpEval(d, T) for d in ds.deformation]
    mapping = _isubd_mapping(ds, cellmap, cellcoords, warps, T)
    base = IsubdBase(_corners_as(corners, T), mapping, adjacency)
    epoch = Ref(0)
    bump(_) = (epoch[] += 1; nothing)
    Makie.on(bump, ds.u)
    for w in warps
        w.u === ds.u || Makie.on(bump, w.u)
        # the deviation memos also depend on the warp's scale and field (the
        # coefficients do not, but one shared epoch is simpler than two, and
        # an occasional redundant prepare! is cheap)
        Makie.on(bump, w.scale)
        Makie.on(bump, w.name)
    end
    return IsubdSubstrate(base, cellmap, groupmap, qpmap, edgemask, cellcoords, warps,
                          unique(cellmap), _grid_diagonal(grid), Dict{Symbol,FieldEvaluator}(),
                          epoch, Dict{Symbol,Dict{UInt64,Float64}}(), Ref(-1))
end

# The dataset's substrate, built on first use. Deliberately behind a
# `Ref{Any}` on the FEData (its concrete type would otherwise leak into the
# dataset's type parameters); `_adaptive_solutionplot!`/`_adaptive_wireframe!`
# immediately pass the result through a function barrier, so the untypedness
# costs one dynamic dispatch per plot creation.
function _substrate(ds::FEData)
    cached = ds.subd_cache[]
    cached === nothing || return cached::IsubdSubstrate
    sub = _build_substrate(ds, ds.sample_type)
    ds.subd_cache[] = sub
    return sub
end

# The sample type every evaluation of this substrate runs in, recovered from
# the base's reference-corner type.
_sample_type(sub::IsubdSubstrate) = eltype(eltype(eltype(sub.base.corners)))

# The substrate-wide evaluator of one dof field, shared (with its coefficient
# buffers) by every plot sampling that field.
_field_evaluator(sub::IsubdSubstrate, dh, name::Symbol) =
    get!(() -> FieldEvaluator(dh, name, _sample_type(sub)), sub.evaluators, name)

# The deviation memo of one criterion term, valid for the current epoch. The
# staleness check clears *all* terms at once (they share the epoch), lazily at
# the first request of a new epoch — so the memos only ever hold deviations of
# the solution state the requesting node is about to refine against.
function _dev_cache(sub::IsubdSubstrate, term::Symbol)
    if sub.dev_epoch[] != sub.epoch[]
        for c in values(sub.dev_caches)
            empty!(c)
        end
        sub.dev_epoch[] = sub.epoch[]
    end
    return get!(Dict{UInt64,Float64}, sub.dev_caches, term)
end

# Refresh the per-cell coefficients of everything a node is about to sample.
# The epoch guard makes repeated calls within one update free, so this runs
# defensively at the top of every sampling node rather than once somewhere
# central — coefficients that do not match the `u` in hand would evaluate
# silently wrong values.
function _refresh_all!(sub::IsubdSubstrate, ev, u::AbstractVector)
    epoch = sub.epoch[]
    for w in sub.warps
        _refresh_warp!(w, sub.used_cells, epoch)
    end
    ev === nothing || _refresh!(ev, sub.used_cells, u, epoch)
    return nothing
end

# One leaf of a derivation chain, ready to evaluate: the dof field's evaluator
# with the solution it is defined on (which need not be the plotted dataset's
# — cf. WarpEval), plus the deviation-memo key its criterion term uses.
struct SourceEval{EV<:FieldEvaluator,UO<:Makie.Observable}
    ev::EV
    u::UO
    name::Symbol
    term::Symbol
end

# A pointwise-derived color, resolved from a `DerivedPointData` record: `fun`
# re-applies the derivation closures at arbitrary (cell, ξ), and `sources`
# names the dof-field leaves. Following the regularity principle — a derived
# quantity is at most as regular as its sources — the refinement criterion
# samples the *sources* (their spans, their deviations, their memos), never
# the derived quantity itself: the mesh follows u, and the colors follow the
# mesh.
struct ChainEvaluator{S<:Tuple,F}
    sources::S
    fun::F        # (cell, ξ) -> derived value; `nothing` outside a subdomain
end

# FieldEvaluator hands back scalars and flat component vectors; the derivation
# closures expect what the static path's `_wrap_row` gives them (a gradient
# row as a Tensor{2}, ...). Matrixized interpolations already evaluate to
# tensors and pass through.
_wrap_value(v::Number, ::Val) = v
_wrap_value(v::Tensors.Vec, ::Val{sdim}) where {sdim} = _wrap_row(v, sdim)
_wrap_value(v, ::Val) = v

function _chain_node(sub::IsubdSubstrate, ds::FEData, src::FieldSource, ::Val{sdim}) where {sdim}
    local se
    if src.dh === ds.dh
        se = SourceEval(_field_evaluator(sub, ds.dh, src.name), src.u, src.name, src.name)
    else
        # a source rebound away by a later filter: private evaluator, and a
        # memo key that cannot collide with a field of the plotted handler
        term = Symbol(src.name, :_, string(objectid(src.dh); base=16))
        se = SourceEval(FieldEvaluator(src.dh, src.name, _sample_type(sub)), src.u, src.name, term)
        # its solution must invalidate the shared epoch too (duplicate
        # listeners from several plots only advance the counter faster)
        src.u === ds.u || Makie.on(_ -> (sub.epoch[] += 1; nothing), src.u)
    end
    leaf = function (cell::Int, ξ)
        v = evaluate_at(se.ev, cell, ξ, se.u[])
        return v === nothing ? nothing : _wrap_value(v, Val(sdim))
    end
    return (se,), leaf
end

function _chain_node(sub::IsubdSubstrate, ds::FEData, rec::DerivedPointData, v::Val)
    children = map(inp -> _chain_node(sub, ds, inp, v), Tuple(rec.inputs))
    sources = reduce((acc, c) -> (acc..., c[1]...), children; init=())  # flatten the leaf tuples
    funs = map(c -> c[2], children)
    f = rec.f
    node = function (cell::Int, ξ)
        vals = map(g -> g(cell, ξ), funs)
        any(x -> x === nothing, vals) && return nothing
        return f(vals...)
    end
    return sources, node
end

function _chain_evaluator(sub::IsubdSubstrate, ds::FEData{dim}, rec::DerivedPointData) where {dim}
    sources, fun = _chain_node(sub, ds, rec, Val(dim))
    return ChainEvaluator(sources, fun)
end

_refresh_chain!(::IsubdSubstrate, ::Nothing) = nothing
function _refresh_chain!(sub::IsubdSubstrate, ch::ChainEvaluator)
    epoch = sub.epoch[]
    for s in ch.sources
        _refresh!(s.ev, sub.used_cells, s.u[], epoch)
    end
    return nothing
end

# Per-vertex derived values at the decoded reference coordinates — the chain
# counterpart of `_transfer_at!`. The chain's output must be a scalar, which
# `_validate_chain` established at plot creation.
function _transfer_chain_at!(out::Vector{Float32}, ch::ChainEvaluator, cellmap::Vector{Int}, mesh)
    resize!(out, length(mesh.refcoords))
    @inbounds for v in eachindex(mesh.refcoords)
        val = ch.fun(cellmap[mesh.vertex_base[v]], mesh.refcoords[v])
        out[v] = val === nothing ? NaN32 : Float32(val)
    end
    return out
end

# Eager scalar check: evaluate the chain once, at the centroid of the first
# base triangle whose cell carries all sources, and fail at plot creation
# rather than mid-render.
function _validate_chain(ch::ChainEvaluator, sub::IsubdSubstrate, fname::Symbol)
    _refresh_chain!(sub, ch)
    for b in 1:length(sub.base.corners)
        c = sub.base.corners[b]
        val = ch.fun(sub.cellmap[b], (c[1] + c[2] + c[3]) / 3)
        val === nothing && continue
        val isa Number && return nothing
        error("adaptive coloring by the derived :$fname needs a scalar; its record evaluates to " *
              "$(typeof(val)) — reduce it first (e.g. Magnitude, VonMises, ExtractComponent)")
    end
    return nothing   # nowhere defined: rendered as NaN, like the static path
end

# Per-vertex field values at the decoded reference coordinates: the adaptive
# counterpart of transfer_solution, evaluating at the sub-triangle vertices
# instead of the static tessellation vertices.
function _transfer_at!(out::Vector{Float32}, ev::FieldEvaluator, cellmap::Vector{Int},
                       mesh, u::AbstractVector; reduce::Bool)
    resize!(out, length(mesh.refcoords))
    @inbounds for v in eachindex(mesh.refcoords)
        val = evaluate_at(ev, cellmap[mesh.vertex_base[v]], mesh.refcoords[v], u)
        out[v] = val === nothing ? NaN32 :
                 reduce ? Float32(LinearAlgebra.norm(val)) : Float32(val[1])
    end
    return out
end

# The scalar drawn as color, as a (base_id, ξ) probe for the solution-error
# estimator; `u` is bound per key-node invocation.
function _scalar_probe(ev::FieldEvaluator, cellmap::Vector{Int}, u::AbstractVector; reduce::Bool)
    function probe(base_id::Int, ξ)
        cell_id = cellmap[base_id]
        val = evaluate_at(ev, cell_id, ξ, u)
        val === nothing && return 0.0
        return reduce ? Float64(LinearAlgebra.norm(val)) : Float64(val[1])
    end
    return probe
end

# Value span and magnitude of the drawn scalar — the span is the reference
# scale for the relative solution tolerance, the magnitude the reference for
# the sample type's noise floor (a large constant field samples with noise
# proportional to its value while its span is zero) — from the field's dof
# values on the sampled cells. For the nodal (Lagrange) bases supported here
# the dof values are the function's values at the interpolation nodes,
# including the edge/face/interior nodes no base-triangle corner visits.
# Sampling the corners instead once collapsed the reference to noise level
# for a field peaking between them, which then refined a visually flat
# surface to the depth cap.
function _field_span(ev::FieldEvaluator, cells, u::AbstractVector; reduce::Bool)
    lo, hi = Inf, -Inf
    n = ev.ncomps
    for cell in cells
        ev.sdh_of_cell[cell] == 0 && continue
        dofs = ev.celldofs_field[cell]
        if n == 1
            @inbounds for d in dofs
                v = Float64(u[d])
                isfinite(v) || continue
                lo, hi = min(lo, v), max(hi, v)
            end
        else
            # component dofs are consecutive per node (as in _gather_nodal!);
            # the drawn scalar of a vector field is its magnitude
            @inbounds for k in 1:(length(dofs) ÷ n)
                o = (k - 1) * n
                s = 0.0
                for c in 1:n
                    s += abs2(Float64(u[dofs[o + c]]))
                end
                v = reduce ? sqrt(s) : Float64(u[dofs[o + 1]])
                isfinite(v) || continue
                lo, hi = min(lo, v), max(hi, v)
            end
        end
    end
    hi >= lo || return (0.0, 0.0)
    return (hi - lo, max(abs(lo), abs(hi)))
end

# The finest relative tolerance worth resolving in sample type `T`: below
# ~100 ulps a measured "deviation" is the pipeline's own rounding noise, and
# the renderer — which draws `T` — cannot show the difference either.
# Requested tolerances are floored here, so the Float32 default caps
# adaptation at ~1e-5 of the reference scale however small a tolerance is
# asked for. (Grids offset far from the origin relative to their size carry
# coordinate noise beyond any diagonal-scaled floor — those need
# `sample_type = Float64`.)
_tol_floor(::Type{T}) where {T<:AbstractFloat} = 100 * Float64(eps(T))

# The tolerance of one solution term: the requested fraction of the span,
# floored at the noise level of the field's magnitude. A (near-)constant
# field never asks for refinement.
_solution_tol(rel::Float64, (span, mag)::NTuple{2,Float64}, ::Type{T}) where {T} =
    max(rel * span, _tol_floor(T) * mag)

function _grid_diagonal(grid)
    nodes = Ferrite.getnodes(grid)
    lo = hi = Ferrite.get_node_coordinate(first(nodes))
    for n in nodes
        x = Ferrite.get_node_coordinate(n)
        lo, hi = min.(lo, x), max.(hi, x)
    end
    return LinearAlgebra.norm(hi - lo)
end

# The wireframe of an adaptively refined mesh: the element edges, subdivided
# exactly as the surface subdivided them.
#
# An element edge *is* a base triangle's split edge (that is what the fan
# construction guarantees), and a point lies on it precisely when its
# barycentric weight for the apex vanishes — which `key_xform` gives exactly,
# since the bisection weights are dyadic. So a leaf contributes a wireframe
# segment when two of its corners have zero apex weight, and conformity
# guarantees the two triangles sharing an element edge subdivide it
# identically. Each edge is therefore emitted once, by the base triangle with
# the smaller id of the pair.
function _element_edge_segments!(out::Vector{PT}, keys::Vector{UInt64}, base::IsubdBase,
                                 edgemask::Vector{Bool}) where {PT}
    empty!(out)
    for k in keys
        b = key_base(k)
        edgemask[b] || continue                 # a Voronoi rim, not an element edge
        nb = base.adjacency[b][EDGE_S]
        (nb[1] != 0 && nb[1] < b) && continue   # the partner draws this one
        X = key_xform(k)
        on = (X[2, 1] == 0.0, X[2, 2] == 0.0, X[2, 3] == 0.0)
        (count(on) == 2) || continue            # a leaf touches the edge with at most one of its own
        c = key_corners(base, k)
        for (i, j) in ((1, 2), (2, 3), (3, 1))
            (on[i] && on[j]) || continue
            push!(out, PT(base.mapping(b, c[i])...), PT(base.mapping(b, c[j])...))
        end
    end
    return out
end

# The warps' scale and field observables as graph inputs, so a slider-driven
# warp re-refines and re-decodes the way a solution update does. The mapping
# closure reads their current values itself; these inputs exist to make the
# nodes rerun (without them an adaptive plot silently kept the stale geometry
# until the next solution update flushed it through).
function _warp_inputs!(graph, sub::IsubdSubstrate)
    names = Symbol[]
    for (i, w) in enumerate(sub.warps)
        s = Symbol(:subd_warpscale_, i)
        n = Symbol(:subd_warpfield_, i)
        ComputePipeline.add_input!(graph, s, w.scale)
        ComputePipeline.add_input!(graph, n, w.name)
        push!(names, s, n)
    end
    return names
end

# The adaptive branch of meshplot's plot!: same refinement machinery as the
# surface, but only the geometry criterion — a wireframe has no field to
# resolve, only a curve to follow. The outer function exists as a barrier past
# the dataset's untyped substrate cache.
_adaptive_wireframe!(WF, ds::FEData) = _wire_adaptive_wireframe!(WF, ds, _substrate(ds))

function _wire_adaptive_wireframe!(WF, ds::FEData{dim}, sub::IsubdSubstrate) where {dim}
    graph = WF.attributes
    ComputePipeline.add_input!(graph, :subd_u, ds.u)
    warp_inputs = _warp_inputs!(graph, sub)
    base = sub.base
    state = (keys=root_keys(base), prev=UInt64[], scratch=UInt64[])
    diag = sub.diag
    T = _sample_type(sub)
    ComputePipeline.register_computation!(graph, [:subd_u, :geometry_tol, :max_depth, warp_inputs...],
                                          [:subd_keys]) do inputs, changed, cached
        _refresh_all!(sub, nothing, inputs.subd_u)
        lod = DeviationLoD(base.mapping, max(Float64(inputs.geometry_tol), _tol_floor(T)) * diag,
                           _dev_cache(sub, :geometry))
        refine_keys!(state.keys, state.scratch, base, lod; max_depth=Int(inputs.max_depth))
        cached !== nothing && state.keys == state.prev && return nothing
        copy!(state.prev, state.keys)
        return (state.keys,)
    end
    # Makie draws 2D/3D points; pad 1D grids with a zero y-coordinate
    segments = GeometryBasics.Point{max(dim, 2),Float32}[]
    Makie.map!(graph, [:subd_keys, :subd_u, warp_inputs...], :edge_lines) do keys, u, _warps...
        _refresh_all!(sub, nothing, u)
        return _element_edge_segments!(segments, keys, base, sub.edgemask)
    end
    return nothing
end

# The adaptive branch of solutionplot's plot!. The dataset, its visibility and
# the color resolution are read eagerly (as everywhere in the recipes); the
# solution, the warps and the tolerances drive the graph.
function _adaptive_solutionplot!(SP, ds::FEData)
    sub = _substrate(ds)

    # color resolution (eager): a dof field, a recorded derivation of dof
    # fields, or the dataset's quadrature-point partition, evaluated per
    # vertex — or a plain color
    colorval = SP.color[]
    qp = ds.qp_partition
    if colorval isa Symbol && qp !== nothing && colorval === qp.output
        v0 = qp.extract(_qp_at(qp.values[], 1, 1))
        v0 isa Number ||
            error("adaptive quadrature-point coloring needs a scalar; `extract` yields " *
                  "$(typeof(v0)) — reduce in `extract` (e.g. extract = s -> vonmises(s.σ))")
        return _wire_adaptive_qpplot!(SP, ds, sub, qp)
    end
    ev = nothing
    chain = nothing
    reduce = false
    fname = :none
    if colorval isa Symbol && (colorval === :default || _data_association(ds, colorval) !== :none)
        fname = _resolve_name(ds, colorval)
        if fname in Ferrite.getfieldnames(ds.dh)
            ev = _field_evaluator(sub, ds.dh, fname)
            reduce = colorval === :default && ev.ncomps > 1
            reduce || ev.ncomps == 1 ||
                error("field :$fname has $(ev.ncomps) components; adaptive coloring needs a scalar " *
                      "(or :default, which reduces to the magnitude)")
        elseif haskey(ds.point_derivations, fname)
            chain = _chain_evaluator(sub, ds, ds.point_derivations[fname])
            _validate_chain(chain, sub, fname)
        else
            error("adaptive solutionplot colors by evaluating a dof field — or a derivation of dof " *
                  "fields (Derive, VonMises, ...) — at the refined vertices; :$fname is a raw data " *
                  "array on the static tessellation and cannot be resampled. " *
                  "Color by a dof field, a derived quantity or a plain color, or use adaptive=false.")
        end
    end

    # barrier past the dataset's untyped substrate cache: everything below
    # specializes on the substrate's (and the evaluators') concrete types
    return _wire_adaptive_solutionplot!(SP, ds, sub, ev, chain, fname, reduce)
end

function _wire_adaptive_solutionplot!(SP, ds::FEData{dim}, sub::IsubdSubstrate, ev, chain, fname::Symbol,
                                      reduce::Bool) where {dim}
    graph = SP.attributes
    ComputePipeline.add_input!(graph, :subd_u, ds.u)
    warp_inputs = _warp_inputs!(graph, sub)
    base = sub.base
    cellmap = sub.cellmap
    state = (keys=root_keys(base), scratch=UInt64[], prev=UInt64[])
    diag = sub.diag
    T = _sample_type(sub)
    span = Ref((NaN, NaN))
    spans = chain === nothing ? NTuple{2,Float64}[] : fill((NaN, NaN), length(chain.sources))
    keyinputs = [:subd_u, :geometry_tol, :solution_tol, :max_depth, warp_inputs...]
    ComputePipeline.register_computation!(graph, keyinputs, [:subd_keys]) do inputs, changed, cached
        _refresh_all!(sub, ev, inputs.subd_u)
        _refresh_chain!(sub, chain)
        geo = DeviationLoD(base.mapping, max(Float64(inputs.geometry_tol), _tol_floor(T)) * diag,
                           _dev_cache(sub, :geometry))
        lods = (geo,)
        if ev !== nothing
            probe = _scalar_probe(ev, cellmap, inputs.subd_u; reduce)
            if isnan(span[][1]) || changed.subd_u
                span[] = _field_span(ev, sub.used_cells, inputs.subd_u; reduce)
            end
            tol = _solution_tol(Float64(inputs.solution_tol), span[], T)
            lods = (lods..., DeviationLoD(probe, tol, _dev_cache(sub, fname)))
        elseif chain !== nothing
            # the regularity principle: the criterion samples the chain's dof
            # field *sources* (a vector source through its magnitude), not the
            # derived quantity — the mesh follows u, the colors follow the mesh
            if isnan(spans[1][1]) || changed.subd_u
                for (i, s) in enumerate(chain.sources)
                    spans[i] = _field_span(s.ev, sub.used_cells, s.u[]; reduce=s.ev.ncomps > 1)
                end
            end
            src_lods = ntuple(length(chain.sources)) do i
                s = chain.sources[i]
                tol = _solution_tol(Float64(inputs.solution_tol), spans[i], T)
                DeviationLoD(_scalar_probe(s.ev, cellmap, s.u[]; reduce=s.ev.ncomps > 1), tol,
                             _dev_cache(sub, s.term))
            end
            lods = (lods..., src_lods...)
        end
        refine_keys!(state.keys, state.scratch, base, CombinedLoD(lods);
                     max_depth=Int(inputs.max_depth))
        # an update that does not change the key set leaves the output clean,
        # so the decode below is skipped entirely
        cached !== nothing && state.keys == state.prev && return nothing
        copy!(state.prev, state.keys)
        return (state.keys,)
    end

    # geometry, connectivity and colors from ONE edge: they can never be
    # resolved against different key sets. positions/faces enter Makie's mesh
    # internals, which may retain them: hand over fresh arrays. The reference
    # coordinates stay in-graph, so the buffer is shared.
    # Connectivity and values are separate nodes: the vertex layout depends
    # only on the key set, so a solution change on an unchanged mesh skips the
    # vertex-sharing lookups and re-evaluates one point per vertex — the same
    # work the static path does. Values still cannot disagree with the mesh
    # they belong to, because they are computed *from* the layout node's
    # output rather than beside it.
    #
    # Vertices are shared within a cell (where the drawn field is continuous)
    # and duplicated across cells, so element-boundary jumps survive.
    # Every node hands out the same buffer it filled, rather than a copy: the
    # pipeline compares an output against its previous value to decide what is
    # stale, and an array that *is* the previous one (same pointer) counts as
    # changed — so in-place reuse propagates exactly like a fresh array, while
    # a mesh that keeps its size stops allocating per update altogether.
    mesh_buf = IsubdMesh(base)
    faces_out = GeometryBasics.GLTriangleFace[]
    positions_out = GeometryBasics.Point{dim,Float32}[]
    colors_out = Float32[]
    ComputePipeline.register_computation!(graph, [:subd_keys], [:subd_ξ, :subd_faces]) do inputs, changed, cached
        decode_topology!(mesh_buf, inputs.subd_keys, base; groups=sub.groupmap)
        resize!(faces_out, length(mesh_buf.faces))
        @inbounds for i in eachindex(mesh_buf.faces)
            f = mesh_buf.faces[i]
            faces_out[i] = GeometryBasics.GLTriangleFace(f[1], f[2], f[3])
        end
        return (mesh_buf.refcoords, faces_out)
    end
    Makie.map!(graph, [:subd_ξ, :subd_u, warp_inputs...], :subd_positions) do _ξ, uu, _warps...
        _refresh_all!(sub, ev, uu)
        decode_positions!(mesh_buf, base)
        return _render_positions!(positions_out, mesh_buf.positions)
    end
    if ev === nothing && chain === nothing
        colornode = SP.color   # plain color, converted by the mesh child
    else
        Makie.map!(graph, [:subd_ξ, :subd_u], :subd_color) do _ξ, u
            _refresh_all!(sub, ev, u)
            _refresh_chain!(sub, chain)
            return chain === nothing ? _transfer_at!(colors_out, ev, cellmap, mesh_buf, u; reduce) :
                   _transfer_chain_at!(colors_out, chain, cellmap, mesh_buf)
        end
        colornode = SP.subd_color
    end

    # plain mesh child from graph nodes; the ShaderAbstractions.Buffer path of
    # `_mesh!` is tied to the static tessellation and does not apply here
    return Makie.mesh!(SP, SP.attributes, SP.subd_positions, SP.subd_faces, color=colornode)
end

# The quadrature-point branch of solutionplot's plot!: the substrate's base is
# the Voronoi-region fans (see `_isubd_base_qp`), the criterion is geometry
# only — piecewise-constant data has nothing to say about refinement, its
# discontinuities would only chase the depth cap, and the region rims are
# preserved structurally as split edges — and the colors are a flat gather:
# every vertex takes its base triangle's quadrature-point value. The values
# observable is a graph input, so updating the internal variables recolors
# the refined mesh the way `update!` moves a field.
function _wire_adaptive_qpplot!(SP, ds::FEData{dim}, sub::IsubdSubstrate, qp::QPPartition) where {dim}
    graph = SP.attributes
    ComputePipeline.add_input!(graph, :subd_u, ds.u)
    ComputePipeline.add_input!(graph, :subd_qpvalues, qp.values)
    warp_inputs = _warp_inputs!(graph, sub)
    base = sub.base
    state = (keys=root_keys(base), scratch=UInt64[], prev=UInt64[])
    diag = sub.diag
    T = _sample_type(sub)
    ComputePipeline.register_computation!(graph, [:subd_u, :geometry_tol, :max_depth, warp_inputs...],
                                          [:subd_keys]) do inputs, changed, cached
        _refresh_all!(sub, nothing, inputs.subd_u)
        lod = DeviationLoD(base.mapping, max(Float64(inputs.geometry_tol), _tol_floor(T)) * diag,
                           _dev_cache(sub, :geometry))
        refine_keys!(state.keys, state.scratch, base, lod; max_depth=Int(inputs.max_depth))
        cached !== nothing && state.keys == state.prev && return nothing
        copy!(state.prev, state.keys)
        return (state.keys,)
    end
    mesh_buf = IsubdMesh(base)
    faces_out = GeometryBasics.GLTriangleFace[]
    positions_out = GeometryBasics.Point{dim,Float32}[]
    colors_out = Float32[]
    ComputePipeline.register_computation!(graph, [:subd_keys], [:subd_ξ, :subd_faces]) do inputs, changed, cached
        decode_topology!(mesh_buf, inputs.subd_keys, base; groups=sub.groupmap)
        resize!(faces_out, length(mesh_buf.faces))
        @inbounds for i in eachindex(mesh_buf.faces)
            f = mesh_buf.faces[i]
            faces_out[i] = GeometryBasics.GLTriangleFace(f[1], f[2], f[3])
        end
        return (mesh_buf.refcoords, faces_out)
    end
    Makie.map!(graph, [:subd_ξ, :subd_u, warp_inputs...], :subd_positions) do _ξ, uu, _warps...
        _refresh_all!(sub, nothing, uu)
        decode_positions!(mesh_buf, base)
        return _render_positions!(positions_out, mesh_buf.positions)
    end
    ncells = Ferrite.getncells(Ferrite.get_grid(ds.dh))
    Makie.map!(graph, [:subd_ξ, :subd_qpvalues], :subd_color) do _ξ, vals
        _qp_ncells(vals) == ncells ||
            error("quadrature point data must have one entry per cell ($ncells), got $(_qp_ncells(vals))")
        resize!(colors_out, length(mesh_buf.refcoords))
        @inbounds for v in eachindex(mesh_buf.refcoords)
            b = mesh_buf.vertex_base[v]
            colors_out[v] = Float32(qp.extract(_qp_at(vals, sub.cellmap[b], sub.qpmap[b])))
        end
        return colors_out
    end
    return Makie.mesh!(SP, SP.attributes, SP.subd_positions, SP.subd_faces, color=SP.subd_color)
end
