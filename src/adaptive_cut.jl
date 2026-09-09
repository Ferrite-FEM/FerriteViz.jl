# Streaming cut rendering. FEData's coarse arrays are an inspectable snapshot;
# plots retain the source domain and re-cut its continuous geometry in their
# ComputeGraph. Only exterior faces and plane/cell candidates enter recursion.
# No refined volume mesh is stored. A common red-refinement level keeps shared
# tetrahedron faces conforming; interval bounds prune invisible subtrees.

function _cut_chain(ds)
    planes = ClipPlane[]
    root = ds
    while root.cut_domain !== nothing
        dom = root.cut_domain
        dom.filter isa Clip || return nothing
        root.deformation == dom.source.deformation || return nothing
        pushfirst!(planes, dom.filter.plane)
        root = dom.source
    end
    return root, planes
end

function _cut_render_capable(ds)
    ds.cut_domain === nothing && return false
    chain = _cut_chain(ds)
    chain === nothing && return false
    root, _ = chain
    root.cells_intact && root.qp_partition === nothing || return false
    # Raw-array warps have no continuous interpretation.
    all(w -> w.field[] in Ferrite.getfieldnames(w.dh), ds.deformation) || return false
    return all(c -> getrefshape(c) in (Ferrite.RefTetrahedron, Ferrite.RefHexahedron,
                                     Ferrite.RefPrism, Ferrite.RefPyramid),
               Ferrite.getcells(get_grid(root.dh)))
end

struct CutSample
    x::Tensors.Vec{3,Float64}
    ξ::Tensors.Vec{3,Float64}
    value::Float64
end
_cut_lerp(a::CutSample, b::CutSample, t) =
    CutSample(a.x + t * (b.x - a.x), a.ξ + t * (b.ξ - a.ξ), a.value + t * (b.value - a.value))
function _cut_distance(p, plane)
    g = p.x ⋅ plane.normal - plane.distance
    tol = 16eps(Float64) * (sum(abs(p.x[i] * plane.normal[i]) for i in 1:3) + abs(plane.distance))
    return abs(g) <= tol ? 0.0 : g
end

# Snap to existing endpoints: this also avoids duplicate on-plane vertices.
function _cut_intersection(a, b, ga, gb)
    ga == 0 && return a
    gb == 0 && return b
    # Canonical endpoint order makes shared-edge intersections bit-identical.
    isless(NTuple{3,Float64}(b.ξ), NTuple{3,Float64}(a.ξ)) && return _cut_intersection(b, a, gb, ga)
    return _cut_lerp(a, b, ga / (ga - gb))
end

function _cut_polygon(poly, plane)
    out = CutSample[]
    isempty(poly) && return out
    a = last(poly); ga = _cut_distance(a, plane)
    for b in poly
        gb = _cut_distance(b, plane)
        (ga <= 0) != (gb <= 0) && push!(out, _cut_intersection(a, b, ga, gb))
        gb <= 0 && push!(out, b)
        a, ga = b, gb
    end
    return out
end

const _CUT_EDGES = ((1,2), (1,3), (1,4), (2,3), (2,4), (3,4))
const _CUT_FACES = ((1,2,3), (1,4,2), (2,4,3), (3,4,1))
const _CUT_MISSING = (4,3,1,2)
const _CUT_RED = ((1,5,6,7), (5,2,8,9), (6,8,3,10), (7,9,10,4),
                  (5,6,7,10), (5,6,8,10), (5,7,9,10), (5,8,9,10))
const _CUT_SUPPORT = (0x01,0x02,0x04,0x08,0x03,0x05,0x09,0x06,0x0a,0x0c)

function _cut_children(ξs, exterior)
    pts = (ξs..., map(e -> (ξs[e[1]] + ξs[e[2]]) / 2, _CUT_EDGES)...)
    orientation = _signed_tet_volume(ξs...)
    return map(_CUT_RED) do child
        if _signed_tet_volume(map(i -> pts[i], child)...) * orientation < 0
            child = (child[2], child[1], child[3], child[4])
        end
        faces = ntuple(4) do f
            support = foldl(|, (_CUT_SUPPORT[child[i]] for i in _CUT_FACES[f]))
            any(j -> exterior[j] && support & (0x01 << (_CUT_MISSING[j]-1)) == 0, 1:4)
        end
        (map(i -> pts[i], child), faces)
    end
end

# Polynomial range enclosure over a reference AABB. Unlike a nodal AABB it
# cannot miss a high-order mode that bulges beyond all the geometric nodes.
function _power_range(a, b, e)
    e == 0 && return (1.0, 1.0)
    x, y = a^e, b^e
    return (iseven(e) && a <= 0 <= b ? 0.0 : min(x,y), max(x,y))
end
function _poly_plane_range(pf::PolyField{3,N}, cell, lo, hi, normal) where {N}
    lower = upper = 0.0
    for m in 1:N
        a = b = 1.0
        for d in 1:3
            c, e = _power_range(lo[d], hi[d], pf.basis.exponents[m][d])
            q = (a*c, a*e, b*c, b*e)
            a, b = minimum(q), maximum(q)
        end
        c = pf.coeffs[m,cell] ⋅ normal
        lower += min(c*a, c*b)
        upper += max(c*a, c*b)
    end
    # Enclose coefficient accumulation roundoff too.
    pad = 64eps(Float64) * max(abs(lower), abs(upper), 1.0)
    return lower-pad, upper+pad
end
_poly_plane_range(::Nothing, args...) = (-Inf, Inf)

struct CutRenderState{D,P,G,W,C}
    root::D
    planes::P
    geometry::G
    warps::W
    color::C
    inputs::Vector{Makie.Observable}
    sources::Vector{Any}
    diag::Float64
end

function _cut_color(ds, color, inputs, sources)
    color isa Symbol || return (cell, ξ) -> 0.0
    name = _resolve_name(ds, color)
    if name in Ferrite.getfieldnames(ds.dh)
        ev = FieldEvaluator(ds.dh, name, Float64)
        push!(inputs, ds.u)
        push!(sources, (ev,ds.u))
        reduce = color === :default && ev.ncomps > 1
        reduce || ev.ncomps == 1 || error("cut coloring needs a scalar field or :default")
        return (cell, ξ) -> begin
            v = evaluate_at(ev, cell, ξ, ds.u[])
            v === nothing ? NaN : reduce ? LinearAlgebra.norm(v) : Float64(v[1])
        end
    elseif haskey(ds.point_derivations, name)
        return _cut_derived(ds.point_derivations[name], inputs, sources)
    elseif haskey(ds.cell_data, name)
        obs = ds.cell_data[name]; push!(inputs, obs)
        return (cell, ξ) -> Float64(obs[][cell])
    elseif _data_association(ds, name) === :none
        return (cell, ξ) -> 0.0
    end
    error("cut coloring cannot resample raw point data :$name")
end
function _cut_derived(src::FieldSource, inputs, sources)
    ev = FieldEvaluator(src.dh, src.name, Float64)
    push!(inputs, src.u)
    push!(sources, (ev,src.u))
    return (cell, ξ) -> begin
        v = evaluate_at(ev, cell, ξ, src.u[])
        v === nothing ? nothing : _wrap_value(v, Val(3))
    end
end
function _cut_derived(rec::DerivedPointData, inputs, sources)
    fs = map(s -> _cut_derived(s, inputs, sources), Tuple(rec.inputs))
    return (cell, ξ) -> begin
        vs = map(f -> f(cell, ξ), fs)
        any(isnothing, vs) ? NaN : rec.f(vs...)
    end
end

function CutRenderState(ds, color)
    root, ps = _cut_chain(ds)
    planes = unique([ClipPlane(Tensors.Vec{3,Float64}(NTuple{3,Float64}(p.normal / LinearAlgebra.norm(p.normal))), Float64(p.distance / LinearAlgebra.norm(p.normal))) for p in ps])
    grid = get_grid(root.dh)
    cells = Ferrite.getcells(grid)
    gips = [Ferrite.geometric_interpolation(typeof(c)) for c in cells]
    coords = [Ferrite.getcoordinates(grid, i) for i in eachindex(cells)]
    geo = _geometry_poly(gips, coords, collect(eachindex(cells)), Val(3), Float64)
    warps = Tuple(WarpEval(d, Float64) for d in root.deformation)
    inputs = Makie.Observable[ds.u]
    for w in warps
        append!(inputs, (w.u, w.scale, w.name))
    end
    sources = Any[]
    c = _cut_color(ds, color, inputs, sources)
    return CutRenderState(root, Tuple(planes), (poly=geo, ips=gips, coords=coords), warps,
                          c, unique(objectid, inputs), sources, _grid_diagonal(grid))
end

function _cut_mapping(st, cell, ξ)
    g = st.geometry
    x = g.poly !== nothing && g.poly.filled[cell] ? evaluate(g.poly, cell, ξ) :
        geometric_map(g.ips[cell], g.coords[cell], ξ)
    for w in st.warps
        v = evaluate_at(w.ev[], cell, ξ, w.u[])
        v === nothing || (x += w.scale[] * v)
    end
    return Tensors.Vec{3,Float64}(ntuple(i -> x[i], 3))
end

function _cut_range(st, cell, ξs, plane)
    lo = ntuple(d -> minimum(ξ[d] for ξ in ξs), 3)
    hi = ntuple(d -> maximum(ξ[d] for ξ in ξs), 3)
    g = st.geometry.poly
    g !== nothing && g.filled[cell] || return (-Inf, Inf)
    a,b = _poly_plane_range(g, cell, lo, hi, plane.normal)
    for w in st.warps
        w.ev[].sdh_of_cell[cell] == 0 && continue
        p = w.ev[].poly
        p !== nothing && p.filled[cell] || return (-Inf, Inf)
        c,d = _poly_plane_range(p, cell, lo, hi, plane.normal)
        s = Float64(w.scale[])
        a += min(s*c,s*d); b += max(s*c,s*d)
    end
    return a-plane.distance, b-plane.distance
end

mutable struct CutRenderMesh
    positions::Vector{GeometryBasics.Point3f}
    faces::Vector{GeometryBasics.GLTriangleFace}
    colors::Vector{Float32}
    edges::Vector{GeometryBasics.Point3f}
    candidates::Int
    visited::Int
    depth::Int
end
CutRenderMesh() = CutRenderMesh(GeometryBasics.Point3f[], GeometryBasics.GLTriangleFace[], Float32[], GeometryBasics.Point3f[], 0,0,0)
function _emit_cut_polygon!(out, poly)
    for k in 2:length(poly)-1
        a,b,c = poly[1],poly[k],poly[k+1]
        _tri_area(a.x,b.x,c.x) > 0 || continue
        i = length(out.positions)
        for p in (a,b,c)
            push!(out.positions, GeometryBasics.Point3f(p.x...))
            push!(out.colors, Float32(p.value))
        end
        push!(out.faces, GeometryBasics.GLTriangleFace(i+1,i+2,i+3))
    end
end

function _cut_leaf!(out, samples, exterior, planes)
    any(p -> all(s -> _cut_distance(s,p) >= 0, samples), planes) && return
    # Original boundary facets. Buried cell faces never enter the output.
    for f in 1:4
        exterior[f] || continue
        poly = [samples[i] for i in _CUT_FACES[f]]
        for p in planes
            poly = _cut_polygon(poly, p)
        end
        _emit_cut_polygon!(out, poly)
    end
    for (k,p) in enumerate(planes)
        gs = map(s -> _cut_distance(s,p), samples)
        # A zero-volume contact must not duplicate the kept neighbor's face.
        minimum(gs) < 0 && maximum(gs) >= 0 || continue
        any(f -> exterior[f] && all(i -> gs[i] == 0, _CUT_FACES[f]), 1:4) && continue
        poly = [samples[i] for i in 1:4 if gs[i] == 0]
        for (i,j) in _CUT_EDGES
            gs[i]*gs[j] < 0 || continue
            s = _cut_intersection(samples[i],samples[j],gs[i],gs[j])
            any(v -> v.x == s.x, poly) || push!(poly,s)
        end
        length(poly) >= 3 || continue
        center = sum(s.x for s in poly) / length(poly)
        u = first(poly).x-center
        v = p.normal × u
        sort!(poly; by=s -> atan((s.x-center) ⋅ v, (s.x-center) ⋅ u))
        for (j,q) in enumerate(planes)
            j == k && continue
            poly = _cut_polygon(poly,q)
        end
        _emit_cut_polygon!(out,poly)
    end
end

function _visit_cut!(out, st, cell, ξs, exterior, depth, target, tolerances, needs_refinement)
    out.visited += 1
    ranges = map(p -> _cut_range(st,cell,ξs,p), st.planes)
    any(r -> r[1] > 0, ranges) && return
    !any(exterior) && all(r -> r[2] < 0, ranges) && return
    if depth < target
        for (child,flags) in _cut_children(ξs,exterior)
            _visit_cut!(out,st,cell,child,flags,depth+1,target,tolerances,needs_refinement)
        end
        return
    end
    samples = map(ξ -> CutSample(_cut_mapping(st,cell,ξ), ξ, Float64(st.color(cell,ξ))), ξs)
    if !needs_refinement[] && any(isfinite, tolerances)
        for (plane,range) in zip(st.planes,ranges)
            gs = map(s -> _cut_distance(s,plane),samples)
            if range[1] < -tolerances[1] && range[2] > tolerances[1] &&
               (minimum(gs) > 0 || maximum(gs) < 0)
                needs_refinement[] = true
            end
        end
        # Edges and the interior: a displacement bubble can vanish on every
        # boundary edge, which is why surface-only refinement was insufficient.
        for weights in ((1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(1,2,3,4))
            ξ = sum(ξs[i] for i in weights)/length(weights)
            x = sum(samples[i].x for i in weights)/length(weights)
            v = sum(samples[i].value for i in weights)/length(weights)
            if LinearAlgebra.norm(_cut_mapping(st,cell,ξ)-x) > tolerances[1] ||
               abs(st.color(cell,ξ)-v) > tolerances[2]
                needs_refinement[] = true
                break
            end
        end
    end
    _cut_leaf!(out,samples,exterior,st.planes)
end

# Preserve explicit Refine(surface=..., edges=...) as the starting layout.
# The last source vertex is the volume fan's centroid; rebuild that centroid
# through the continuous map instead of retaining its physical affine mean.
function _cut_reference_tessellation(ds, cell)
    offset = ds.cell_vertex_offsets[cell]
    verts = vertices_on_cell(ds,cell)
    coords = [Tensors.Vec{3,Float64}(ntuple(d -> ds.reference_coords[v,d],3))
              for v in first(verts):last(verts)-1]
    triangles = [ntuple(j -> convert(Int,ds.all_triangles[t][j])-offset,3)
                 for t in triangles_on_cell(ds,cell)]
    edges = [map(v -> v-offset,ds.all_edges[e]) for e in edges_on_cell(ds,cell)]
    return ReferenceTessellation(coords,triangles,edges)
end

function _render_cut(st::CutRenderState, cfg)
    cells = Ferrite.getcells(get_grid(st.root.dh))
    ids = findall(st.root.solid)
    for w in st.warps
        _refresh_warp!(w, ids, w.ev[].poly === nothing ? 0 : w.ev[].poly.epoch[]+1)
    end
    for (ev,u) in st.sources
        prepare!(ev,ids,u[])
    end
    # Candidate discovery precedes any tetrahedral decomposition.
    selected = Int[]
    tessellations = Dict{DataType,ReferenceTessellation{3,Float64}}()
    spanlo,spanhi = Inf,-Inf
    for cell in ids
        tess = get!(() -> _cut_reference_tessellation(st.root,cell), tessellations, typeof(cells[cell]))
        ranges = map(p -> _cut_range(st,cell,tess.coords,p), st.planes)
        any(r -> r[1] > 0,ranges) && continue
        surface = any(f -> _is_surface_facet(st.root,cell,f), eachindex(Ferrite.reference_faces(getrefshape(cells[cell]))))
        !surface && all(r -> r[2] < 0,ranges) && continue
        push!(selected,cell)
        for ξ in tess.coords
            v = st.color(cell,ξ)
            isfinite(v) || continue
            spanlo,spanhi = min(spanlo,v),max(spanhi,v)
        end
    end
    span = spanlo <= spanhi ? spanhi-spanlo : 0.0
    tol = cfg === nothing ? (Inf,Inf) :
        (max(cfg.geometry_tol[],_tol_floor(cfg.sample_type))*st.diag,
         _solution_tol(cfg.solution_tol[],(span,max(abs(spanlo),abs(spanhi))),cfg.sample_type))
    maxrounds = cfg === nothing ? 0 : max(0,fld(cfg.max_depth[],3))
    out = CutRenderMesh(); out.candidates = length(selected)
    for round in 0:maxrounds
        empty!(out.positions); empty!(out.faces); empty!(out.colors); empty!(out.edges)
        needs_refinement = Ref(false)
        for cell in selected
            tess = tessellations[typeof(cells[cell])]
            center = sum(tess.coords)/length(tess.coords)
            tri = 0
            faces = Ferrite.reference_faces(getrefshape(cells[cell]))
            subdivisions = length(tess.triangles) ÷ sum(face -> length(face)==3 ? 1 : 4, faces)
            for (f,face) in enumerate(Ferrite.reference_faces(getrefshape(cells[cell])))
                exterior = _is_surface_facet(st.root,cell,f)
                for _ in 1:(subdivisions * (length(face)==3 ? 1 : 4))
                    tri += 1
                    ξs = (map(i -> tess.coords[i],tess.triangles[tri])...,center)
                    _visit_cut!(out,st,cell,ξs,(exterior,false,false,false),0,round,tol,needs_refinement)
                end
            end
        end
        out.depth = 3round
        (!needs_refinement[] || round == maxrounds) && break
    end
    # Wireframe uses the same dyadic subdivisions as the volume boundary.
    n = 2^div(out.depth,3)
    for cell in selected
        tess = tessellations[typeof(cells[cell])]
        for (i,j) in tess.edges, k in 0:n-1
            aξ = tess.coords[i] + (k/n)*(tess.coords[j]-tess.coords[i])
            bξ = tess.coords[i] + ((k+1)/n)*(tess.coords[j]-tess.coords[i])
            a = CutSample(_cut_mapping(st,cell,aξ),aξ,0.0)
            b = CutSample(_cut_mapping(st,cell,bξ),bξ,0.0)
            kept = true
            for p in st.planes
                ga,gb = _cut_distance(a,p),_cut_distance(b,p)
                if ga > 0 && gb > 0
                    kept = false; break
                elseif ga > 0
                    a = _cut_intersection(a,b,ga,gb)
                elseif gb > 0
                    b = _cut_intersection(a,b,ga,gb)
                end
            end
            kept && a.x != b.x && append!(out.edges,(GeometryBasics.Point3f(a.x...),GeometryBasics.Point3f(b.x...)))
        end
    end
    return out
end

function _wire_cut_plot!(plot, ds; wireframe=false)
    graph = plot.attributes
    color = wireframe ? :black : plot.color[]
    st = CutRenderState(ds,color)
    inputs = Symbol[]
    for (i,obs) in enumerate(st.inputs)
        key = Symbol(:cut_input_,i)
        ComputePipeline.add_input!(graph,key,obs); push!(inputs,key)
    end
    if ds.adaptivity !== nothing
        _adaptivity_inputs!(graph,ds.adaptivity; solution=true)
        append!(inputs,(:geometry_tol,:solution_tol,:max_depth))
    end
    ComputePipeline.register_computation!(graph,inputs,[:cut_positions,:cut_faces,:cut_colors,:cut_edges]) do _,_,_
        out = _render_cut(st,ds.adaptivity)
        return (out.positions,out.faces,out.colors,out.edges)
    end
    if wireframe
        Makie.map!(identity,graph,:cut_edges,:edge_lines)
    else
        named = color isa Symbol && (color === :default || _data_association(ds,color) !== :none)
        Makie.mesh!(plot,plot.attributes,plot.cut_positions,plot.cut_faces;
                    color=named ? plot.cut_colors : plot.color)
    end
    return nothing
end
