# Gradient interpolation machinery, consumed by the Gradient filter.

# FIXME this is not correct in general and may lead to loss of details
get_gradient_interpolation(::Lagrange{shape,order}) where {sdim,shape<:Ferrite.AbstractRefShape{sdim},order} =
    VectorizedInterpolation{sdim}(DiscontinuousLagrange{shape,order-1}())
get_gradient_interpolation(::VectorizedInterpolation{vdim,shape,order,<:Lagrange{shape,order}}) where {sdim,vdim,shape<:Ferrite.AbstractRefShape{sdim},order} =
    MatrixizedInterpolation{vdim,sdim}(DiscontinuousLagrange{shape,order-1}())

"""
    _tensorsjl_gradient_accessor(v::Tensors.Vec, field_dim_idx::Int, spatial_dim_idx::Int)

This is a helper to access the correct value in Tensors.jl entities, because the gradient index is the outermost one.
"""
@inline _tensorsjl_gradient_accessor(v::Tensors.Vec{dim}, field_dim_idx::Int, spatial_dim_idx::Int) where {dim} = v[spatial_dim_idx]
@inline _tensorsjl_gradient_accessor(m::Tensors.Tensor{2,dim}, field_dim_idx::Int, spatial_dim_idx::Int) where {dim} = m[field_dim_idx, spatial_dim_idx]

function _check_full_domain(dh, what::String)
    length(dh.subdofhandlers) == 1 ||
        error("$what supports only DofHandlers with a single subdofhandler (single subdomain)")
    length(dh.subdofhandlers[1].cellset) == Ferrite.getncells(Ferrite.get_grid(dh)) ||
        error("$what supports only DofHandlers covering the full grid (the subdofhandler covers a subset of the cells)")
    return nothing
end

function _gradient_dofhandler(dh::DofHandler, field_name::Symbol, copy_fields::Vector{Symbol})
    field_idx = Ferrite.find_field(dh, field_name)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    dh_gradient = Ferrite.DofHandler(Ferrite.get_grid(dh))
    add!(dh_gradient, :gradient, get_gradient_interpolation(ip)) # field dim × spatial dim components
    for fieldname in copy_fields
        _field_idx = Ferrite.find_field(dh, fieldname)
        _ip = Ferrite.getfieldinterpolation(dh, _field_idx)
        add!(dh_gradient, fieldname, _ip)
    end
    Ferrite.close!(dh_gradient)
    return dh_gradient
end

function _compute_gradient_values(dh::DofHandler{spatial_dim}, dh_gradient::DofHandler, u::AbstractVector,
                                  field_name::Symbol, copy_fields::Vector{Symbol}) where {spatial_dim}
    field_idx = Ferrite.find_field(dh, field_name)
    ip = Ferrite.getfieldinterpolation(dh, field_idx)
    field_dim = Ferrite.n_components(dh, field_idx)
    ip_gradient = Ferrite.getfieldinterpolation(dh_gradient, Ferrite.find_field(dh_gradient, :gradient))

    # FIXME this does not work for mixed grids
    ip_geom = Ferrite.geometric_interpolation(typeof(Ferrite.getcells(Ferrite.get_grid(dh), 1)))
    ref_coords_gradient = Ferrite.reference_coordinates(ip_gradient)
    qr_gradient = QuadratureRule{getrefshape(Ferrite.getcells(get_grid(dh), 1))}(ones(length(ref_coords_gradient)), ref_coords_gradient)
    cv = CellValues(qr_gradient, ip, ip_geom)

    cell_dofs = zeros(Int, Ferrite.ndofs_per_cell(dh))
    cell_dofs_gradient = zeros(Int, Ferrite.ndofs_per_cell(dh_gradient))

    u_gradient = zeros(Ferrite.ndofs(dh_gradient))
    # In general uᵉ_gradient is an order 3 tensor [field_dim, spatial_dim, nqp]
    uᵉ_gradient = zeros(length(Ferrite.dof_range(dh_gradient.subdofhandlers[1], :gradient)))
    uᵉ_gradient_view = reshape(uᵉ_gradient, (spatial_dim, field_dim, getnquadpoints(qr_gradient)))

    for (cell_num, cell) in enumerate(Ferrite.CellIterator(dh))
        Ferrite.celldofs!(cell_dofs, dh, cell_num)
        uᵉ = @views u[cell_dofs[Ferrite.dof_range(dh, field_name)]]

        Ferrite.reinit!(cv, cell)

        # Evaluate the gradient at the basis function locations of the gradient field
        for i ∈ 1:Ferrite.getnquadpoints(qr_gradient)
            uᵉgradi = Ferrite.function_gradient(cv, i, uᵉ)
            for ds in 1:spatial_dim, df in 1:field_dim
                uᵉ_gradient_view[ds, df, i] = _tensorsjl_gradient_accessor(uᵉgradi, df, ds)
            end
        end

        Ferrite.celldofs!(cell_dofs_gradient, dh_gradient, cell_num)
        u_gradient[cell_dofs_gradient[Ferrite.dof_range(dh_gradient, :gradient)]] .+= uᵉ_gradient

        for fieldname in copy_fields
            u_gradient[cell_dofs_gradient[Ferrite.dof_range(dh_gradient, fieldname)]] .= u[cell_dofs[Ferrite.dof_range(dh, fieldname)]]
        end
    end
    return u_gradient
end

"""
    interpolate_gradient_field(dh::DofHandler, u::AbstractVector, field_name::Symbol; copy_fields::Vector{Symbol})

Compute the piecewise discontinuous gradient field for `field_name`. Returns the flux dof handler and the corresponding flux dof values.
If the additional keyword argument `copy_fields` is provided with a non empty `Vector{Symbol}`, the corresponding fields of `dh` will be
copied into the returned flux dof handler and flux dof value vector.
"""
function interpolate_gradient_field(dh::DofHandler, u::AbstractVector, field_name::Symbol; copy_fields::Vector{Symbol}=Symbol[])
    _check_full_domain(dh, "interpolate_gradient_field")
    dh_gradient = _gradient_dofhandler(dh, field_name, copy_fields)
    return dh_gradient, _compute_gradient_values(dh, dh_gradient, u, field_name, copy_fields)
end

# Foundation block for https://github.com/Ferrite-FEM/Ferrite.jl/issues/398 - Remove this below after the issue is resolved.
##################################################
# MatrixizedInterpolation{<:ScalarInterpolation} #
##################################################
abstract type MatrixInterpolation{vdim1, vdim2, refshape, order} <: Ferrite.Interpolation{refshape, order} end

struct MatrixizedInterpolation{vdim1, vdim2, refshape, order, SI <: ScalarInterpolation{refshape, order}} <: MatrixInterpolation{vdim1, vdim2, refshape,order}
    ip::SI
    function MatrixizedInterpolation{vdim1, vdim2}(ip::SI) where {vdim1, vdim2, refshape, order, SI <: Ferrite.ScalarInterpolation{refshape, order}}
        return new{vdim1, vdim2, refshape, order, SI}(ip)
    end
end

Ferrite.mapping_type(::MatrixizedInterpolation) = Ferrite.IdentityMapping()

Ferrite.typeof_N(   ::Type{T}, ::MatrixizedInterpolation{dim,dim}, ::VectorizedInterpolation{dim, <: Ferrite.AbstractRefShape{dim}}) where {T, dim} = Tensor{2, dim, T}
Ferrite.typeof_dNdx(::Type{T}, ::MatrixizedInterpolation{dim,dim}, ::VectorizedInterpolation{dim, <: Ferrite.AbstractRefShape{dim}}) where {T, dim} = Tensor{3, dim, T}
Ferrite.typeof_dNdξ(::Type{T}, ::MatrixizedInterpolation{dim,dim}, ::VectorizedInterpolation{dim, <: Ferrite.AbstractRefShape{dim}}) where {T, dim} = Tensor{3, dim, T}

Ferrite.n_components(::MatrixizedInterpolation{vdim1, vdim2}) where {vdim1, vdim2} = vdim1*vdim2
Ferrite.adjust_dofs_during_distribution(ip::MatrixizedInterpolation) = Ferrite.adjust_dofs_during_distribution(ip.ip)

# Matrixize to reference dimension by default
function MatrixizedInterpolation(ip::ScalarInterpolation{shape}) where {refdim, shape <: Ferrite.AbstractRefShape{refdim}}
    return MatrixizedInterpolation{refdim,refdim}(ip)
end

Base.:(^)(ip::VectorizedInterpolation{vdim1}, vdim2::Int) where {vdim1} = MatrixizedInterpolation{vdim1, vdim2}(ip)
function Base.literal_pow(::typeof(^), ip::VectorizedInterpolation{vdim1}, ::Val{vdim2}) where {vdim1,vdim2}
    return MatrixizedInterpolation{vdim1, vdim2}(ip.ip)
end

function Base.show(io::IO, mime::MIME"text/plain", ip::MatrixizedInterpolation{vdim1, vdim2}) where {vdim1, vdim2}
    show(io, mime, ip.ip)
    print(io, "^", vdim1 , "×", vdim2)
end

# Helper to get number of copies for DoF distribution
Ferrite.get_n_copies(::MatrixizedInterpolation{vdim1, vdim2}) where {vdim1, vdim2} = vdim1*vdim2

function Ferrite.getnbasefunctions(ipv::MatrixizedInterpolation{vdim1, vdim2}) where {vdim1, vdim2}
    return vdim1 * vdim2 * getnbasefunctions(ipv.ip)
end
function Ferrite.reference_shape_value(ipv::MatrixizedInterpolation{vdim, vdim, shape}, ξ::Tensors.Vec{refdim, T}, I::Int) where {vdim, refdim, shape <: Ferrite.AbstractRefShape{refdim}, T}
    # First flatten to vector
    i0, c0 = divrem(I - 1, vdim^2)
    i = i0 + 1
    v = Ferrite.reference_shape_value(ipv.ip, ξ, i)

    # Then compute matrix index
    ci0, cj0 = divrem(c0, vdim)
    ci = ci0 + 1
    cj = cj0 + 1
    return Ferrite.Tensor{2, vdim, T}((k, l) -> k == ci && l == cj ? v : zero(v))
end

Ferrite.reference_coordinates(ip::MatrixizedInterpolation) = Ferrite.reference_coordinates(ip.ip)

Ferrite.conformity(ip::MatrixizedInterpolation) = Ferrite.conformity(ip.ip)

Ferrite.InterpolationInfo(ip::MatrixizedInterpolation) = Ferrite.InterpolationInfo(ip.ip, Ferrite.get_n_copies(ip))
