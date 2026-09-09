# Per-cell polynomial coefficients, the fast path behind every pointwise
# evaluation the adaptive tessellation does.
#
# Evaluating an FE field the usual way costs one shape function per basis
# function per point — for a quadratic quad, nine tensor-product Lagrange
# polynomials to produce one number. But on a fixed cell the field *is* a
# polynomial in the reference coordinate, so it can be written once in a
# monomial basis
#
#     u(ξ) = Σ_m c_m ξ₁^{e₁ᵐ} ξ₂^{e₂ᵐ} …
#
# and then evaluated with a handful of multiply-adds. The change of basis is
# a small matrix solve per cell (the nodal values are the field at the
# interpolation's own reference nodes, so `V c = u` with `V` the monomials
# sampled there), done once per solution update and amortized over every
# sample point the error estimators ask for — which outnumber the cells by
# orders of magnitude.
#
# This is also the representation the fragment-shader track needs: a shader
# gets a per-cell coefficient buffer and the interpolated reference
# coordinate, and evaluates the same sum per pixel (#161).
#
# Not every interpolation is covered — serendipity spaces, anything whose
# basis count matches no tensor or total-degree monomial set, and any
# mapped-value (Piola) space are refused, and callers fall back to summing
# shape functions. `PolyBasis` verifies itself numerically before it is used,
# so a silently wrong fit cannot slip through.

# `P` is the number of powers a coordinate needs (maxdeg + 1) and lives in
# the type so the power table is a fixed, unrolled tuple rather than a loop
# over a runtime degree.
struct PolyBasis{refdim,N,P}
    exponents::NTuple{N,NTuple{refdim,Int}}
    Vinv::Matrix{Float64}       # nodal values -> coefficients
end

# Monomials for `nbase` basis functions in `refdim` dimensions: the tensor
# set for quad/hex Lagrange, the total-degree set for triangle/tet. Both are
# tried in that order and only an exact count is accepted.
function _monomial_exponents(refdim::Int, nbase::Int)
    for p in 0:8
        if (p + 1)^refdim == nbase
            exps = vec(collect(Iterators.product(ntuple(_ -> 0:p, refdim)...)))
            return exps, p
        end
    end
    for p in 0:8
        if binomial(p + refdim, refdim) == nbase
            exps = [e for e in Iterators.product(ntuple(_ -> 0:p, refdim)...) if sum(e) <= p]
            return exps, p
        end
    end
    return nothing, 0
end

"""
    PolyBasis(ip) -> PolyBasis or nothing

The change of basis from an interpolation's nodal values to monomial
coefficients, or `nothing` when the interpolation is not a nodal polynomial
space this can represent. Built once per interpolation type.
"""
function PolyBasis(ip::Ferrite.ScalarInterpolation)
    refdim = Ferrite.getrefdim(ip)
    nbase = Ferrite.getnbasefunctions(ip)
    exps, maxdeg = _monomial_exponents(refdim, nbase)
    exps === nothing && return nothing
    nodes = Ferrite.reference_coordinates(ip)
    length(nodes) == nbase || return nothing
    V = [prod(ntuple(d -> nodes[k][d]^exps[m][d], refdim)) for k in 1:nbase, m in 1:nbase]
    # a badly conditioned fit would evaluate to noise; refuse it instead
    c = LinearAlgebra.cond(V)
    (isfinite(c) && c < 1e10) || return nothing
    basis = PolyBasis{refdim,nbase,maxdeg + 1}(ntuple(m -> exps[m], nbase), inv(V))
    return _verifies(basis, ip) ? basis : nothing
end

PolyBasis(ip::Ferrite.VectorizedInterpolation) = PolyBasis(ip.ip)
PolyBasis(::Ferrite.Interpolation) = nothing

# The fit must reproduce every shape function at points that are not nodes.
function _verifies(basis::PolyBasis{refdim,N,P}, ip) where {refdim,N,P}
    probes = (ntuple(d -> 0.137 + 0.041d, refdim), ntuple(d -> -0.219 - 0.017d, refdim),
              ntuple(d -> 0.0, refdim))
    for p in probes
        ξ = Tensors.Vec(p)
        Ferrite.getrefshape(ip) <: Ferrite.RefHypercube || all(x -> x >= 0, p) || continue
        mono = monomials(basis, ξ)
        for i in 1:N
            # coefficients of shape function i are column i of Vinv
            approx = sum(basis.Vinv[m, i] * mono[m] for m in 1:N)
            exact = Ferrite.reference_shape_value(ip, ξ, i)
            isapprox(approx, exact; atol = 1e-9, rtol = 1e-7) || return false
        end
    end
    return true
end

# The monomials at ξ, from powers computed once per dimension. Computed in
# ξ's number type: the adaptive pipeline samples in its (typically Float32)
# type end to end, the basis verification in Float64.
@inline function monomials(basis::PolyBasis{refdim,N,P}, ξ) where {refdim,N,P}
    pw = ntuple(d -> _powers(ξ[d], Val(P)), refdim)
    return ntuple(m -> begin
                      e = basis.exponents[m]
                      v = one(ξ[1])
                      @inbounds for d in 1:refdim
                          v *= pw[d][e[d] + 1]
                      end
                      v
                  end, Val(N))
end

# Unrolled at the type level, so each entry is a literal power the compiler
# turns into multiplications.
@inline _powers(x::Real, ::Val{P}) where {P} = ntuple(i -> x^(i - 1), Val(P))

"""
    PolyField(basis, ncells, T) -> PolyField

Per-cell monomial coefficients of one field. `refresh!` recomputes them from
the current nodal values; `evaluate` then costs `N` multiply-adds.
"""
struct PolyField{refdim,N,P,T}
    basis::PolyBasis{refdim,N,P}
    coeffs::Matrix{T}           # N × ncells
    filled::Vector{Bool}        # per cell, whether coefficients are present
    # Which solution epoch the coefficients were computed at (see
    # `IsubdSubstrate`): a refresh for the same epoch is a no-op, which is what
    # lets several plots share one evaluator without refreshing it per node.
    # `-1` means never; the geometry field, refreshed once from the node
    # coordinates, stays there.
    epoch::Base.RefValue{Int}
end

PolyField(basis::PolyBasis{refdim,N,P}, ncells::Int, ::Type{T}) where {refdim,N,P,T} =
    PolyField{refdim,N,P,T}(basis, zeros(T, N, ncells), zeros(Bool, ncells), Ref(-1))

# c = V⁻¹ u, one small matvec per cell.
function refresh_cell!(pf::PolyField{refdim,N,P,T}, cell::Int, nodal) where {refdim,N,P,T}
    Vinv = pf.basis.Vinv
    @inbounds for m in 1:N
        acc = zero(T)
        for k in 1:N
            acc += Vinv[m, k] * nodal[k]
        end
        pf.coeffs[m, cell] = acc
    end
    pf.filled[cell] = true
    return pf
end

@inline function evaluate(pf::PolyField{refdim,N,P,T}, cell::Int, ξ) where {refdim,N,P,T}
    mono = monomials(pf.basis, ξ)
    acc = zero(T)
    @inbounds for m in 1:N
        acc += pf.coeffs[m, cell] * mono[m]
    end
    return acc
end
