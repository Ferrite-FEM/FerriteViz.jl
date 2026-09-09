# Run with: jld --idle-timeout=2h run benchmarks/cutting.jl
# Warm construction, clipping, ten-level extraction, and live cut traversal.
# Run the same file in a checkout of the parent commit for comparison.
# Example warm run on Julia 1.12.7 / Ferrite 1.7.0 / Makie 0.24.14, n=24:
#                        parent 9ee6a70       this change
# FEData allocations     125.2 MB             110.3 MB
# Clip allocations       112.1 MB              91.6 MB
# 10 isosurfaces          474 ms / 453.6 MB    163 ms / 110.9 MB
# Live cut traversal: ~66 ms, 31.8 MB; 2,256 of 13,824 cells tessellated.
# Times vary; these include CPU allocations, not GPU uploads or compilation.
using FerriteViz, Ferrite
import FerriteViz as FV
function cutbench(n)
    grid = generate_grid(Hexahedron,(n,n,n))
    dh = DofHandler(grid); add!(dh,:p,Lagrange{RefHexahedron,1}()); close!(dh)
    u = zeros(ndofs(dh)); Ferrite.apply_analytical!(u,dh,:p,x->x[1])
    FEData(dh,u;adaptivity=false)
    GC.gc(); src = @timed FEData(dh,u;adaptivity=false)
    ds = src.value
    f = Clip(ClipPlane(Vec(1.,0.,0.),0.123))
    f(ds)
    GC.gc(); cut = @timed f(ds)
    iso = ExtractIsosurfaces(collect(range(-0.8,0.8;length=10)))
    iso(ds)
    GC.gc(); ext = @timed iso(ds)
    println((n=n, source_seconds=src.time, source_bytes=src.bytes,
             clip_seconds=cut.time, clip_bytes=cut.bytes,
             iso_seconds=ext.time, iso_bytes=ext.bytes))
    if isdefined(FV,:CutRenderState)
        st = FV.CutRenderState(cut.value,:p)
        FV._render_cut(st,nothing)
        GC.gc(); r = @timed FV._render_cut(st,nothing)
        println((render_seconds=r.time,render_bytes=r.bytes,candidates=r.value.candidates,
                 cells=n^3,triangles=length(r.value.faces)))
    end
end
cutbench(12)
cutbench(24)
