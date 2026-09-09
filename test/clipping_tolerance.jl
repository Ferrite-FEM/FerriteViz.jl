@testset "point-local clipping tolerance on graded grids" begin
    # A distant cell must not determine the tolerance of a micron-sized one.
    small = generate_grid(Hexahedron,(1,1,1),Vec(-1e-7,-1e-7,-1e-7),Vec(1e-7,1e-7,1e-7))
    far = generate_grid(Hexahedron,(1,1,1),Vec(9.,9.,9.),Vec(11.,11.,11.))
    grid=Grid([only(getcells(small)),Hexahedron(ntuple(i -> i+8,8))],vcat(getnodes(small),getnodes(far)))
    dh=DofHandler(grid);add!(dh,:u,Lagrange{RefHexahedron,1}());close!(dh)
    ds=FEData(dh,zeros(ndofs(dh)))
    cut=ds |> Clip(ClipPlane(Vec(1.,0.,0.),0.0))
    @test cut.solid == [true,false]
    @test _kept_volume(cut) ≈ 4e-21 rtol=1e-6
end
