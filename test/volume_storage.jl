@testset "lazy volume connectivity" begin
    large,_,_ = _hexds(12)
    cut = large |> Clip(ClipPlane(Vec(1.,0.,0.),0.1))
    @test large.simplices isa FerriteViz.FanSimplices
    @test cut.simplices isa FerriteViz.ClippedSimplices
    @test length(cut.simplices.cuts) < length(cut.simplices) ÷ 2
end
