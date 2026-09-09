using Ferrite

elements = [Triangle((4, 17, 9)), Triangle((11, 16, 5)), Triangle((5, 16, 12)), Triangle((13, 17, 4)), Triangle((9, 18, 3)), Triangle((3, 18, 11)), Triangle((12, 15, 6)), Triangle((6, 15, 13)), Triangle((17, 18, 9)), Triangle((12, 16, 15)), Triangle((16, 17, 15)), Triangle((11, 18, 16)), Triangle((16, 18, 17)), Triangle((15, 17, 13)), Quadrilateral((1, 7, 14, 10)), Quadrilateral((10, 14, 9, 3)), Quadrilateral((7, 2, 8, 14)), Quadrilateral((14, 8, 4, 9))]
nodes = Node.(Vec{2, Float64}.([[-0.5, -1.0], [0.5, -1.0], [-0.5, 0.0], [0.5, 0.0], [-0.5, 1.0], [0.5, 1.0], [-1.3751222383007189e-12, -1.0], [0.5, -0.5000000000020595], [1.3751222383007189e-12, 0.0], [-0.5, -0.5000000000020595], [-0.5, 0.4999999999986921], [-1.3751222383007189e-12, 1.0], [0.5, 0.4999999999986921], [-5.664046543035832e-24, -0.5000000000020595], [0.206249999999433, 0.7062499999994473], [-0.12500000000067074, 0.6249999999999999], [0.14375000000014615, 0.35208333333295083], [-0.21874999999926012, 0.28124999999956607]]))
facesets = Dict{String, Set{FaceIndex}}("bottom" => Set([FaceIndex((15, 1)), FaceIndex((17, 1))]), "top" => Set([FaceIndex((7, 3)), FaceIndex((3, 3))]))
cellsets = Dict{String, Set{Int64}}("quad" => Set([15, 16, 18, 17]), "triangle" => Set([5, 12, 8, 1, 6, 11, 9, 14, 3, 7, 4, 13, 2, 10]))

grid = Grid(elements,nodes,facesets=facesets,cellsets=cellsets)

dh = DofHandler(grid)
sdh_tri = SubDofHandler(dh, getcellset(grid, "triangle"))
add!(sdh_tri, :u, Lagrange{RefTriangle,1}()^2)
add!(sdh_tri, :p, Lagrange{RefTriangle,1}())
sdh_quad = SubDofHandler(dh, getcellset(grid, "quad"))
add!(sdh_quad, :u, Lagrange{RefQuadrilateral,1}()^2)
close!(dh)

u = zeros(ndofs(dh))

for cell in CellIterator(dh,collect(dh.subdofhandlers[1].cellset))
    celldofs_ = celldofs(cell)
    u[celldofs_] .= 1
end
for cell in CellIterator(dh,collect(dh.subdofhandlers[2].cellset))
    celldofs_ = celldofs(cell)
    dof_range_ = Ferrite.dof_range(dh.subdofhandlers[2],:u)
    u[celldofs_[dof_range_]] .= 0.5
end
