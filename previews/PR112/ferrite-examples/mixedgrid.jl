using FerriteGmsh
using Ferrite

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("demo")

lc = 0.2
gmsh.model.geo.addPoint(-0.5, -1, 0, lc, 1)
gmsh.model.geo.addPoint(0.5, -1, 0, lc, 2)
gmsh.model.geo.addPoint(-0.5, 0, 0, lc, 3)
gmsh.model.geo.addPoint(0.5, 0, 0, lc, 4)
gmsh.model.geo.addPoint(-0.5, 1, 0, lc, 5)
gmsh.model.geo.addPoint(0.5, 1, 0, lc, 6)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 4, 2)
gmsh.model.geo.addLine(4, 3, 3)
gmsh.model.geo.addLine(1, 3, 4)
gmsh.model.geo.addLine(3, 5, 5)
gmsh.model.geo.addLine(5, 6, 6)
gmsh.model.geo.addLine(4, 6, 7)

gmsh.model.geo.addCurveLoop([1, 2, 3, -4], 1)
gmsh.model.geo.addCurveLoop([-3, 7, -6, -5], 2)
gmsh.model.geo.addPlaneSurface([1], 1)
gmsh.model.geo.addPlaneSurface([2], 2)
gmsh.model.geo.mesh.setTransfiniteCurve(1, 3)
gmsh.model.geo.mesh.setTransfiniteCurve(2, 3)
gmsh.model.geo.mesh.setTransfiniteCurve(3, 3)
gmsh.model.geo.mesh.setTransfiniteCurve(4, 3)
gmsh.model.geo.mesh.setTransfiniteCurve(5, 3)
gmsh.model.geo.mesh.setTransfiniteCurve(6, 3)
gmsh.model.geo.mesh.setTransfiniteCurve(7, 3)
gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.geo.mesh.setRecombine(2, 1)

gmsh.model.addPhysicalGroup(2, [1], 1)
gmsh.model.setPhysicalName(2, 1, "quad")

gmsh.model.addPhysicalGroup(2, [2], 2)
gmsh.model.setPhysicalName(2, 2, "triangle")

gmsh.model.addPhysicalGroup(1, [6], 3)
gmsh.model.setPhysicalName(1, 3, "top")

gmsh.model.addPhysicalGroup(1, [1], 4)
gmsh.model.setPhysicalName(1, 4, "bottom")

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)

nodes = tonodes()
elements, gmsh_eleidx = toelements(2)
boundarydict = toboundary(1)
facesets = tofacesets(boundarydict,elements)
cellsets = tocellsets(2,gmsh_eleidx)
grid = Grid(elements,nodes,facesets=facesets,cellsets=cellsets)

dh = MixedDofHandler(grid)
push!(dh,FieldHandler([Field(:p,Lagrange{2,RefTetrahedron,1}(),1),Field(:u,Lagrange{2,RefTetrahedron,1}(),2)], getcellset(grid,"triangle")))
push!(dh,FieldHandler([Field(:u,Lagrange{2,RefCube,1}(),2)], getcellset(grid,"quad")))
close!(dh)

u = zeros(ndofs(dh))

for cell in CellIterator(dh,collect(dh.fieldhandlers[1].cellset))
    celldofs_ = celldofs(cell)
    u[celldofs_] .= 1
end
for cell in CellIterator(dh,collect(dh.fieldhandlers[2].cellset))
    celldofs_ = celldofs(cell)
    dof_range_ = Ferrite.dof_range(dh.fieldhandlers[2],:u)
    u[celldofs_[dof_range_]] .= 0.5
end
