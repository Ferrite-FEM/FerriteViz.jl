# Developer Documentation

Note that these functions could be removed or change in behavior between minor version changes! Use and dispatch on these with care!

```@docs
FerriteViz.num_vertices
FerriteViz.vertices
FerriteViz.transfer_quadrature_face_to_cell
FerriteViz.decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Union{FerriteViz.Ferrite.AbstractCell{2,N,3}, FerriteViz.Ferrite.AbstractCell{3,3,1}}) where {N}
FerriteViz.decompose!(coord_offset, coord_matrix::Vector{Point{space_dim,T}}, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::Union{FerriteViz.Ferrite.AbstractCell{2,N,4}, FerriteViz.Ferrite.AbstractCell{3,4,1}}) where {N,space_dim,T}
FerriteViz.decompose!(coord_offset, coord_matrix, ref_coord_matrix, triangle_offset, triangle_matrix, grid, cell::FerriteViz.Ferrite.AbstractCell{3,N,M}) where {N,M}
FerriteViz.transfer_solution
FerriteViz.postprocess
FerriteViz._tensorsjl_gradient_accessor
FerriteViz.linear_face_cell
```
