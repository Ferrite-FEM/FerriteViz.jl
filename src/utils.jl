#Visualization is just fancy triangle plotting, every element needs to be translatable to a triangle *sadnoises*
to_triangle(cell::Ferrite.AbstractCell{2,N,3}) where N = [Ferrite.vertices(cell)[1:3]]
to_triangle(cell::Ferrite.AbstractCell{3,N,4}) where N = [Ferrite.vertices(cell)[[1,2,3]], Ferrite.vertices(cell)[[1,2,4]], Ferrite.vertices(cell)[[2,3,4]], Ferrite.vertices(cell)[[1,4,3]]]
to_triangle(cell::Ferrite.AbstractCell{2,N,4}) where N = [Ferrite.vertices(cell)[[1,2,3]], Ferrite.vertices(cell)[[3,4,1]]]
to_triangle(cell::Ferrite.AbstractCell{3,N,6}) where N = [Ferrite.vertices(cell)[[1,2,3]], Ferrite.vertices(cell)[[3,4,1]],
                                                          Ferrite.vertices(cell)[[1,5,6]], Ferrite.vertices(cell)[[6,2,1]],
                                                          Ferrite.vertices(cell)[[2,6,7]], Ferrite.vertices(cell)[[7,3,2]],
                                                          Ferrite.vertices(cell)[[3,7,8]], Ferrite.vertices(cell)[[8,4,3]],
                                                          Ferrite.vertices(cell)[[1,4,8]], Ferrite.vertices(cell)[[8,5,1]],
                                                          Ferrite.vertices(cell)[[5,8,7]], Ferrite.vertices(cell)[[7,6,5]]]

function postprocess(node_values)
    dim = length(node_values)
    if dim == 1
        return node_values
    else 
        return sqrt(sum(node_values.^2))
    end 
end

function dof_to_node(dh::Ferrite.AbstractDofHandler, u::Array{T,1}; field::Int=1, process::Function=postprocess) where T
    fieldnames = Ferrite.getfieldnames(dh)  
    field_dim = Ferrite.getfielddim(dh, field)
    data = fill(NaN, Ferrite.getnnodes(dh.grid), field_dim) 
    offset = Ferrite.field_offset(dh, fieldnames[field])

    for cell in Ferrite.CellIterator(dh)
        _celldofs = Ferrite.celldofs(cell)
        counter = 1
        for node in cell.nodes
            for d in 1:field_dim
                data[node, d] = u[_celldofs[counter + offset]]
                counter += 1
            end
        end
    end
    return mapslices(process, data, dims=[2])
end
