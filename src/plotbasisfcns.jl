using Ferrite, GLMakie

function get_domain(ip::Interpolation{ref_shape,order}, meshsize::Int) where {ref_shape<:Union{RefQuadrilateral,RefTriangle},order}
    rcs = Ferrite.reference_coordinates(ip)
    lwr = min(vcat(collect.(rcs)...)...)
    upr = max(vcat(collect.(rcs)...)...)
    y = range(upr,lwr,length=meshsize)
    x = range(lwr,upr,length=meshsize)
    ref_z = get_ref_z(ip,meshsize)
    return x,y,ref_z
end

get_ref_z(ip::Interpolation{ref_shape,order},meshsize::Int) where {ref_shape,order} = zeros(meshsize,meshsize)
function get_ref_z(ip::Interpolation{ref_shape,order},meshsize::Int) where {ref_shape<:RefTriangle,order}
    z = zeros(meshsize,meshsize)
    for ix in 1:size(z,2), iy in 1:size(z,1)
        (ix>iy) ? z[iy,ix]=NaN : nothing
    end
    println(typeof(ip))
    if typeof(ip)==CrouzeixRaviart{RefTriangle,1,Nothing}
        z = collect(z')
        println("i dont work")
    end
    display(z)
    return z
end

function initialize_figure(ip::Interpolation{ref_shape,N}) where {ref_shape<:Union{RefQuadrilateral,RefTriangle},N}
    rcs = Ferrite.reference_coordinates(ip)
    min_val = Int(min(vcat(rcs...)...))
    rcs = [[c[1]-min_val,c[2]-min_val] for c in rcs]
    scale = round(inv(min(filter(!iszero,vcat(rcs...))...));digits=15)
    max_id = Int(round(max(vcat(rcs...)...)*scale))
    refpos = [Int.([max_id-coord[2]+1, coord[1]+1]) for coord in rcs.*scale]

    fig = Figure()
    ax = [Axis3(fig[pos...];title="node = "*string(round.(rcs[i];digits=2))) for (i,pos) in enumerate(refpos)]
    return fig, ax
end

#function initialize_figure(ip::Interpolation{ref_shape,N,unused}) where {ref_shape<:RefTriangle,N,unused}
#    rcs = Ferrite.reference_coordinates(ip)
#    scale = abs(inv(min(filter(!iszero,vcat(rcs...))...)))
#    max_id = Int(max(vcat(rcs...)...)*scale)
#    refpos = [Int.([max_id-coord[2]+1, coord[1]+1]) for coord in rcs.*scale]
#  
#    fig = Figure()
#    ax = [Axis3(fig[pos...];title="node = "*string(round.(rcs[i];digits=2))) for (i,pos) in enumerate(refpos)]
#    return fig,ax
#end

function show_basis_function(ip::Interpolation{ref_shape,N}) where {ref_shape<:Union{RefTriangle,RefQuadrilateral},N}
    x,y,ref_z = get_domain(ip,20)

    # initialize axes
    fig,ax = initialize_figure(ip)

    for i in 1:getnbasefunctions(ip)
        local z = copy(ref_z)
        for (iy,_y) in enumerate(y), (ix,_x) in enumerate(x)
            !isnan(z[iy,ix]) && (z[iy,ix]=Ferrite.value(ip,i,Ferrite.Vec{2,Float64}((_x,_y))))
        end
        vertices, clist = get_triangulation(x,y,z)
        mesh!(ax[i],vertices,clist; shading=false, fxaa=true, transparency=false, color=[vertices[i,3] for i in 1:size(vertices)[1]], colormap=:viridis)
    end
    current_figure()
end

function show_basis_function(ip::Lagrange{RefLine,N}) where {N}
    x = range(-1,1,length=100)
    fig = Figure()

    # initialize axes
    _refpos = [coord.*[1] for coord in Ferrite.reference_coordinates(ip).*(N*0.5)] # get and scale reference coordinates to unit intervals. Also mirror y coordinate
    shift = max(vcat(_refpos...)...)
    _refpos .+= [Tensors.Vec((shift+1))]
    refpos = [Int.(rp) for rp in _refpos]
    @show refpos
    ax = [Axis(fig[1,i]) for i in 1:getnbasefunctions(ip)]
    for i in 1:getnbasefunctions(ip)
        local z = [Ferrite.value(ip,i,Ferrite.Vec{1,Float64}((_x,))) for _x in x]
        lines!(ax[refpos[i]...],x,z)
    end
    current_figure()
end

function get_triangulation(x,y,z::Matrix)
any(isnan.(z[end,:])) ? error("case is not implemented") : nothing
    sz = size(z)
    nvertices_percol = (.!isnan.(z)')*ones(Int,sz[2])
    nvertices = sum(nvertices_percol)
# ===============================
# ===== compute vertice list =====
# ================================
#   vertice numbering according to numbering:
#   ┌                   ┐      ┌                   ┐                ┌                   ┐
#   │ 1 NaN NaN NaN NaN │      │ 1 NaN NaN NaN NaN │  not allowed:  │ 1 NaN NaN NaN NaN │
#   │ 2  6  NaN NaN NaN │      │ 2  6  NaN 13  17  │ -connectivity  │ 2  6  NaN 13  16  │
#   │ 3  7  10  NaN NaN │  or  │ 3  7  10  14  18  │  list is not   │ 3  7  10  14  17  │
#   │ 4  8  11  13  NaN │      │ 4  8  11  15  19  │  implemented   │ 4  8  11  15  18  │
#   │ 5  9  12  14  15  │      │ 5  9  12  16  20  │  for this      │ 5  9  12  NaN NaN │
#   └                          └                                    └                    
    vertices = zeros(nvertices,3)
    cnt = 1
    for col in 1:sz[2], row in 1:sz[1]
        if !isnan(z[row,col]) 
            vertices[cnt,:] = [x[col], y[row], z[row,col]]
            cnt += 1
        end
    end
# =====================================    
# ===== compute connectivity list =====
# =====================================
    ntriangles = 0 # number of triangles necessary
    for (i,vertincol) in enumerate(nvertices_percol)
        if i!=length(nvertices_percol)
            ntriangles+=(vertincol-1)*2 # if current collum has same length as next one add 2 triangles per "gap" between two nodes
            if vertincol==nvertices_percol[i+1]-1 # if next col is larger add 1
                ntriangles+=1
            elseif vertincol==nvertices_percol[i+1]+1 # if next col is smaller subtract 1
                ntriangles+=-1
            else
                !(vertincol==nvertices_percol[i+1]) ? error("triangulation not possible") : nothing
            end
        end
    end
    clist = zeros(Int,ntriangles,3)

    cnt = 1
    lastvertidpercol = [sum(nvertices_percol[1:c]) for c in 1:sz[2]]
    for vert in 1:lastvertidpercol[end-1]
        newcol = vert in vcat(1,lastvertidpercol[1:end-1].+1) # vertex first id in a column??
        col = findfirst(i->vert<=i,lastvertidpercol)
        vert==sum(nvertices_percol[1:col]) && continue # skip last vertex per column
        # fill clist
        offset = nvertices_percol[col]==nvertices_percol[col+1]-1 ? nvertices_percol[col]+1 : nvertices_percol[col] #set ofset according to next collumn size
        clist[cnt,:] = [vert, vert+1, vert+offset]
        cnt+=1
        if ((nvertices_percol[col]==nvertices_percol[col+1]+1)&&(!newcol)) || ((nvertices_percol[col]==nvertices_percol[col+1]-1)&&(newcol)) # if next column is shorter && not first vertex of collumn or if next column is longer and first vertex in collumn
            clist[cnt,:] = [vert, vert+offset-1, vert+offset]
            cnt+=1
        end
        if (nvertices_percol[col]==nvertices_percol[col+1]) || (nvertices_percol[col]==nvertices_percol[col+1]-1) # if next row is longer or shorter
            clist[cnt,:] = [vert+1, vert+offset+1, vert+offset]
            cnt+=1
        end
    end
    return vertices, clist
end
