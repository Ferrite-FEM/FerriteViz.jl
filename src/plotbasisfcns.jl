using Ferrite, GLMakie

function show_basis_function(ip::Union{Lagrange{RefQuadrilateral,N},Serendipity{RefQuadrilateral,N}}) where {N}
    x = range(-1,1,length=50)
    y = range(-1,1,length=50)
    #grid = generate_grid(Quadrilateral,(1,1))

    fig = Figure()

    # initialize axes
    rcs = Ferrite.reference_coordinates(ip)
    _refpos = rcs.*(N*0.5) # get and scale reference coordinates to unit intervals. Also mirror y coordinate
    shift = max(vcat(_refpos...)...)
    _refpos .+= [Tensors.Vec((shift+1,shift+1))]
    lngth = max(vcat(_refpos...)...)
    refpos = [Int.(rp) for rp in [[lngth-_rp[2]+1,_rp[1]] for _rp in _refpos]]
    ax = [Axis3(fig[pos...];title="node = "*string(round.(rcs[i];digits=2))) for (i,pos) in enumerate(refpos)]
    for i in 1:getnbasefunctions(ip)
        local z = [Ferrite.value(ip,i,Ferrite.Vec{2,Float64}((_x,_y))) for _x in x, _y in y]
        surface!(ax[i],x,y,z)
        #FerriteViz.wireframe!(grid,markersize=10,strokewidth=2)
    end
    current_figure()
end

function show_basis_function(ip::Union{Lagrange{RefTriangle,N,unused},BubbleEnrichedLagrange{RefTriangle,N,unused}}) where {N,unused}
    rcs = Ferrite.reference_coordinates(ip)
    lwr = min(vcat(collect.(rcs)...)...)
    upr = max(vcat(collect.(rcs)...)...)
    y = range(upr,lwr,length=20)
    x = range(lwr,upr,length=20)

    # initialize axes
    if typeof(ip)==Lagrange{RefTriangle,N,unused}
        __refpos = [coord.*[1, 1] for coord in rcs.*(N)] # get and scale reference coordinates to unit intervals. Also mirror y coordinate
        shift = min(vcat(__refpos...)...)
        __refpos .+= [Tensors.Vec((shift+1,shift+1))]
        _refpos = [Int.(reverse(rp)) for rp in __refpos]
        lngth = max(vcat(_refpos...)...)
        refpos = [[lngth-rp[1] ,rp[2]] for rp in _refpos]
    elseif typeof(ip)==BubbleEnrichedLagrange{RefTriangle,1,unused}
        refpos = [[3,3],[1,1],[3,1],[2,2]]
    else
        throw(ArgumentError("method not implemented for type $(typeof(ip))"))
    end    
    fig = Figure()
    ax = [Axis3(fig[pos...];title="node = "*string(round.(rcs[i];digits=2))) for (i,pos) in enumerate(refpos)]
    for i in 1:getnbasefunctions(ip)
        local z = [Ferrite.value(ip,i,Ferrite.Vec{2,Float64}((_x,_y))) for _y in y, _x in x]
        for i in 1:length(x), j in 1:length(y)
            i>j ? z[j,i]=NaN : nothing
        end
        vertices, clist = get_triangulation(x,y,z)
        #surface!(ax[i],x,y,z)
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
