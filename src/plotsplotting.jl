function solutionplot_pgfplots(plotter;
                               axis_options=@pgf{},
                               plot_options=@pgf{"shader" = "faceted"})
    solution = reshape(transfer_solution(plotter,plotter.u[]; field_idx=1, process=norm), num_vertices(plotter))
    u_matrix = transfer_solution(plotter,plotter.u[]; field_idx=1, process=identity)
    coords = plotter.physical_coords .+ u_matrix
    @pgf PGFPlotsX.Axis({axis_options..., "scale mode" = "scale uniformly"},
    Plot3(
    {
        plot_options...,
        patch,
        "axis equal image",
        "table/row sep" = "\\\\",
        patch_table = TableData(plotter.triangles .- 1),
    },
    Table(
        {
            point_meta = raw"\thisrow{c}"
        },
        :x => coords[:,1],
        :y => coords[:,2],
        :z => coords[:,3],
        :c => solution)))
end

function cellplot_pgfplots(plotter,celldata;
                           axis_options=@pgf{},
                           plot_options=@pgf{"shader" = "faceted"})
    solution = reshape(transfer_solution(plotter,plotter.u[]; field_idx=1, process=norm), num_vertices(plotter))
    u_matrix = transfer_solution(plotter,plotter.u[]; field_idx=1, process=identity)
    solution =  reshape(transfer_scalar_celldata(plotter, celldata; process=identity), num_vertices(plotter))
    coords = plotter.physical_coords .+ u_matrix
    @pgf PGFPlotsX.Axis({axis_options...,"scale mode" = "scale uniformly"},
    Plot3(
    {
        plot_options...,
        patch,
        "axis equal image",
        "table/row sep" = "\\\\",
        patch_table = TableData(plotter.triangles .- 1),
    },
    Table(
        {
            point_meta = raw"\thisrow{c}"
        },
        :x => coords[:,1],
        :y => coords[:,2],
        :z => coords[:,3],
        :c => solution)))
end

function wireframe_pgfplots(plotter::Plotter{dim};
                            axis_options=@pgf{},
                            node_options=@pgf{mark_size = "0.5pt"},
                            line_options=@pgf{opacity=0.3,"line width" = "0.2pt"}) where dim
    u_matrix = dof_to_node(plotter.dh, plotter.u[]; field=Ferrite.find_field(plotter.dh,:u), process=identity)
    coords = plotter.gridnodes .+ u_matrix
    figure = @pgf PGFPlotsX.Axis({axis_options...,"scale mode" = "scale uniformly"},Plot3(
    {
        node_options...,
        only_marks,
        "axis equal image",
    },
    Table(
        :x => coords[:,1],
        :y => coords[:,2],
        :z => coords[:,3])))
    for cell in Ferrite.getcells(plotter.dh.grid)
        boundaryentities = dim < 3 ? Ferrite.faces(cell) : Ferrite.edges(cell)
        cell_boundary_coords = [coords[e,:] for boundary in boundaryentities for e in boundary] |> x->reduce(hcat,x)' #TODO unnecessary complicated
        cellboundary = Table(:x=>cell_boundary_coords[:,1],:y=>cell_boundary_coords[:,2],:z=>cell_boundary_coords[:,3])
        boundaryplot = @pgf Plot3({line_options...,"mark=none"},cellboundary)
        push!(figure, boundaryplot)
    end
    return figure
end
