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

function wireframe_pgfplots(plotter,args...;kwargs...) 
    u_matrix = transfer_solution(plotter,plotter.u[]; field_idx=1, process=identity)
    coords = plotter.physical_coords .+ u_matrix
    @pgf PGFPlotsX.Axis(Plot3(
    {
        only_marks,
    },
    Table(
        :x => coords[:,1],
        :y => coords[:,2],
        :z => coords[:,3])))
end
