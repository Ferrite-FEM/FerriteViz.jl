function solutionplot_pgfplots(plotter,args...;kwargs...)
    solution = reshape(transfer_solution(plotter,plotter.u[]; field_idx=1, process=norm), num_vertices(plotter))
    u_matrix = transfer_solution(plotter,plotter.u[]; field_idx=1, process=identity)
    coords = plotter.physical_coords .+ u_matrix
    @pgf PGFPlotsX.Axis({
        "scale mode" = "scale uniformly",
    }
    Plot3(
    {
        patch,
        "table/row sep" = "\\\\",
        "shader" = "flat",
        patch_table = TableData(plotter.triangles .- 1)
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
