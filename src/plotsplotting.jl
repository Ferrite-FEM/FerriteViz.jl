function solutionplot_pgfplots!(figure,plotter::Plotter{dim};
                      field_idx=1,
                      process=norm,
                      scale=1.0,
                      deformation_field=:default,
                      plot_options=@pgf{"shader" = "interp"}) where dim
    solution = reshape(transfer_solution(plotter,plotter.u[]; field_idx=field_idx, process=process), num_vertices(plotter))
    u_matrix = deformation_field != :default ? transfer_solution(plotter,plotter.u[]; field_idx=Ferrite.find_field(plotter.dh,deformation_field), process=identity) : zeros(0,3)
    coords = deformation_field != :default ? plotter.physical_coords .+ scale .* u_matrix : plotter.physical_coords
    if dim == 3
        plot = @pgf Plot3(
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
            :c => solution))
        push!(figure,plot)
    elseif dim == 2
        plot = @pgf Plot(
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
            :c => solution))
        push!(figure,plot)
    end
end

function solutionplot_pgfplots(plotter::Plotter;
                  field_idx=1,
                  process=norm,
                  scale=1.0,
                  deformation_field=:default,
                  axis_options=@pgf{},
                  plot_options=@pgf{"shader" = "interp"})
    figure = @pgf PGFPlotsX.Axis({axis_options...,"scale mode" = "scale uniformly"})
    solutionplot_pgfplots!(figure,plotter;field_idx,process,scale,deformation_field,plot_options)
end

function cellplot_pgfplots!(figure,plotter::Plotter{dim},celldata;
                  process=identity,
                  scale=1.0,
                  deformation_field=:default,
                  plot_options=@pgf{"shader" = "interp"}) where dim
    u_matrix = deformation_field != :default ? transfer_solution(plotter,plotter.u[]; field_idx=Ferrite.find_field(plotter.dh,deformation_field), process=identity) : zeros(0,3)
    coords = deformation_field != :default ? plotter.physical_coords .+ scale .* u_matrix : plotter.physical_coords
    solution = reshape(transfer_scalar_celldata(plotter, celldata; process=process), num_vertices(plotter))
    if dim == 3
        plot = @pgf Plot3(
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
            :c => solution))
        push!(figure,plot)
    elseif dim == 2
        plot = @pgf Plot(
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
            :c => solution))
        push!(figure,plot)
    end
end

function cellplot_pgfplots(plotter::Plotter,celldata;
                  process=identity,
                  scale=1.0,
                  deformation_field=:default,
                  axis_options=@pgf{},
                  plot_options=@pgf{"shader" = "interp"})
    figure = @pgf PGFPlotsX.Axis({axis_options...,"scale mode" = "scale uniformly"})
    cellplot_pgfplots!(figure,plotter,celldata;process,scale,deformation_field,plot_options)
end

function wireframe_pgfplots!(figure, plotter::Plotter{dim};
                   scale=1.0,
                   deformation_field=:default,
                   node_options=@pgf{mark_size = "0.5pt"},
                   line_options=@pgf{opacity=0.3,"line width" = "0.2pt"}) where dim
    u_matrix = deformation_field != :default ? dof_to_node(plotter.dh, plotter.u[]; field=Ferrite.find_field(plotter.dh,deformation_field), process=identity) : zeros(0,3)
    coords = deformation_field != :default ? plotter.gridnodes .+ scale .* u_matrix : plotter.gridnodes
    if dim == 3
        nodeplot = @pgf Plot3(
        {
            node_options...,
            only_marks,
            "axis equal image",
        },
        Table(
            :x => coords[:,1],
            :y => coords[:,2],
            :z => coords[:,3]))
        push!(figure,nodeplot)
        for cell in Ferrite.getcells(plotter.dh.grid)
            boundaryentities = Ferrite.edges(cell)
            cell_boundary_coords = [coords[e,:] for boundary in boundaryentities for e in boundary] |> x->reduce(hcat,x)' #TODO unnecessary complicated
            cellboundary = Table(:x=>cell_boundary_coords[:,1],:y=>cell_boundary_coords[:,2],:z=>cell_boundary_coords[:,3])
            boundaryplot = @pgf Plot3({line_options...,"mark=none"},cellboundary)
            push!(figure, boundaryplot)
        end
        return figure
    elseif dim == 2
        nodeplot = @pgf Plot(
        {
            node_options...,
            only_marks,
            "axis equal image",
        },
        Table(
            :x => coords[:,1],
            :y => coords[:,2]))
        push!(figure,nodeplot)
        for cell in Ferrite.getcells(plotter.dh.grid)
            boundaryentities = Ferrite.faces(cell)
            cell_boundary_coords = [coords[e,:] for boundary in boundaryentities for e in boundary] |> x->reduce(hcat,x)' #TODO unnecessary complicated
            cellboundary = Table(:x=>cell_boundary_coords[:,1],:y=>cell_boundary_coords[:,2])
            boundaryplot = @pgf Plot({line_options...,"mark=none"},cellboundary)
            push!(figure, boundaryplot)
        end
        return figure
    end
end

function wireframe_pgfplots(plotter::Plotter;
                   scale=1.0,
                   deformation_field=:default,
                   axis_options = @pgf{},
                   node_options=@pgf{mark_size = "0.5pt"},
                   line_options=@pgf{opacity=0.3,"line width" = "0.2pt"})
    figure = @pgf PGFPlotsX.Axis({axis_options...,"scale mode" = "scale uniformly"})
    wireframe_pgfplots!(figure,plotter;scale,deformation_field,node_options,line_options) 
end
