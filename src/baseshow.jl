using Makie
using MakiePlotter

function my_plot_layout(scene::AbstractPlotting.SceneLike)
    xaxis = scene.axes[1].axis
    yaxis = scene.axes[2].axis

    xticks = Makie.xticks(scene)
    yticks = Makie.yticks(scene)

    layout = @layout [
        a[Plot(scene)],
        b[EmptyBox()],
        c[XAxis(xaxis, ticks=xticks)],
        d[YAxis(yaxis, ticks=yticks)],
        e[EmptyBox()]
    ]

    return layout
end
