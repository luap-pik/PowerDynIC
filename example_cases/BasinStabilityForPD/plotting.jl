using LaTeXStrings
using Plots
using GraphPlot, LightGraphs
using Cairo, Compose
using Makie
using Colors
gr()

function plot_res(result, powergrid,disturbed_node)
    ω_indices = findall(:ω .∈ symbolsof.(powergrid.nodes))

    #ω_colors = reshape(Plots.get_color_palette(cgrad(:inferno), 25)[ω_indices], (1,length(ω_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices], (1, length(ω_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))
    v_labels = reshape([latexstring(string(raw"v", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))

    pl_ω = Plots.plot(result, ω_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels)
    pl_p = Plots.plot(result, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_v = Plots.plot(result, :, :v, legend = (0.8, 0.7), ylabel=L"v [p.u.]", label=v_labels)

    pl = Plots.plot(
        pl_p, pl_ω,pl_v;
        layout=(3,1),
        size = (500, 500),
        lw=3,
        xaxis = (L"t[s]")
    )
    display(pl)
end

function distribution_Makie(f::ODEFunction, surface_func::Function, symbol_list, initial_con,
                            dist_args, dist::UnionAll, method::String, method_arg;
                            throws = 1000)

    x_vec = Array{Float64}(undef, throws)
    y_vec = Array{Float64}(undef, throws)
    z_vec = Array{Float64}(undef, throws)

    if method == "AmbientForcing"
        for i in 1:throws
            new_cond = AmbientForcing(f, initial_con, dist_args, dist, method_arg)
            x_vec[i] = new_cond[1]
            y_vec[i] = new_cond[2]
            z_vec[i] = new_cond[3]
        end
    end
    if method == "RandomWalk"
        for i in 1:throws
            new_cond = RandomWalkManifold(f, initial_con, dist_args, dist, method_arg)
            x_vec[i] = new_cond[1]
            y_vec[i] = new_cond[2]
            z_vec[i] = new_cond[3]
        end
    end

    x_ext = extrema(x_vec)
    y_ext = extrema(y_vec)
    z_ext = extrema(z_vec)

    scene = Scene(resolution = (1000, 1000));

    my_blue =  RGB{Float64}(102/255, 102/255, 255/255)
    my_red =   RGB{Float64}(255/255, 0/255, 127/255)
    const_color = cgrad([RGB{Float64}(153/255,204/255,255/255) for _ in 1:2 ])

    m_size = y_ext[2] * 40

    Makie.scatter!(scene, x_vec, y_vec, z_vec, color = my_blue, markersize = m_size, strokewidth = 0.0)
    Makie.scatter!(scene, [initial_con[1]],[initial_con[2]],[initial_con[3]], color = my_red, markersize = m_size * 2, strokewidth = 0.0)

    x = x_ext[1]:x_ext[2];
    y = y_ext[1]:y_ext[2];

    Makie.surface!(scene, x, y, surface_func, transparency = true, colormap = const_color, shading = false)


    axis = scene[Axis] # get axis
    axis[:names, :axisnames] = ("\\bf{y1}", "\\bf{y2}", "\\bf{y3}")
    tstyle = axis[:names] # or just get the nested attributes and work directly with them

    tstyle[:textsize] = 10
    tstyle[:font] =  "helvetica"

    wh = widths(scene)

    display(scene)
    return scene
end


"""
Produces a plot of the calulated basin stability against the nodes index.
Inputs:
    nodesarray: Array containg the index of all nodes that contain varible
    basinstability: Array containing the basin stability of each node
    standarddev: Array of the standard deviations of the basin stability
"""
function PlotBasinStability(nodesarray, basinstability, standarddev; labtext = "Basin Stability")
    Plots.plot(
        nodesarray,
        basinstability,
        xaxis = ("Node Index"),
        ylabel = "Basin Stability",
        marker = (:circle, 8, "red"),
        line = (:path, 2,"gray"),
        label = labtext,
        grid = false,
        show = true
    )
end

function PlotBasinStability!(nodesarray, basinstability, standarddev; labtext = "Basin Stability")
    Plots.plot!(
        nodesarray,
        basinstability,
        xaxis = ("Node Index"),
        ylabel = "Basin Stability",
        marker = (:circle, 8, "blue"),
        line = (:path, 2,"gray"),
        label = labtext,
        grid = false,
        show = true
    )
end

function plot_powergrid(pg)
    vec = pg.nodes |> collect
    types = typeof.(vec) # fetching only the node types
    membership = []

    for i in 1:length(types)
        type_str = string(types[i])
        if i == 1
            global type_dict = Dict([(type_str, 1)])
        end
        if get(type_dict, type_str, 0) == 0
            idx, last_type = findmax(type_dict)
            type_dict[string(types[i])] = idx +1
        end
        append!(membership, type_dict[type_str])
    end
    idx, type = findmax(type_dict)
    nodecolor = distinguishable_colors(idx, colorant"red");
    nodefillc = nodecolor[membership];

    nodelabel= string.(types)

    gplt = gplot(powergrid.graph, nodelabel = nodelabel, layout=circular_layout, nodefillc=nodefillc)
end
