using LaTeXStrings
using Plots
gr()

function plot_res(result, powergrid,disturbed_node)
    ω_indices = findall(n -> isa(n, SwingEqLVS), powergrid.nodes)
    append!(ω_indices,findall(n -> isa(n, VSIMinimal), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, PVInverterWithFrequencyControl), powergrid.nodes))
    #append!(ω_indices,findall(n -> isa(n, WindTurbineGenType4_RotorControl), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, FourthOrderEq), powergrid.nodes))
    append!(ω_indices,findall(n -> isa(n, ThirdOrderEq),powergrid.nodes))


    #ω_colors = reshape(Plots.get_color_palette(cgrad(:inferno), 25)[ω_indices], (1,length(ω_indices)))
    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices], (1, length(ω_indices)))
    p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))
    v_labels = reshape([latexstring(string(raw"v", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))

    pl_p = plot(result, :, :p, legend = (0.8, 0.95), ylabel=L"p [p.u.]", label=p_labels)
    pl_v = plot(result, :, :v, legend = (0.8, 0.7), ylabel=L"v [p.u.]", label=v_labels)
    pl_ω = plot(result, ω_indices, :ω, legend = (0.8, 0.7), ylabel=L"\omega \left[rad/s\right]", label=ω_labels)
    pl = plot(
        pl_p, pl_ω,pl_v;
        layout=(3,1),
        size = (500, 500),
        lw=3,
        xguide=L"t[s]"
    )
    display(pl)
end

function plot_distribution(g::Function, symbol_list,fixpoint, op::State,
                            Distribution::UnionAll, Proj::String , node::Int, throws = 1000)
    x_vec = Array{Float64}(undef, throws)
    y_vec = Array{Float64}(undef, throws)
    z_vec = Array{Float64}(undef, throws)
    for i in 1:throws
        dz = RandPertWithConstrains(g, symbol_list,
                                    fixpoint, op, node,
                                    [-π, π],[-100,100], Distribution;
                                    ProjectionMethod = Proj, nsteps = 1)
            x_vec[i] = dz[1]
            y_vec[i] = dz[2]
            z_vec[i] = dz[3]
    end

    plot(x_vec, y_vec, z_vec,
         zcolor = reverse(z_vec), seriestype = :scatter,
         m = (5, 0.8, :blues, Plots.stroke(0)),
         label = "Pertubations", cbar = true, w = 5,
         title ="Distribution: " * string(Distribution) *
         ", Projection Method: " * string(Proj) *
         ", Node: " * string(node))

    plot!([fixpoint[1]],[fixpoint[2]], [fixpoint[3]],
          seriestype = :scatter,
          m = (10, 0.8, :red, Plots.stroke(0)),
          label = "Fix Point")

    xlabel!(string(symbol_list[1]))
    plt = ylabel!(string(symbol_list[2]))
    display(plt)
    extrem_y = y_vec |> extrema
    extrem_z = z_vec |> extrema

    return extrem_y, extrem_z
end


function plot_distance(g::Function, symbol_list, fixpoint, op::State,
                            Distribution::UnionAll, Proj::String , node::Int, throws = 1000)
    x_vec = Array{Float64}(undef, throws)
    y_vec = Array{Float64}(undef, throws)
    z_vec = Array{Float64}(undef, throws)
    for i in 1:throws
            dz = RandPertWithConstrains(g, symbol_list,
                                        fixpoint, op, node,
                                        [-π, π],[-100,100], Distribution;
                                        ProjectionMethod = Proj, nsteps = 1)

            x_vec[i] = dz[1] - fixpoint[1]
            y_vec[i] = dz[2] - fixpoint[2]
            z_vec[i] = dz[3] - fixpoint[3]
    end

    plot(x_vec, y_vec, z_vec,
         zcolor = reverse(z_vec), seriestype = :scatter,
         m = (5, 0.8, :blues, Plots.stroke(0)),
         label = "Pertubations", cbar = true, w = 5,
         title ="Distribution: " * string(Distribution) *
         ", Projection Method: " * string(Proj) *
         ", Node: " * string(node))

    plot!([0],[0],[0],
          seriestype = :scatter,
          m = (10, 0.8, :red, Plots.stroke(0)),
          label = "Fix Point")

    xlabel!(string(symbol_list[1]))
    ylabel!(string(symbol_list[2]))

    return max(y_vec), min(y_vec), max(z_vec), min(z_vec)
end
