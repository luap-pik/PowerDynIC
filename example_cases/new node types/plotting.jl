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
    append!(ω_indices,findall(n -> isa(n, ThirdOrderEqGovernorIEEEGone),powergrid.nodes))


    ω_colors = reshape(Plots.get_color_palette(cgrad(:inferno), 25)[ω_indices], (1,length(ω_indices)))
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
