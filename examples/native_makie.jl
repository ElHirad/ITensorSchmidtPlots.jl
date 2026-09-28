using ITensors, ITensorMPS, ITensorSchmidtPlots
using CairoMakie, LaTeXStrings, Random

sites = siteinds("Qubit", 10)
psi = random_mps(MersenneTwister(42), ComplexF64, sites; linkdims=16)
data = schmidt_spectrum(psi)

fig, ax, hm = schmidtplot(data;
    colormap=:inferno, colorrange=(1e-5, 1.0),
    axis=(xticks=1:9, yticks=([1,2,4,8,16], [L"1",L"2",L"4",L"8",L"16"]),
        xtickalign=1, ytickalign=1),
    figure=(size=(850,430), fontsize=18, fonts=theme_latexfonts().fonts))
limits!(ax, 0.5, 9.5, 0.5, 16)

Colorbar(fig[1,2], hm; label=L"\sigma_{j,k}", width=16,
    ticks=([1e-5,1e-3,1e-1,1.0], [L"10^{-5}",L"10^{-3}",L"10^{-1}",L"1"]))

# Each curve is an ordinary Makie plot with independent styling.
styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
for (cutoff, style) in zip(10.0 .^ (-1:-1:-5), styles)
    schmidttruncation!(ax, data; cutoff, linestyle=style,
        color=:black, linewidth=2)
end
Legend(fig[2,1:2], ax; orientation=:horizontal, framevisible=false,
    patchsize=(35,14))

out = joinpath(@__DIR__, "..", "output")
mkpath(out)
save(joinpath(out,"native_makie.pdf"), fig; pt_per_unit=1)
