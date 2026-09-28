const _LINESTYLES = [:dash, :dashdot, :dot, :solid]

function _rankticks(rmax, base)
    exps = collect(0:floor(Int, log(base, max(1, rmax))))
    step = max(1, cld(length(exps), 7))
    exps = exps[1:step:end]
    return (Float64.(base .^ exps), [e == 0 ? L"1" : LaTeXString("$(base)^{$e}") for e in exps])
end

function _cutofflabel(c; cap_limited=false)
    e = iszero(c) ? Inf : -log10(c)
    expression = isfinite(e) && isapprox(e, round(e); atol=1e-10) ?
        "C = $(Int(round(e)))" : "\\epsilon = $(c)"
    cap_limited && (expression *= raw"\ \mathrm{(cap\ limited)}")
    return LaTeXString(expression)
end

"""
    plot_schmidt!(layout, spectrum; style=:paper2026, cutoffs=10.0.^(-2:-1:-5), ...)

Add a spectrum plot to a Makie GridLayout. The sole preset, `:paper2026`, uses
horizontal bonds, log₂ rank and inferno colors. Colors represent singular values,
not their squares. Stepped curves show independently computed retained-rank boundaries.

Keywords: `bond_labels`, `title`, `color_limits` (positive singular-value bounds),
`maxdim`, `legend`, `rasterize` (heatmap only; default false, vector PDF).
Values below `color_limits[1]` saturate at the lowest color; exact zeros are blank.
Each rank-j cell ends at j, so a contour at χ encloses precisely χ values.
The first cell begins at 1/2 because zero cannot appear on a log axis.
Text uses bundled Computer Modern fonts and mathematical labels use LaTeXStrings.
Titles are empty by default; no explanatory captions are added.
"""
function plot_schmidt!(layout, data::SchmidtSpectrum;
        style::Symbol=:paper2026, cutoffs=10.0 .^ (-2:-1:-5),
        maxdim::Integer=typemax(Int), title="", bond_labels=nothing,
        color_limits=nothing, legend::Bool=true, rasterize=false)
    style == :paper2026 || throw(ArgumentError("The only supported preset is :paper2026."))
    n = length(data.values)
    rmax = maximum(length, data.values)
    labels = isnothing(bond_labels) ? LaTeXString.(string.(1:n)) : collect(bond_labels)
    length(labels) == n || throw(ArgumentError("Provide one bond label per internal bond ($n)."))
    cutoffs = collect(cutoffs)
    length(cutoffs) <= 4 || throw(ArgumentError("At most four cutoff curves fit the reference preset."))
    profiles = [truncation_profile(data; cutoff=c, maxdim) for c in cutoffs]
    limits = isnothing(color_limits) ? (1e-6, 1.0) : color_limits
    length(limits) == 2 && all(isfinite, limits) && 0 < limits[1] < limits[2] ||
        throw(ArgumentError("color_limits must contain two finite, increasing positive values."))
    lo, hi = log10.(limits)
    displayrange = (lo,hi)
    cmap = :inferno
    rankticks = _rankticks(rmax, 2)
    fonts = theme_latexfonts().fonts
    regular = fonts.regular[]
    attrs = (; title, xgridvisible=false, ygridvisible=false, xtickalign=1, ytickalign=1,
        xticksmirrored=true, yticksmirrored=true, spinewidth=1.2, titlesize=19,
        xticklabelsize=15, yticklabelsize=15, titlefont=fonts.bold[],
        xlabelfont=regular, ylabelfont=regular, xticklabelfont=regular, yticklabelfont=regular)
    ax = Axis(layout[1,1]; attrs..., xlabel=L"\mathrm{Bond}\ k", ylabel=L"\mathrm{Schmidt\ index}\ j",
        yscale=log2, xticks=(collect(1:n), labels), yticks=rankticks,
        xticklabelrotation=isnothing(bond_labels) ? 0 : pi/2)
    limits!(ax, 0.5, n+0.5, 0.5, max(1, rmax))
    for (k,s) in enumerate(data.values)
        # Plot ragged columns directly, avoiding a padded N × max-rank allocation.
        rankedges = vcat(0.5, Float64.(1:length(s)))
        logvalues = [iszero(v) ? NaN : clamp(log10(v), displayrange...) for v in s]
        attrs_hm = (; colormap=cmap, colorrange=displayrange, nan_color=:transparent,
            interpolate=false, rasterize)
        heatmap!(ax, [k-0.5,k+0.5], rankedges, reshape(logvalues,1,:); attrs_hm...)
    end
    contours = []
    for (i,p) in enumerate(profiles)
        b = [x for k in 1:n for x in (k-0.5,k+0.5)]
        r = repeat(Float64.(p.ranks); inner=2)
        label = _cutofflabel(p.cutoff; cap_limited=any(!,p.meets_cutoff))
        line = lines!(ax,b,r; color=:black, linewidth=2, linestyle=_LINESTYLES[i], label)
        push!(contours, line)
    end
    tickstep = max(1, ceil(Int, (hi-lo)/4))
    cticks = sort!(unique(vcat(lo, collect(ceil(lo):tickstep:floor(hi)), hi)))
    ticklabels = [LaTeXString("10^{$(isinteger(t) ? string(Int(t)) : string(round(t;digits=2)))}") for t in cticks]
    colorlabel = data.normalized ? L"\frac{\sigma_{j,k}}{\sqrt{\sum_i\sigma_{i,k}^2}}" : L"\sigma_{j,k}"
    cbattrs = (; colormap=cmap, limits=displayrange, ticks=(cticks,ticklabels), label=colorlabel,
        ticklabelsize=14, labelsize=16, tickalign=1, spinewidth=1,
        ticklabelfont=regular, labelfont=regular)
    cb = Colorbar(layout[1,2]; cbattrs..., width=16)
    if legend && !isempty(contours)
        Legend(layout[2,1:2], contours,
            [_cutofflabel(p.cutoff; cap_limited=any(!,p.meets_cutoff)) for p in profiles];
            orientation=:horizontal, framevisible=false,
            labelsize=14, labelfont=regular, patchsize=(30,14), colgap=16, rowgap=3)
    end
    return (; axis=ax, colorbar=cb, profiles, contours)
end

"""Create a Makie Figure; `size` is the figure's width and height in points at PDF export."""
function plot_schmidt(data::SchmidtSpectrum; style::Symbol=:paper2026, size=nothing, kwargs...)
    style == :paper2026 || throw(ArgumentError("The only supported preset is :paper2026."))
    dims = isnothing(size) ? (850,430) : size
    fig = Figure(; size=dims, fontsize=18, fonts=theme_latexfonts().fonts)
    plot_schmidt!(fig[1,1],data; style,kwargs...)
    return fig
end
plot_schmidt(psi::MPS; kwargs...) = plot_schmidt(schmidt_spectrum(psi); kwargs...)

"""
    save_schmidt(path, spectrum_or_mps; kwargs...)

Save a publication-style PDF and return its absolute path. Parent directories
must exist. Only `.pdf` is accepted. Plot keywords are passed to `plot_schmidt`.
"""
function save_schmidt(path::AbstractString, data; kwargs...)
    lowercase(splitext(path)[2]) == ".pdf" || throw(ArgumentError("Output filename must end in .pdf."))
    fig = plot_schmidt(data; kwargs...)
    CairoMakie.save(path,fig; pt_per_unit=1)
    return abspath(path)
end
