"""
    schmidtplot(spectrum_or_mps; axis=(...), figure=(...), kwargs...)
    schmidtplot!(axis, spectrum_or_mps; kwargs...)

Native Makie Schmidt-spectrum heatmap. Returns the usual FigureAxisPlot or Plot
object. Use `Colorbar(fig[position], plot)` and `schmidttruncation!` to add a
colorbar and cutoff curves. Colors are actual singular values, with `colorscale=log10`
by default; colorbar tick positions and `colorrange` are in singular-value units.

Supports ordinary Makie colormap, colorrange, colorscale, alpha, clipping, visibility,
rasterization, axis, figure, layout, and theme controls. Zero coefficients are blank.
`orientation` is `:horizontal` (bonds on x) or `:vertical` (bonds on y).
Use `theme_latexfonts()` for Computer Modern fonts, as in the supplied examples.

Observable data updates may change ranks but must preserve the number of bonds.
Create a new plot when changing the number of sites. Mutating a stored vector in
place requires notifying its Observable, following ordinary Makie conventions.
"""
Makie.@recipe SchmidtPlot (spectrum,) begin
    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_colormap_attributes()...
    "Whether bonds are plotted horizontally or vertically."
    orientation = :horizontal
    colormap = :inferno
    colorrange = (1e-6, 1.0)
    colorscale = log10
    "Rasterize the heatmap in vector output; labels remain vectors."
    rasterize = false
end

Makie.convert_arguments(::Type{<:SchmidtPlot}, data::SchmidtSpectrum) = (data,)
Makie.convert_arguments(::Type{<:SchmidtPlot}, psi::MPS) = (schmidt_spectrum(psi),)

function _horizontal(orientation)
    orientation in (:horizontal, :vertical) ||
        throw(ArgumentError("orientation must be :horizontal or :vertical."))
    return orientation == :horizontal
end

function Makie.plot!(p::SchmidtPlot)
    n = length(p.spectrum[].values)
    # Keep each ragged bond as a separate heatmap instead of allocating N × maxrank.
    # Compute all three coordinates together so rank changes update atomically.
    for k in 1:n
        names = [Symbol("bond_x_",k), Symbol("bond_y_",k), Symbol("bond_z_",k)]
        Makie.map!(p, [:spectrum, :orientation], names) do data, orientation
            length(data.values) == n || throw(ArgumentError("Data updates must preserve the number of bonds; create a new plot for a different chain length."))
            horizontal = _horizontal(orientation)
            s = data.values[k]
            edges = vcat(0.5, Float64.(1:length(s)))
            values = [iszero(v) ? NaN : v for v in s]
            bonds = [k-0.5,k+0.5]
            return horizontal ? (bonds, edges, reshape(values,1,:)) :
                (edges, bonds, reshape(values,:,1))
        end
        heatmap!(p, p[names[1]], p[names[2]], p[names[3]];
            colormap=p.colormap, colorrange=p.colorrange, colorscale=p.colorscale,
            lowclip=p.lowclip, highclip=p.highclip, nan_color=p.nan_color,
            alpha=p.alpha, visible=p.visible, transparency=p.transparency,
            inspectable=p.inspectable, rasterize=p.rasterize, interpolate=false)
    end
    return p
end

# All columns share the same color attributes, so one column defines the colorbar.
Makie.extract_colormap(p::SchmidtPlot) = Makie.extract_colormap(first(p.plots))

function Makie.preferred_axis_attributes(::Type{<:Axis}, p::SchmidtPlot)
    horizontal = _horizontal(p.orientation[])
    n = length(p.spectrum[].values)
    ranks = _rankticks(maximum(length,p.spectrum[].values), horizontal ? 2 : 4)
    bonds = (collect(1:n), LaTeXString.(string.(1:n)))
    fonts = theme_latexfonts().fonts
    common = (; xgridvisible=false, ygridvisible=false,
        xlabelfont=fonts.regular[], ylabelfont=fonts.regular[],
        xticklabelfont=fonts.regular[], yticklabelfont=fonts.regular[])
    return horizontal ?
        (; common..., xlabel=L"\mathrm{Bond}\ k", ylabel=L"\mathrm{Schmidt\ index}\ j",
            yscale=log2, xticks=bonds, yticks=ranks) :
        (; common..., xlabel=L"\mathrm{Schmidt\ index}\ \alpha", ylabel=L"\mathrm{Bond}\ n",
            xscale=log2, yreversed=true, xaxisposition=:top, xticks=ranks, yticks=bonds)
end

"""
    schmidttruncation!(axis, spectrum_or_mps; cutoff=1e-4, maxdim=typemax(Int), kwargs...)

Native Makie line recipe for the minimum retained rank at each bond.
Accepts standard `color`, `linewidth`, `linestyle`, `label`, and other line
attributes. Add as many curves as needed and use `axislegend(axis)` or
`Legend(fig[position], axis)` to place their legend.
The default label shows the cutoff and indicates a restrictive `maxdim` cap.
Data, cutoff, and maxdim can be updated using normal Makie Observables/attributes.
"""
Makie.@recipe SchmidtTruncation (spectrum,) begin
    Makie.filtered_attributes(Makie.Lines, exclude=(:cycle,))...
    cutoff = 1e-4
    maxdim = typemax(Int)
    orientation = :horizontal
    color = :black
    linewidth = 2
    linestyle = :dash
    label = Makie.automatic
end

Makie.convert_arguments(::Type{<:SchmidtTruncation}, data::SchmidtSpectrum) = (data,)
Makie.convert_arguments(::Type{<:SchmidtTruncation}, psi::MPS) = (schmidt_spectrum(psi),)

function Makie.plot!(p::SchmidtTruncation)
    Makie.map!(p, [:spectrum, :cutoff, :maxdim], :profile) do data, cutoff, maxdim
        return truncation_profile(data; cutoff, maxdim)
    end
    Makie.map!(p, [:profile, :orientation], :positions) do profile, orientation
        horizontal = _horizontal(orientation)
        return [horizontal ? Point2d(b,r) : Point2d(r,b)
            for (k,r) in enumerate(profile.ranks) for b in (k-0.5,k+0.5)]
    end
    Makie.map!(p, [:profile, :label], :curve_label) do profile, label
        return label === Makie.automatic ?
            _cutofflabel(profile.cutoff; cap_limited=any(!,profile.meets_cutoff)) : label
    end
    lines!(p, p.attributes, p.positions; label=p.curve_label)
    return p
end

# Expose the labeled line to Makie's ordinary axislegend/Legend machinery.
Makie.get_plots(p::SchmidtTruncation) = p.plots
