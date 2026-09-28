import Makie

@testset "Native Makie composition and reactive customization" begin
    data = SchmidtSpectrum([[1.], [0.8,0.6,0.], [1.]]; normalized=true)
    observable = Observable(data)
    fig, ax, hm = schmidtplot(observable;
        colormap=:inferno, colorrange=(1e-5,1.0),
        axis=(xticks=[1,3], ylabel="custom rank"),
        figure=(size=(700,400), fonts=theme_latexfonts().fonts))
    @test hm isa SchmidtPlot
    @test ax.xticks[] == [1,3]
    @test ax.ylabel[] == "custom rank"
    @test ax.yscale[] === log2
    @test isempty(ax.title[])
    @test length(fig.content) == 1 # No implicit caption, legend, or colorbar.

    cb = Colorbar(fig[1,2],hm; ticks=([1e-5,1.0],["low","high"]))
    @test cb.limits[] ≈ [1e-5,1.0] rtol=eps(Float32)
    @test cb.scale[] === log10
    hm.colormap = :viridis
    hm.colorrange = (1e-4,1.0)
    @test Makie.to_colormap(cb.colormap[]) == Makie.to_colormap(:viridis)
    @test cb.limits[] ≈ [1e-4,1.0] rtol=eps(Float32)
    @test cb.ticks[] == ([1e-5,1.0],["low","high"])
    @test isnan(hm.plots[2][3][][1,3])
    @test hm.plots[2][3][][1,1:2] ≈ [0.8,0.6] # True values, not log values.

    curve = schmidttruncation!(ax,observable;
        cutoff=0.1, linewidth=3, color=:cyan, label="custom cutoff")
    @test curve isa SchmidtTruncation
    @test curve.profile[].ranks == [1,2,1]
    @test only(Makie.get_plots(curve)).label[] == "custom cutoff"
    curve.cutoff = 0.4
    @test curve.profile[].ranks == [1,1,1]
    curve.linewidth = 4
    @test only(curve.plots).linewidth[] == 4
    @test axislegend(ax; position=:lt) isa Legend

    # A rank-changing update must keep heatmap coordinates, cells, and curves aligned.
    observable[] = SchmidtSpectrum([[1.,1.], [1.,1.,1.,1.], [1.]];normalized=true)
    @test size(hm.plots[2][3][]) == (1,4)
    @test hm.plots[2][2][] == [0.5,1,2,3,4]
    @test hm.plots[2][3][] ≈ fill(0.5,1,4)
    @test curve.profile[].ranks == [2,3,1]
    curve.maxdim = 1
    @test curve.profile[].ranks == [1,1,1]
    @test curve.profile[].meets_cutoff == [false,false,true]

    vertical_fig = Figure(fonts=theme_latexfonts().fonts)
    vax, vertical = schmidtplot(vertical_fig[1,1], data;orientation=:vertical)
    @test vax.xscale[] === log2
    @test vax.yreversed[]
    @test vax.xaxisposition[] == :top
    @test size(vertical.plots[2][3][]) == (3,1)
    vc = schmidttruncation!(vax,data;orientation=:vertical,cutoff=0.1)
    @test vc.positions[][3:4] == [Point2d(2,1.5),Point2d(2,2.5)]
    existing = Axis(vertical_fig[1,2];xlabel="preserved",yscale=identity)
    @test schmidtplot!(existing,data) isa SchmidtPlot
    @test existing.xlabel[] == "preserved"
    @test existing.yscale[] === identity

    # No four-curve preset limit in the native API; legends use ordinary line plots.
    for cutoff in (0.1,0.01,0.001,0.0001,0.00001)
        schmidttruncation!(existing,data;cutoff)
    end
    @test Legend(vertical_fig[2,1:2],existing;orientation=:horizontal) isa Legend
    psi = MPS(siteinds("Qubit",3),"0")
    @test schmidtplot!(existing,psi).spectrum[].values == [[1.],[1.]]
    @test schmidttruncation!(existing,psi).profile[].ranks == [1,1]

    mktempdir() do dir
        for (name,f) in (("native",fig),("composed",vertical_fig))
            path = joinpath(dir,name*".pdf")
            save(path,f)
            @test startswith(read(path,String),"%PDF")
            @test filesize(path) > 1000
        end
    end
end
