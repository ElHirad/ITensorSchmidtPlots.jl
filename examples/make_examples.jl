using ITensorSchmidtPlots, CairoMakie, LaTeXStrings, LinearAlgebra
include("states.jl")

BLAS.set_num_threads(1)
out = joinpath(@__DIR__, "..", "output", "pdf")
mkpath(out)
report = ["# Numerical validation", "", "Generated with Julia $(VERSION), ITensors $(pkgversion(ITensors)), ITensorMPS $(pkgversion(ITensorMPS)), CairoMakie $(pkgversion(CairoMakie)).",
    "Every internal bond is compared with a direct dense SVD of the original vector.",
    "Truncation ranks are independently checked against dense-SVD squared tails.", "",
    "| Example | Max absolute SVD error | Analytic error | Largest tail error | Rank checks |",
    "| --- | ---: | ---: | ---: | ---: |"]

for case in (pair_example(), field_example(), random_example())
    psi = state_mps(case.v,case.n)
    data = schmidt_spectrum(psi)
    dense = dense_spectra(case.v,case.n)
    padded = [vcat(s,zeros(length(dense[k])-length(s))) for (k,s) in enumerate(data.values)]
    error = maximum(maximum(abs.(padded[k]-dense[k])) for k in eachindex(dense))
    analytic_error = isnothing(case.analytic) ? nothing :
        maximum(maximum(abs.(padded[k]-case.analytic[k])) for k in eachindex(dense))
    @assert error < 5e-13
    @assert isnothing(analytic_error) || analytic_error < 5e-13
    tail_error, checks = 0.0, 0
    profiles = [truncation_profile(data;cutoff=10.0^-c) for c in 2:5]
    for p in profiles, k in eachindex(dense)
        r=p.ranks[k]
        independent_tail=sum(abs2,dense[k][r+1:end])
        independent_rank=findfirst(j -> sum(abs2,dense[k][j+1:end]) <= p.cutoff,1:length(dense[k]))
        @assert r == independent_rank
        @assert independent_tail <= p.cutoff + 1e-13
        @assert r == 1 || sum(abs2,dense[k][r:end]) > p.cutoff - 1e-13
        tail_error=max(tail_error,abs(independent_tail-p.discarded_weights[k]))
        checks+=1
    end
    # Store the plotted values and exact truncation choices for reproducibility.
    open(joinpath(out,case.name*"_values.csv"),"w") do io
        println(io,"bond,schmidt_index,singular_value,dense_svd_value")
        for k in eachindex(dense),j in eachindex(padded[k])
            println(io,"$k,$j,$(padded[k][j]),$(dense[k][j])")
        end
    end
    open(joinpath(out,case.name*"_truncation.csv"),"w") do io
        println(io,"bond,C,retained_rank,discarded_squared_weight,cutoff_satisfied")
        for (c,p) in zip(2:5,profiles),k in eachindex(p.ranks)
            println(io,"$k,$c,$(p.ranks[k]),$(p.discarded_weights[k]),$(p.meets_cutoff[k])")
        end
    end
    fig=Figure(size=(1120,680),fontsize=18,
        fonts=(;regular="Times New Roman",bold="Times New Roman Bold"),figure_padding=24)
    Label(fig[0,1:2],case.title;fontsize=24,font=:bold,tellwidth=false)
    left=GridLayout(fig[1,1])
    right=GridLayout(fig[1,2])
    plot_schmidt!(left,data;style=:paper2026,bond_labels=case.labels,
        title="(a)  Figure 6 layout",color_limits=(1e-6,1.))
    plot_schmidt!(right,data;style=:paper2022,cutoffs=[1e-4],
        title="(b)  Supplementary Figure 1 layout",color_limits=(1e-16,1.))
    colsize!(fig.layout,1,Relative(0.66))
    colgap!(fig.layout,34)
    Label(fig[2,1:2],case.description;fontsize=15,tellwidth=false)
    Label(fig[3,1:2],L"\mathrm{Contours:}\ \sum_{j>\chi_k}\sigma_{j,k}^2/\sum_j\sigma_{j,k}^2\leq10^{-C}\quad\mathrm{(independent\ bond\ truncations)}";
        fontsize=17,tellwidth=false)
    Label(fig[4,1:2],"Same data in both panels; separate color ranges. Positive values below each color floor saturate; exact zeros are blank.";
        fontsize=13,color=:gray30,tellwidth=false)
    save(joinpath(out,case.name*".pdf"),fig;pt_per_unit=1)
    push!(report,"| $(case.title) | $error | $(isnothing(analytic_error) ? "—" : analytic_error) | $tail_error | $checks / $checks |")
    println(case.name,": max |Δσ| = ",error,"; tail error = ",tail_error,"; ",checks," rank checks passed")
end
append!(report,["", "The first example additionally agrees with the analytic product of the seven pair spectra.",
    "Values are normalized singular values (not squared probabilities). Tests also cover complex amplitudes, QN blocks, gauge changes, input preservation, zero states, rank caps, and PDF export.",
    "These are synthetic examples, not reconstructions of the papers' DNS data."])
write(joinpath(out,"VALIDATION.md"),join(report,"\n")*"\n")
