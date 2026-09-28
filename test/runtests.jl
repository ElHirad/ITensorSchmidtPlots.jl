using Test, ITensorSchmidtPlots, ITensors, ITensorMPS, LinearAlgebra, Random, CairoMakie
include(joinpath(@__DIR__, "..", "examples", "states.jl"))
BLAS.set_num_threads(1)

@testset "Independent dense and analytic spectra" begin
    for case in (pair_example(), field_example(), random_example())
        psi = state_mps(case.v, case.n)
        old = deepcopy(psi)
        data = schmidt_spectrum(psi)
        exact = dense_spectra(case.v,case.n)
        @test length(data.values) == case.n-1
        for k in eachindex(exact)
            s = data.values[k]
            padded = vcat(s,zeros(length(exact[k])-length(s)))
            @test padded ≈ exact[k] atol=5e-13 rtol=5e-12
            @test sum(abs2,s) ≈ 1 atol=5e-14
            @test psi[k] == old[k]
            if !isnothing(case.analytic)
                @test padded ≈ case.analytic[k] atol=5e-13 rtol=5e-12
            end
        end
        @test psi[end] == old[end]
        @test (ITensorMPS.leftlim(psi), ITensorMPS.rightlim(psi)) ==
            (ITensorMPS.leftlim(old), ITensorMPS.rightlim(old))
        for cutoff in (1e-2,1e-3,1e-4,1e-5)
            profile = truncation_profile(data; cutoff)
            for (k,r) in enumerate(profile.ranks)
                tail = sum(abs2, exact[k][r+1:end])
                @test r == findfirst(j -> sum(abs2,exact[k][j+1:end]) <= cutoff,1:length(exact[k]))
                @test tail <= cutoff + 1e-13
                @test profile.discarded_weights[k] ≈ tail atol=1e-13
                @test r == 1 || sum(abs2,exact[k][r:end]) > cutoff - 1e-13
            end
        end
    end
end

include("native_makie.jl")

@testset "Gauge independence, raw norm, and QN blocks" begin
    sites = siteinds("Qubit",6)
    psi = random_mps(MersenneTwister(8), ComplexF64, sites; linkdims=5)
    original = schmidt_spectrum(psi)
    orthogonalize!(psi,3)
    # In-place storage modification deliberately leaves orthogonality metadata stale.
    psi[1] .*= 3.0
    @test all(isapprox.(schmidt_spectrum(psi).values,original.values))
    @test all(isapprox.(schmidt_spectrum(psi;normalize=false).values,3 .* original.values))
    # Nonunitary gauge change across a link preserves the full state.
    link = linkind(psi,3)
    G = ITensor(link,prime(link))
    Gi = ITensor(link,prime(link))
    for j in 1:dim(link)
        G[link=>j,prime(link)=>j] = 0.5 + j
        Gi[link=>j,prime(link)=>j] = inv(0.5+j)
    end
    psi[3] = noprime(psi[3]*G)
    psi[4] = noprime(Gi*psi[4])
    @test all(isapprox.(schmidt_spectrum(psi).values,original.values))
    qsites = siteinds("S=1/2",6; conserve_qns=true)
    qpsi = random_mps(MersenneTwister(12),qsites,[isodd(k) ? "Up" : "Dn" for k in 1:6];linkdims=4)
    qdata = schmidt_spectrum(qpsi)
    full = prod(qpsi)
    a = vec(Array(full,qsites...))
    exact = dense_spectra(a,6)
    for k in 1:5
        s=qdata.values[k]
        @test vcat(s,zeros(length(exact[k])-length(s))) ≈ exact[k] atol=1e-13
    end
end

@testset "Truncation boundary and invalid inputs" begin
    data=SchmidtSpectrum([[sqrt(0.75),0.5]])
    # Use an exactly representable squared-weight boundary for the equality test.
    equal=SchmidtSpectrum([[1.,1.]])
    @test truncation_profile(equal;cutoff=0.5).ranks == [1]
    @test truncation_profile(equal;cutoff=prevfloat(0.5)).ranks == [2]
    @test truncation_profile(data;cutoff=0.249).ranks == [2]
    cap=truncation_profile(data;cutoff=0.01,maxdim=1)
    @test cap.ranks == [1]
    @test cap.discarded_weights ≈ [0.25]
    @test cap.meets_cutoff == [false]
    @test truncation_profile(SchmidtSpectrum([[1.,1e-200,0.]]);cutoff=0).ranks == [2]
    @test truncation_profile(SchmidtSpectrum([[1.,0.]]);cutoff=0).ranks == [1]
    @test SchmidtSpectrum([[1e300,1e300]];normalized=true).values[1] ≈ [inv(sqrt(2)),inv(sqrt(2))]
    for values in ([],[Float64[]],[[0.]],[[NaN]],[[Inf]],[[-1.]],[[0.2,0.8]])
        @test_throws ArgumentError SchmidtSpectrum(values)
    end
    for c in (-1.,1.,NaN,Inf)
        @test_throws ArgumentError truncation_profile(data;cutoff=c)
    end
    @test_throws ArgumentError truncation_profile(data;maxdim=0)
    @test_throws ArgumentError schmidt_spectrum(MPS(siteinds("Qubit",1),"0"))
    zero=0.0*MPS(siteinds("Qubit",2),"0")
    @test_throws ArgumentError schmidt_spectrum(zero)
end

@testset "PDF export and plotting edge cases" begin
    data=SchmidtSpectrum([[1.],[0.9,sqrt(0.19)],[1.]];normalized=true)
    mktempdir() do dir
        path=joinpath(dir,"paper2026.pdf")
        @test save_schmidt(path,data;style=:paper2026) == abspath(path)
        @test read(path,String)[1:4] == "%PDF"
        @test filesize(path) > 1000
        fig=Figure()
        plotted=plot_schmidt!(fig[1,1],data;style=:paper2026)
        @test !plotted.axis.yreversed[]
        @test_throws ArgumentError save_schmidt(joinpath(dir,"removed.pdf"),data;style=:paper2022)
        @test !isfile(joinpath(dir,"removed.pdf"))
        @test_throws ArgumentError save_schmidt(joinpath(dir,"bad.svg"),data)
    end
    @test plot_schmidt(SchmidtSpectrum([[1.]]);cutoffs=Float64[]) !== nothing
    @test plot_schmidt(data;maxdim=1) !== nothing
    @test_throws ArgumentError plot_schmidt(data;style=:wrong)
    @test_throws ArgumentError plot_schmidt(data;bond_labels=["one"])
    @test_throws ArgumentError plot_schmidt(data;color_limits=(0.,1.))
end
