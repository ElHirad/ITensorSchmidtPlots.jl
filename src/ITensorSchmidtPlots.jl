module ITensorSchmidtPlots

using ITensors
using ITensorMPS
using LinearAlgebra: norm
using CairoMakie
using LaTeXStrings

export SchmidtSpectrum, TruncationProfile, schmidt_spectrum, truncation_profile,
       plot_schmidt, plot_schmidt!, save_schmidt

include("spectrum.jl")
include("plotting.jl")

end
