module ITensorSchmidtPlots

using ITensors
using ITensorMPS
using LinearAlgebra: norm
using CairoMakie
using LaTeXStrings
import Makie

export SchmidtSpectrum, TruncationProfile, schmidt_spectrum, truncation_profile,
       plot_schmidt, plot_schmidt!, save_schmidt,
       SchmidtPlot, schmidtplot, schmidtplot!,
       SchmidtTruncation, schmidttruncation, schmidttruncation!

include("spectrum.jl")
include("plotting.jl")
include("recipes.jl")

end
