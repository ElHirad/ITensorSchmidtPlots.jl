# Small, reproducible full states used only for examples and independent validation.
# Production schmidt_spectrum does not form a full wavefunction.
using ITensors, ITensorMPS, LinearAlgebra, Random, LaTeXStrings

function state_mps(v, n)
    sites = siteinds("Qubit", n)
    tensor = ITensor(reshape(v, ntuple(_ -> 2, n)), sites...)
    return ITensors.@disable_warn_order MPS(tensor, sites; cutoff=0.0)
end

function pair_example()
    probabilities = [0.0003, 0.002, 0.02, 0.07, 0.15, 0.25, 0.45]
    m = length(probabilities)
    n = 2m
    v = zeros(2^n)
    for bits in 0:(2^m-1)
        index, amplitude = 0, 1.0
        for j in 1:m
            bit = (bits >> (j-1)) & 1
            index |= (bit << (j-1)) | (bit << (n-j))
            amplitude *= sqrt(bit == 0 ? 1-probabilities[j] : probabilities[j])
        end
        v[index+1] = amplitude
    end
    analytic = [sort!([prod(sqrt(((bits >> (j-1)) & 1) == 0 ? 1-probabilities[j] : probabilities[j])
                    for j in 1:min(k,n-k)) for bits in 0:(2^min(k,n-k)-1)]; rev=true) for k in 1:n-1]
    return (; name="01_entangled_pairs", title="Nonuniform entangled pairs", n, v,
        labels=nothing, analytic,
        description="14 qubits; seven mirrored pairs with known analytic Schmidt coefficients.")
end

function field_example()
    levels = 7
    n = 2levels
    v = zeros(2^n)
    rng = MersenneTwister(210605782)
    # A synthetic multiscale field, not DNS or data from the reference paper.
    modes = [(rand(rng,1:24),rand(rng,-24:24),rand(rng)*2pi,randn(rng)) for _ in 1:70]
    for bits in 0:(2^n-1)
        ix = sum(((bits >> (2j-2)) & 1) << (levels-j) for j in 1:levels)
        iy = sum(((bits >> (2j-1)) & 1) << (levels-j) for j in 1:levels)
        x, y = ix/2^levels, iy/2^levels
        coarse = sin(2pi*x)*cos(2pi*y) + 0.35cos(4pi*x+0.5)*sin(6pi*y)
        vortex = 0.8exp(-((x-0.38)^2+(y-0.57)^2)/0.012)*sin(30pi*(x+y))
        fine = sum(a*sin(2pi*(kx*x+ky*y)+phase)/(kx^2+ky^2)^0.8 for (kx,ky,phase,a) in modes)
        v[bits+1] = coarse + vortex + fine
    end
    v ./= norm(v)
    labels = [LaTeXString(isodd(k) ? "x_{$(div(k+1,2))}" : "y_{$(div(k,2))}") for k in 1:n-1]
    return (; name="02_multiscale_field", title="Synthetic multiscale field", n, v, labels,
        analytic=nothing, description="128 × 128 field; coarse-to-fine interleaved x/y bits; synthetic data.")
end

function random_example()
    n = 14
    v = randn(MersenneTwister(260828869), ComplexF64, 2^n)
    v ./= norm(v)
    return (; name="03_complex_random", title="Complex random state", n, v, labels=nothing,
        analytic=nothing, description="14 qubits; seeded complex Gaussian amplitudes; broad entanglement spectrum.")
end

function dense_spectra(v,n)
    return [svdvals(reshape(v, 2^k, :)) ./ norm(v) for k in 1:n-1]
end
