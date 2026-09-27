"""
    SchmidtSpectrum(values; normalized=false)

Schmidt singular values in descending order, one vector per internal bond.
The constructor validates and copies the vectors. `normalized=true` normalizes
each vector to unit Euclidean norm. Zero states and nonfinite values are rejected.
"""
struct SchmidtSpectrum
    values::Vector{Vector{Float64}}
    normalized::Bool
    function SchmidtSpectrum(values; normalized::Bool=false)
        isempty(values) && throw(ArgumentError("At least one internal bond is required."))
        result = Vector{Float64}[]
        for v in values
            a = Float64.(collect(v))
            isempty(a) && throw(ArgumentError("A spectrum cannot be empty."))
            all(x -> isfinite(x) && x >= 0, a) ||
                throw(ArgumentError("Singular values must be finite and nonnegative."))
            issorted(a; rev=true) || throw(ArgumentError("Sort singular values in descending order."))
            a[1] > 0 || throw(ArgumentError("The zero state has no normalized Schmidt spectrum."))
            if normalized
                # Scale first so very large/small state norms do not overflow/underflow.
                a ./= a[1]
                a ./= norm(a)
            end
            push!(result, a)
        end
        new(result, normalized)
    end
end

"""
    schmidt_spectrum(psi::MPS; normalize=true, svd_alg="qr_iteration")

Compute the true bipartition Schmidt spectrum of the state represented by `psi`
at every internal bond. Recanonicalize a private copy, then SVD the orthogonality
center with `cutoff=0`. No dense wavefunction is formed and the input is unchanged.
Works with real, complex, and QN-conserving MPSs. Output values are Float64.

Already-discarded singular values cannot be reconstructed. Supply the original
untruncated MPS to visualize a proposed truncation. This is not a record of the
SVDs used during an earlier, sequentially truncated MPS construction.
"""
function schmidt_spectrum(psi::MPS; normalize::Bool=true, svd_alg="qr_iteration")
    length(psi) >= 2 || throw(ArgumentError("An MPS needs at least two sites to have an internal bond."))
    # A new MPS resets the orthogonality metadata, including stale metadata from
    # in-place tensor edits. Deep copying guarantees no tensor storage is shared.
    work = MPS([deepcopy(psi[k]) for k in 1:length(psi)])
    for k in 1:(length(work)-1)
        length(commoninds(work[k], work[k+1])) == 1 ||
            throw(ArgumentError("Expected exactly one link between adjacent sites $k and $(k+1)."))
    end
    orthogonalize!(work, 1)
    nrm = norm(work[1])
    isfinite(nrm) && nrm > 0 || throw(ArgumentError("MPS norm must be finite and nonzero."))
    normalize && (work[1] /= nrm)
    values = Vector{Float64}[]
    for k in 1:(length(work)-1)
        orthogonalize!(work, k)
        left = uniqueinds(work[k], work[k+1])
        _, S, _ = svd(work[k], left; cutoff=0.0, alg=svd_alg)
        s = sort!([Float64(abs(S[j,j])) for j in 1:dim(inds(S)[1])]; rev=true)
        push!(values, s)
    end
    return SchmidtSpectrum(values; normalized=normalize)
end

"""Retained ranks, actual discarded squared weights, and cutoff status for each bond."""
struct TruncationProfile
    cutoff::Float64
    maxdim::Int
    ranks::Vector{Int}
    discarded_weights::Vector{Float64}
    meets_cutoff::BitVector
end

function _tails(s)
    # Reverse accumulation avoids cancellation in 1 - cumulative retained weight.
    weights = abs2.(s ./ first(s))
    weights ./= sum(weights)
    tail = zeros(length(s)+1)
    for j in length(s):-1:1
        tail[j] = tail[j+1] + weights[j]
    end
    return tail
end

"""
    truncation_profile(spectrum; cutoff=1e-4, maxdim=typemax(Int))

For each bond choose the smallest χ ≥ 1 with
`sum(abs2, s[χ+1:end]) / sum(abs2, s) <= cutoff`, capped by `maxdim`.
The cutoff is a relative discarded **squared** weight, not an amplitude threshold.
If the cap prevents meeting the cutoff, `meets_cutoff[k]` is false.
These are independent single-bond truncations of the supplied state: they do not
predict the ranks or total error of a sequential truncation of all MPS bonds.
"""
function truncation_profile(data::SchmidtSpectrum; cutoff::Real=1e-4, maxdim::Integer=typemax(Int))
    isfinite(cutoff) && 0 <= cutoff < 1 || throw(ArgumentError("cutoff must be finite and in [0, 1)."))
    1 <= maxdim <= typemax(Int) || throw(ArgumentError("maxdim must be a positive Int-sized integer."))
    ranks, weights, meets = Int[], Float64[], BitVector()
    for s in data.values
        tail = _tails(s)
        # For cutoff=0 preserve even positive values whose square underflows.
        rank = iszero(cutoff) ? findlast(>(0), s) : findfirst(r -> tail[r+1] <= cutoff, 1:length(s))
        rank = min(rank, maxdim)
        push!(ranks, rank)
        push!(weights, tail[rank+1])
        push!(meets, tail[rank+1] <= cutoff && (!iszero(cutoff) || all(iszero, s[rank+1:end])))
    end
    return TruncationProfile(Float64(cutoff), Int(maxdim), ranks, weights, meets)
end
