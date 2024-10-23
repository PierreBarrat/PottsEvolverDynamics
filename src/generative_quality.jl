function _fitting_quality(nat::Alignment, sample::Alignment; kwargs...)

    f1_nat = site_specific_frequencies(nat; as_vec=true)
    C_nat = pairwise_correlations(nat; as_mat=true)

    f1_sample = site_specific_frequencies(sample; as_vec=true)
    C_sample = pairwise_correlations(sample; as_mat=true)

    pC = plot()
    plot!([-0.2,0.2], [-0.2,0.2], line = (:black, :dash), alpha = 0.5, label="")
    scatter!(
        vec(C_nat)[1:10:end], vec(C_sample)[1:10:end];
        marker=(4, .5, stroke(0)), label="", color=1,
    )
    lims = (
        min(minimum(vec(C_nat)[1:10:end]), minimum(vec(C_sample)[1:10:end]) * 1.25),
        max(maximum(vec(C_nat)[1:10:end]), maximum(vec(C_sample)[1:10:end]) * 1.25),
    )
    plot!(xlabel = "Alignment 1", ylabel = "Alignment 2", title = "Connected correlations")
    plot!(; lims)
    plot!(; kwargs...)

    pf1 = plot()
    plot!([0,1], [0,1], line = (:black, :dash), label="", alpha = 0.5)
    scatter!(
        vec(f1_nat), vec(f1_sample);
        marker=(4, stroke(0.5)), label="", color=1
    )
    plot!(xlabel = "Alignment 1", ylabel = "Alignment 2", title = "Single site frequencies")
    plot!(; kwargs...)

    return pf1, pC
end

"""
    fitting_quality(
        X::Alignment, Y::Alignment;
        weights_X = nothing, weights_Y = nothing, compute_weights=false
    )
    fitting_quality(
        X::AbstractString, Y::AbstractString; kwargs...
    )
"""
function fitting_quality(
    X::Alignment, Y::Alignment;
    weights_X = nothing, weights_Y = nothing, compute_weights=false
)
    !isnothing(weights_X) && (X.weights = _normalize(read_weights(weights_X)))
    !isnothing(weights_Y) && (Y.weights = _normalize(read_weights(weights_Y)))

    if compute_weights
        !isnothing(weights_X) && compute_weights!(X)
        !isnothing(weights_Y) && compute_weights!(Y)
    end

    return _fitting_quality(X, Y)
end
function fitting_quality(X::AbstractString, Y::AbstractString; kwargs...)
    return fitting_quality(read_fasta(X), read_fasta(Y); kwargs...)
end
_normalize(X) = X / sum(X)


read_weights(w::AbstractVector) = w
read_weights(w::AbstractString) = vec(readdlm(w))


function pairwise_hamming_histogram(
    X::Alignment, Y::Alignment;
    max_seqs = 2_500, # number of hamming distances: max_seqs^2
    step_X = Int(ceil(size(X, 2) / max_seqs)),
    step_Y = Int(ceil(size(Y, 2) / max_seqs)),
    label_X = "", label_Y = "",
    kwargs...
)
    Hx = pairwise_hamming(X; step=step_X, kwargs...)
    Hy = pairwise_hamming(Y; step=step_Y, kwargs...)

    L = size(X, 1)
    p = plot()

    hvals, w1, w2 = _my_histogram(Hx, Hy, 2/L)
    plot!(p, hvals, w1, label=label_X)
    plot!(p, hvals, w2, label=label_Y)

    plot!(xlabel = "Hamming distance", title = "Pairwise hamming distance", yscale = :log10)
    plot!(legend = :bottomleft)
    plot!(; kwargs...)

    return p
end

function _my_histogram(H1, H2, binwidth)
    hmin = min(minimum(H1), minimum(H2))
    hmax = max(maximum(H1), maximum(H2))
    hedges = collect(range(hmin, hmax, step=binwidth))
    hvals = (hedges[1:end-1] + hedges[2:end])/2

    w1 = fit(Histogram, H1, hedges).weights
    w1 = w1 / sum(w1) / binwidth
    w1 = convert(Vector{Union{Float64, Missing}}, w1)

    w2 = fit(Histogram, H2, hedges).weights
    w2 = w2 / sum(w2) / binwidth
    w2 = convert(Vector{Union{Float64, Missing}}, w2)

    for i in 1:length(hvals)
        w1[i] = (w1[i] == 0 ? missing : w1[i])
        w2[i] = (w2[i] == 0 ? missing : w2[i])
    end

    return hvals, w1, w2
end
