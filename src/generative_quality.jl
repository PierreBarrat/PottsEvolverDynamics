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
    hvals, w1, w2 = _shared_histogram(Hx, Hy, 2/L)

    p = plot()
    plot!(p, hvals, w1, label=label_X)
    plot!(p, hvals, w2, label=label_Y)

    plot!(xlabel = "Hamming distance", title = "Pairwise hamming distance", yscale = :log10)
    plot!(legend = :bottomleft)
    plot!(; kwargs...)

    return p
end

function energy_histogram(
    X::Alignment, Y::Alignment, potts::PottsGraph;
    label_X = "", label_Y = "", kwargs...
)
    potts = deepcopy(potts)
    PottsEvolver.set_gauge!(potts, :zero_sum)
    Ex = map(s -> energy(s, potts), X)
    Ey = map(s -> energy(s, potts), Y)

    Evals, w1, w2 = _shared_histogram(Ex, Ey)
    p = plot()
    plot!(p, Evals, w1, label=label_X)
    plot!(p, Evals, w2, label=label_Y)

    plot!(xlabel = "Energy", title = "Energy distribution", yscale = :log10)
    plot!(legend = :topright)
    plot!(; kwargs...)

    return p
end

function frobnorm_histogram(potts::PottsGraph; kwargs...)
    frob(J) = sum(x -> x^2, J)
    (; L, q) = size(potts)
    F = zeros(Float64, Int(L*(L-1)/2))
    k = 1
    for i in 1:L, j in (i+1):L
        F[k] = frob(potts.J[:, :, j, i])
        k += 1
    end

    fvals, w = _histogram(F)
    p = plot(fvals, w, label="")
    plot!(p, xlabel="Frobenius norm", title="Distribution of |J|^2")
    plot!(p; kwargs...)
    return p
end

function pca(X::Alignment, Y::Alignment)
    q = length(X.alphabet)
    L = size(X,1)

    X_01 = reshape(BioSequenceMappings.onehot(X).data, q*L, size(X,2))
    Y_01 = reshape(BioSequenceMappings.onehot(Y).data, q*L, size(Y,2))

    pca_components = fit(PCA, X_01, maxoutdim=2)

    proj_X = predict(pca_components, X_01)
    proj_Y = predict(pca_components, Y_01)

    p = plot()

    # Y with markers
    x, y = proj_Y[1,:], proj_Y[2,:]
    scatter!(p, x, y, marker=(2, 0.3, stroke(0), :blue), label="Potts")

    # X with a contour plot
    x, y = proj_X[1,:], proj_X[2,:]
    k = kde((x, y))
    contour!(p, k, c = :Reds, linewidth = 3, levels=6, label="Natural")

    plot!(xlabel = "PC1", ylabel = "PC2")
    return p
end

function _histogram(H, binwidth = auto_bindwidth(H))
    low, high = extrema(H)
    edges = collect(range(low, high, step=binwidth))
    xvals = (edges[1:end-1] + edges[2:end])/2

    w = fit(Histogram, H, edges).weights
    w = w / sum(w) / binwidth

    return xvals, w
end

function _shared_histogram(H1, H2, binwidth=auto_bindwidth(H1, H2))
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

function auto_bindwidth(X::AbstractVector, Y::AbstractVector, n=50)
    low, high = (min(minimum(X), minimum(Y)), max(maximum(X), maximum(Y)))
    return (high - low) / n
end
auto_bindwidth(X, n=50) = (maximum(X) - minimum(X))/n
