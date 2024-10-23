"""
    infer_profile_model
"""
function infer_profile_model(family::AbstractString; kwargs...)
    alignment = read_fasta(PED.ref_families[family]["alignment"])
    weights = readdlm(PED.ref_families[family]["weights"]) |> vec
    return infer_profile_model(alignment, weights; kwargs...)
end
function infer_profile_model(
    A::AbstractAlignment, weights::AbstractVector{Float64}; pc = 1e-3,
)
    f1 = site_specific_frequencies(A, weights/sum(weights))
    q, L = size(f1)
    # add pc
    f1 .= (1-pc)*f1 .+ pc/q

    # construct h and J
    h = log.(f1)
    J = zeros(eltype(h), q, q, L, L)
    return PottsGraph(; h, J, alphabet = Alphabet(A))
end
