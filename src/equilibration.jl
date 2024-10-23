# represents a set of MCMC chains who share time steps and parameters
@kwdef mutable struct MCMCEqChains{T<:AbstractSequence}
    init_sequences :: Vector{T}
    chains :: Vector{Vector{T}}
    tvals :: Vector{Int}
    info :: Vector{Vector{Any}}
    params :: Dict
    alphabet :: Alphabet
    avruntime :: Float64
end

"""
    equilibration_dynamics(
        g::PottsGraph, M, sampling_parameters::SamplingParameters;
        n_chains = 20,
        Teq_long = 1000*size(g).L,
        distance = hamming,
        codons = true, init_sequences = nothing,
)

Start `n_chains` MCMC chains from already equilibrated sequences.
Return the sequences and time steps for the chains.

## Kwargs
- `n_chains`: number of mcmc chains
- `init_sequences`: vector of sequences or `nothing`. Length should be `n_chains`.
- `Teq_long`: used if `init_sequences` not provided. Equilibration time between sequences at the start of chains.

The number of configurations in a chain is set by parameter `M`, and the number of steps
between two configurations is set by `sampling_parameters.Teq`.

The equilibration time for the initial start of chains is set by `Teq_long`.

The `distance` kwarg will be called like so `distance(x, y; normalize=true)` where x and y are two integer vectors
"""
function equilibration_chains(
    g::PottsGraph, M, sampling_parameters::SamplingParameters;
    init_sequences = :random_codon, n_chains = 20, Teq_long = 1000*size(g).L, kwargs...
)
    sampling_parameters.burnin != 0 && @warn "Burnin should probably be set to 0"

    init_sampling_parameters = @set sampling_parameters.Teq = Teq_long
    init_sampling_parameters = @set init_sampling_parameters.burnin = 5*Teq_long
    init_sequences = map(1:n_chains) do _
        mcmc_sample(
            g, 1, init_sampling_parameters;
            init=init_sequences, alignment_output=false, kwargs...
        ).sequences[end]
    end

    return equilibration_chains(g, M, sampling_parameters, init_sequences)
end

function equilibration_chains(
    g::PottsGraph, M, sampling_parameters, init_sequences;
    kwargs...
)
    sampling_parameters.burnin != 0 && @warn "Burnin should probably be set to 0"
    runtime = @elapsed chains = map(init_sequences) do s0
        mcmc_sample(g, M, sampling_parameters; init=s0, alignment_output=false, kwargs...)
    end
    tvals = chains[1].tvals
    params = chains[1].params
    info = map(x -> x.info, chains)
    chains = map(x -> x.sequences, chains)
    avruntime = runtime / length(init_sequences)
    return MCMCEqChains(;
        init_sequences, chains, tvals, info, params, g.alphabet, avruntime,
    )
end

