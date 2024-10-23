using DrWatson
@quickactivate "PottsEvolverDynamics"

using Accessors
using DataFrames # for collect_results
using Dates
using PottsEvolver
using PottsEvolverDynamics
using ProgressMeter
using YAML

# Global parameters
global_params = Dict(
    "step_type" => :gibbs,
    "step_meaning" => :accepted,
    "random_init" => :random_codon,
    "nchains" => 100,
    "Tshort_multiplier" => 0.5, # Teq_short = this * L ~ it's in sweep
    "Tlong_multiplier" => 1000, # in sweeps
    "chain_sweeps" => 200, # MEASURED IN SWEEPS
    "fgap" => [.1, .25, .5, .75, .9, .95],
)
@tag!(global_params)
id = hash(global_params) # id will depend on above parameters and git commit / script line
global_params["time"] = now()

params_template = SamplingParameters(;
    step_type = global_params["step_type"],
    step_meaning = global_params["step_meaning"],
    Teq = 1,
)

function get_init_sequences(potts, nchains, Teq_long)
    init_sampling_parameters = @set params_template.burnin = Teq_long
    @reset init_sampling_parameters.fraction_gap_step = .9 # should be decent
    init_sequences = map(1:nchains) do _
        mcmc_sample(
            potts, nchains, init_sampling_parameters;
            init = global_params["random_init"], alignment_output = false,
        ).sequences[1]
    end
    return init_sequences
end

function _custom_collect_results(savedir)
    collect_results!(
        joinpath(dirname(savedir), "results_$(basename(savedir))_id$(id).jld2"), savedir;
        verbose=true, update=true, rinclude=[Regex("id=$id")],
    )
end

# per simulation parameters

# families = ["PF00076", "PF00072"]
families = collect(keys(PED.ref_families))

for fam in families
    @info fam
    savedir = datadir("equilibration_fgap", fam, "chains")
    # if simulation was already there, no need to run it again
    params_file = joinpath(savedir, "parameters_id=$(id).yaml")
    if isfile(params_file)
        @warn """
            Simulation already performed for these parameters.
            See $(params_file)
            Updating results dataframe (in case) and skipping simulation.
        """
        _custom_collect_results(savedir)
        continue
    else
        safesave(params_file, YamlWrapper(global_params))
    end

    potts_file = PED.ref_families[fam]["potts"]["file"]
    potts = read_graph(potts_file)

    # setting times
    L = size(potts).L
    Teq_short = round(Int, L * global_params["Tshort_multiplier"])
    Teq_long = round(Int, L * global_params["Tlong_multiplier"])
    chainsteps = round(Int, global_params["chain_sweeps"] * L / Teq_short) + 1

    # extra parameters / identifiers
    @unpack fgap, nchains = global_params

    # init sequences
    @info "Sampling initial sequences"
    init_sequences = get_init_sequences(potts, nchains, Teq_long)
    sim_params = @strdict nchains fgap Teq_long Teq_short chainsteps id
    @info "Starting chains"
    @showprogress dt=1 showspeed=true for sparams in dict_list(sim_params)
        # Simulating
        sparams = convert(Dict{String, Any}, sparams)
        sampling_parameters = @set params_template.Teq = Teq_short
        @reset sampling_parameters.burnin = 0
        @reset sampling_parameters.fraction_gap_step = sparams["fgap"]
        chains = PED.equilibration_chains(potts, chainsteps, sampling_parameters, init_sequences);

        # Saving structure
        sparams["chains"] = chains
        sparams["runtime"] = chains.avruntime
        sparams["stepmeaning"] = global_params["step_meaning"]
        sparams["time"] = now()
        # saving
        safesave(
            joinpath(savedir, savename(sparams, "jld2"; ignores=["runtime", "time"])),
            sparams
        )
    end
    _custom_collect_results(savedir)
end

