using DrWatson
@quickactivate "PottsEvolverDynamics"

using Accessors
using Chain
using DataFrames # for collect_results
using Dates
using PottsEvolver
using PottsEvolverDynamics
using ProgressMeter
using YAML

# Global parameters
global_params = Dict(
    "model" => ["potts", "profile"],
    "step_type" => :gibbs,
    "stepmeaning" => :accepted,
    "init" => [:random_codon, :random_aa],
    "nchains" => 50,
    "Tshort_multiplier" => 5, # Teq_short = this * L ~ it's in sweep
    "Tlong_multiplier" => 10000, # in sweeps - for initial sequences
    "chain_sweeps" => 10000, # MEASURED IN SWEEPS
)
@tag!(global_params)
id = hash(global_params) # id will depend on above parameters and git commit / script line
global_params["time"] = now()

params_template = SamplingParameters(;
    step_type = global_params["step_type"],
    step_meaning = global_params["stepmeaning"],
    Teq = 1,
)

function get_init_sequences(potts, nchains, Teq_long)
    init_sampling_parameters = @set params_template.burnin = Teq_long
    @reset init_sampling_parameters.fraction_gap_step = .9 # should be decent
    init_sequences = map(1:nchains) do _
        mcmc_sample(
            potts, nchains, init_sampling_parameters;
            init = :random_codon, alignment_output = false,
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

slow_families = ["PF01817"]
# families = collect(keys(PED.ref_families))
families = slow_families
profile_ratio = 20

for fam in families
    @info fam
    #=====================================================================#
    ######################## Simulation parameters ########################
    #=====================================================================#
    # choosing fgap
    fgap = @chain PED.ref_families[fam] begin
        getindex("potts")
        getindex("fgap")
    end

    # setting times
    L = PED.ref_families[fam]["L"]
    Teq_short = round(Int, L * global_params["Tshort_multiplier"])
    Teq_long = round(Int, L * global_params["Tlong_multiplier"])
    chainsteps = round(Int, global_params["chain_sweeps"] * L / Teq_short) + 1

    # extra parameters / identifiers
    @unpack model, init, stepmeaning, nchains  = global_params
    sim_params = @strdict nchains model init fgap Teq_long Teq_short chainsteps stepmeaning id

    #=============================================================================================#
    ################################ Already performed for params? ################################
    #=============================================================================================#

    # if simulation was already there, no need to run it again
    savedir = datadir("equilibration", fam, "chains")
    params_file = joinpath(savedir, "parameters_id=$(id).yaml")
    data_file(sparams) = joinpath(
        savedir,
        savename(sparams, "jld2"; ignores=["runtime", "time", "fgap", "Teq_long"])
    )
    if isfile(params_file) && all(p -> isfile(data_file(p)), dict_list(sim_params))
        @warn """
            Simulation already performed for these parameters.
            See $(params_file) and $data_file
            Updating results dataframe and skipping simulation.
        """
        _custom_collect_results(savedir)
        continue
    else
        safesave(params_file, YamlWrapper(global_params))
    end

    #================================================================================================#
    ################################# Initial sequences and sampling #################################
    #================================================================================================#
    # init sequences - pick with potts model
    # init sequences are always CodonSequence at this point -- translate later if needed
    @info "Sampling initial sequences"
    potts_file = PED.ref_families[fam]["potts"]["file"]
    potts = read_graph(potts_file)
    init_sequences = get_init_sequences(potts, nchains, Teq_long)

    # Simulate chains
    @info "Starting chains"
    @showprogress dt=1 showspeed=true for sparams in dict_list(sim_params)
        # picking model
        model_file = PED.ref_families[fam][sparams["model"]]["file"]
        model = read_graph(model_file)

        # Setting initial sequences if aa
        local_init_sequences = if sparams["init"] == :random_aa
            map(PottsEvolver.translate, init_sequences)
        elseif sparams["init"] == :random_codon
            init_sequences
        else
            error()
        end

        # Simulating
        sparams = convert(Dict{String, Any}, sparams)
        chains = if sparams["model"] == "profile"
            sampling_parameters = @set params_template.Teq = Teq_short/profile_ratio
            @reset sampling_parameters.burnin = 0
            @reset sampling_parameters.fraction_gap_step = fgap
            chainsteps_profile = round(Int, chainsteps/profile_ratio)
            PED.equilibration_chains(
                model, chainsteps_profile, sampling_parameters, local_init_sequences,
            )
        else
            sampling_parameters = @set params_template.Teq = Teq_short
            @reset sampling_parameters.burnin = 0
            @reset sampling_parameters.fraction_gap_step = fgap
            chainsteps_profile = round(Int, chainsteps)
            PED.equilibration_chains(
                model, chainsteps, sampling_parameters, local_init_sequences,
            )
        end

        # Saving
        sparams["chains"] = chains
        sparams["runtime"] = chains.avruntime
        sparams["time"] = now()
        safesave(data_file(sparams), sparams)
    end
    # Collect chains into dataframe and save it
    _custom_collect_results(savedir)
end

