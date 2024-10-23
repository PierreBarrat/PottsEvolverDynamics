using DrWatson
@quickactivate "PottsEvolverDynamics"

using Accessors
using Chain
using Dates
using PottsEvolver
using PottsEvolverDynamics
using ProgressMeter
using YAML

# Global parameters
global_params = Dict(
    "model" => "potts",
    "step_type" => :gibbs,
    "stepmeaning" => "accepted",
    "init" => :random_codon,
    "nsequences" => 10_000,
    "burnin_factor" => 5,
)
@tag!(global_params)
id = hash(global_params) # id will depend on above parameters and git commit / script line
global_params["time"] = now()

params_template = SamplingParameters(;
    step_type = global_params["step_type"],
    step_meaning = global_params["stepmeaning"],
    Teq = 0,
)

families = collect(keys(PED.ref_families))

for fam in families
    @info fam
    savedir = datadir("equilibrium_samples", fam)
    params_file = joinpath(savedir, "parameters_id=$(id).yaml")
    sample_file = joinpath(savedir, "$(fam)_eqsample_id=$(id).fasta")
    sample_file_jld2 = joinpath(savedir, "$(fam)_eqsample_id=$(id).jld2")
    if isfile(params_file) && isfile(sample_file)
        @warn """
            Simulation already performed for these parameters.
            See $(params_file) and $(sample_file)
        """
        continue
    end

    # familiy specific settings
    fam_params = deepcopy(global_params)

    # extracting params
    L = PED.ref_families[fam]["L"]
    @unpack model, init, stepmeaning, burnin_factor, nsequences  = fam_params

    # Sampling parameters
    Teq = PED.ref_families[fam][model]["equilibration_codons"][string(stepmeaning)].Teq_steps
    burnin = burnin_factor * Teq
    fgap = @chain PED.ref_families[fam] begin
        getindex("potts")
        getindex("fgap")
    end
    fam_params["fgap"] = fgap
    fam_params["Teq"] = Teq
    fam_params["burnin"] = burnin
    sampling_parameters = @set params_template.Teq = Teq
    @reset sampling_parameters.burnin = burnin
    @reset sampling_parameters.fraction_gap_step = fgap

    # Reading model
    model_file = PED.ref_families[fam][model]["file"]
    model = read_graph(model_file)

    @info "Sampling..."
    runtime = @elapsed begin
        sample_result = mcmc_sample(
            model, nsequences, sampling_parameters;
            init, verbose=0, alignment_output=true, translate_output=true
        )
    end
    # Saving
    wsave(params_file, YamlWrapper(global_params))

    aln = sample_result[1]
    write(sample_file, aln)

    fam_params["runtime"] = runtime
    fam_params["time"] = now()
    result = Dict("sample" => sample_result, "params" => fam_params)
    wsave(sample_file_jld2, result)

end


