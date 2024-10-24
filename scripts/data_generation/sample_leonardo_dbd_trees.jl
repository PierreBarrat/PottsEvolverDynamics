using DrWatson

@quickactivate "PottsEvolverDynamics"

using Chain
using PottsEvolver
using PottsEvolverDynamics
using TreeTools

# branch lengths of both these trees are O(1)
working_dir = datadir("leonardo_dbd_tree")
trees = Dict(
    "star" => read_tree(joinpath(working_dir, "star_tree.nwk")),
    "dbd" => read_tree(joinpath(working_dir, "DBD_tree_noX.nwk")),
)

# Family specific parameters
fam = "PF00105"
potts = read_graph(PED.ref_families[fam]["potts"]["file"])
Teq_ref = PED.ref_families[fam]["potts"]["equilibration_codons"]["accepted"].Teq_steps
fgap = PED.ref_families[fam]["potts"]["fgap"]
sampling_params = SamplingParameters(
    Teq=0,
    burnin=5*Teq_ref, # to equilibrate root sequence
    step_meaning=:accepted,
    fraction_gap_step=fgap,
)

# Simulation parameters
μ = @chain logrange(Teq_ref/100, Teq_ref, length=5) map(x -> round(x, sigdigits=3), _)
tree = ["dbd", "star"]
allparams = @strdict μ tree
params_list = dict_list(allparams)

function scale_branches!(tree, μ)
    for node in nodes(tree; skiproot=true)
        t = branch_length(node)
        branch_length!(node, round(Int, μ*t))
    end
    return tree
end

function sample_tree_and_write(params; kwargs...)
    @info params
    @unpack μ, tree = params
    # scale tree branches
    tree = copy(trees[tree])
    scale_branches!(tree, μ)

    # sample
    @info "Sampling..."
    time = @timed begin
        (; leaf_sequences) = mcmc_sample(
            potts, tree, sampling_params; init=:random_codon, kwargs...
        )
    end
    @info "done in $(time.time) seconds"
    # storing runtime
    params["runtime"] = round(time.time; sigdigits=3)
    alloc_fields = [:malloc, :realloc, :poolalloc, :bigalloc]
    params["nallocated"] = sum(p -> getproperty(time.gcstats, p), alloc_fields)
    params["allocated"] = round(time.bytes; sigdigits=2)

    outfile = joinpath(
        working_dir,
        savename("leaves", params, "fasta")
    )
    @info "Writing output alignment to $outfile"
    write(outfile, leaf_sequences)
end

for params in params_list
    sample_tree_and_write(params)
end
