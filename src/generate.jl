function generate_phyloformer_training(;
    # model
    family = nothing,
    model = :potts,
    model_file = nothing,
    # tree
    nleaves = 20,
    tree_height = 1, # measured in steps/L
    coalescent = :yule,
    # dynamics
    step_meaning = :accepted,
    # init sequence
    init = :random_codon,
)



end


function generate_phyloformer_training(params::Dict{Symbol, Any})
    return generate_phyloformer_training(; params...)
end


"""
"""
function generate_tree(coalescent, nleaves, tree_height)
    coa = if coalescent == :yule
        YuleCoalescent(n=nleaves, b=1)
    elseif coalescent == :kingman
        KingmanCoalescent(n=nleaves, N=1)
    else
        error("Unrecognized coalescent $coalescent. Choose among `[:yule, :kingman]`")
    end
    tree = genealogy(coa)
    # rescaling branches
    H = distance(first(leaves(tree)), root(tree)) # current height
    for node in nodes(tree; skiproot=true)
        τ = branch_length(node)
        branch_length!(node, round(Int, τ / H * tree_height))
    end

    return tree
end

function load_model(family, model, model_file)
    return if !isnothing(model_file)
        read_graph(model_file)
    elseif !isnothing(family)
        model_file = PED.ref_families[string(family)][string(model)]["file"]
        read_graph(model_file)
    else
        error("Need input `family` or `model_file` to load a model. Got `nothing` and `nothing`.")
    end
end
