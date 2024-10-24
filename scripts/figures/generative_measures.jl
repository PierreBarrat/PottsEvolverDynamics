using DrWatson
@quickactivate "PottsEvolverDynamics"

using BioSequenceMappings
using Chain
using DelimitedFiles
using Infiltrator
using Measures
using Plots
using PottsEvolver
using PottsEvolverDynamics


include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))
Plots.default(; pubfig(24)...)

families = keys(PED.ref_families)
# families = ["PF00105"]
# families = ["PF00014"]

for fam in families
    @info fam
    #==========================================#
    ############### Prior checks ###############
    #==========================================#

    # check files are all here
    aln_naturals = PED.ref_families[fam]["alignment"]
    weights_naturals = PED.ref_families[fam]["weights"]
    aln_naturals, weights_natural = if !isfile(aln_naturals) || !isfile(weights_naturals)
        @error "Could not find $aln_naturals or $weights_naturals"
        continue
    else
        read_fasta(aln_naturals), vec(readdlm(weights_naturals))
    end

    potts = PED.ref_families[fam]["potts"]["file"]
    potts = if !isfile(potts)
        @error "Could not find potts model $potts for $fam"
        continue
    else
        read_graph(potts)
    end

    if !isdir(datadir("equilibrium_samples", fam))
        @error"No sample for $fam - looked for $(datadir(\"equilibrium_samples\", fam))"
        continue
    end
    expr = Regex("$(fam)_eqsample_id=([0-9]+)\\.fasta")
    aln_sampled = @chain datadir("equilibrium_samples", fam) begin
        readdir(; join=true)
        filter(file -> occursin(expr, file), _)
    end
    aln_sampled, id = if length(aln_sampled) > 1
        @error """
            Found samples $(aln_sampled), but expected only one.
            Fix this script, or remove extra samples.
        """
        continue
    elseif length(aln_sampled) == 0
        @error"No sample for $fam - looked in $(datadir(\"equilibrium_samples\", fam))"
        continue
    else
        (
            read_fasta(first(aln_sampled)),
            parse(UInt, match(expr, first(aln_sampled)).captures[1])
        )
    end

    #====================================#
    ############# For saving #############
    #====================================#

    savedir = plotsdir("generative", fam)
    mkpath(savedir)

    #=====================#
    ######## Plots ########
    #=====================#

    shared_title = "$fam - $(PED.ref_families[fam]["long_name"])"

    # f1 and C
    plt_f1, plt_corr = PED.fitting_quality(
        aln_naturals, aln_sampled; weights_X = weights_natural
    )
    plot!(plt_f1, xlabel = "Naturals", ylabel = "Potts")
    plot!(plt_f1, title = "$(shared_title)\n Single site frequencies")

    plot!(plt_corr, xlabel = "Naturals", ylabel = "Potts")
    plot!(plt_corr, title = "$(shared_title)\n Connected correlations")

    savefig(plt_f1, joinpath(savedir, "$(fam)_fit_f1_id$(id).png"))
    savefig(plt_corr, joinpath(savedir, "$(fam)_fit_corr_id$(id).png"))

    # pairwise hamming
    plt_hamming = PED.pairwise_hamming_histogram(
        aln_naturals, aln_sampled; label_X = "Natural", label_Y = "Potts",
    )
    plot!(plt_hamming, title = "$(shared_title)\n Pairwise hamming distances")
    savefig(plt_hamming, joinpath(savedir, "$(fam)_pw_hamming_id$(id).png"))

    # Energy distribution
    plt_energy = PED.energy_histogram(
        aln_naturals, aln_sampled, potts; label_X = "Natural", label_Y = "Potts",
    )
    plot!(plt_energy, title = "$(shared_title)\n Energy distribution")
    savefig(plt_energy, joinpath(savedir, "$(fam)_energy_id$(id).png"))

    # Frobenius norm distribution
    plt_frob = PED.frobnorm_histogram(potts)
    plot!(plt_frob, title = "$(shared_title)\n Frobenius norm distribution")
    savefig(plt_frob, joinpath(savedir, "$(fam)_frobenius_id$(id).png"))

    # PCA
    plt_pca = PED.pca(aln_naturals, aln_sampled)
    plot!(plt_pca, title = "$(shared_title)\n PCA", margin=10mm)
    savefig(plt_pca, joinpath(savedir, "$(fam)_pca_id$(id).png"))

    ## Extra plots if useful
    # #= pairwise hamming excluding gaps =#
    # alphabet = aln_naturals.alphabet
    # gap_integer = findfirst(i -> alphabet(i) == '-', 1:length(alphabet))
    # plt_hamming_nogaps = PED.pairwise_hamming_histogram(
    #     aln_naturals, aln_sampled;
    #     label_X = "Natural", label_Y = "Potts", exclude_state = gap_integer,
    # )
    # plot!(plt_hamming_nogaps, title = "$(shared_title)\n Pairwise hamming distances - ignore gaps")
    # savefig(plt_hamming_nogaps, joinpath(savedir, "$(fam)_pw_hamming_nogaps_id$(id).png"))
end
