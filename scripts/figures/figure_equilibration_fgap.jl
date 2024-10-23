### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 17ca728e-871d-11ef-0de1-5d7831969c4d
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "PottsEvolverDynamics")
end

# ╔═╡ 108a02bf-35a3-4475-9721-a8cf2f17211c
begin
	using DataFrames
	using DataFramesMeta
	using Measures
	using PottsEvolverDynamics
	using PlutoUI
	using StatsPlots
	using StatsBase
end

# ╔═╡ 36ee37e6-2d5f-4bc8-a12b-4579ff22cab1
include(joinpath(homedir(), ".julia/config/plot_defaults.jl"))

# ╔═╡ 7ea22e29-f130-4f6b-bc66-14a521b8510b
Plots.default(; pubfig(24)...)

# ╔═╡ 4b0d145c-6c8a-4ae2-b409-7af8830408f3
function datafiles(fam)
    famdir = datadir("equilibration_fgap", fam)
    return if isdir(famdir)
        filter(readdir(famdir; join=true)) do f
            occursin(r"hamming_inter_intra_chains_id[0-9]+\.jld2", f)
        end
    else
        ()
    end
end



# ╔═╡ 399cc06d-07da-4f5e-9a29-abb36abd4e4e
function make_save_figures(fam)
    L = PED.ref_families[fam]["L"]
    if isempty(datafiles(fam))
        @error """
        No data found for $fam
        Maybe run 'data_analysis/equilibration_fgap_extract_hamming_distances_codons.jl'
        """
    end

    # Plot settings
    dpi = 200
    left_margin = 15mm
    bottom_margin = 15mm
    layout_size = (2100, 1400)
    xscale = :log10
    fmt(x) = round(x, sigdigits=2)

    for datafile in datafiles(fam)
        @info datafile
        # Reading and parsing datas
        data = let
            @unpack data_measures = wload(datafile)
            sort!(data_measures, :fgap)
        end
        stepmeaning, run_id = let
            @assert length(unique(data.stepmeaning)) == 1
            @assert length(unique(data.id)) == 1
            first(data.stepmeaning), first(data.id)
        end

        # Normal hamming
        plts_hamming = map(eachrow(data)) do r
            idx = (xscale == :linear) ? (1:length(r.tvals)) : (2:length(r.tvals))

            plot(
                r.tvals[idx], r.av_hamming_inter_aa[idx], ribbon=r.std_hamming_inter_aa[idx];
                label = "Inter-chain distance",
            )
            plot!(
                r.tvals[idx], r.av_hamming_init_aa[idx];
                label = "Distance to init",
            )
            plot!(;
                xscale,
                ylim = (0, L),
                title = "fgap = $(fmt(r.fgap)) - runtime = $(fmt(r.runtime))",
                frame = :box,
                xlabel = "MCMC steps",
                ylabel = "Hamming",
                legend = :bottomright,
            )
        end
        fig_all = plot(
            plts_hamming...;
            layout = grid(2,3), size = layout_size, bottom_margin, left_margin, dpi,
            plot_title = "$fam - stepmeaning: $stepmeaning - Hamming distance, all",
            plot_titlefontsize=30,
        )

        # Hamming: only gaps
        plts_hamming_onlygaps = map(eachrow(data)) do r
            idx = (xscale == :linear) ? (1:length(r.tvals)) : (2:length(r.tvals))

            plot(
                r.tvals[idx], r.av_hamming_inter_aa_onlygaps[idx];
                ribbon=r.std_hamming_inter_aa_onlygaps[idx], label = "Inter-chain distance",
            )
            plot!(
                r.tvals[idx], r.av_hamming_init_aa_onlygaps[idx];
                label = "Distance to init",
            )
            fmt(x) = round(x, sigdigits=2)
            plot!(;
                xscale,
                ylim = (0, 10),
                title = "fgap = $(fmt(r.fgap)) - runtime = $(fmt(r.runtime))",
                frame = :box,
                xlabel = "MCMC steps",
                ylabel = "Hamming - only gaps"
            )
        end
        fig_onlygaps = plot(
            plts_hamming_onlygaps...;
            layout = grid(2,3), size = layout_size, bottom_margin, left_margin, dpi,
            plot_title = "$fam - stepmeaning: $stepmeaning - Hamming distance, only gaps",
            plot_titlefontsize=30,
        )

        # Hamming: no gaps
        plts_hamming_nogaps = map(eachrow(data)) do r
            idx = (xscale == :linear) ? (1:length(r.tvals)) : (2:length(r.tvals))

            plot(
                r.tvals[idx], r.av_hamming_inter_aa_nogaps[idx];
                ribbon=r.std_hamming_inter_aa_nogaps[idx], label = "Inter-chain distance",
            )
            plot!(
                r.tvals[idx], r.av_hamming_init_aa_nogaps[idx];
                label = "Distance to init",
            )
            plot!(;
                xscale,
                ylim = (0, L),
                title = "fgap = $(fmt(r.fgap)) - runtime = $(fmt(r.runtime))",
                frame = :box,
                xlabel = "MCMC steps",
                ylabel = "Hamming - ignore gaps",
                legend = :bottomright,
            )
        end
        fig_ignoregaps = plot(
            plts_hamming_nogaps...;
            layout = grid(2,3), size = layout_size, bottom_margin, left_margin, dpi,
            plot_title = "$fam - stepmeaning: $stepmeaning - Hamming distance, ignore gaps",
            plot_titlefontsize=30,
            # tickfontsize = 16,
        )

        # Saving

        suffix = "$(fam)_fgaps_hamming_inter_vs_intra_$(stepmeaning)_id$(run_id)"
        figs = Dict(
            suffix * "_all.png" => fig_all,
            suffix * "_ignoregaps.png" => fig_ignoregaps,
            suffix * "_onlygaps.png" => fig_onlygaps,
        )
        for (name, fig) in figs
            mkpath(plotsdir("equilibration_fgap", fam))
            savefig(fig, plotsdir("equilibration_fgap", fam, name))
        end
    end # file loop
end

# ╔═╡ 5079b59b-ebdb-4cb1-8f74-09ef1185b6e1
for fam in keys(PED.ref_families)
    @info fam
    make_save_figures(fam)
end





# ╔═╡ 17dc09a4-3140-46ac-8523-3e19a6c28db0

# ╔═╡ 972ce37c-2d29-407b-9680-38bc7fbf951d

# ╔═╡ 2d7579cb-3ca8-4250-b150-5f2c25dc904b

# ╔═╡ 56bd74f5-4857-4b46-bd0e-f381dcd734c7

# ╔═╡ 396e8348-61ca-40dc-a721-12f9b4f7cb4e


# ╔═╡ 46d8f56c-1040-4902-a815-3794d29fd314


# ╔═╡ 08f1c0ed-6cdd-4db2-b0b8-f609d64d504f


# ╔═╡ 39abfc57-7aa5-43d4-b6aa-c3e7ae710172


# ╔═╡ 7663500e-75a5-4de4-bac9-0f7dba02f699


# ╔═╡ abda1b27-eff1-44b7-8a63-cf4d861c8b2c

# ╔═╡ 89973699-a206-4a42-8f6c-5d4c066169f1

# ╔═╡ 9170d9bb-aba5-436e-b03d-614d23c22f22


# ╔═╡ 31b22d53-0633-4004-951f-3e02b53ba229


# ╔═╡ c9522000-8a7e-4c55-8176-c9bf0c59c30c

# ╔═╡ Cell order:
# ╠═17ca728e-871d-11ef-0de1-5d7831969c4d
# ╠═108a02bf-35a3-4475-9721-a8cf2f17211c
# ╠═36ee37e6-2d5f-4bc8-a12b-4579ff22cab1
# ╠═7ea22e29-f130-4f6b-bc66-14a521b8510b
# ╠═399cc06d-07da-4f5e-9a29-abb36abd4e4e
# ╠═4b0d145c-6c8a-4ae2-b409-7af8830408f3
# ╠═5079b59b-ebdb-4cb1-8f74-09ef1185b6e1
# ╟─17dc09a4-3140-46ac-8523-3e19a6c28db0
# ╠═396e8348-61ca-40dc-a721-12f9b4f7cb4e
# ╠═46d8f56c-1040-4902-a815-3794d29fd314
# ╠═e87d9c25-6d0e-470a-9770-5b8c6465b2fd
# ╠═08f1c0ed-6cdd-4db2-b0b8-f609d64d504f
# ╠═39abfc57-7aa5-43d4-b6aa-c3e7ae710172
# ╠═7663500e-75a5-4de4-bac9-0f7dba02f699
# ╠═972ce37c-2d29-407b-9680-38bc7fbf951d
# ╠═2d7579cb-3ca8-4250-b150-5f2c25dc904b
# ╠═56bd74f5-4857-4b46-bd0e-f381dcd734c7
# ╠═abda1b27-eff1-44b7-8a63-cf4d861c8b2c
# ╟─89973699-a206-4a42-8f6c-5d4c066169f1
# ╟─9170d9bb-aba5-436e-b03d-614d23c22f22
# ╟─31b22d53-0633-4004-951f-3e02b53ba229
# ╠═c9522000-8a7e-4c55-8176-c9bf0c59c30c
