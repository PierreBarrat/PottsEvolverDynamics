### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ 0c5fe0c0-86ee-11ef-05c2-f7cbf0e82146
begin
	using Revise
	using DrWatson
	quickactivate(@__DIR__, "PottsEvolverDynamics")
end

# ╔═╡ aa55d93a-1a7b-43aa-b95e-083c87eaacb5
begin
	using BioSequenceMappings
	using DataFrames
	using DataFramesMeta
	using Dates
	using PottsEvolver
	using PottsEvolverDynamics
	using Random
	using StatsBase
end

# ╔═╡ 507356ee-8e8a-4021-bed5-71e86bbf0f56
begin
    function input_data(fam)
       # return a vector of ("result_chain_idxxx.jld2", xxx) for every matching folder in datadir("equilibration_fgap", fam)
        return if isdir(datadir("equilibration_fgap", fam))
            @chain datadir("equilibration_fgap", fam) begin
                readdir
                map(s -> match(r"results_chains_id([0-9]+)\.jld2", s), _)
                filter(!isnothing, _)
                map(_) do m
                    (datadir("equilibration_fgap", fam, m.match), parse(UInt, m.captures[1]))
                end
            end
        else
            []
        end
    end
    # datadir("equilibration_fgap", fam, "results_chains.jld2")
	output_data(fam, id) = datadir(
	    "equilibration_fgap", fam, "hamming_inter_intra_chains_id$(id).jld2"
	)
end

# ╔═╡ a4da1f14-46da-44c2-9803-cc40906375c4
begin
	# settings
	normalize = false
end

# ╔═╡ a9c0f474-48c0-4f5e-b027-28863ea2f25e
md"## Custom hamming distance"

# ╔═╡ 666b3076-3907-4dba-9e33-9aab2e9d4948
begin
	function hamming_onlygaps(
		X::AbstractVector{<:Integer}, Y::AbstractVector{<:Integer}, alphabet; 
		normalize=false, 
	)
		H = sum(zip(X, Y)) do (x, y)
			x_is_gap = PottsEvolver.isgap(x, alphabet)
			y_is_gap = PottsEvolver.isgap(y, alphabet)
			(x_is_gap || y_is_gap) && x != y ? 1 : 0
		end
		return normalize ? H/length(X) : H
	end
	function hamming_onlygaps(
		x::CodonSequence, y::CodonSequence; source=:codon, kwargs...
	)
		return if source == :codon
			hamming_onlygaps(x.seq, y.seq, codon_alphabet; kwargs...)
		elseif source == :aa
			hamming_onlygaps(x.aaseq, y.aaseq, aa_alphabet; kwargs...)
		else
			error("")
		end
	end
end

# ╔═╡ f40ea29f-86f6-4a13-a807-37dbc6e884f7
begin
	function hamming_nogaps(
		X::AbstractVector{<:Integer}, Y::AbstractVector{<:Integer}, alphabet; 
		normalize=false, 
	)
		H = sum(zip(X, Y)) do (x, y)
			x_is_gap = PottsEvolver.isgap(x, alphabet)
			y_is_gap = PottsEvolver.isgap(y, alphabet)
			!x_is_gap && !y_is_gap && x != y ? 1 : 0
		end
		return normalize ? H/length(X) : H
	end
	function hamming_nogaps(
		x::CodonSequence, y::CodonSequence; source=:codon, kwargs...
	)
		return if source == :codon
			hamming_nogaps(x.seq, y.seq, codon_alphabet; kwargs...)
		elseif source == :aa
			hamming_nogaps(x.aaseq, y.aaseq, aa_alphabet; kwargs...)
		else
			error("")
		end
	end
end

# ╔═╡ 6b8602e4-c0bd-4413-b9f1-ab43bef880d3
md"## Measures"

# ╔═╡ 2173ceb7-a959-407a-a7a6-01fb2301004c
md"""
A measure `m` should work like: `m(chain, x) :: Vector{Float64}` with a length equal to the length  of the chain. 
`chain` is a vector of `AbstractSequence`.

The second argument `x` is the initial sequence if `intra_chain_measures` or another chain of the same length (and implicitely the same step values) if `inter_chain_measures`
"""

# ╔═╡ 688e7917-7ba7-4753-8fe8-29b8494c9a8b
intra_chain_measures = [
	"hamming_init_aa" => (chain, init) -> map(chain) do x
		hamming(x, init; source=:aa, normalize)
	end,
	"hamming_init_aa_nogaps" => (chain, init) -> map(chain) do x 
		hamming_nogaps(x, init; source=:aa, normalize)
	end,
	"hamming_init_aa_onlygaps" => (chain, init) -> map(chain) do x 
		hamming_onlygaps(x, init; source=:aa, normalize)
	end,
]

# ╔═╡ 432b8918-eb21-48e6-93cc-a64ef822a607
function do_intra_chain_measures(chains)
	# Assume N chains of length L
	measures = Dict()
	for (name, m) in intra_chain_measures
		# N vectors of length L: hamming distance to init
		X = map(zip(chains.init_sequences, chains.chains)) do (init, chain)
			m(chain, init)
		end
		measures["av_"*name] = mean(X)
		measures["std_"*name] = std(X; mean=measures["av_"*name])
	end
	measures["tvals"] = chains.tvals
	return measures
end

# ╔═╡ f4015b51-0643-4fc7-b84e-579a06ee019e
inter_chain_measures = [
	"hamming_inter_aa" => (C1, C2) -> map(zip(C1, C2)) do (x,y)
		hamming(x, y; source=:aa, normalize)
	end,
	"hamming_inter_aa_nogaps" => (C1, C2) -> map(zip(C1, C2)) do (x,y) 
		hamming_nogaps(x, y; source=:aa, normalize)
	end,
	"hamming_inter_aa_onlygaps" => (C1, C2) -> map(zip(C1, C2)) do (x,y)
		hamming_onlygaps(x, y; source=:aa, normalize)
	end,
]

# ╔═╡ 15d637c1-5afd-4ef4-ac5d-cc375c80056a
function do_inter_chain_measures(chains)
	# Assume N chains of length L
	# take at most 100 random pairs from those
	max_pairs = 500
	N = length(chains.chains)
	index_pairs = [(i, j) for i in 1:N for j in (i+1):N] |> shuffle
	n = length(index_pairs)
	index_pairs = n > max_pairs ? index_pairs[1:max_pairs] : index_pairs
	# compute measures for these pairs
	measures = Dict()
	for (name, m) in inter_chain_measures
		# N vectors of length L: hamming distance to init
		X = map(index_pairs) do (i, j)
			m(chains.chains[i], chains.chains[j])
		end
		measures["av_"*name] = mean(X)
		measures["std_"*name] = std(X; mean=measures["av_"*name])
	end
	measures["tvals"] = chains.tvals
	return measures
end



# ╔═╡ 732d8a9d-22f9-4071-b6c3-7d3ad0ae582d
function do_measures(results_file)
    # ╔═╡ 4f845544-3b4d-4c80-82a3-931d7f0aeeca
    all_data = try
        @unpack df = wload(results_file)
        df
    catch err
        return nothing
    end

    # ╔═╡ bcbdf01c-0a4b-4018-a20e-a75ead42a6e8
    intra_chain_data = let
        df = @select all_data Not(:Teq_long, :time, :chainsteps)
        # X ~ [Dict(measure_name => vector)], with one element per df row
        X = map(r -> do_intra_chain_measures(r.chains), eachrow(df))
        for key in keys(X[1])
            @transform! df $key = map(d -> d[key], X)
        end
        @select df Not(:chains)
    end;

    # ╔═╡ d8715b97-40de-40d6-a8db-4dae15c3d5de
    inter_chain_data = let
        df = @select all_data Not(:Teq_long, :time, :chainsteps)
        # X ~ [Dict(measure_name => vector)], with one element per df row
        X = map(r -> do_inter_chain_measures(r.chains), eachrow(df))
        for key in keys(X[1])
            @transform! df $key = map(d -> d[key], X)
        end
        @select df Not(:chains)
    end;

    # ╔═╡ f0cc2068-c57b-4067-8478-c9255662ad31
    data_measures = begin
        shared_cols = intersect(names(inter_chain_data), names(intra_chain_data))
        innerjoin(inter_chain_data, intra_chain_data, on=shared_cols)
    end;

    # ╔═╡ 07ecf0ea-185f-48e2-acaa-f0d9cb9c7359
    to_save = let
        time = now()
        d = @strdict data_measures time
        @tag! d
    end;

    # ╔═╡ c2fb6b78-4fb8-48ef-9981-acaefa1bbd2f

    return to_save
end

# ╔═╡ 6a48a189-6b2a-40eb-8499-badf8d9b758b
for fam in keys(PED.ref_families)
    @info input_data(fam)
    for (results_file, id) in input_data(fam)
        data = do_measures(results_file)
        !isnothing(data) && wsave(output_data(fam, id), data)
    end
end

# ╔═╡ Cell order:
# ╠═0c5fe0c0-86ee-11ef-05c2-f7cbf0e82146
# ╠═aa55d93a-1a7b-43aa-b95e-083c87eaacb5
# ╠═507356ee-8e8a-4021-bed5-71e86bbf0f56
# ╠═a4da1f14-46da-44c2-9803-cc40906375c4
# ╟─a9c0f474-48c0-4f5e-b027-28863ea2f25e
# ╟─666b3076-3907-4dba-9e33-9aab2e9d4948
# ╟─f40ea29f-86f6-4a13-a807-37dbc6e884f7
# ╟─6b8602e4-c0bd-4413-b9f1-ab43bef880d3
# ╟─2173ceb7-a959-407a-a7a6-01fb2301004c
# ╠═688e7917-7ba7-4753-8fe8-29b8494c9a8b
# ╠═432b8918-eb21-48e6-93cc-a64ef822a607
# ╠═f4015b51-0643-4fc7-b84e-579a06ee019e
# ╠═15d637c1-5afd-4ef4-ac5d-cc375c80056a
# ╠═732d8a9d-22f9-4071-b6c3-7d3ad0ae582d
# ╟─6a48a189-6b2a-40eb-8499-badf8d9b758b
# ╠═4f845544-3b4d-4c80-82a3-931d7f0aeeca
# ╠═bcbdf01c-0a4b-4018-a20e-a75ead42a6e8
# ╠═d8715b97-40de-40d6-a8db-4dae15c3d5de
# ╠═f0cc2068-c57b-4067-8478-c9255662ad31
# ╠═07ecf0ea-185f-48e2-acaa-f0d9cb9c7359
# ╠═c2fb6b78-4fb8-48ef-9981-acaefa1bbd2f
