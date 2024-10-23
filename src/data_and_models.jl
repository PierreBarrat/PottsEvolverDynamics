# From the equilibration_fgaps notebook: decide on fgap for each type of step
# I choose these parameters because they seem to equilibrate gaps as fast as the other states.
# Teq is kind of a "conservative decision" here: wait for intra chain hamming to reach inter chain hamming
const potts_equilibration_codons = SortedDict(
    "PF00014" => SortedDict(
        "accepted" => (Teq_steps = 5000, Hmax = 34),
        "changed" => nothing,
        "proposed" => nothing,
    ),
    "PF00072" => SortedDict(
        "accepted" => (Teq_steps = 50_000, Hmax = 83),
        "changed" => nothing,
        "proposed" => nothing,
    ),
    "PF00076" => SortedDict(
        "accepted" => (Teq_steps = 5_000, Hmax = 54),
        "changed" => nothing,
        "proposed" => nothing, # did not run long enough simulations yet, and doubt I'll need it
    ),
    "PF00105" => SortedDict(
        "accepted" => (Teq_steps = 50_000, Hmax = 48),
        "changed" => nothing,
        "proposed" => nothing,
    ),
    "PF00583" => SortedDict(
        "accepted" => (Teq_steps = 100_000, Hmax = 94),
        "changed" => nothing,
        "proposed" => nothing,
    ),
    "PF00595" => SortedDict(
        "accepted" => (Teq_steps = 3_000, Hmax = 60),
        "changed" => nothing,
        "proposed" => nothing,
    ),
    "PF01817" => SortedDict(
        "accepted" => (Teq_steps = 1_000_000, Hmax = 75),
        "changed" => nothing,
        "proposed" => nothing,
    ),
)

const potts_models = SortedDict(
    "PF00072" => OrderedDict(
        "file" => projectdir("models/PF00072/parameters_dcatools_PF00072.dat"),
        "sample" => nothing,
        "method" => "dcatools",
        "equilibration_codons" => potts_equilibration_codons["PF00072"],
        "fgap" => 0.5,
    ),
    "PF00076" => OrderedDict(
        "file" => projectdir("models/PF00076/parameters_adaBM_PF00076.dat"),
        "sample" => nothing,
        "method" => "adabmdca",
        "equilibration_codons" => potts_equilibration_codons["PF00076"],
        "fgap" => 0.9,
    ),
    "PF00595" => OrderedDict(
        "file" => projectdir("models/PF00595/parameters_adaBM_PF00595.dat"),
        "sample" => nothing,
        "method" => "adabmdca",
        "equilibration_codons" => get(potts_equilibration_codons, "PF00595", nothing),
        "fgap" => 0.9,
    ),
    "PF00105" => OrderedDict(
        "file" => projectdir("models/PF00105/Parameters_conv_Thr-PCD40.dat"),
        "sample" => nothing,
        "method" => nothing,
        "equilibration_codons" => get(potts_equilibration_codons, "PF00105", nothing),
        "fgap" => 0.95,
    ),
    "PF00583" => OrderedDict(
        "file" => projectdir("models/PF00583/Parameters_conv_AAC6.dat"),
        "sample" => nothing,
        "method" => nothing,
        "equilibration_codons" => get(potts_equilibration_codons, "PF00583", nothing),
        "fgap" => 0.95,
    ),
    "PF01817" => OrderedDict(
        "file" => projectdir("models/PF01817/Parameters_conv_CM-PCD40.dat"),
        "sample" => nothing,
        "method" => nothing,
        "equilibration_codons" => get(potts_equilibration_codons, "PF01817", nothing),
        "fgap" => 0.95,
    ),
    "PF00014" => OrderedDict(
        "file" => projectdir("models/PF00014/parameters_adaBM_PF00014.dat"),
        "sample" => nothing,
        "method" => nothing,
        "equilibration_codons" => get(potts_equilibration_codons, "PF00014", nothing),
        "fgap" => 0.75,
    )
)

const profile_models = SortedDict(
    "PF00014" => OrderedDict(
        "file" => projectdir("models/PF00014/parameters_profile_PF00014.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
    "PF00072" => OrderedDict(
        "file" => projectdir("models/PF00072/parameters_profile_PF00072.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
    "PF00076" => OrderedDict(
        "file" => projectdir("models/PF00076/parameters_profile_PF00076.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
    "PF00105" => OrderedDict(
        "file" => projectdir("models/PF00105/parameters_profile_PF00105.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
    "PF00583" => OrderedDict(
        "file" => projectdir("models/PF00583/parameters_profile_PF00583.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
    "PF00595" => OrderedDict(
        "file" => projectdir("models/PF00595/parameters_profile_PF00595.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
    "PF01817" => OrderedDict(
        "file" => projectdir("models/PF01817/parameters_profile_PF01817.dat"),
        "sample" => nothing,
        "pseudocount" => 1e-3,
    ),
)

const ref_families = SortedDict(
    "PF00072" => OrderedDict(
        "name" => "PF00072",
        "long_name" => "Response regulator receiver domain",
        "alignment" => projectdir("models/PF00072/PF00072_mgap6_subsample.fasta"),
        "L" => 112,
        "weights" => projectdir("models/PF00072/PF00072_mgap6_subsample_weights.csv"),
        "notes" => "ArDCA paper - subsample because original alignment very large",
        "potts" => get(potts_models, "PF00072", nothing),
        "profile" => get(profile_models, "PF00072", nothing),
    ),
    "PF00076" => OrderedDict(
        "name" => "PF00076",
        "long_name" => "RNA recognition motif",
        "alignment" => projectdir("models/PF00076/PF00076_mgap6.fasta"),
        "L" => 70,
        "weights" => projectdir("models/PF00076/PF00076_mgap6_weights.csv"),
        "notes" => "ArDCA paper",
        "potts" => get(potts_models, "PF00076", nothing),
        "profile" => get(profile_models, "PF00076", nothing),
    ),
    "PF00595" => OrderedDict(
        "name" => "PF00595",
        "long_name" => "PDZ domain",
        "alignment" => projectdir("models/PF00595/PF00595_mgap6.fasta"),
        "notes" => "ArDCA paper",
        "L" => 82,
        "weights" => projectdir("models/PF00595/PF00595_mgap6_weights.csv"),
        "potts" => get(potts_models, "PF00595", nothing),
        "profile" => get(profile_models, "PF00595", nothing),
    ),
    "PF00105" => OrderedDict(
        "name" => "PF00105",
        "long_name" => "DNA binding domain",
        "long_name_2" => "Zinc finger C4 type",
        "alignment" => projectdir("models/PF00105/DBD_alignment.uniref90.cov80.noX.a2m"),
        "L" => 76,
        "weights" => projectdir("models/PF00105/DBD_alignment_weights.uniref90.cov80.noX.csv"),
        "notes" => "EvoLeonardo",
        "potts" => get(potts_models, "PF00105", nothing),
        "profile" => get(profile_models, "PF00105", nothing),
    ),
    "PF00583" => OrderedDict(
        "name" => "PF00583",
        "long_name" => "Acetyltransferase",
        "alignment" => projectdir("models/PF00583/PF00583_noinsert_max11gaps_nodupl_noclose_noX.faa"),
        "L" => 117,
        "weights" => projectdir("models/PF00583/PF00583_noinsert_max11gaps_nodupl_noclose_noX_weights.csv"),
        "notes" => "EvoLeonardo",
        "potts" => get(potts_models, "PF00583", nothing),
        "profile" => get(profile_models, "PF00583", nothing),
    ),
    "PF01817" => OrderedDict(
        "name" => "PF01817",
        "long_name" => "Chorismate mutase type II",
        "alignment" => projectdir("models/PF01817/CM_130530_MC_noX.fasta"),
        "L" => 96,
        "weights" => projectdir("models/PF01817/CM_130530_MC_noX_weights.csv"),
        "notes" => "EvoLeonardo - the HMM on interpro is of length 80, why?",
        "potts" => get(potts_models, "PF01817", nothing),
        "profile" => get(profile_models, "PF01817", nothing),
    ),
    "PF00014" => OrderedDict(
        "name" => "PF00014",
        "long_name" => "Trypsin inhibitor Kunitz domain",
        "alignment" => projectdir("models/PF00014/PF00014_mgap6.fasta"),
        "L" => 53,
        "weights" => projectdir("models/PF00014/PF00014_mgap6_weights.csv"),
        "notes" => "ArDCA paper",
        "potts" => get(potts_models, "PF00014", nothing),
        "profile" => get(profile_models, "PF00014", nothing),
    ),
)

_make_paths_local(data) = data
function _make_paths_local(data::AbstractDict)
    for (k, v) in data
        if v isa AbstractString && ispath(v)
            data[k] = relpath(v)
        else
            vnew = _make_paths_local(v)
            data[k] = vnew
        end
    end
    return data
end

open(projectdir("models/families.yml"), "w") do io
    D = deepcopy(ref_families)
    _make_paths_local(D)
    YAML.write(io, D)
end


