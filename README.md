# PottsEvolverDynamics

## Structure

Contains
- alignments of natural sequences
- potts and profile models learned on these alignments
- simple figures to estimate the equilibration time of said potts models and their generative capacities

To unzip the models, run `python unzip_all.py`. The `pyyaml` package should be installed

The files are organized like so
- `models/families.yml` contains information about all families.  
  To read it in Julia, `using YAML; families = YAML.load_file("models/families.yml")`. 
  This returns a dictionary, with pfam identifiers as keys. 
- `models/PFXXX` contains the alignment and models for the corresponding family. 
  How the alignment was obtained is in `families[PFXXX][notes]`. 
  If I have the information on how the Potts model was inferred, it's in `families[PFXXX]["potts"]["method"]`. 
- `plots/equilibration` contains figures about equilibration of the models. 
  For each model (potts/profile) and sampling technique (aa/codons) I start `n=50` chains at random equilibrated sequences, and look at the average Hamming distance to the start of the chain, and the average inter-chain Hamming distance. 
- `plots/generative` contains figures about generative properties of the Potts models: 
   - fitting quality for the first two moments (single site frequencies/connected correlations)
   - pairwise Hamming distance between sequences of a sample, vs natural sequences


## Instructions to run scripts

*Note*: First, unzip the models by running `python unzip_all.py` at the root of the repo. 
The `pyyaml` package should be installed.

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> PottsEvolverDynamics

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "PottsEvolverDynamics"
```
which auto-activate the project and enable local path handling from DrWatson.
