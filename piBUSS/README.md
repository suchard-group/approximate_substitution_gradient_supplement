# Posterior predictive simulations using piBUSS

## Dependencies
The posterior-predictive simulations require:
- piBUSS (which is already available to users who have followed the instructions for installing BEAST in the main README)
- R, including the package `phangorn` (and all its dependencies). The versions need not be particularly recent.

## File overview

- R scripts
  - `piBUSS/src/model_adequacy.R` Contains helper functions for computing summaries.
  - `piBUSS/src/simulate_and_summarize.R` As named, simulates posterior-predictive alignments and summarizes them.
  - `piBUSS/src/compute_pppvals.R` As named, computes the posterior-predictive p-values.
- Template XMLs used to generate piBUSS instructions.
  - `piBUSS/src/piBUSS_GTR.xml`
  - `piBUSS/src/piBUSS_randomEffects.xml`

## Running the pipeline

### Assumptions
The pipeline makes a number of relatively strong assumptions about the existence of files and their placement.

There must be a file `piBUSS/data/aln.fasta` which contains the alignment for the sequences in fasta format.
This is needed to compute the posterior predictive p-values, and to ensure simulated alignments have ambiguities that match the real one.

Additionally, it is assumed that at least one complete analysis each of the SARS-CoV-2 analysis with both the random-effects model (using HMC) and the GTR model are completed, with the run lengths and thinning settings unchanged.
(The number of replicates is set on line 9 of `simulate_and_summarize.R`, while changes to the length and thinning will necessitate changing lines 12 and 16 to match.)
The output of these completed runs must be placed as follows:
- GTR analyses (log and tree files) go in `approximate_substitution_gradient_supplement/output/SC2_GTR/job_<i>` where `<i>` goes from 1 to the number of replicate chains run (minimum 1)
- HMC-inferred random-effects analyses (log and tree files) go in `approximate_substitution_gradient_supplement/output/SC2_HMC/job_<i>`
If you wish to change this directory structure, amend lines 18, 41-49, and 55-56 of `simulate_and_summarize.R` as needed.

BEAGLE is assumed to be in `/usr/local/lib` as in the main README, and BEAST is assumed to exist in `~/git_repos/beast-mcmc/build/dist/beast.jar` (otherwise amend lines 125-126 of `simulate_and_summarize.R`).

### Simulating
After the appropriate files been added/moved to the appropriate locations, open R, set the working directory to the top level of `approximate_substitution_gradient_supplement`, and run the file `simulate_and_summarize.R` then `compute_pppvals.R`.

From R,
```
setwd(path/to/approximate_substitution_gradient_supplement)
source(piBUSS/src/simulate_and_summarize.R)
source(piBUSS/src/compute_pppvals.R)
```
Note that the simulation and summarization step (even though piBUSS is impressively fast for the number of genomes being simulated) takes some time.

Temporary files, named `sequences.fasta` and `tmp_pibuss.xml` will be generated in the top level of `approximate_substitution_gradient_supplement` and may be deleted afterwards.