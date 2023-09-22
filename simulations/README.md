# Generating XMLs to analyze simulated data

## Dependencies
To generate the XMLs required to run the analyses of simulated data requires:
- piBUSS (which is already available to users who have followed the instructions for installing BEAST in the main README).
- R, including the package `phangorn` (and all its dependencies). The versions need not be particularly recent.

## File overview

- R scripts
  - `simulations/src/helper.R` Contains helper function for formatting XML alignment.
  - `simulations/src/simulate_and_make_xml.R` As named, simulates alignments and makes them into runnable XMLs for analyses.
- Template XMLs used to generate XMLs for piBUSS and for analysis.
  - `piBUSS/src/piBUSS_randomEffects.xml`
  - `simulations/src/mcmc_template.xml`

## Running the pipeline

### Assumptions
Like the posterior predictive simulation pipeline, this pipeline makes a number of relatively strong assumptions about the existence of files and their placement.

It is assumed that at least one complete analysis of the SARS-CoV-2 analysis with the random-effects model (using HMC) has been completed, with the run lengths, thinning settings, and file output names unchanged.
(The number of replicates is set on line 9 of `simulate_and_make_xml.R`, while changes to the length and thinning will necessitate changing lines 14 and 18 to match.)

The output of all completed HMC-inferred random-effects analyses (log and tree files) must go in `approximate_substitution_gradient_supplement/output/SC2_HMC/job_<i>`.

If you wish to change this directory structure, amend lines 21, 23-26, and 33-34 of `simulate_and_make_xml.R` as needed.

BEAGLE is assumed to be in `/usr/local/lib` as in the main README, and BEAST is assumed to exist in `~/git_repos/beast-mcmc/build/dist/beast.jar` (otherwise amend lines 125-126 of `simulate_and_make_xml.R`).

### Simulating
After the appropriate files been added/moved to the appropriate locations, open R, set the working directory to the top level of `approximate_substitution_gradient_supplement`, and run the file `simulate_and_make_xml.R`.

This will produce a folder `approximate_substitution_gradient_supplement/simulations/xml` and a file `approximate_substitution_gradient_supplement/simulations/simulating_values.csv`.
Each simulation replicate's true parameter values are recorded in the CSV file.
Each simulation replicate produces 4 XMLs in the XML directory, which analyze the simulated data under different settings of the Bayesian bridge exponent parameter.
These analysis XMLs may be run in BEAST in the usual way.

From R,
```
setwd(path/to/approximate_substitution_gradient_supplement)
source(piBUSS/src/simulate_and_make_xml.R)
```
Note that the simulation and XML generation step (even though piBUSS is impressively fast for the number of genomes being simulated) takes some time.

Temporary files, named `sequences.fasta` and `tmp_pibuss.xml` will be generated in the top level of `approximate_substitution_gradient_supplement` and may be deleted afterwards.