# Random-effects substitution models for phylogenetics via scalable gradient approximations
This repository contains the scripts and XML files required to reproduce the analyses performed in the "Random-effects substitution models for phylogenetics via scalable gradient approximations" paper by Magee et al.

## Installing BEAST/BEAGLE
These analyses require installation of the `hmc-clock` branches of both BEAST and BEAGLE.

## Setting up BEAGLE
Please follow the [BEAGLE installation instructions](https://github.com/beagle-dev/beagle-lib), but be sure to get the `hmc-clock` branch.

For Mac users, the following commands will compile the CPU version of BEAGLE.
Follow the [instructions](https://github.com/beagle-dev/beagle-lib) if you need to install any other dependent software, ignore the first 2 lines if you already have all requisite dependencies installed.

```
xcode-select --install
brew install libtool autoconf automake
git clone https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib
git checkout hmc-clock
mkdir build
cd build
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
sudo make install
```


For Linux users, the commands are similar.

```
sudo apt-get install build-essential autoconf automake libtool git pkg-config openjdk-9-jdk
git clone https://github.com/beagle-dev/beagle-lib.git
cd beagle-lib
git checkout hmc-clock
mkdir build
cd build
cmake -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF ..
sudo make install
```

The libraries are installed into `/usr/local/lib`.

### Setting up BEAST

The following commands will compile the `hmc-clock` branch of BEAST.

```
git clone https://github.com/beast-dev/beast-mcmc.git
cd beast-mcmc
git checkout hmc-clock
ant
```

For Mac users, you may need to install ant using `brew install ant` via [Homebrew](https://brew.sh/).

For Linux users, you can install ant by `sudo apt-get install ant`.

This will compile the `jar` file `beast.jar` in `beast-mcmc/build/dist/beast.jar` which can be called to execute BEAST.

## Reproducing the analyses

You may use the following commands to run the analyses described in the manuscript.

From the top level of `approximate_substitution_gradient_supplement`, create and change directories to `output`.
```
mkdir output
cd output
```

As written, run commands will lead to the output files being placed in `approximate_substitution_gradient_supplement/output`.
If you wish them to be elsewhere, you may set the working directory to any directory where you want to store the resulting log files first, and amend the path to the xmls in the commands as needed.

Where four replicate chains were run, the first seed is provided. The other seeds are 6662, 6663, and 6664.

### SARS-CoV-2
Due to data sharing limitations, the provided XMLs run the SARS-CoV-2 analyses in the paper under the prior.

The SARS-CoV-2 dataset is all genomes from China collected through 16 July 2020 and queried from GISAID (refer to [Pekar et al., 2021,](https://www.science.org/doi/full/10.1126/science.abf8003) for further details).
The GISAID accession IDs are available in `acknowledgements_table.xslx`.

The analyses on this dataset are:
- Random-effects with HMC
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/SC2_HMC.xml
  ```
- Random-effects with MH-MCMC
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/SC2_MH.xml
  ```
- GTR analysis (standard MH-MCMC moves, for comparison to random-effects with posterior-predictive simulations)
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/SC2_GTR.xml
  ```

The posterior-predictive simulations using `piBUSS` are explained in the README in that folder.

### Influenza A

This analysis requires a set of empirical trees, which have been compressed.
First these must be extracted (and placed in the `xmls` folder).
On Mac or Linux this can be accomplished with `tar`:
```
tar -xzf ../xmls/airCommunitiesMM_500.trees.tar.gz -C ../xmls/
```

With this accomplished, the following commands will run the MH and HMC analyses.
- MH
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/airCommunitiesMM_MH.xml
  ```
- HMC
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/airCommunitiesMM_HMC.xml
  ```

### Metazoa

These XMLs are for parameter optimzation and such routines output directly to screen unless it is captured.
The key lines of the output are:
- `X: [...]`: the location of the mode (comma-separated)
- `Gradient: [...]`: the gradient at the mode (comma-separated)
- `Fx:`: the (negative) value of the posterior density at the mode


The following commands will run MAP optimization with numeric and approximate gradients.
- Numeric
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/optimize_aa_num.xml > aa_numerical_output.txt
  ```
- Approximate
  ```
  java -Djava.library.path=/usr/local/lib \
  -jar /path/to/beast-mcmc/build/dist/beast.jar \
  -seed 6661 -overwrite \
  ../xmls/optimize_aa_approx.xml > aa_approximate_output.txt
  ```
- Laplace approximation. 
    - The approximate inference run is also set up to make a Laplace approximation to the posterior at the mode. 
    - In a Laplace approximation, the posterior is approximated as a multivariate Normal. From this we can also obtain a univariate Normal approximation to the marginal distribution of each variable.
    - This adds several lines to the output. The key lines are:
        - `Mean vector: [...]`: the mean vector of the multivariate Normal approximation (in the unbounded transformed space). For this analysis this is the same as `X`.
        - `Stdev vector: [...]`: for each variable, this is the standard deviation of the Normal approximation to its marginal posterior distribution.
