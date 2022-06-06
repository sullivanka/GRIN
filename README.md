# GRIN
Geneset Refinement using Interacting Networks

A command line script to filter gene sets using biological networks. This tool builds off of the [RandomWalkRestartMH](https://github.com/alberto-valdeolivas/RandomWalkRestartMH) R package, and can be combined with the [RWRtoolkit](https://github.com/dkainer/RWRtools) R package using the [RWR_make_multiplex.R](https://github.com/dkainer/RWRtools/packageAttempt1/R/RWR_make_multiplex.R) function to build a user-defined multiplex network.

## Installation using Conda (Conda version 4.12.0)
1. Verify Conda version (version 4.12.0).
  ```
  (base) $ conda -V
  conda 4.12.0
  ```
  
  If conda version is not updated, update conda:
`(base) $ conda update -n base -c defaults conda`

1. Clone GRIN:
    `git clone https://github.com/sullivanka/GRIN.git`
2. Install the conda env (~5 mins):
    `(base) $ conda create --name GRIN --channel conda-forge r-base=4.0.2 r-devtools r-doparallel r-essentials r-igraph r-optparse r-signal`
3. Activate conda environment:
    `(base) $ conda activate GRIN`
4. Activate R and install additional packages. If prompted update all packages as necessary. Note that this process will take some time (~20 mins) from a clean environment.
    ```
    (GRIN) $ R
    > library(devtools)
    > devtools::install_github("dkainer/RandomWalkRestartMH")
    > devtools::install_github("agentlans/KneeArrower")
    ```
## Testing Installation
1. Activate GRIN conda environment:
  `(base) $ conda activate GRIN`
2. Run `testGRIN.sh`:
  `(GRIN) $ bash testGRIN.sh`
  NOTE: If running on a non-Linux system, you will need to ensure `gawk` is installed at the command line.  
  Mac users: Install `gawk` by installing [Homebrew](https://brew.sh/), followed by  
  ` $ brew install gawk`  
  Windows users: Install [Gawk for Windows](http://gnuwin32.sourceforge.net/packages/gawk.htm)  
3. Upon successful running of GRIN, the test script will run GRIN followed by the following output printed to console:
```
GRIN Installation successful!
Removing test files...
Installation test complete.
```
    
## Running GRIN
1. Ensure your input files are in 2 tab-separated columns. The first column corresponds to the user-defined gene set name, and the second column corresponds to the gene ID for the gene of interest. NOTE: this gene ID MUST match the gene IDs used in the multiplex network. See example file in `test/TestGenes.txt`.
2. Activate GRIN conda environment using `conda activate GRIN` as described above.
3. Use pre-assembled multiplex network RData object (example multiplex is `test/suicide_weighted_Multiplex_0.5Delta.RData`) and command line parameters to identify which genes are retained or removed by GRIN.
Example command line arguments:
` (GRIN) $ Rscript $GRIN_DIR/R/GRIN.R -d $GRIN_DIR/test/suicide_weighted_Multiplex_0.5Delta.RData -g $GRIN_DIR/test/TestGenes.txt -r 0.7 -m Test_Install_User --tau 1,1,1,1,1,1,1,1,1,1 -o $GRIN_DIR/test/test_output `
4. Flags to include for command line script:
`-d or --data (Required):` Path to multiplex RData object, built using RandomWalkRestartMH R package. This multiplex network contains a delta value indicating the probability of staying in a given network layer or jumping between layers (example multiplex network uses delta = 0.5).  
`-g or --geneset (Required):` Path to tab-separated user geneset. See step 1 for instructions about gene set formatting.  
`-o or --outdir (Required):` Path for GRIN output files to be written.  
`-r or --restart (Required, default value = 0.7):` Probability for random walk algorithm to restart at a input gene in the multiplex network at any given point.  
`-t or --tau (Required, default value = 1 for each layer):` Comma-separated values indicating relative probability for random walk to restart in a given network layer.  
`-m or --modname (Optional):` User-defined name for an individual GRIN run in order to distinguish outputs from independent runs.  
`-p or --plot (Optional):` If flag used, outputs Mann-Whitney U sliding window plot and associated matrix of p-values between null distribution and user geneset at each window.  
`-s or --simple-filenames (Optional):` If flag used, shortens the names of output files written by GRIN.  
`--threads (Optional, default value = all system cores - 1):` If present, defines the number of threads to use for GRIN based on system capabilities.  
`-v or --verbose (Optional):` If flag used, allows more output for user to see as GRIN is running.  

5. Retained genes and removed genes are written to file at the directory indicated by `--outdir`. Any duplicate genes or genes not in the multiplex network from the initial user gene set will also be written to file in this directory, along with the sliding window matrix and Mann-Whitney U plot if `--plot` is used.
  
