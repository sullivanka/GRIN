# GRIN
Geneset Refinement using Interacting Networks


A command line script to filter gene sets using biological networks. This tool builds off of the [RandomWalkRestartMH](https://github.com/alberto-valdeolivas/RandomWalkRestartMH) R package, and can be combined with the [RWRtoolkit](https://github.com/dkainer/RWRtoolkit) R package using the [RWR_make_multiplex.R](https://github.com/dkainer/RWRtoolkit/R/RWR_make_multiplex.R) function to build a user-defined multiplex network.

To understand more about GRIN and how we validated this tool, [check out our preprint on medRxiv](https://www.medrxiv.org/content/10.1101/2022.04.20.22273895v2).

#### Developed by:
* David Kainer (@github/dkainer)
* Kyle Sullivan (@github/sullivanka)

#### Overview of GRIN

Biological studies such as GWAS or differential expression experiments produce sets of genes that pass a significance threshold. Such a set of genes is likely to contain genes that truly affect the trait or condition, but also false positives, and you really can't tell which ones are which! One can assume, however, that true positive genes associated with a trait are probably functionally related in some manner within the inner workings of molecular biology. False positive genes, on the other hand, are probably functionally random with respect to each other and the true positive genes. GRIN operates on this concept. It takes functional networks that define the functional relationships (e.g. PPI, coexpression, knockout effect, metabolic pathways) between all genes in the genome, and figures out which genes in the gene set appear to be connected more strongly than random genes using rankings generated by the Random Walk with Restart algorithm. The output is two gene sets: Retained (the functional genes), and Removed (the apparently random false positive genes).

#### Using GRIN

To use our implementation of GRIN you need to preform a few steps:

0. [Install GRIN](https://github.com/sullivanka/GRIN/edit/main/README.md#installation-using-conda-conda-version-4110)
1. Get a set of genes of interest, such as candidate genes from a GWAS
2. Generate a multiplex network that captures biological relationships between genes
3. Run GRIN using the inputs from (1) and (2) 

For step (2), the multiplex network can be as produced from just a single network layer, such as co-expression edges. However, the more lines of evidence used in your network, the better the outcomes of GRIN. The **RWR_make_multiplex.R** command in the [RWRtoolkit](https://github.com/dkainer/RWRtoolkit) R package is the method to take one or more functional networks and turn them into a multiplex network containing multiple lines of evidence. The multiplex network will be saved in the form of an RData object file, which you then use as an input to GRIN. Comand line instructions are [here](https://github.com/sullivanka/GRIN#running-grin).



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
3. Install Git LFS (~1 min):
    `(base) $ git lfs install`
    For more information on Git LFS, see `https://git-lfs.github.com/`.
4. Activate conda environment:
    `(base) $ conda create --name GRIN --channel conda-forge r-base=4.0.2 r-devtools r-doparallel r-essentials r-igraph r-optparse r-signal`
5. Activate conda environment:
    `(base) $ conda activate GRIN`
6. Activate R and install additional packages. If prompted update all packages as necessary. Note that this process will take some time (~20 mins) from a clean environment.
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
  `(GRIN) $ source GRIN.sh`
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
  
