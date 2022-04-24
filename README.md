# gerhardSigtyp2022



This repository contains the code for my submission to the [**SIGTYP 2022 Shared Task: Prediction of Cognate Reflexes**](https://github.com/sigtyp/ST2022). 

## Software

The following software has to be installed to run the code:

- julia v.1.7.2
- R v.3.5.3 
- R-package *ape* (to be installed with `install.packages("ape")`)
- python 3.7.6
- Python package *lingpy* (to be installed from the shell via `pip install lingpy`)
- (ideally mpi-enabled version of) *MrBayes* 

### Using Conda

Instead of installing the right versions of *python*, *R* and *ape*, you can create *conda* environment. Edit the `prefix` line in the *yml*-files according to your local configuration and run the command in the repository’s root directory:

`conda env create -f sigtyp2022.yml && conda activate sigtyp2022`

For the *MrBayes* part, you can do

`conda env create -f mrbayes.yml && conda activate mrbayes`

### Julia project directory

Julia’s project directory has to be set to the subdirectory `code`. You can achieve this either with setting the environment variable

`export JULIA_PROJECT=code`

or to run `julia` with the flag `--project=code`.

### Using Singularity

The least painful way to reproduce my workflow is to use a Singularity container. `gjSigtyp2022.def` is a Singularity definition file for a Singularity image. You can also download the Singularity image from 'https://unitc-my.sharepoint.com/:u:/g/personal/nwsja01_cloud_uni-tuebingen_de/EbCMvXu0LSBLlqqv84H56qUBGmeTv90-dl0k9q4C-KtyMA?e=cu6jMt&download=1' and put it into the repository’s root directory.

Be warned, it’s 3 GB.

Then you can run all commands via the image without further worries. Save the image in the project’s root directory and prefix each command with

`singularity exec gjSigtyp2022.sif`

When using the Singularity container, you don’t have to alter the Julia project directory.



## Workflow

- place a copy of https://github.com/sigtyp/ST2022 into the project’s root directory, e.g., by running
  `git clone https://github.com/sigtyp/ST2022`

- `<trainingFile>` is the path to the training file you want to work with, e.g. 

  `ST2022/data/allenbai/test-0.10.tsv`

- run the following commands in that order:

  ```shell
  export dir=`echo <trainingFile> | sed 's/.training//' | sed 's/.tsv//' | sed 's/ST2022/tmp/'` 
  julia code/prepareWordlists.jl <trainingFile>
  julia code/baumWelchTraining.jl <trainingFile>
  julia code/upgmaTree.jl <trainingFile>
  julia code/tcoffeeAlignment.jl <trainingFile>
  julia code/createMBfile.jl <trainingFile>
  mpirun -np 4 mb $dir/mrbayes/mb.nex
  bash code/extractPosteriorTrees.sh <trainingFile>
  julia code/estimateRates.jl <trainingFile>
  julia code/reconstruction.jl <trainingFile>
  ```

  

During computation, a directory `tmp` will be created, containing lots of files with intermediate stages of the workflow. The results of the computation are saved as result tsv files inside the `ST2022` directory.
