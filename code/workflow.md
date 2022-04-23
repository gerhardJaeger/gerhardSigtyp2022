# Workflow

- `wget -O gjSigtyp2022.sif  'https://unitc-my.sharepoint.com/:u:/g/personal/nwsja01_cloud_uni-tuebingen_de/EbCMvXu0LSBLlqqv84H56qUBGmeTv90-dl0k9q4C-KtyMA?e=cu6jMt&download=1'`
- `julia code/prepareWordlists.jl <project>` 
- `julia code/baumWelchTraining.jl <project>`
- `julia code/upgmaTree.jl <project>`
- `julia code/tcoffeeAlignment.jl <project>`
- `julia code/createMBfile.jl <project>`
- `mpirun -np 4 mb <project>/mrbayes/mb.nex`
- `bash code/extractPosteriorTrees.sh <project>`
- `julia code/estimateRates.jl <project>`
- `julia code/reconstruction.jl <project>`



