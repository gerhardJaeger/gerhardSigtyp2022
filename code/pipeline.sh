julia code/prepareWordlists.jl $1 $2
julia code/baumWelchTraining.jl $1 $2
julia code/upgmaTree.jl $1 $2
julia code/tcoffeeAlignment.jl $1 $2
julia code/createMBfile.jl $1 $2
mpirun -np 4 mb $1/mrbayes/mb.nex
bash code/extractPosteriorTrees.sh $1 $2
julia code/estimateRates.jl $1 $2
julia code/reconstruction.jl $1 $2