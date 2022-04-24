export dir=`echo $1 | sed 's/.training//' | sed 's/.tsv//' | sed 's/ST2022/tmp/'` 
julia code/prepareWordlists.jl $1
julia code/baumWelchTraining.jl $1
julia code/upgmaTree.jl $1
julia code/tcoffeeAlignment.jl $1
julia code/createMBfile.jl $1
mpirun -np 4 mb $dir/mrbayes/mb.nex
bash code/extractPosteriorTrees.sh $dir
julia code/estimateRates.jl $1
julia code/reconstruction.jl $1

