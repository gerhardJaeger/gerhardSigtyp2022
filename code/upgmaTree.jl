project = ARGS[1]
dataset, proportion = split(project,"-")


cd(project)

##

SC = :ASJP

##
using CSV, DataFrames
using SequenceAlignment
using Random
Random.seed!(12345)
using ProgressMeter
using Statistics
using MCPhylo
##


tr = CSV.File("data/training.tsv", types=String) |> DataFrame


##


ld_ = []

cClasses = unique(tr.COGID)

for cc in cClasses
    ccData = filter(x -> x.COGID == cc, tr)
    for i in 1:nrow(ccData), j in 1:nrow(ccData)
        l1, l2 = ccData.TAXON[[i, j]]
        w1, w2 = ccData[[i, j], SC]
        d = ldn(w1, w2)
        push!(
            ld_,
            [l2, l1, d]
        )
    end
end



##

ld = DataFrame(
    l1=[x[1] for x in ld_],
    l2=[x[2] for x in ld_],
    ldn=[x[3] for x in ld_]
)

##

ldnDF = combine(
    groupby(ld, [:l1, :l2]),
    :ldn => mean => :dist
)

##

dstMtx_ = unstack(ldnDF, :l1, :l2, :dist)

languages = dstMtx_[:,:l1]

dstMtx = dstMtx_[:,languages] |> Array#

mn = mean(ldnDF.dist)
dstMtx[ismissing.(dstMtx)] .= mn

dstMtx = convert(Matrix{Float64}, dstMtx)

tree = upgma(dstMtx, languages)

open("data/upgmaTree.tre", "w") do file
    write(
        file,
        newick(tree)
    )
end
