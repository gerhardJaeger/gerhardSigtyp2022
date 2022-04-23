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

##


tr = CSV.File("data/training.tsv", types=String) |> DataFrame

ts = CSV.File("data/test.tsv", types=String) |> DataFrame

##

trPairs = []
cClasses = tr.COGID |> unique

for cc in cClasses
    ccData = filter(x -> x.COGID == cc, tr)[:,SC]
    for i in 2:length(ccData), j in 1:i
        w1 = convert(Vector{String}, split(ccData[i]))
        w2 = convert(Vector{String}, split(ccData[j]))
        push!(trPairs, [w1, w2])
    end
end

shuffle!(trPairs)

##


alphabet = unique(vcat(vcat(trPairs...)...))

##

EPOCHS = 2
BATCHES = 20
##

batches = []
counter = 1
batchLength = length(trPairs) ÷ BATCHES
for i in 1:BATCHES
    e = minimum([counter + batchLength, length(trPairs)])
    push!(
        batches,
        trPairs[counter:e],
    )
    global counter = e + 1
end


##

function updateBatch(ph, batch)
    phe = phmmExpectations(ph0)
    for pr in batch[1:end]
        phe += baumWelch(pr..., ph)[2]
    end
    Phmm(phe)
end

##

ph0 = SequenceAlignment.levPhmm(alphabet)

##

phmmSequence = [ph0]
@showprogress for e in 1:EPOCHS, batch in batches
    phNew = updateBatch(phmmSequence[end], batch)
    push!(phmmSequence, phNew)
end

##

write("phmmTraining.json", phmmSequence[end])