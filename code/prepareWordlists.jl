

##
flag = string(strip(ARGS[1], '/'))

_, type, dataset, fn = split(flag, "/")

proportion = replace(fn, r"[^-]*-" => "", ".tsv" => "")

project = "tmp/$type/$dataset-$proportion"

##

try
    mkdir("tmp")
catch e
end

try
    mkdir("tmp/$type")
catch e
end

try
    mkdir("tmp/$type/$dataset-$proportion")
catch e
end


##


cd(project)

##

using DataFrames
using CSV
using DataStructures
using PyCall

##

lp = pyimport("lingpy")

##

trainingF = "../../../ST2022/$type/$dataset/training-$proportion.tsv"
testF = "../../../ST2022/$type/$dataset/test-$proportion.tsv"

tr = CSV.File(trainingF, types=String) |> DataFrame
ts = CSV.File(testF, types=String) |>  DataFrame

##

languages = names(tr)[2:end]

##

tr_stacked = dropmissing(
    stack(
        tr,
        languages,
        variable_name=:TAXON,
        value_name=:TOKENS,
    ),
    :TOKENS,
)
insertcols!(tr_stacked, 1, :GLOSS => tr_stacked.COGID)
insertcols!(tr_stacked, :IPA => replace.(tr_stacked.TOKENS, " " => ""));


##

ts_stacked = dropmissing(
    stack(
        ts,
        languages,
        variable_name=:TAXON,
        value_name=:TOKENS,
    ),
    :TOKENS,
)


insertcols!(ts_stacked, 1, :GLOSS => ts_stacked.COGID)
insertcols!(ts_stacked, :IPA => replace.(ts_stacked.TOKENS, " " => ""))
insertcols!(ts_stacked, 1, :ID => 1:nrow(ts_stacked))

##

function ipa2sca(w)
    if w=="?"
        return w
    end
    v = lp.tokens2class(split(w), "sca")
    filter!(x -> x != "+", v)
    join(v, " ")
end

function ipa2dolgo(w)
    if w == "?"
        return w
    end
    v = lp.tokens2class(split(w), "dolgo")
    filter!(x -> x != "+", v)
    join(v, " ")
end

function ipa2asjp(w)
    if w == "?"
        return w
    end
    v = lp.tokens2class(split(w), "asjp")
    filter!(x -> x != "+", v)
    join(v, " ")
end



##

insertcols!(
    tr_stacked,
    :SCA => ipa2sca.(tr_stacked.TOKENS)
)

insertcols!(
    tr_stacked,
    :DOLGO => ipa2dolgo.(tr_stacked.TOKENS)
)

insertcols!(
    tr_stacked,
    :ASJP => ipa2asjp.(tr_stacked.TOKENS)
)


##


insertcols!(
    ts_stacked,
    :SCA => ipa2sca.(ts_stacked.TOKENS)
)

insertcols!(
    ts_stacked,
    :DOLGO => ipa2dolgo.(ts_stacked.TOKENS)
)

insertcols!(
    ts_stacked,
    :ASJP => ipa2asjp.(ts_stacked.TOKENS)
)


##

try
    mkdir("data")
catch e
end

CSV.write("data/training.tsv", tr_stacked, delim="\t")

CSV.write("data/test.tsv", ts_stacked, delim="\t")
