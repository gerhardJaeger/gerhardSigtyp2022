
##
flag = string(strip(ARGS[1], '/'))

_, type, dataset, fn = split(flag, "/")

proportion = replace(fn, r"[^-]*-" => "", ".tsv" => "")

project = "tmp/$type/$dataset-$proportion"


cd(project)


##

using MCPhylo, MCPhyloTree
using CSV, DataFrames
using Pipe
using StatsFuns
using ProgressMeter
using Random
Random.seed!(12345)

##



import MCPhylo: datafortree


function datafortree(
    df::Matrix{A},
    leave_names::Vector{B},
    tree::T,
    symbols::Vector{A},
    gap::A,
    miss::A,
)::Array{Float64} where {
    T<:GeneralNode,
    A<:AbstractString,
    B<:AbstractString
}
    n_nodes = length(post_order(tree))
    n_states = length(symbols)
    n_sites = size(df, 2)
    my_df = Array{Float64}(undef, n_states, n_sites, n_nodes)


    # iterate through the data frame and get the node information
    for (ind, row) in enumerate(eachrow(df))
        #data_vec = zeros(Float64, (2, n_sites))
        mn = find_by_name(tree, leave_names[ind])
        mind = mn.num
        for (ind, i) in enumerate(row)
            ent = string(i)
            index = findfirst(x -> x == ent, symbols)
            if index === nothing
                if ent == gap || ent == miss
                    my_df[:, ind, mind] .= 1.0
                else
                    throw("unknown symbol $ent, $symbols")
                end
            else
                my_df[:, ind, mind] .= 0.0
                my_df[index, ind, mind] = 1.0
            end # if
        end # for
    end # for
    my_df
end

##

function getLikelihood(
    tree::GeneralNode{Float64,Int64},
    taxa::Vector{String},
    char::Vector{String},
    r::Float64,
)
    dt = datafortree(
        hcat(char),
        taxa,
        tree,
        unique(char),
        "?",
        "."
    )
    nSymbols = length(unique(char))
    eq_freq = ones(nSymbols) ./ nSymbols
    ds = PhyloDist(
        tree,
        eq_freq,
        [1.0],
        [r],
        JC
    )
    logpdf(ds, dt)
end
##

# missing value imputation
function mvi(
    char::Vector{String},
    trees::Vector{GeneralNode},
    languages::Vector{String},
    rates::Vector{Float64},
    targetIndex::Int,
    n::Int=500,
)
    charSymbols = filter(x -> x ∉ [".", "?"], unique(char))

    if length(charSymbols) == 1
        return first(charSymbols)
    else
        posteriors = []
        for _ in 1:n
            treeIndex = rand(1:length(trees))
            ltree = deepcopy(trees[treeIndex])
            rateIndex = rand(1:length(rates))
            r = rates[rateIndex]
            ll = []
            for s in charSymbols
                char_s = copy(char)
                char_s[targetIndex] = s
                push!(
                    ll,
                    getLikelihood(
                        ltree,
                        languages,
                        char_s,
                        r,
                    )
                )
            end
            ps = exp.(ll .- logsumexp(ll))
            if !any(isnan.(ps))
                push!(
                    posteriors,
                    ps,
                )
            end
        end
        posterior = mapslices(
            mean,
            hcat(posteriors...),
            dims=2,
        ) |> vec
        return charSymbols[argmax(posterior)]
    end
end


##
trees = ParseNewick("mrbayes/posterior.tree")

rates = CSV.File("data/posteriorRates.csv", header=false) |>
        DataFrame |> Array |> vec

##
d = CSV.File("data/testAligned.tsv", types=String) |> DataFrame

languages = unique(d.TAXON)
nLanguages = length(languages)

##

cClasses = unique(d.COGID)

##

recs = []
@showprogress for cc in cClasses
    ccData = filter(x -> x.COGID == cc, d)
    ccLanguages = ccData.TAXON
    missingLanguages = setdiff(languages, ccLanguages)
    targetLanguage = first(ccData.TAXON[ccData.IPA.=="?"])
    predLanguages = @pipe ccData |>
                          filter(x -> x.IPA != "?", _) |>
                          _.TAXON |>
                          Array{String,1}(_)
    charDict = Dict{String,Vector{String}}()
    for (l, a) in zip(ccData.TAXON, ccData.A_IPA)
        charDict[l] = string.(split(a))
    end
    nChars = maximum(length.(values(charDict)))
    charDict[targetLanguage] = fill("?", nChars)
    for l in missingLanguages
        charDict[l] = fill(".", nChars)
    end
    charMtx = fill(".", nLanguages, nChars)
    for (i, l) in enumerate(languages)
        charMtx[i, :] = charDict[l]
    end
    targetIndex = indexin([targetLanguage], languages)[1]
    rec_ = []
    for i in 1:size(charMtx, 2)
        char = charMtx[:, i]
        push!(
            rec_,
            mvi(char, trees, languages, rates, targetIndex)
        )
    end
    rec = join(rec_, " ")
    push!(
        recs,
        rec
    )
end


##

rDict = Dict(zip(cClasses, recs))

##
testF = "../../../ST2022/$type/$dataset/test-$proportion.tsv"


ts = CSV.File(testF, types=String) |> DataFrame

results = deepcopy(ts)

for i in 1:size(results,1), j in 2:size(results,2)
    if !ismissing(results[i,j]) && results[i,j] == "?"
        cc = results.COGID[i]
        w = rDict[cc]
        results[i,j] = join(filter(x -> x != "-", split(w)), " ")
    else
        results[i,j] = ""
    end
end

##

resultF = ""
if type=="data"
    global resultF = "../../../ST2022/systems/jaeger-prediction/training/$dataset/result-$proportion.tsv"
else
    global resultF = "../../../ST2022/systems/jaeger-prediction/surprise/$dataset/result-$proportion.tsv"
end

CSV.write(resultF, results, delim="\t")

##
