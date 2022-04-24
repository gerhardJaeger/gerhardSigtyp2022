
##
flag = string(strip(ARGS[1], '/'))
##
_, type, dataset, fn = split(flag, "/")

proportion = replace(fn, r"[^-]*-" => "", ".tsv" => "")

project = "tmp/$type/$dataset-$proportion"


cd(project)


##

using CSV, DataFrames, MCPhyloTree, Pipe
using Distributions, Distributed
using MCPhylo
using DelimitedFiles
using Random
Random.seed!(12345)


##

import MCPhylo: SamplerTune, SymDistributionType, Sampler, ElementOrVector
import MCPhylo: sample!

include("../code/prior_sampler.jl")

##

nchains = 2

##

function charCode(x::Vector{String})::Vector{String}
    xSymbols = [s for s in x if s != "."] |> unique
    symCode = Dict{String,String}()
    symCode["."] = "."
    for (i, s) in enumerate(xSymbols)
        symCode[s] = string(i)
    end
    [symCode[s] for s in x]
end

##



import MCPhylo: datafortree


function datafortree(
    df::Matrix{A},
    leave_names::Vector{A},
    tree::T,
    symbols::Vector{A},
    gap::A,
    miss::A,
)::Array{Float64} where {T<:GeneralNode,A<:AbstractString}
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





trees = ParseNewick("mrbayes/posterior.tree")


##

tra = CSV.File("data/trainingAligned.tsv", types=String) |> DataFrame

cClasses = unique(tra.COGID)

doculects = convert(Vector{String},unique(tra.TAXON))
nDoculects = length(doculects)
##
chMtx_ = []
for cc in cClasses
    ccData = filter(x -> x.COGID == cc, tra)
    doc2alg = Dict()
    for rw in eachrow(ccData)
        doc2alg[rw.TAXON] = rw.A_IPA
    end
    nColumns = length(split(first(values(doc2alg))))
    ccMtx = fill(".", nDoculects, nColumns)
    for (i,l) in enumerate(doculects)
        if l in ccData.TAXON
            ccMtx[i,:] = split(doc2alg[l])
        end
    end
    push!(
        chMtx_,
        ccMtx
    )
end

chMtx = hcat(chMtx_...)

##


#charMatrix = convert.(String, Array(characters[:, 2:end]))
#languages = convert.(String, characters.DOCULECT)
charMtxCoded = hcat(
    [charCode(chMtx[:, i]) for i in 1:size(chMtx, 2)]...
)

##


totalSymbols = filter(x -> x ∉ [".", "-"], charMtxCoded) |> unique
nSymbols = length(totalSymbols)
data = datafortree(
    charMtxCoded,
    doculects,
    trees[1],
    totalSymbols,
    ".",
    ".",
);
eq_freq = ones(nSymbols) ./ nSymbols
ntrees = length(trees)

##

model = Model(
    data=Stochastic(3,
        (idx, rates) -> PhyloDist(
            trees[Int(ceil(idx))],
            eq_freq, [1.0], rates, JC),
        false,
    ),
    idx=Stochastic(
        () -> Uniform(0, ntrees),
        false,
    ),
    rates=Logical(1,
        (a) -> discrete_gamma_rates(a, a, 4),
        false,
    ),
    a=Stochastic(
        () -> Exponential(),
        true,
    )
)
scheme = [Slice([:a], 1.0), Prior([:idx])]

##

setsamplers!(model, scheme);
data_dictionary = Dict{Symbol,Any}(
    :data => data,
);
inits = [
    Dict{Symbol,Union{Any,Real}}(
        :idx => rand(Uniform(0, ntrees)),
        :data => data_dictionary[:data],
        :a => rand(Exponential()),
    ) for i in 1:nchains
];

##
sim = mcmc(
    model,
    data_dictionary,
    inits,
    1000,
    burnin=0,
    thin=10,
    chains=nchains,
    trees=false,
    verbose=true,
)

##

bi = size(sim,1) ÷ 2

gd = gelmandiag(sim[(bi+1):end, :a, :]).value[1,1,1]

while gd > 1.1
    global sim = mcmc(sim, 1000)
    global bi = size(sim,1) ÷ 2
    global gd = gelmandiag(sim[(bi+1):end, :a, :]).value[1,1,1]
end

# open("data/gelmandiag.tsv", "w") do file
#     show(file, gelmandiag(sim[:, :a, :]))
# end

##

open("data/posteriorRates.csv", "w") do file
    writedlm(file, vec(sim[(bi+1):end, :a, :].value))
end
