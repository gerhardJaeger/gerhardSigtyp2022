project = ARGS[1]
dataset, proportion = split(project,"-")


cd(project)

##

using CSV, DataFrames
using SequenceAlignment
using MCPhylo
using MCPhyloTree
using JSON
using ProgressMeter


##


SC = :ASJP

##
s = read("phmmTraining.json", String)
j = JSON.Parser.parse(s)
alphabet = string.(j["alphabet"])
lp = Matrix{Float64}(hcat(j["lp"]...))
lq = Vector{Float64}(j["lq"])
lt = convert(Matrix{Float64}, hcat(replace.(j["lt"], nothing => -Inf)...))
ph = Phmm{eltype(alphabet)}(alphabet, lt, lp, lq);

##

tr = CSV.File("data/training.tsv", types=String) |> DataFrame

ts = CSV.File("data/test.tsv", types=String) |> DataFrame;

# filter!(x -> x.TOKENS != "?", ts)

##
insertcols!(tr, 1, :ID => 1:size(tr, 1))

#insertcols!(ts, 1, :ID => 1:size(ts, 1))

##

tree = parsing_newick_string(
    read("data/upgmaTree.tre", String)
)
MCPhyloTree.initialize_tree!(tree)


languages = [x.name for x in get_leaves(tree)]

## create library


"""
Takes two gapped strings and returns the hamming distance between
them. Positions containing a gap in at least one string are ignored.
"""
function sHamming(al)
    ds = []
    for i in 1:size(al, 1)
        s1, s2 = al[i, :]
        if !ismissing(s1) && !ismissing(s2)
            push!(
                ds,
                Int(s1 != s2)
            )
        end
    end
    mean(ds)
end


##

"""
Takes a pairwise alignment (i.e. a pair of gapped strings with identical length)
as input and returns a matrix representation M as output.
The matrix M is defined as M[i,j] = 1 if x[i] is matched with y[j]
in the alignment, 0 else (where x,y are the two ungapped strings to be aligned).
"""

function algnMatrix(al)
    w1 = filter(x -> !ismissing(x), al[:, 1])
    w2 = filter(x -> !ismissing(x), al[:, 2])
    dm = zeros(Int, length(w1), length(w2))
    i, j = 1, 1
    for k in 1:size(al, 1)
        s1, s2 = al[k, :]
        if !ismissing(s1)
            if !ismissing(s2)
                dm[i, j] += 1
                i += 1
                j += 1
            else
                i += 1
            end
        else
            j += 1
        end
    end
    dm
end

##

"""
Takes a list of sequences and returns a library in the sense of the
T-Coffee algorithm. A library is a dictionary with sequence pairs
as keys and pairwise alignments in matrix format as columns.
"""

function createLibrary(words, ph)
    library = Dict()
    for w1 in words, w2 in words
        if (w2, w1) in keys(library)
            x = library[w2, w1]
            library[w1, w2] = Matrix{Int}(x[1]'), x[2]
        else
            v1 = string.(split(w1))
            v2 = string.(split(w2))
            al = viterbi(v1, v2, ph)[1]
            library[w1, w2] = algnMatrix(al), 1 - sHamming(al)
        end
    end
    library
end

##

"""
Takes a list of sequences and returns an extended library in the
sense of the T-Coffee algorithm. An extended library is a dictionary with
sequence pairs as keys and a score matrix as values.
For a pair of sequences x,y and a corresponding score matrix M,
M[i,j] is the score for aligning x[i] with y[j].
"""

function createExtendedLibrary(words, ph)
    library = createLibrary(words, ph)
    extLibrary = Dict()
    for w1 in words, w2 in words
        n, m = length.(split.([w1, w2]))
        dm = zeros(n, m)
        for w3 in words
            a1, s1 = library[w1, w3]
            a2, s2 = library[w3, w2]
            dm += (s1 + s2) * (a1 * a2)
        end
        extLibrary[w1, w2] = dm
    end
    extLibrary
end

##

"""
Returns the index of gappedString[i] in the
ungapped version thereof.
If gappedString[i] is a gap, returns -1
"""

function pos(alVector, i)
    if ismissing(alVector[i])
        return -1
    end
    return i - sum(ismissing.(alVector[1:i]))
end

##


"""
Needleman-Wunsch alignment of two aligned blocks b1 and b2,
using the scores in the extended library lib.
"""



function nwBlock(al1, al2, extLibrary)
    words1 = []
    for i in 1:size(al1, 2)
        push!(
            words1,
            join(filter(x -> !ismissing(x), al1[:, i]), " ")
        )
    end
    words2 = []
    for i in 1:size(al2, 2)
        push!(
            words2,
            join(filter(x -> !ismissing(x), al2[:, i]), " ")
        )
    end
    n, m = size(al1, 1), size(al2, 1)
    dp = zeros(n + 1, m + 1)
    pointers = zeros(Int, n + 1, m + 1)
    pointers[1, 2:end] .= 3
    pointers[2:end, 1] .= 2
    for i in 2:(n+1), j in 2:(m+1)
        insert = dp[i-1, j]
        delet = dp[i, j-1]
        match = dp[i-1, j-1]
        for k in 1:length(words1), l in 1:length(words2)
            if !ismissing(al1[i-1, k]) && !ismissing(al2[j-1, l])
                pos1 = pos(al1[:, k], i - 1)
                pos2 = pos(al2[:, l], j - 1)
                w1, w2 = words1[k], words2[l]
                match += extLibrary[w1, w2][pos1, pos2]
            end
        end
        dp[i, j] = maximum([insert, delet, match])
        pointers[i, j] = argmax(([match, insert, delet]))
    end
    i, j = size(dp)
    indices = []
    while maximum([i, j]) > 1
        p = pointers[i, j]
        if p == 1
            i -= 1
            j -= 1
            pushfirst!(indices, (i, j))
        elseif p == 2
            i -= 1
            pushfirst!(indices, (i, -1))
        else
            j -= 1
            pushfirst!(indices, (-1, j))
        end
    end
    alNew = Array{Union{String,Missing}}(
        missing,
        length(indices),
        size(al1, 2) + size(al2, 2))

    for (k, (i, j)) in enumerate(indices)
        if i == -1
            x1 = fill(missing, size(al1, 2))
        else
            x1 = al1[i, :]
        end
        if j == -1
            x2 = fill(missing, size(al2, 2))
        else
            x2 = al2[j, :]
        end
        alNew[k, :] = vcat(x1, x2)
    end
    alNew
end


##



function tCoffee(words, ccLangs, tree, ph)
    extLibrary = createExtendedLibrary(words, ph)
    alHistory = Dict()
    nums = Dict()
    for nd in post_order(tree)
        if length(nd.children) == 0
            l = nd.name
            if l ∈ ccLangs
                w = words[l.==ccLangs][1]
                alHistory[nd.num] = reshape(string.(split(w)), :, 1)
                nums[nd.num] = indexin([l], ccLangs)
            else
                alHistory[nd.num] = []
                nums[nd.num] = []
            end
        else
            ch1, ch2 = nd.children
            al1 = alHistory[ch1.num]
            al2 = alHistory[ch2.num]
            nums1 = nums[ch1.num]
            nums2 = nums[ch2.num]
            if length(al1) == 0
                alHistory[nd.num] = al2
                nums[nd.num] = nums2
            elseif length(al2) == 0
                alHistory[nd.num] = al1
                nums[nd.num] = nums1
            else
                alHistory[nd.num] = nwBlock(al1, al2, extLibrary)
                nums[nd.num] = vcat(nums1, nums2)
            end
        end
    end
    return alHistory[tree.num][:, sortperm(nums[tree.num])]
end

##

algDict = Dict()

cClasses = unique(tr.COGID)

for cc in cClasses
    ccData = filter(x -> x.COGID == cc, tr)
    ccLangs = ccData.TAXON
    words = ccData[:, SC]
    alg = tCoffee(words, ccLangs, tree, ph)
    alg[ismissing.(alg)] .= "-"
    for (i, j) in enumerate(ccData.ID)
        algDict[j] = join(alg[:, i], " ")
    end
end

insertcols!(
    tr,
    Symbol("A_" * string(SC)) => [algDict[i] for i in tr.ID]
)


##
function newPhmm(ph, words)
    newSymbols = setdiff(vcat(split.(words)...), ph.alphabet)
    newAlphabet = string.(vcat(ph.alphabet, newSymbols))
    p = zeros(length(newAlphabet), length(newAlphabet))
    p = p .+ exp(minimum(ph.lp))
    p[1:length(ph.alphabet), 1:length(ph.alphabet)] = exp.(ph.lp)
    p /= sum(p)
    new_lp = log.(p)
    q = zeros(length(newAlphabet))
    q = q .+ exp(minimum(ph.lq))
    q[1:length(alphabet)] = exp.(ph.lq)
    q /= sum(q)
    new_lq = log.(q)
    new_ph = Phmm{String}(newAlphabet,
        ph.lt, new_lp, new_lq)
end

##

ccTs = unique(ts.COGID)
algDict = Dict()
for (i,s) in zip(ts.ID, ts.TOKENS)
    if s == "?"
        algDict[i] = "?"
    end
end

##

for cc in ccTs
    ccData = filter(x -> x.COGID == cc, ts)
    filter!(x -> x.TOKENS != "?", ccData)
    words = ccData[:, SC]
    new_ph = newPhmm(ph, words)
    alg = tCoffee(words, ccData.TAXON, tree, new_ph)
    alg[ismissing.(alg)] .= "-"
    for (i, j) in enumerate(ccData.ID)
        algDict[j] = join(alg[:, i], " ")
    end
end

##

aColumn = Symbol("A_" * string(SC))

insertcols!(
    ts,
    aColumn => [algDict[i] for i in ts.ID]
)

##

function copyAlignment(v, s)
    sa = []
    for x in v
        if x ∈ ["-", "?"]
            push!(sa, x)

        else
            push!(sa, popfirst!(s))
        end
    end
    return join(sa, " ")
end

##

if SC != :TOKENS
    aTokens = []
    for (v, s) in zip(tr[:, aColumn], tr.TOKENS)
        push!(
            aTokens,
            copyAlignment(split(v), split(s))
        )
    end
    insertcols!(
        tr,
        :A_IPA => aTokens
    )
    aTokens = []
    for (v, s) in zip(ts[:, aColumn], ts.TOKENS)
        push!(
            aTokens,
            copyAlignment(split(v), split(s))
        )
    end
    insertcols!(
        ts,
        :A_IPA => aTokens
    )
end

##


if SC != :ASJP
    aTokens = []
    for (v, s) in zip(tr[:, aColumn], tr.ASJP)
        push!(
            aTokens,
            copyAlignment(split(v), split(s))
        )
    end
    insertcols!(
        tr,
        :A_ASJP => aTokens
    )
    aTokens = []
    for (v, s) in zip(ts[:, aColumn], ts.ASJP)
        push!(
            aTokens,
            copyAlignment(split(v), split(s))
        )
    end
    insertcols!(
        ts,
        :A_ASJP => aTokens
    )
end

##


if SC != :DOLGO
    aTokens = []
    for (v, s) in zip(tr[:, aColumn], tr.DOLGO)
        push!(
            aTokens,
            copyAlignment(split(v), split(s))
        )
    end
    insertcols!(
        tr,
        :A_DOLGO => aTokens
    )
    aTokens = []
    for (v, s) in zip(ts[:, aColumn], ts.DOLGO)
        push!(
            aTokens,
            copyAlignment(split(v), split(s))
        )
    end
    insertcols!(
        ts,
        :A_DOLGO => aTokens
    )
end

##





CSV.write("data/trainingAligned.tsv", tr, delim="\t")

CSV.write("data/testAligned.tsv", ts, delim="\t")