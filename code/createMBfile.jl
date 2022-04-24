flag = string(strip(ARGS[1], '/'))

_, type, dataset, fn = split(flag, "/")

proportion = replace(fn, r"[^-]*-" => "", ".tsv" => "")

project = "tmp/$type/$dataset-$proportion"


cd(project)

##

using CSV, DataFrames
using Pipe

##

tra = CSV.File("data/trainingAligned.tsv", types=String) |> DataFrame;

##

try
    mkdir("mrbayes")
catch e
end

##

doculects = unique(tra.TAXON)

nDoculects = length(doculects)

soundClasses = unique(vcat(split.(tra.A_DOLGO)...))

##

soundFreqs = @pipe tra |>
                   select(_, :A_DOLGO) |>
                   Array |>
                   split.(_) |>
                   vcat(_...) |>
                   DataFrame(sca=_) |>
                   filter(x -> x.sca != "-", _) |>
                   groupby(_, :sca) |>
                   combine(_, nrow => :freq) |>
                   sort(_, :freq, rev=true)
##

sc2char = Dict()
for (i, s) in enumerate(soundFreqs.sca)
    if i <= 10
        sc2char[s] = string(i - 1)
    else
        sc2char[s] = "-"
    end
end
sc2char["-"] = "-"

##

cClasses = unique(tra.COGID)

charMtx_ = []

for cc in cClasses
    ccData = filter(x -> x.COGID == cc, tra)
    doc2a = Dict(zip(ccData.TAXON, ccData.:A_DOLGO))
    nColumns = length(split(first(values(doc2a))))
    ccMtx = fill("-", nDoculects, nColumns)
    for (i, l) in enumerate(doculects)
        if l ∈ ccData.TAXON
            chars = split(doc2a[l])
            for (j, s) in enumerate(chars)
                ccMtx[i, j] = sc2char[s]
            end
        end
    end
    push!(
        charMtx_,
        ccMtx
    )
end

charMtx = hcat(charMtx_...)

##

pad = maximum(length.(doculects)) + 5

nex = """#Nexus
BEGIN DATA;
DIMENSIONS ntax=$(length(doculects)) nchar = $(size(charMtx,2));
FORMAT DATATYPE=Standard GAP=? MISSING=- interleave=yes;
MATRIX

"""

for (i, l) in enumerate(doculects)
    global nex *= rpad(l, pad)
    nex *= join(charMtx[i, :]) * "\n"
end

nex *= ";\nEND"

open("mrbayes/data.nex", "w") do file
    write(file, nex)
end

##

mb = """
#Nexus
Begin MrBayes;
    set seed=34567;
    execute $project/mrbayes/data.nex;
    lset rates=gamma coding=all;
    prset brlenspr = clock:uniform;
    mcmcp stoprule=yes stopval=0.01 filename=$project/mrbayes/posterior samplefreq=1000;
    mcmc ngen=100000000;
    sumt;
    q;
end;
"""

##

open("mrbayes/mb.nex", "w") do file
    write(file, mb)
end
