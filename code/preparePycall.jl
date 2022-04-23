cd(@__DIR__)
cd("..")

println(pwd())
using Pkg

Pkg.activate("code")
Pkg.instantiate()

pycommand = read(`which python`, String)

ENV["PYTHON"] = strip(pycommand)

Pkg.add("PyCall")

Pkg.build("PyCall")

run(`pip install lingpy`)
