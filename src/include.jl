using JuMP
using AmplNLWriter, Bonmin_jll, Couenne_jll
using Ipopt
using Plots
using MathOptInterface

const MOI = MathOptInterface

include("types.jl")
include("utils.jl")
include("input.jl")
include("model.jl")
include("enumeration.jl")
include("output.jl")

const density = 1000
const g = 9.81
const deepest = 200.0