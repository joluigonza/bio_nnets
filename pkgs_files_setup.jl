#Pkg.instantiate()
import Pkg;

using LinearAlgebra

using IntervalArithmetic
using GenericLinearAlgebra
using GenericSchur

using RCall

using Statistics
using Random, Distributions

using PyPlot

using TickTock

using Serialization

using IntervalUnionArithmetic
#################################################

include("int_nnets_functs.jl")

include("act_nnets_functs.jl")

R"""
source(file="functions.R")

"""

include("nnet_setup.jl")

#########################################################