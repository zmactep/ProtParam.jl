using ProtParam
using Base.Test, FactCheck

include("functions_tests.jl")
#add acidparam_tests.jl - to check if aminoacid params are saved in constructor
FactCheck.exitstatus()
