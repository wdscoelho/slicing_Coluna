##############################################################################
# Created by: Wesley DA SILVA COELHO
##############################################################################
#This file contains the main script for NSDP model execution
#: https://hal.archives-ouvertes.fr/hal-02448028v3


##############################################################################
#                   INCLUDE PACKAGES
##############################################################################
using NBInclude
using JuMP,CPLEX,MathProgBase,GLPK,BlockDecomposition, Coluna
using LightGraphs, MetaGraphs,LightGraphsFlows  

# path to some folders
input_folder = ""
include("my_functions.jl") # our data structures 

instance_number = 1
instance = get_Instance(input_folder,instance_number)  

optimal_solution_value = create_my_model(instance)
best_Coluna_solution_value = create_COLUNA(instance)

println("optimal_solution_value = $(optimal_solution_value)")
println("best_Coluna_solution_value = $(best_Coluna_solution_value)")
