#
#  main.jl
#  PDE_method_of_lines
#
#  Created by Dionn Hargreaves
#
#

# Load simulate function
"./" in LOAD_PATH ? nothing : push!(LOAD_PATH,"./")
using Simulate

# Notes so file is identifiable - parameters will be saved as a text file
Notes = "Test_nondim_FIG7a"

# Run simulation
simulate(Notes)
