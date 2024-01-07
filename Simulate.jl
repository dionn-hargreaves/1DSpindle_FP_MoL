#
#  Simulate.jl
#  PDE_method_of_lines
#
#  Created by Dionn Hargreaves
#
#

module Simulate

# Import Julia packages
using DifferentialEquations
using LinearAlgebra
using SparseArrays
using JLD2

# Import local modules
using Model
using AMatrix

@inline @views function simulate(Notes)

    # Define system parameters
    tMax       = 8000.0
    savePoints = 1.0   # Output interval
    Y          = 600   # Spatial resolution
    yMax       = 6.0
    y          = LinRange(0,yMax,Y)

    Γ          = 0.02

    γ          = 2
    ω_off      = 0.001*exp.(γ.*y)  # off rate is extension sensitive
    ω_on       = 0.003 # on rate is a constant,
    Db          = 0.008
    Du          = 0.004
    K          = 0.05
    N          = 45.0 # number of force generators PER SIDE
    ξ          = 0.625 # friction coeff
    h          = y[2]-y[1]

    # Create output files
    folderNotFound = 0
    folderCounter = 0
    folderName = "output/MoLRunData_$Notes$folderCounter"
    while folderNotFound==0
        if isdir("output/MoLRunData_$Notes$folderCounter")
            folderCounter += 1
        else
            folderName = "output/MoLRunData_$Notes$folderCounter"
            folderNotFound = 1
        end
    end

    mkpath(folderName)
    save_params = open("$folderName/run_parameters.txt","w")
    write(save_params, "Γ γ ω_off ω_on Db Du K N ξ h ymax Y tMax savePoints \n")
    write(save_params, string(Γ,", ", γ,", ", last(ω_off),", ", ω_on,", ", Db,", ", Du,", ", K,", ", N,", ", ξ,", ", h,", ", last(y), ", ", Y, ", ", tMax, ", ", savePoints ))
    close(save_params)

    # Set up ODE matrix
    A = aMatrix(Y,h, Γ,ξ,N,K,y,ω_on,ω_off)

    # Set up vectors to be reused in model
    ypt = zeros(Y)
    ymt = zeros(Y)
    F   = zeros(8*Y+1)

    # Fill parameter list
    p  = (A, y, K, Γ, ξ, N, h, Db, Du, Y, ypt, ymt, F)

    # Set initial conditions
    Up_0 = sqrt(1/(π*0.02))*exp.(-y.^2.0./0.02)
    Bp_0 = sqrt(1/(π*0.08))*exp.(-(y.-1).^2.0./0.02)
    Um_0 = sqrt(1/(π*0.02))*exp.(-y.^2.0./0.02)
    Bm_0 = sqrt(1/(π*0.08))*exp.(-(y.-1).^2.0./0.02)
    Z_0  = 10
    # Combine into initial state
    u0   = vcat(Up_0, Bp_0, Um_0, Bm_0, Z_0)

    # Set time range
    tspan = (0, tMax)

    # begin solving section
    println("let's go...")
    prob  = ODEProblem(model!,u0,tspan,p)
    @time sol = solve(prob, Tsit5(), reltol=1e-6, saveat=savePoints, maxiters = 1e12)
    #@time sol = solve(prob, lsoda()), saveat=savePoints, maxiters = 1e8)
    # Output solution to file
    JLD2.@save "$folderName/results.jld2" sol

    # Output parameters to file
    params = (p, savePoints, tMax)
    JLD2.@save "$folderName/parameters.jld2" params

    return 1

end

export simulate

end
