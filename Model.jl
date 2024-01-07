#
#  Model.jl
#  PDE_method_of_lines
#
#  Created by Dionn Hargreaves on DD/MM/YYYY.
#
#
# define IVP to solve

module Model

using SparseArrays

@inline @views function model!(du, u, p, t)

    # unpack parameters
    A, y, K, Γ, ξ, N, h, Db, Du, Y, ypt, ymt, F = p

    ypt .= 1.0 .- y .- du[4*Y+1]
    ymt .= 2.0.*(1.0 .- y) .- ypt

    # fill flux, note flux = 0 at boundaries i = 1, Y
    for i in 2:Y-1
        F[i]     = Γ*(-y[i]*u[i]-(Du/(2*h))*(u[i+1]-u[i-1]))
        F[2*Y+i] = ypt[i]*u[Y+i]-(Db/(2*h))*(u[Y+i+1]-u[Y+i-1])
        F[4*Y+i] = Γ*(-y[i]*u[2*Y+i]-(Du/(2*h))*(u[2*Y+i+1]-u[2*Y+i-1]))
        F[6*Y+i] = ymt[i]*u[3*Y+i]-(Db/(2*h))*(u[3*Y+i+1]-u[3*Y+i-1])
    end

    F[Y+1:2*Y]   .= u[1:Y]
    F[3*Y+1:4*Y] .= u[Y+1:2*Y]
    F[5*Y+1:6*Y] .= u[2*Y+1:3*Y]
    F[7*Y+1:8*Y] .= u[3*Y+1:4*Y]
    F[8*Y+1]      = u[4*Y+1]

    du .= A*F

end

export model!

end
