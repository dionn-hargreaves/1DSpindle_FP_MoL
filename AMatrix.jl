#
#  AMatrix.jl
#  PDE_method_of_lines
#
#  Created by Dionn Hargreaves on DD/MM/YYYY.
#
#
# Return A matrix

module AMatrix

using LinearAlgebra
using SparseArrays

@inline @views function aMatrix(Y,h, Γ,ξ,N,K,y,ω_on,ω_off)
    scale = 1
    ### Creating matrix A
    A = spzeros(4*Y+1, 8*Y+1)

    for i in 1:Y
        if i!=1 && i!=Y

            ### Up
            # J contributions, in first section of F
            A[i,i-1]    = 1/(2*h)
            A[i,i+1]    = -1/(2*h)
            # upper unbound contribution, 2nd section of F
            A[i, Y+i]   = -ω_on
            # upper bound contribution, 4th section of F
            A[i, 3*Y+i] = ω_off[i]

            ### Bp
            # J contributions, in third section of F
            A[Y+i,2*Y+i-1] = 1/(2*h)
            A[Y+i,2*Y+i+1] = -1/(2*h)
            # upper unbound contribution, 2nd section of F
            A[Y+i, Y+i]    = ω_on
            # upper bound contribution, 4th section of F
            A[Y+i, 3*Y+i]  = -ω_off[i]

            ### Um
            # J contributions, in 5th section of F
            A[2*Y+i,4*Y+i-1] = 1/(2*h)
            A[2*Y+i,4*Y+i+1] = -1/(2*h)
            # lower unbound contribution, 6th section of F
            A[2*Y+i, 5*Y+i]  = -ω_on
            # lower bound contribution, 8th section of F
            A[2*Y+i, 7*Y+i]  = ω_off[i]*scale

            ### Bm
            # J contributions, in 7th section of F
            A[3*Y+i,6*Y+i-1] = 1/(2*h)
            A[3*Y+i,6*Y+i+1] = -1/(2*h)
            # lower unbound contribution, 6th section of F
            A[3*Y+i, 5*Y+i]  = ω_on
            # lower bound contribution, 8th section of F
            A[3*Y+i, 7*Y+i]  = -ω_off[i]*scale

            ### Z
            A[4*Y+1,3*Y+i]  = (1/(ξ))*N*h*y[i]
            A[4*Y+1, 7*Y+i] = -(1/(ξ))*N*h*y[i]
        end
        if i == 1
            ### Uu
            A[1,2]      = -1/h;
            A[1, Y+1]   = - ω_on
            A[1, 3*Y+1] = ω_off[1]

            ### Bu
            A[Y+1, 2*Y+2] = -1/h
            A[Y+1, Y+1]   = ω_on
            A[Y+1, 3*Y+1] = -ω_off[1]

            ### Ul
            A[2*Y+1,4*Y+2]  = -1/h
            A[2*Y+1, 5*Y+1] = - ω_on
            A[2*Y+1, 7*Y+1] = ω_off[1]*scale

            ### Bl
            A[3*Y+1, 6*Y+2] = -1/h
            A[3*Y+1, 5*Y+1] = ω_on
            A[3*Y+1, 7*Y+1] = -ω_off[1]*scale

            ### Z
            A[4*Y+1,3*Y+1]  = -(1/(ξ))*N*h*y[1]*0.5
            A[4*Y+1, 7*Y+1] = (1/(ξ))*N*h*y[1]*0.5
        end
        if i == Y
            ### Uu
            A[Y,Y-1]      = 1/h
            A[Y,2*Y]      = -ω_on
            A[Y, 4*Y]     = last(ω_off)

            ### Bu
            A[2*Y,3*Y-1]  = 1/h
            A[2*Y,2*Y]    = ω_on
            A[2*Y, 4*Y]   = -last(ω_off)

            ### Ul
            A[3*Y,5*Y-1]  = 1/h
            A[3*Y,6*Y]    = -ω_on
            A[3*Y, 8*Y]   = last(ω_off)*scale

            ### Bl
            A[4*Y,7*Y-1]  = 1/h
            A[4*Y,6*Y]    = ω_on
            A[4*Y, 8*Y]   = -last(ω_off)*scale

            ### Z
            A[4*Y+1,4*Y]  = -(1/(ξ))*N*h*last(y)*0.5
            A[4*Y+1, 8*Y] = (1/(ξ))*N*h*last(y)*0.5
        end
    end

    A[4*Y+1, 8*Y+1] = -(K/(ξ))

    return A

end

export aMatrix

end
