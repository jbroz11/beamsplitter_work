# found here: https://groups.google.com/g/julia-users/c/BswpTOoHdfA
# lexicographic premutations generation, By Donald Knuth
"""
    lpermutations(a::Vector)
Returns all unique permutations of `a`.
""" 
function lpermutations(a)
    b = Vector()
    sort!(a)
    n = length(a)
    while(true)
        push!(b, copy(a))    
        j = n-1
        while(a[j] >= a[j+1])
            j -= 1
            j == 0 && return(b)
        end
        l=n
        while(a[j] >= a[l])
            l -= 1
        end
        tmp = a[l]
        a[l] = a[j]
        a[j] = tmp
        k = j+1
        l = n
        while(k < l)
            tmp = a[k]
            a[k] = a[l]
            a[l] = tmp
            k += 1
            l -= 1
        end
    end
end

"""
    generate_basis_state(N::Int, j::Int, k::Int)

Generate basis column vector of length `N` with `j` excitations (assuming at most one excitation per channel).
`1 ≤ k ≤ binomial(N, j)` is a unique basis state in the `j`-excitation manifold.
"""
function generate_basis_state(N::Int, j::Int, k::Int)
    state = zeros(Int, N)
    length = binomial(N, j)
    @assert j ≤ N "j > N (we assume one excitation per input at most)"
    @assert 1 ≤ k ≤ length "k<1 or k>binomial(N, j)"
    for i in 1:j
        state[i] = 1
    end
    # return state
    return convert(Vector{ComplexF64}, lpermutations(state)[k])
end

# https://stackoverflow.com/questions/37101324/how-to-zero-out-small-values-in-an-array#answer-37101609
"""
    sparsify(x::Array, eps::Real)
Sets `x[i, j] < eps = 0` for all `i,j`. If `eps` is unspecified, defaults to `eps=1e-8`.
"""
function sparsify(x, eps) 
    if abs(x) < eps
        return 0.0
    elseif abs(real(x)) < eps
        return imag(x) * 1im
    elseif abs(imag(x)) < eps
        return real(x)
    end
    return x
end

function sparsify(x) 
    eps = 1e-8
    if abs(x) < eps
        return 0.0
    elseif abs(real(x)) < eps
        return imag(x) * 1im
    elseif abs(imag(x)) < eps
        return real(x)
    end
    return x
end

##################################################################################################################################
# Convenience functions for computing Tₘₙ(θ,ϕ) and Tₘₙ⁻¹(θ,ϕ) as defined in http://dx.doi.org/10.1364/OPTICA.3.001460 for n=m+1
##################################################################################################################################
function Tdiagonal(N, m, θ, ϕ; inv=false)
    values = ones(ComplexF64, N)
    values[m] = exp(1im * (-1)^Int(inv) * ϕ) * cos(θ)
    values[m+1] = cos(θ)
    return values
end

function Tlowerdiagonal(N, m, θ, ϕ; inv=false)
    values = zeros(ComplexF64, N-1)
    if inv
        values[m] = -sin(θ)
    else
        values[m] = exp(1im * ϕ) * sin(θ)
    end
    return values
end

function Tupperdiagonal(N, m, θ, ϕ; inv=false)
    values = zeros(ComplexF64, N-1)
    if inv
        values[m] = exp(-1im * ϕ) * sin(θ)
    else   
        values[m] = -sin(θ)
    end
    return values
end

"""
    Tmn(N::Int, m::Int, θ::Real, ϕ::Real; inv=false, , convention="clements")

Compute Tₘₙ(θ,ϕ) of dimension N×N (or Tₘₙ⁻¹(θ,ϕ) when `inv=true`) as described in [this paper](http://dx.doi.org/10.1364/OPTICA.3.001460) when `n=m+1`.
If `convention="clements"`, the `N=2` version corresponds to the convention in the reference. If `convention="asymmetric"`,
the `N=2` version corresponds to the standard asymmetric beam splitter convention with phase shift on upper output:
[exp(iϕ)cos(θ)   exp(iϕ)sin(θ)]
[sin(θ)          -cos(θ)      ]
"""
function Tmn(N::Int, m::Int, θ::Real, ϕ::Real; inv=false, convention="clements")
    if convention == "clements"
        diagm(0 => Tdiagonal(N, m, θ, ϕ, inv=inv), 
             -1 => Tlowerdiagonal(N, m, θ, ϕ, inv=inv), 
              1 => Tupperdiagonal(N, m, θ, ϕ, inv=inv)
        )
    elseif convention == "asymmetric"
         ϕ += π
        diagm(0 => Tdiagonal(N, m, θ, ϕ, inv=inv), 
             -1 => Tlowerdiagonal(N, m, θ, ϕ, inv=inv), 
              1 => Tupperdiagonal(N, m, θ, ϕ, inv=inv)
        )
    end
end

# function Tmn(N::Int, m::Int, n::Int, θ::Real, ϕ::Real; inv::Bool=false, convention="clements")
#     a = 1
#     if convention == "clements"
#         nothing
#     elseif convention == "asymmetric"
#         a = -1
#         ϕ += π
#     end
#     T = diagm(0 => a * Complex[1. for _ in 1:N])
#     T[m, m] = exp(1im * (-1)^Int(inv) * ϕ) * cos(θ)
#     T[n, n] = cos(θ)
#     if inv
#         T[m, n] = exp(-1im * ϕ) * sin(θ)
#         T[n, m] = -sin(θ)
#     else
#         T[m, n] = -sin(θ)
#         T[n, m] = exp(1im * ϕ) * sin(θ)
#     end
#     return T
# end

"""
    swap(N, m, n)
Returns NxN unitary description of operation that swaps input m with n on outputs.
"""
function swap(N, m, n)
    U = diagm(0 => Complex[i == m || i == n ? 0 : 1 for i in 1:N])
    U[m, n] = 1.
    U[n, m] = 1.
    return U
end