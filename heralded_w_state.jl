"""
    commonelements(A, B)

Given two arrays `A` and `B` returns `true` if they have at least one element in common (false, otherwise).
"""
function commonelements(A, B)
    for c1 in A, c2 in B
        c1 == c2 && return true
    end
    return false
end

"""
    W(x...)

`x` is a vector of phases. Returns a generalized W-state of dimension length(x) with phases as given by elements of `x`\n
eg. x = 1, 1, -1 returns (|001⟩ + |010⟩ - |100⟩)/√3
"""
function W(x...) 
    N = length(x)
    lp = lpermutations(vcat([0 for _ in 1:N-1], [1]))
    return [(x[i] / √N, lp[i]) for i in 1:length(x)]
end

"""
    M(x...)

Same as `W(x...)` but with 1 ⟷ 0
"""
function M(x...) 
    N = length(x)
    lp = lpermutations(vcat([1 for _ in 1:N-1], [0]))
    return [(x[i] / √N, lp[i]) for i in 1:length(x)]
end

"""
    GHZ(x; N=2)

Returns (|00...0⟩ + x|11...1⟩)/√2 with dimension `N`.
"""
GHZ(x; N=2) = [(1/√2, [0 for _ in 1:N]), ((-1)^((x%2)+1)/√2, [1 for _ in 1:N])]

"""
    K(x)

Returns the 6 permutations of the product state |0011⟩.
"""
K(x) = [(1, lpermutations([0, 0, 1, 1])[x])]

function statetostring(state)
    s = ""
    for (i, el) in enumerate(state)
        el == 0 && continue
        if i%2 == 1
            s *= "H$(i÷2 + 1)_"
        else
            s *= "V$(i÷2)_"
        end
    end
        s = strip(s, '_')
    return join(sort(split(s, "_")), " ")
end

"""
    singlephotonclicks(U::Array{Complex}, input_state)

We assume the input state corresponds to exactly one photon incident at each input port, which can either be
horizontally (0) or vertically (1) polarized, though this input state will generally be a superposition.

For a 2-photon GHZ state (|HH⟩ - |VV⟩)√2, we would set `input_state = [(1, [0, 0]), (-1, [1, 1])]`
Each element of the vector describes one branch of a superposition state where the first element of the tuple
is the complex amplitude of the branch and the second element is a vector describing the configuration of H, V.

U is the unitary that will act on the input state to yield the output, which is an array with elements of the form: 
`["H2 H3 H4 V1", 0.0625, [0, 1, 1, 0, 1, 0, 1, 0]]`
where the first element is the click pattern (we only keep track of click patterns where there is at most one photon
per channel (mode + polarization), the second element is the probability of this click pattern and the last element
is vector of the form [H1, V1, H2, V2, ..., HN, VN].

"""
function singlephotonclicks(U, input_state)
    output_operator = 0
    N = length(input_state[1][2])
    @syms aH[1:N] aV[1:N]
    
    for branch in input_state
        amplitude, product_state = branch
        branch_output_operator = 1
        
        for (i, polarization) in enumerate(product_state)
            # polarization == 0 means H-polarized light
            # polarization == 1 means V-polarized light
            product_state_output_operator = 0
            
            for (j, U_ji) in enumerate(view(U, :, i))
                if polarization == 0
                    product_state_output_operator += U_ji * aH[j]
                else
                    product_state_output_operator += U_ji * aV[j]
                end
            end
            
            branch_output_operator *= product_state_output_operator
        end
        
        output_operator += amplitude * branch_output_operator
    end
    
    output_operator = expand(output_operator)
    permutations = lpermutations(vcat([0 for _ in 1:N], [1 for _ in 1:N]))
    coincidence_counts = []
    #=
    This can be done much more efficiently, the number of coefficients scales as N^N presently
    but we only really care about the single photon clicks 
    =#
    for permutation in permutations
        operator = 1
        for i in 1:N
            if permutation[2i - 1] == 1
                operator *= aH[i]
            end
            if permutation[2i] == 1
               operator *= aV[i] 
            end
        end
        if haskey(output_operator.dict, operator)
            prob = abs(output_operator.dict[operator])^2
            if prob > 1e-8
                push!(coincidence_counts, [statetostring(permutation), prob, permutation])
            end
        end
    end
    return coincidence_counts
end

"""
    actual_output_state(U, input_state)

Computes the full quantum output state after the beam splitter network described by `U`. `input_state` should be
in the same form as in `singlephotonclicks`.
"""
function full_output_state(U, input_state)
    N = length(input_state[1][2])
    b = FockBasis(N);
    I = identityoperator(b)
    a = [reduce(⊗, [i == j ? create(b) : I for i in 1:2N]) for j in 1:2N]
    ground_state = reduce(⊗, [fockstate(b, 0) for _ in 1:2N])
    output_state = ground_state * 0
    zero_array = a[1] * 0
    id = identityoperator(a[1].basis_l)
    
    for branch in input_state
        amplitude, product_state = branch
        total_operator = copy(id)
        
        for (i, polarization) in enumerate(product_state)
            operator = copy(zero_array)
            
            for (j, U_ji) in enumerate(view(U, :, i))
                operator += U_ji * a[2j - 1 + polarization]
            end
            
            total_operator *= operator
        end
        
        output_state += amplitude * total_operator * ground_state
    end
    
    return output_state
end

"""
    clicksfromfullstate(state)

Returns the same thing as `singlephotonclicks` but takes as input the full quantum state.
"""
function clicksfromfullstate(state)
    N = length(state.basis.bases) ÷ 2
    b = FockBasis(N)
    results = []
    permutations = lpermutations(vcat([0 for _ in 1:N], [1 for _ in 1:N]))
    for permutation in permutations
        mstate = mapreduce(x -> x == 0 ? fockstate(b, 0) : fockstate(b, 1), ⊗, permutation)
        prob = round(abs(mstate' * state)^2, digits=4)
        prob > 1e-6 && push!(results, [statetostring(permutation), prob, permutation])
    end
    return results
end