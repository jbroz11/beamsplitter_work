"""
    reference_plot(N; colors=false, save_plot_as=nothing, figsize=5, io_labels=true, grid=false)

Plot the beam splitter geometry for a Clements decomposition of an arbitrary N×N unitary.
"""
function reference_plot(N; colors=false, save_plot_as=nothing, figsize=5, io_labels=true, grid=false)
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    # function annotate(i, j, x, y)
    #     text = plt.annotate(
    #         L"$T_{$i$j}$", (x, y), color="red", zorder=10,
    #         bbox=Dict([("boxstyle", "rarrow,pad=0.3"), ("fc", "white"), ("ec", "r"), ("lw", 2)])
    #     )
    #     text.set_rotation(45)
    # end
    
    # construct line segments
    lines = []
    for i in 1:N
        line = [(0, -i*1.), (1/2, -i)]
        if i%2 == 0 || i == N
            if i == N
                if N%2 == 1
                    push!(line, (3/2, -i), (N+1/2, -1))
                else
                    push!(line, (N-1/2, -1))
                end
            else
                push!(line, (0, -i*1.), (1/2, -i), (i, -1/2), (N+1/2, -N+i-1))
            end
        else
            if i == 1
                push!(line, (N-1/2, -N))
            else
                push!(line, (N-i+1, -N-1/2), (N+1/2, -N+i-1))
            end
        end
        push!(line, (N+1, -N+i-1))
        push!(lines, line)
    end
    
    # setup plot configs
    mc = plt.matplotlib.collections
    if colors
        colors = [plt.matplotlib.colors.to_rgba(c) for c in plt.matplotlib.rcParams["axes.prop_cycle"].by_key()["color"]]
        lc = mc.LineCollection(lines, color=colors, linewidths=2)
    else
        lc = mc.LineCollection(lines, color="k", linewidths=2)
    end
    ax.add_collection(lc)
    plt.gca().set_aspect("equal")
    ax.set_xticklabels([])
    ax.spines["top"].set_visible(false)
    ax.spines["right"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.spines["left"].set_visible(false)
    plt.tick_params(bottom=false, left=false)
    if io_labels
        ax.set_yticks(collect(-1:-1:-N))
        ax.set_yticklabels(collect(N:-1:1))
        ax.tick_params(labelright=true)
        pad = -20 / 5 * figsize
        ax.yaxis.set_tick_params(pad=pad)
    else
        ax.set_yticklabels([])
    end
    ax.margins(0.1)
    if grid
        plt.grid(grid)
        ax.set_xticks(collect(1:N))
        plt.ylim(-1/2, -(N + 1/2))
    end
    if !isnothing(save_plot_as)
        plt.savefig(save_plot_as)
    end
end

##################################################################################################################################
# Convenience functions for computing appropriate nulling θ, ϕ for Tₘₙ(θ,ϕ) and Tₘₙ⁻¹(θ,ϕ) as defined in 
# http://dx.doi.org/10.1364/OPTICA.3.001460
##################################################################################################################################

# F[1] (F[2]) are the real (imaginary) expressions for [Tₘₙ(θ,ϕ)×U]ₖₗ
function f!(F, x, k::Int, l::Int, u1, u2)
    θ = x[1]
    ϕ = x[2]
    F[1] = real(u1) * cos(ϕ) * sin(θ) - imag(u1) * sin(ϕ) * sin(θ) + real(u2) * cos(θ)
    F[2] = real(u1) * sin(ϕ) * sin(θ) + imag(u1) * cos(ϕ) * sin(θ) + imag(u2) * cos(θ) 
end

# The Jacobian for the system of equations described by F[1], F[2] in f! as a function of θ, ϕ
function j!(J, x, k::Int, l::Int, u1, u2)
    θ = x[1]
    ϕ = x[2]
    J[1, 1] =  real(u1) * cos(ϕ) * cos(θ) - imag(u1) * sin(ϕ) * cos(θ) - real(u2) * sin(θ)
    J[1, 2] = -real(u1) * sin(ϕ) * sin(θ) - imag(u1) * cos(ϕ) * sin(θ)
    J[2, 1] =  real(u1) * sin(ϕ) * cos(θ) + imag(u1) * cos(ϕ) * cos(θ) - imag(u2) * sin(θ)
    J[2, 2] =  real(u1) * cos(ϕ) * sin(θ) - imag(u1) * sin(ϕ) * sin(θ)
end

# F[1] (F[2]) are the real (imaginary) expressions for [U×Tₘₙ⁻¹(θ,ϕ)]ₖₗ
function finv!(F, x, k::Int, l::Int, u1, u2)
    θ = x[1]
    ϕ = x[2]
    F[1] =  real(u1) * cos(ϕ) * cos(θ) + imag(u1) * sin(ϕ) * cos(θ) - real(u2) * sin(θ)
    F[2] = -real(u1) * sin(ϕ) * cos(θ) + imag(u1) * cos(ϕ) * cos(θ) - imag(u2) * sin(θ)
end

# The Jacobian for the system of equations described by F[1], F[2] in finv! as a function of θ, ϕ
function jinv!(J, x, k::Int, l::Int, u1, u2)
    θ = x[1]
    ϕ = x[2]
    J[1, 1] = -real(u1) * cos(ϕ) * sin(θ) - imag(u1) * sin(ϕ) * sin(θ) - real(u2) * cos(θ)
    J[1, 2] = -real(u1) * sin(ϕ) * cos(θ) + imag(u1) * cos(ϕ) * cos(θ)
    J[2, 1] =  real(u1) * sin(ϕ) * sin(θ) - imag(u1) * cos(ϕ) * sin(θ) - imag(u2) * cos(θ)
    J[2, 2] = -real(u1) * cos(ϕ) * cos(θ) - imag(u1) * sin(ϕ) * cos(θ)
end

"""
    clements_decomposition(U::Array{<:Complex})

Takes an arbitrary unitary and decomposes it into elementary beam splitter operations according to the 
[Clements algorithm](http://dx.doi.org/10.1364/OPTICA.3.001460). Returns a tuple with the first element
a diagonal matrix `D` describing some combination of single mode phase shifts and  a list of tuples
specificying elementary beam splitter operations `(m, n=m+1, θ, ϕ)` (see ref) such that one may construct
`U = D × .. × Tₘₙ(θ, ϕ).
"""
function clements_decomposition(U; p0=[0, π])
    N, M = size(U)
    right_multipliers = []
    left_multipliers = []
    reverse_left_multipliers = []
    @assert N == M && U * U' ≈ diagm(0 => [1 for i in 1:N]) "only implemented for unitary arrays"
    Uc = copy(U)
    # T = diagm(0 => Complex[1. for _ in 1:N])
    for i in 1:N-1, j in 1:i
        θ, ϕ = 0., 0.
        if i%2 == 1
            m = i - j + 1
            n = m + 1
            k = N - j + 1
            l = m
            # Solve for Tₘₙ⁻¹(θ,ϕ) such that [U×Tₘₙ⁻¹(θ,ϕ)]ₖₗ = 0
            if abs(Uc[k, l]) > 1e-6
                u1 = Uc[k, l]
                u2 = Uc[k, l+1]
                s = nlsolve((F, x) -> finv!(F, x, k, l, u1, u2), (J, x) -> jinv!(J, x, k, l, u1, u2), p0)
                @assert s.f_converged "Solver failed to converge when nulling U_{$k, $l}."
                θ, ϕ = s.zero
                θ = θ #% (π/2)
                ϕ = ϕ #% (2π)
            end
            push!(right_multipliers, (m, n, θ, ϕ))
            Uc = Uc * Tmn(N, m, θ, ϕ, inv=true)
        else
            m = N - i + j - 1
            n = m + 1
            k = n
            l = j
            # Solve for Tₘₙ(θ,ϕ) such that [Tₘₙ(θ,ϕ)×U]ₖₗ = 0
            if abs(Uc[k, l]) > 1e-6
                u1 = Uc[k-1, l]
                u2 = Uc[k, l]
                s = nlsolve((F, x) -> f!(F, x, k, l, u1, u2), (J, x) -> j!(J, x, k, l, u1, u2), p0)
                @assert s.f_converged "Solver failed to converge when nulling U_{$k, $l}."
                θ, ϕ = s.zero
                θ = θ #% (π/2)
                ϕ = ϕ #% (2π)
            end
            push!(left_multipliers, (m, n, θ, ϕ))
            Uc = Tmn(N, m, θ, ϕ, inv=false) * Uc
        end
    end
    for multiplier in reverse(left_multipliers)
        m, n, θ, ϕ1 = multiplier
        φ1 = atan(imag(Uc[m, m]), real(Uc[m, m]))
        φ2 = atan(imag(Uc[n, n]), real(Uc[n, n]))
        Uc[m, m] = exp(1im * (φ2 - ϕ1 + π))
        ϕ2 = φ1 - π - φ2
        push!(reverse_left_multipliers, (m, n, θ, ϕ2))
    end
    return Uc, vcat(right_multipliers, reverse_left_multipliers)
end