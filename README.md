# Installation

Clone the repo with `git clone https://github.com/jbroz11/beamsplitter_work.git`.

From within the local copy of the repo open up a Julia REPL open the package manager with `]` and input:
```
pkg> activate .
pkg> instantiate
```

From here you're ready to open up the Jupyter notebooks.

# clements_algorithm.ipynb

`clements_decomposition(U::Array{<:Complex}; p0=[0, π])`

Perfroms Clements decomposition of arbitrary N-port transfer matrix to network of 2-port beam splitters +
phase shifters as described here: [https://doi.org/10.1364/OPTICA.3.001460](https://doi.org/10.1364/OPTICA.3.001460).

`U` is an arbitrary NxN unitary that describes an N-port beam splitter network. `p0` are the initial conditions
for the solver.

This returns a tuple with:

+ The first element a diagonal matrix corresponingto a final phase shift on the 
ouputs. 
+ The second element an array of tuples `(m, n, θ, ϕ)` which describe the parameters of a 
`Tₘₙ(θ, ϕ)` matrix as defined here [https://doi.org/10.1364/OPTICA.3.001460](https://doi.org/10.1364/OPTICA.3.001460).
`Tₘₙ(θ, ϕ)` itself describes a 2-port beam splitter + phase shift acting on inputs `m, n`. The array is ordered in
accordance with the action on the input state.

Also useful might be:

`Tmn(N::Int, m::Int, θ::Real, ϕ::Real; inv::Bool=false)`

which constructs `Tₘₙ(θ, ϕ)` (`Tₘₙ⁻¹(θ, ϕ)`) when `inv=false` (`inv`=true) and assuming `n=m+1` as is always the
case for the Clements decomposition.

And:

```
reference_plot(N; colors=false, save_plot_as=nothing, figsize=5, io_labels=true, grid=false)

Plot the beam splitter geometry for a Clements decomposition of an arbitrary N×N unitary.
```

# heralded_w_state.ipynb


