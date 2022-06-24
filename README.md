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
`noiserun(σreflection::Real, σphase::Real)` does the following:

1) Add noise to the ideal beam splitter parameters and compute the resulting noisy unitary where `σreflection` is the standard deviation of the Gaussian noise aded on top of the reflection coefficient parameter and `σphase` is the standard deviation of the noise injected independently for each beam splitter.
2) Compute the click patterns corresponding to each W/M state in the noisy case.
3) Iterate through all noiseless W/M click patterns (computed previously) and then, for each choice, iterate through all noisy click patterns. If there’s a match and the noiseless/noisy click patterns both result from the same input state, we call this a true positive. Otherwise, we call this a false positive. We sum the probabilities for all true positive cases (same for false positives).

The following for loops, executes this procedure many times in parallel for different noise realizations and averages the results. The following plots illustrate these averages as a function of the noise strength parameters.
