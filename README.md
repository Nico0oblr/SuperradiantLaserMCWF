# SuperradiantLasing

Julia tools for simulating permutation-symmetric open spin systems in the Dicke basis, with complementary stochastic and deterministic workflows for populations, coherence-sector correlation functions, and spectra.

The core idea is to reduce the many-body dynamics to structured motion on Dicke sectors `(S, M)`. In this representation, the code supports:

- population-sector Monte Carlo on the Dicke lattice
- deterministic rate-matrix evolution for the same population dynamics
- coherence-sector dynamics for first-order correlators and spectra
- comparison utilities against a separate bad-cavity mean-field model

## Repository structure

### Main entry point

- `SuperradiantLasing.jl`  
  Top-level module that includes all submodules and re-exports the main public interface.

### Core modules

- `PermRates.jl`  
  Dicke-basis branching factors and channel-resolved rate helpers for collective decay/pump, local decay/pump, and local dephasing.

- `PermBasis.jl`  
  Basis helpers for the Dicke-sector representation.

- `PopulationMC.jl`  
  Gillespie-type Monte Carlo simulation of population dynamics on Dicke states `(S, M)`.

- `CoherenceMC.jl`  
  Reduced coherence-sector dynamics for objects labeled by `(S, M) ↔ |S, M-1⟩⟨S, M|`, used to compute first-order correlators and spectra.

- `ODE.jl`  
  Sparse rate-matrix construction and deterministic propagation for both population and coherence sectors.

- `Spectrum.jl`  
  Frequency-domain post-processing from time-domain correlators, including FFT-based spectra and linewidth extraction.

### Comparison / analysis modules

- `MeanField.jl`  
  Bad-cavity mean-field model for cavity-field, dipole, and inversion dynamics.

- `Comparisons.jl`  
  Helpers for comparing the Dicke-basis simulations against the mean-field model.

### Example / validation scripts

- `compare_population_mc_vs_ode.jl`
- `compare_coherence_mc_vs_ode.jl`
- `compare_mc_vs_meanfield.jl`
- `check_stationary_sampling.jl`

### Notebooks

- `CavityMonteCarlo.ipynb`
- `SuperradiantVsMeanField2.ipynb`

## Conceptual structure

The code naturally splits into three layers.

### 1. Population sector

The many-body open-system problem is reduced to a classical jump process over Dicke states `(S, M)`. The allowed jumps are generated channel by channel:

- collective decay / pump: `(S, M) -> (S, M±1)`
- local decay / pump: `(S, M) -> (S±1, M∓1)` and `(S, M∓1)`
- local dephasing: `(S, M) -> (S±1, M)`

This sector can be simulated in two complementary ways:

- **stochastically**, via a Gillespie-type trajectory simulation
- **deterministically**, via a sparse rate matrix

The stochastic route is convenient for trajectories and large systems. The deterministic route is convenient for validation, steady states, and direct propagation of observables.

### 2. Coherence sector

For first-order correlators, the code introduces a reduced coherence sector labeled by `(S, M)` with the convention

`(S, M) ↔ |S, M-1⟩⟨S, M|`.

A population-sector branch rate `r(S, M)` and its neighboring rate `r(S, M-1)` are converted into:

- an internal coherence-sector rate `sqrt(r(S, M) r(S, M-1))`
- a diagonal loss `(r(S, M) + r(S, M-1)) / 2`
- a sink term equal to the remaining probability flow

This gives an efficient reduced description for time-domain correlators and spectra.

### 3. Mean-field comparison layer

The cavity mean-field model evolves a complex cavity field `α`, dipole `s`, and inversion `s^z` using an ODE system. It is conceptually separate from the Dicke-basis machinery and is best viewed as a comparison model.

## Dependencies

From the current files, the code uses at least:

- `Random`
- `Statistics`
- `PyPlot`
- `OrdinaryDiffEq`
- `LinearAlgebra`

There is currently no `Project.toml` in the uploaded files, so dependency installation appears to be manual.

## Loading the code

At the moment the project is structured as an include-based Julia module rather than a registered package. A typical workflow is:

```julia
include("SuperradiantLasing.jl")
using .SuperradiantLasing
```

The example scripts use this pattern in commented form.

## Quick start

### Population-sector Monte Carlo

```julia
include("SuperradiantLasing.jl")
using .SuperradiantLasing, Random

obs_exc(S, M, model) = M + model.N / 2

model = PopulationModel(
    200,      # N
    100,      # Jmax = N ÷ 2
    1.0,      # global_decay
    0.0,      # global_pump
    0.0,      # local_decay
    0.3,      # local_pump
    0.0,      # local_dephasing
)

t, avg_exc = simulate_population_ensemble(
    model,
    obs_exc,
    2000;
    t_max = 5.0,
    n_grid = 1000,
    S0 = model.Jmax,
    M0 = -model.Jmax,
    rng = Xoshiro(1234),
)
```

### Deterministic population evolution

```julia
using OrdinaryDiffEq

params = Dict(
    "global_decay" => model.global_decay,
    "global_pump" => model.global_pump,
    "local_pump" => model.local_pump,
    "local_decay" => model.local_decay,
    "local_dephasing" => model.local_dephasing,
)

basis, idx, R, nstates = set_up_rate_matrix(model.Jmax, params)
u0 = initial_state(idx, (model.Jmax, -model.Jmax), nstates)

f!(du, u, p, t) = mul!(du, R, u)
prob = ODEProblem(f!, u0, (0.0, 5.0))
sol = solve(prob, Tsit5())
```

### Stationary sampling

There are two stationary population samplers:

- `get_stationary_population_samples`: many independent runs up to `t_ss`
- `get_stationary_population_samples_ergodic`: one long run sampled at random times after burn-in

The ergodic version is often more efficient for steady-state observables when the process is ergodic.

```julia
S_samples, M_samples = get_stationary_population_samples_ergodic(
    model,
    10000;
    t_ss = 20.0,
    t_sample = 20.0,
    S0 = model.Jmax,
    M0 = model.Jmax,
    rng = Xoshiro(1234),
)
```

### Correlators and spectra

A typical coherence-sector workflow is:

1. obtain stationary population samples
2. convert the population model into a coherence model
3. simulate the coherence-sector correlator
4. Fourier transform to obtain the spectrum

```julia
coh_model = CoherenceModel(model)

τ_grid, Cτ = simulate_coherence_correlator_from_samples(
    S_samples,
    M_samples,
    coh_model;
    τ_max = 10.0,
    n_grid = 1200,
    rng = Xoshiro(1234),
)

ω_grid, Sω = spectrum_even_fft(τ_grid, Cτ)
```

## Main exported interface

### Models

- `PopulationModel`
- `CoherenceModel`

### Population simulation

- `initial_state`
- `run_population_trajectory`
- `simulate_population_ensemble`
- `get_stationary_population_samples`
- `get_stationary_population_samples_ergodic`

### Coherence / correlators

- `simulate_coherence_correlator_from_samples`
- `simulate_correlator`
- `simulate_correlator_ergodic`

### Deterministic generators

- `set_up_rate_matrix`
- `set_up_coherence_rate_matrix`
- `steady_state_from_rate_matrix`
- `coherence_initial_from_population_ss`

### Spectrum utilities

- `spectrum_from_correlator`
- `spectrum_even_fft`
- `spectrum_fwhm`
- `fit_linewidth_vs_gamma`

### Mean-field comparison

- `compare_mc_to_meanfield`
- `compatible_mf_params`
- the helper functions exported from `MeanField.jl`

## Suggested workflow

For most problems, a natural path is:

1. define a `PopulationModel`
2. validate transient observables with population MC and/or rate-matrix evolution
3. compute stationary populations using either independent runs or ergodic stationary sampling
4. initialize the coherence sector from the stationary population state
5. propagate coherence dynamics to obtain `C(τ)`
6. Fourier transform `C(τ)` to obtain spectra and linewidths

This keeps the population and coherence sectors conceptually separate while allowing direct cross-checks between stochastic and deterministic calculations.

## Notes

- The example scripts are currently the best usage documentation.
- The code mixes reusable module code with analysis scripts and notebooks in a lightweight way.
- If the project grows further, it would likely benefit from a `Project.toml`, clearer examples, and a small `examples/` or `scripts/` directory structure.

## One-sentence summary

This codebase implements Dicke-basis simulation tools for superradiant open-spin systems, combining stochastic and deterministic population dynamics, reduced coherence-sector correlators, and spectrum extraction.
