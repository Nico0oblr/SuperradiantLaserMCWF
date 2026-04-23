module SuperradiantLasing

using Random
using PyPlot
using Statistics
using OrdinaryDiffEq


include("PermRates.jl")
include("PermBasis.jl")
include("PopulationMC.jl")
include("CoherenceMC.jl")
include("ODE.jl")
include("MeanField.jl")
include("Spectrum.jl")
include("Comparisons.jl")
include("Util.jl")

export PopulationModel, CoherenceModel
export initial_state

export run_population_trajectory,
       simulate_population_ensemble,
       get_stationary_population_samples,
       get_stationary_population_samples_ergodic

export simulate_coherence_correlator_from_samples,
       simulate_correlator,
       simulate_correlator_ergodic

export set_up_rate_matrix,
       set_up_coherence_rate_matrix,
       steady_state_from_rate_matrix,
       coherence_initial_from_population_ss

export cavity_cooperativity,
       pump_thresholds,
       steady_state_sz,
       steady_state_alpha2,
       steady_state_s2,
       has_lasing_steady_state,
       solve_mf_superradiant_laser,
       alpha_of_u,
       s_of_u,
       sz_of_u,
       photon_number,
       dipole_coherence,
       inversion

export spectrum_from_correlator,
       spectrum_even_fft,
       spectrum_fwhm,
       fit_linewidth_vs_gamma

export compare_mc_to_meanfield, compatible_mf_params
export scan_gamma_linewidth, coherence_readout_vector, in_lasing_regime

end