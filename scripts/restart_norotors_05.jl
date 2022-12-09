include("convergence.jl")

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.005, # approximately
    add_rotors=false,
    n=20, n_ccb=20,
    nrevs=1/5*60,
    p_per_step=5
)
