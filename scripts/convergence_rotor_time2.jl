include("convergence.jl")

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.0005, # approximately
    add_rotors=true,
    n=20, n_ccb=10,
    nrevs=1/5*60,
    p_per_step=1
)
