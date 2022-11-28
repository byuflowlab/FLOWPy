include("convergence.jl")

# p per step
run_ecrm2(;
        alpha_deg = 6.0,
        dt = 0.001, # approximately
        add_rotors=true,
        n=20, n_ccb=20,
        nrevs=1/5*45,
        p_per_step=2
    )
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=20, n_ccb=20,
    nrevs=1/5*45,
    p_per_step=3
)
