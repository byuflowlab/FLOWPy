include("convergence.jl")

# base case
run_ecrm2(;
        alpha_deg = 6.0,
        dt = 0.0005, # approximately
        add_rotors=false,
        n=20, n_ccb=20,
        nrevs=1/5*60
    )
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.002, # approximately
    add_rotors=true,
    n=20, n_ccb=20,
    nrevs=1/5*60
)
