include("convergence.jl")

# base case
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.0005, # approximately
    add_rotors=true,
    n=80, n_ccb=10,
    nrevs=1/5*80
)

# n_ccb
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.0005, # approximately
    add_rotors=true,
    n=80, n_ccb=20,
    nrevs=1/5*80
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.0005, # approximately
    add_rotors=true,
    n=80, n_ccb=40,
    nrevs=1/5*80
)
