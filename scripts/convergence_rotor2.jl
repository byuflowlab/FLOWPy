include("convergence.jl")

# n
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=20, n_ccb=10,
    nrevs=1/5*60
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=40, n_ccb=10,
    nrevs=1/5*60
)
