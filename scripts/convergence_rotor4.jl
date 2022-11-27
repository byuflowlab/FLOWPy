include("convergence.jl")

# nr
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=80, n_ccb=40,
    nrevs=1/5*60
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=80, n_ccb=20,
    nrevs=1/5*60
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=160, n_ccb=40,
    nrevs=1/5*60
)
