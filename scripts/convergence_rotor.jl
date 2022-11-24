include("convergence.jl")

# base case
run_ecrm2(;
        alpha_deg = 6.0,
        dt = 0.001, # approximately
        add_rotors=true,
        n=40, n_ccb=10,
        nrevs=1/5*80
    )

# dt
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.0005, # approximately
    add_rotors=true,
    n=40, n_ccb=10,
    nrevs=1/5*80
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.00025, # approximately
    add_rotors=true,
    n=40, n_ccb=10,
    nrevs=1/5*80
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.002, # approximately
    add_rotors=true,
    n=40, n_ccb=10,
    nrevs=1/5*80
)

# n
run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=20, n_ccb=10,
    nrevs=1/5*80
)

run_ecrm2(;
    alpha_deg = 6.0,
    dt = 0.001, # approximately
    add_rotors=true,
    n=80, n_ccb=10,
    nrevs=1/5*80
)
