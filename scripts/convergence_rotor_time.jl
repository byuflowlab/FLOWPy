include("convergence.jl")

run_restart(;
        alpha_deg = 6.0,
        dt = 0.00025, # approximately
        add_rotors=true,
        n=20, n_ccb=10,
        nrevs=1/5*60,
        p_per_step=1,
        restart_vpmfile="restart/a6_dt0.001_n20_nr20_rotorstrue_pfield.300.h5",
        restart_nsteps=300
    )

run_restart(;
        alpha_deg = 6.0,
        dt = 0.0001, # approximately
        add_rotors=true,
        n=20, n_ccb=10,
        nrevs=1/5*60,
        p_per_step=1,
        restart_vpmfile="restart/a6_dt0.001_n20_nr20_rotorstrue_pfield.300.h5",
        restart_nsteps=300
    )
