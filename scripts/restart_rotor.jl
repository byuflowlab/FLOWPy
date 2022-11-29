include("restart.jl")

run_restart(;
        alpha_deg = 6.0,
        dt = 0.001, # approximately
        add_rotors=true,
        n=20, n_ccb=20,
        nrevs=1/5*2,
        p_per_step=5,
        restart_vpmfile="restart/a6_dt0.001_n20_nr20_rotorstrue_pfield.300.h5",
        restart_nsteps=300
    )
