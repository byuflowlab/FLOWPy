import FLOWPy as fp

function run_restart(;
        alpha_deg = 6.0,
        dt = 0.001, # approximately
        add_rotors=true,
        n=20, n_ccb=20,
        nrevs=1/5*2,
        p_per_step=5,
        restart_vpmfile="../restart/a6_dt0.001_n20_nr20_rotorstrue_pfield.300.h5",
        restart_nsteps=300
    )
    # prep restart file
    rpm = 1800.0
    duration = nrevs / (rpm/60)
    nsteps = Int(ceil(duration/dt))

    vinf = 67.0 # m/s
    alpha_rad = alpha_deg * 3.1415926/180.0 # aoa in radians
    run_name="restart_a$(Int(ceil(alpha_deg)))_dt$(round(dt; digits=3))_n$(n)_nr$(n_ccb)_rotors$(add_rotors)_pps$(p_per_step)"
    vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile, restart_nsteps, add_rotors, run_name, n, n_ccb, p_per_step)

    # fp.update_ecrm_002(vpmdata, data_list)
    println("Running $(run_name):")

    for i in 0:nsteps
        save_file = i == nsteps ? true : false
        fp.step_vpm(vpmdata, vinf, alpha_rad; save_file)
    end
end
