import FLOWPy as fp

function run_ecrm2(;
        alpha_deg = 6.0,
        dt = 0.001, # approximately
        add_rotors=true,
        n=10, n_ccb=10,
        nrevs=1/5*80
    )
    # prep restart file
    rpm = 1800.0
    duration = nrevs / (rpm/60)
    nsteps = Int(ceil(duration/dt))

    vinf = 67.0 # m/s
    alpha_rad = alpha_deg * 3.1415926/180.0 # aoa in radians
    restart_vpmfile = nothing
    run_name="a$(Int(ceil(alpha_deg)))_dt$(round(dt; digits=3))_n$(n)_nr$(n_ccb)_rotors$(add_rotors)"
    vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile, add_rotors, run_name, n, n_ccb)

    # fp.update_ecrm_002(vpmdata, data_list)

    for i in 0:nsteps
        save_file = i == nsteps ? true : false
        fp.step_vpm(vpmdata, vinf, alpha_rad; save_file)
    end
end
