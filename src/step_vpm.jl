function step_vpm(vpmdata, vinf, alpha_rad; save_file=false)
    # unpack simulation struct
    nsteps = vpmdata.nsteps
    simulation = vpmdata.simulation
    pfield = vpmdata.pfield
    static_particles_function = vpmdata.static_particles_function
    runtime_function = vpmdata.runtime_function
    run_name = vpmdata.run_name
    save_path = vpmdata.save_path
    nsteps_save = vpmdata.nsteps_save
    save_time = vpmdata.save_time
    v_lvl = vpmdata.v_lvl
    verbose_nsteps = vpmdata.verbose_nsteps
    data_history = vpmdata.data_history
    
    i_step = simulation.nt
    dt = simulation.ttot / nsteps

    # update Vinf
    Uinf(t) = vinf .* [cos(alpha_rad), 0, sin(alpha_rad)]
    pfield.Uinf = Uinf

    if i_step%verbose_nsteps==0
        println("\t"^(v_lvl+1),"Time step $i_step out of $nsteps\tParticles: $(vpm.get_np(pfield))")
    end

    # Relaxation step
    relax = pfield.relaxation != vpm.relaxation_none &&
            pfield.relaxation.nsteps_relax >= 1 &&
            i_step>0 && (i_step%pfield.relaxation.nsteps_relax == 0)

    org_np = vpm.get_np(pfield)

    # Time step
    if i_step>-1
        # Add static particles
        remove = static_particles_function(pfield, pfield.t, dt)

        # Step in time solving governing equations
        vpm.nextstep(pfield, dt; relax=relax)

        # Remove static particles (assumes particles remained sorted)
        if remove==nothing || remove
            for pi in vpm.get_np(pfield):-1:(org_np+1)
                vpm.remove_particle(pfield, pi)
            end
        end
    end

    # Calls user-defined runtime function
    breakflag = runtime_function(pfield, pfield.t, dt;
                                vprintln= (str)-> i_step%verbose_nsteps==0 ?
                                        vprintln(str, v_lvl+2) : nothing)

    # Save particle field
    if save_path!=nothing && (i_step%nsteps_save==0 || i_step==nsteps || breakflag)
        overwrite_time = save_time ? nothing : pfield.nt
        vpm.save(pfield, run_name*"_pfield"; path=save_path, add_num=true,
                                    overwrite_time=overwrite_time)
    end

    # save internal variables
    if pfield.nt > 2
        data_history["cl"][:,i_step-1] .= simulation.vehicle.system.wings[1].sol["cl"]
        data_history["cd"][:,i_step-1] .= simulation.vehicle.system.wings[1].sol["cd"]
        data_history["CL"][i_step-1] = simulation.vehicle.system.wings[1].sol["CLwing"]
        data_history["CD"][i_step-1] = simulation.vehicle.system.wings[1].sol["CDwing"]
    end

    # save to file
    if save_file
        println("Saving file...")
        write_history_h5(run_name*"_history.h5", vpmdata)
    end
end