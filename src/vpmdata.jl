struct VPMData{TF}
    nsteps::Int64
    simulation::uns.Simulation
    pfield::vpm.ParticleField
    static_particles_function::Function
    runtime_function::Function
    S::TF
    b::TF
    y2b::Vector{TF}
    run_name::String
    save_path::Union{Nothing,String}
    nsteps_save::Int64
    save_time::Bool
    v_lvl::Int64
    verbose_nsteps::Int64
    data_history::Dict
end

function VPMData(nsteps, duration; # Note: dt = duration/nsteps
        # data_path="/home/ryan/.julia/dev/FLOWUnsteady/data/", # contains two subdirectories:
        #                                                        # rotors/ and airfoils/
        data_path=joinpath(dirname(pathof(uns)), "..", "data/"),
        # SIMULATION OPTIONS
        Vref = 67.0, # initial value
        rpm = 1800.0, # initial value
        sound_spd=343,             # (m/s) speed of sound
        add_rotors=true,
        VehicleType=uns.VLMVehicle,
        n=20, n_ccb=20,
        rho=1.225,                 # (kg/m^3) air density
        # SOLVERS OPTIONS
        p_per_step=5,              # Particle sheds per time step
        vpm_formulation=vpm.formulation_rVPM, # VPM formulation
        vpm_kernel=vpm.gaussianerf,# VPM kernel
        vpm_UJ=vpm.UJ_fmm,         # VPM particle-to-particle interaction calculation
        vpm_integration=vpm.rungekutta3, # VPM time integration scheme
        vpm_transposed=true,       # VPM transposed stretching scheme
        vpm_viscous=vpm.Inviscid(),# VPM viscous diffusion scheme
        vpm_fmm=vpm.FMM(; p=4, ncrit=50, theta=0.4, phi=0.5), # VPM's FMM options
        vpm_relaxation=vpm.pedrizzetti, # VPM relaxation scheme
        vpm_surface=true,          # Whether to include surfaces in the VPM
        vlm_vortexsheet=false,     # Whether to spread the surface circulation of the VLM as a vortex sheet in the VPM
        vlm_vortexsheet_overlap=2.125, # Overlap of particles that make the vortex sheet
        vlm_vortexsheet_distribution=uns.g_pressure, # Circulation distribution of vortex sheet
        vlm_rlx=0.2,                # VLM relaxation
        vlm_init=false,            # Initialize the first step with the VLM semi-infinite wake solution
        hubtiploss_correction=vlm.hubtiploss_nocorrection, # Hub and tip loss correction of rotors (ignored by quasi-steady solver)
        wake_coupled=true,         # Couple VPM wake on VLM solution
        shed_unsteady=true,        # Whether to shed unsteady-loading wake
        unsteady_shedcrit=0.001,    # Criterion for unsteady-loading shedding
        shed_starting=false,       # Whether to shed starting vortex (only with shed_unsteady==true)
        shed_boundarylayer=false,  # Whether to shed vorticity from boundary layer of surfaces
        boundarylayer_prescribedCd=0.1, # Prescribed Cd for boundary layer shedding used for wings
        boundarylayer_d=0.0,       # Dipole width for boundary layer shedding
        omit_shedding=[],          # Indices of elements in `sim.vehicle.wake_system` on which omit shedding VPM particles
        extra_runtime_function=(PFIELD,T,DT; optargs...)->false,
        # REGULARIZATION OPTIONS
        sigma_vlm_solver=-1,       # Regularization of VLM solver (internal)
        sigmafactor_vpm=1.0,       # Core overlap of wake particles
        sigmafactor_vpmonvlm=1,    # Shrinks the particles by this factor when calculating VPM-on-VLM/Rotor induced velocities
        # RESTART OPTIONS
        restart_vpmfile=nothing,   # VPM particle field restart file
        restart_dt=0.001,
        restart_n=20,
        restart_nr=20,
        restart_pps=5,
        # OUTPUT OPTIONS
        run_name="FLOWUnsteadyPy",
        verbose=true, v_lvl=0, verbose_nsteps=10,
        debug=false,               # Output extra states for debugging
        nsteps_save=1,             # Save vtks every this many steps
        save_horseshoes=false,     # Save VLM horseshoes
        save_static_particles=false,# Whether to save particles to represent the VLM
        save_wopwopin=false,        # Generate inputs for PSU-WOPWOP
        L_dir=[0,0,1],      # Direction of lift component
        D_dir=[-1,0,0],      # Direction of drag component
            # ----- Force calculation
        KJforce_type                = "regular",     # KJ force evaluated at middle of BV
        # KJforce_type                = "averaged"  # KJ force evaluated at average vortex sheet (if vlm_vortexsheet also true)
        # KJforce_type                = "weighted"  # KJ force evaluated at strength-weighted vortex sheet (if vlm_vortexsheet also true)

        include_trailingboundvortex = false,     # Include bound vortices in force calculations

        include_unsteadyforce       = true,      # Include unsteady force
        add_unsteadyforce           = true,     # Whether to add the unsteady force to Ftot or to simply output it

        include_parasiticdrag       = true,      # Include parasitic-drag force
        add_skinfriction            = true,      # If false, the parasitic drag is purely form, meaning no skin friction
        calc_cd_from_cl             = false,     # Whether to calculate cd from cl or effective AOA

        wing_polar_file  = "xf-n0012-il-500000-n5.csv",   # Polar file from airfoiltools.com for wing parasitic-drag force.

        debug_forces                = false,      # Force calculations will output intermediate fields if true
    )


    #####
    ##### build simulation
    #####
    vehicle = build_ecrm_002(; data_path, add_rotors, VehicleType, n, n_ccb)
    wing = vehicle.system.wings[1]
    b = wing._ywingdcr[end] - wing._ywingdcr[1]
    S_ref = 10.48

    y2b = 2*[vlm.getControlPoint(wing, i)[2] for i in 1:vlm.get_m(wing)]/b

    maneuver = ecrm_002_maneuver(; add_rotors)

    Vinf(x,t) = Vref.*[1.0,0,0] # update this when running
    dt = duration/nsteps

    Vinit = Vref*maneuver.Vvehicle(0)       # (m/s) initial vehicle velocity
    Winit = pi/180 * (maneuver.anglevehicle(0+1e-12)-
        maneuver.anglevehicle(0))/(duration*1e-12) # (rad/s) initial vehicle angular velocity

    sim = uns.Simulation(vehicle, maneuver, Vref, rpm, duration;
        Vinit=Vinit, Winit=Winit)

    # options
    if add_rotors
        R = vehicle.rotor_systems[1][1]._wingsystem.wings[1]._yn[end]
    else
        R = 1.0
    end
    lambda = 1.2*2.125
    sigma_vpm_overwrite = lambda * (2*pi*rpm/60*R + Vref)*dt / p_per_step   # Overwrite cores of wake to this value (ignoring sigmafactor_vpm)
    sigma_vlm_surf=R/50         # Size of embedded particles representing VLM surfaces (for VLM-on-VPM and VLM-on-Rotor)
    sigma_rotor_surf=R/20       # Size of embedded particles representing rotor blade surfaces (for Rotor-on-VPM, Rotor-on-VLM, and Rotor-on-Rotor)
    vlm_vortexsheet_sigma_tbv=sigma_vpm_overwrite # Size of particles in trailing bound vortices (defaults to sigma_vlm_surf if not given)
    vpm_SFS = vpm.DynamicSFS(vpm.Estr_fmm, vpm.pseudo3level_positive; alpha=0.999,
                                clippings=[vpm.clipping_backscatter],
                                controls=[vpm.control_directional, vpm.control_magnitude])

    ############################################################################
    # SOLVERS SETUP
    ############################################################################
    if sigma_vlm_solver<=0
        vlm.VLMSolver._blobify(false)
    else
        vlm.VLMSolver._blobify(true)
        vlm.VLMSolver._smoothing_rad(sigma_vlm_solver)
    end

    # ---------------- SCHEMES -------------------------------------------------
    # Set up viscous scheme
    if vpm.isinviscid(vpm_viscous) == false
        vpm_viscous.nu = mu/rho
        if vpm.iscorespreading(vpm_viscous) && sigma_vpm_overwrite!=nothing
            vpm_viscous.sgm0 = sigma_vpm_overwrite
        end
    end

    # static particles
    blade_n = add_rotors ? n_ccb : 0
    n_blades = add_rotors ? length(vehicle.rotor_systems[1][1]._wingsystem.wings) : 0
    n_rotors = add_rotors ? length(vehicle.rotor_systems[1]) : 0
    wing_n = n
    n_wings = 2
    shed_locations = blade_n * n_blades * n_rotors * 3 + wing_n * n_wings * 3
    max_static_particles = shed_locations

    # Initiate particle field

    b = 10.64
    gamma = 9 * pi / 180
    ctip = 0.819
    croot = 1.149
    dx_cutoff_f = add_rotors ? 1.06 : 1.2 # vinf = 67, RPM = 1800, eCRM-002, main wing
    swept_b2 = (1-0.099) * b/2
    xrotor = ctip + swept_b2 * tan(gamma)
    Rrotor = 1.525950
    dx_cutoff = dx_cutoff_f*b
    width_rect = dx_cutoff + croot + 2*ctip + swept_b2 * tan(gamma)
    height_rect = 2*Rrotor + b/2
    Rcutoff = sqrt((width_rect/2)^2 + height_rect^2)
    xcenter = width_rect/2 - (2*ctip + swept_b2 * tan(gamma))

    nsteps_crit = Int(round(dx_cutoff/67.0/dt * 1.35))
    max_particles = shed_locations * nsteps_crit * p_per_step + max_static_particles

    # for restart file
    if !(restart_vpmfile == nothing)
        nsteps_crit_restart = Int(round(dx_cutoff/67.0/restart_dt * 1.35))
        shed_locations_restart = restart_n * 2 * 3 + restart_nr * 5 * 2 * 3
        restart_max_particles = shed_locations_restart * nsteps_crit_restart * restart_pps + shed_locations_restart
        @show restart_max_particles max_particles
        max_particles = maximum(restart_max_particles, max_particles)
    end
    @show max_static_particles max_particles
    vpm_solver = [
                    (:formulation, vpm_formulation),
                    (:viscous, vpm_viscous),
                    (:kernel, vpm_kernel),
                    (:UJ, vpm_UJ),
                    (:SFS, vpm_SFS),
                    (:integration, vpm_integration),
                    (:transposed, vpm_transposed),
                    (:relaxation, vpm_relaxation),
                    (:fmm, vpm_fmm),
                 ]
    Xdummy = zeros(3)
    pfield = vpm.ParticleField(max_particles; Uinf=t->Vinf(Xdummy, t),
                vpm_solver...)

    staticpfield = vpm.ParticleField(max_static_particles; Uinf=t->Vinf(Xdummy, t),
                vpm_solver...)

    # Save settings
    save_path=run_name
    vpm.save_settings(pfield, run_name; path=save_path)
    if restart_vpmfile!=nothing
        vpm.read!(pfield, restart_vpmfile; overwrite=true, load_time=false)
    end

    if vpm_surface
        static_particles_function = uns.generate_static_particle_fun(pfield, staticpfield,
                                    sim.vehicle, sigma_vlm_surf, sigma_rotor_surf;
                                    vlm_vortexsheet=vlm_vortexsheet,
                                    vlm_vortexsheet_overlap=vlm_vortexsheet_overlap,
                                    vlm_vortexsheet_distribution=vlm_vortexsheet_distribution,
                                    vlm_vortexsheet_sigma_tbv=vlm_vortexsheet_sigma_tbv,
                                    save_path=save_static_particles ? save_path : nothing,
                                    run_name=run_name, nsteps_save=nsteps_save)
    else
        static_particles_function(pfield, t, dt) = nothing
    end

    # ------------ ROUTINE FOR CALCULATION OF FORCES ---------------------------

    forces = []

    # Calculate Kutta-Joukowski force
    kuttajoukowski = uns.generate_calc_aerodynamicforce_kuttajoukowski(KJforce_type,
                                    sigma_vlm_surf, sigma_rotor_surf,
                                    vlm_vortexsheet, vlm_vortexsheet_overlap,
                                    vlm_vortexsheet_distribution,
                                    vlm_vortexsheet_sigma_tbv;
                                    vehicle=vehicle)
    push!(forces, kuttajoukowski)

    # Force due to unsteady circulation
    if include_unsteadyforce
        unsteady(args...; optargs...) = uns.calc_aerodynamicforce_unsteady(args...; add_to_Ftot=add_unsteadyforce, optargs...)
        push!(forces, unsteady)
    end

    # Parasatic-drag force (form drag and skin friction)
    if include_parasiticdrag
        parasiticdrag = uns.generate_aerodynamicforce_parasiticdrag(wing_polar_file;
                                                                read_path=joinpath(data_path, "airfoils"),
                                                                calc_cd_from_cl=calc_cd_from_cl,
                                                                add_skinfriction=add_skinfriction)
        push!(forces, parasiticdrag)
    end

    # Calculation and addition of all forces
    function calc_forces(vlm_system, args...; per_unit_span=false, optargs...)

        # Delete any previous force field
        fieldname = per_unit_span ? "ftot" : "Ftot"
        if fieldname in keys(vlm_system.sol)
            pop!(vlm_system.sol, fieldname)
        end

        Ftot = nothing

        for (fi, force) in enumerate(forces)
            Ftot = force(vlm_system, args...; per_unit_span=per_unit_span, optargs...)
        end

        return Ftot
    end

    # critical length to ignore horseshoe force
    lencrit_f = 0.5
    meanchord = S_ref / b # Sref / b
    lencrit = lencrit_f*meanchord/vlm.get_m(sim.vehicle.system.wings[1])

    # ----------------- WAKE TREATMENT FUNCTION --------------------------------
    rmv_strngth = 2*2/p_per_step * dt/(30/(4*5400))
    minmaxGamma = rmv_strngth*[0.00001, 0.5]
    minmaxsigma = sigma_vpm_overwrite*[0.01, 10]
    # @show minmaxGamma
    # @show minmaxsigma
    wake_treatment_strength = uns.remove_particles_strength( minmaxGamma[1]^2, minmaxGamma[2]^2; every_nsteps=1)
    wake_treatment_sigma = uns.remove_particles_sigma( minmaxsigma[1], minmaxsigma[2]; every_nsteps=1)

    b = 10.64
    gamma = 9 * pi / 180
    ctip = 0.819
    croot = 1.149
    dx_cutoff_f = add_rotors ? 1.06 : 1.2 # vinf = 67, RPM = 1800, eCRM-002, main wing
    swept_b2 = (1-0.099) * b/2
    xrotor = ctip + swept_b2 * tan(gamma)
    Rrotor = 1.525950
    dx_cutoff = dx_cutoff_f*b
    width_rect = dx_cutoff + croot + 2*ctip + swept_b2 * tan(gamma)
    height_rect = 2*Rrotor + b/2
    Rcutoff = sqrt((width_rect/2)^2 + height_rect^2)
    xcenter = width_rect/2 - (2*ctip + swept_b2 * tan(gamma))

    wake_treatment_sphere = uns.remove_particles_sphere(Rcutoff^2, 1; Xoff=[xcenter, 0, 0]) # radius is the span
    remove_particles(args...; optargs...) = (wake_treatment_strength(args...; optargs...) ||
                                             wake_treatment_sigma(args...; optargs...) ||
                                             wake_treatment_sphere(args...; optargs...))

    # remove_particles(args...; optargs...) = false


    function runtime_function(PFIELD, T, DT; vprintln=(args...)->nothing)

        # Move tilting systems, and translate and rotate vehicle
        uns.nextstep_kinematic(sim, dt)

        Vinf = (x,t) -> PFIELD.Uinf(t)

        # Solver-specific pre-calculations
        uns.precalculations(sim, Vinf, PFIELD, T, DT)

        # Shed semi-infinite wake
        uns.shed_wake(sim.vehicle, Vinf, PFIELD, DT, sim.nt; t=T,
                            unsteady_shedcrit=-1,
                            p_per_step=p_per_step, sigmafactor=sigmafactor_vpm,
                            overwrite_sigma=sigma_vpm_overwrite,
                            omit_shedding=omit_shedding)

        # Solve aerodynamics of the vehicle
        uns.solve(sim, Vinf, PFIELD, wake_coupled, DT, vlm_rlx,
                sigma_vlm_surf, sigma_rotor_surf, rho, sound_spd,
                staticpfield, hubtiploss_correction;
                init_sol=vlm_init, sigmafactor_vpmonvlm=sigmafactor_vpmonvlm,
                debug=debug)

        # Shed unsteady-loading wake with new solution
        if shed_unsteady
            uns.shed_wake(sim.vehicle, Vinf, PFIELD, DT, sim.nt; t=T,
                        unsteady_shedcrit=unsteady_shedcrit,
                        shed_starting=shed_starting,
                        p_per_step=p_per_step, sigmafactor=sigmafactor_vpm,
                        overwrite_sigma=sigma_vpm_overwrite,
                        omit_shedding=omit_shedding)
        end

        if shed_boundarylayer
            uns.shed_wake(sim.vehicle, Vinf, PFIELD, DT, sim.nt; t=T,
                        unsteady_shedcrit=-1,
                        p_per_step=p_per_step, sigmafactor=sigmafactor_vpm,
                        overwrite_sigma=sigma_vpm_overwrite,
                        omit_shedding=omit_shedding,
                        shed_boundarylayer=true,
                        prescribed_Cd=boundarylayer_prescribedCd,
                        dipole_d=boundarylayer_d)
        end

        # calculate force
        if PFIELD.nt>2
            # Force at each VLM element
            wing = sim.vehicle.system.wings[1]
            prev_wing = sim.vehicle.prev_data[1].wings[1]
            Ftot = calc_forces(wing, prev_wing, PFIELD, Vinf, DT,
                                                            rho; t=PFIELD.t,
                                                            spandir=uns.cross(L_dir, D_dir),
                                                            lencrit=lencrit,
                                                            include_trailingboundvortex=include_trailingboundvortex,
                                                            debug=debug)
            L, D, S = uns.decompose(Ftot, L_dir, D_dir)
            vlm._addsolution(wing, "L", L)
            vlm._addsolution(wing, "D", D)
            vlm._addsolution(wing, "S", S)


            # Force per unit span at each VLM element
            ftot = calc_forces(wing, prev_wing, PFIELD, Vinf, DT,
                                        rho; t=PFIELD.t, per_unit_span=true,
                                        spandir=uns.cross(L_dir, D_dir),
                                        lencrit=lencrit,
                                        include_trailingboundvortex=include_trailingboundvortex,
                                        debug=debug)
            l, d, s = uns.decompose(ftot, L_dir, D_dir)

            vinf = Vinf([0.0,0,0],0)
            vinf_mag2 = vinf'*vinf
            qinf_ref = 0.5 * rho * vinf_mag2

            # Lift of the wing
            Lwing = sum(L)
            Lwing = sign(uns.dot(Lwing, L_dir))*uns.norm(Lwing)
            CLwing = Lwing/(qinf_ref*S_ref)
            cl = broadcast(x -> sign(uns.dot(x, L_dir)), l) .* uns.norm.(l) / (qinf_ref*S_ref/b)

            # Drag of the wing
            Dwing = sum(D)
            Dwing = sign(uns.dot(Dwing, D_dir))*uns.norm(Dwing)
            CDwing = Dwing/(qinf_ref*S_ref)
            cd = broadcast(x -> sign(uns.dot(x, D_dir)), d) .* uns.norm.(d) / (qinf_ref*S_ref/b)

            vlm._addsolution(wing,"Lwing",Lwing)
            vlm._addsolution(wing,"CLwing",CLwing)
            vlm._addsolution(wing,"cl",cl)
            vlm._addsolution(wing,"Dwing",Dwing)
            vlm._addsolution(wing,"CDwing",CDwing)
            vlm._addsolution(wing,"cd",cd)
            # =#
        end

        # Simulation-specific postprocessing
        breakflag = extra_runtime_function(PFIELD, T, DT; vprintln=vprintln) || remove_particles(sim, PFIELD, T, DT)

        # Output vtks
        if save_path!=nothing && PFIELD.nt%nsteps_save==0
            strn = uns.save_vtk(sim, run_name; path=save_path,
                            save_horseshoes=save_horseshoes,
                            save_wopwopin=save_wopwopin)
        end

        return breakflag
    end

    save_time = false

    cl_history = zeros(length(sim.vehicle.system.wings[1]._xm), nsteps-2)
    cd_history = zeros(length(sim.vehicle.system.wings[1]._xm), nsteps-2)
    CL_history = zeros(nsteps-2)
    CD_history = zeros(nsteps-2)
    data_history = Dict(
        "cl" => cl_history,
        "cd" => cd_history,
        "CL" => CL_history,
        "CD" => CD_history,
        "y2b" => y2b
    )

    vpmdata = VPMData(nsteps, sim, pfield, static_particles_function, runtime_function, S_ref, b, y2b, run_name, save_path, nsteps_save, save_time, v_lvl, verbose_nsteps, data_history)

    return vpmdata
end
