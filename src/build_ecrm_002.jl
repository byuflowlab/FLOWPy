function build_ecrm_002(;
    # maneuver parameters
    vinf=67.0,
    rpm=1800,
    # atmospheric parameters
    rho=1.225,
    mu=1.789e-5,
    # wing parameters
    add_wings=true,
    b = 10.64,
    AR = b^2/(10.48 + 10.64*0.099*1.149),
    n = 10,
    pos = [0,0.099,1.0],
    clen = [1.149, 1.149, 0.819],
    twist = [0.0,0,1],
    sweep = [0,-9.0],
    dihed = [0,6.0],
    # rotor parameters
    add_rotors=true,
    CW_starboard = false,
    n_ccb=10,
    blade_r=b/2 * tand(-sweep[1]),
    pitch=0,
    ReD=rho*sqrt(vinf^2 + (rpm*2*pi/60*blade_r)^2)*blade_r/10 / mu,
    rotor_file="uber_ecrm002.csv",
    # rotor_file="/home/ryan/.julia/dev/FLOWUnsteady/data/rotors/uber_ecrm002.csv"
    # vehicle parameters
    VehicleType=uns.VLMVehicle,
    # FLOWUnsteady parameters
    v_lvl=1,
    verbose=true,
    xfoil=false,
    data_path=joinpath(dirname(pathof(uns)), "..", "data/"),
)


    # build wing
    wing = vlm.complexWing(b, AR, n, pos, clen, twist, sweep, dihed)

    # build rotors
    if add_rotors
        rotor_starboard = uns.generate_rotor(rotor_file; pitch=pitch,
            n=n_ccb, blade_r=blade_r, CW=CW_starboard, ReD=ReD,
            verbose=verbose, v_lvl=v_lvl+2, xfoil=xfoil,
            data_path=data_path, plot_disc=false)
        rotor_port = uns.generate_rotor(rotor_file; pitch=pitch,
            n=n_ccb, blade_r=blade_r, CW=!CW_starboard, ReD=ReD,
            verbose=verbose, v_lvl=v_lvl+2, xfoil=xfoil,
            data_path=data_path, plot_disc=false)
        swept_b2 = (pos[end]-pos[end-1]) * b/2
        this_O = [swept_b2 * tand(sweep[end]) - clen[end], b/2, swept_b2 * tand(dihed[end])]
        this_Oaxis = [
            1.0 0 0
            0 1 0
            0 0 1
        ]
        vlm.setcoordsystem(rotor_starboard, this_O, this_Oaxis; user=true)
        vlm.setcoordsystem(rotor_port, this_O .* [1,-1,1], this_Oaxis; user=true)
        rotors = [rotor_starboard,rotor_port]
    end

    # # Position of main wing
    # O_w = [l_pos_w + x_off_w, 0, h_pos_w]
    # Oaxis_w = gt.rotation_matrix(0.0, 0.0, 0.0)
    # vlm.setcoordsystem(main_wing, O_w, Oaxis_w)

    # system object
    system = vlm.WingSystem()
    vlm.addwing(system, "MainWing", wing)
    if add_rotors
        for (i,rotor) in enumerate(rotors)
            vlm.addwing(system, "rotor$i", rotor)
        end
    end

    # Tilting systems
    tilting_systems = ()

    # Rotors grouped by systems of the same RPM
    if add_rotors
        rotor_systems = (rotors,)
    else
        rotor_systems = ()
    end

    # System to solve through the VLM solver
    vlm_system = vlm.WingSystem()
    if add_wings
        vlm.addwing(vlm_system, "MWing", wing)
    end

    # Wake-shedding system (`vlm_system`+`rotors`)
    wake_system = vlm.WingSystem()
    if add_wings
        vlm.addwing(wake_system, "SolveVLM", vlm_system)
    end

    if add_rotors
        if VehicleType==uns.VLMVehicle
            for (i, rotor) in enumerate(rotors)
                vlm.addwing(wake_system, "Rotor$i", rotor)
            end
        else
            # Mute colinear warnings. This is needed since the quasi-steady solver
            #   will probe induced velocities at the lifting line of the blade
            uns.vlm.VLMSolver._mute_warning(true)
        end
    end

    vehicle = VehicleType(system;
            tilting_systems=tilting_systems,
            rotor_systems=rotor_systems,
            vlm_system=vlm_system,
            # wake_system=vlm_system,
            wake_system=wake_system,
            # grids=grids
        )
end

function ecrm_002_maneuver(; add_rotors=true)
    angle(t) = [0.0,0,0]
    RPM(t) = 1.0
    Vvehicle(t) = [0.0,0,0]
    anglevehicle(t) = [0.0,0,0]
    rpm_sys = add_rotors ? (RPM,) : ()
    uns.KinematicManeuver((), rpm_sys, Vvehicle, anglevehicle)
end
