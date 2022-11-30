import julia
from julia import FLOWPy as fp

nsteps = 100 # timesteps
duration = 0.1 # seconds
vinf = 67.0 # m/s
alpha_rad = 6.0 * 3.1415926/180.0 # aoa in radians
restart_vpmfile = "restart/a6_dt0.001_n20_nr20_rotorstrue_pfield.300.h5""
add_rotors = True
n=20
n_ccb=20
vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile, add_rotors, n, n_ccb)

# duration / nsteps = 0.001

for i in range(1,nsteps):
    # fp.update_ecrm_002(vpmdata, data_list)
	fp.step_vpm(vpmdata, vinf, alpha_rad)
    if i>3: # note: since the vpm is an unsteady method, it takes a couple of steps
            # before forces can be evaluated
        forces = fp.forces(vpmdata) # list of 3-length vectors corresponding
                                    # to each control point; length 2*n
