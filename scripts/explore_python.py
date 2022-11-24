import julia
from julia import FLOWPy as fp

nsteps = 100 # timesteps
duration = 0.1 # seconds
vinf = 67.0 # m/s
alpha_rad = 4.0 * 3.1415926/180.0 # aoa in radians
restart_vpmfile = "pathtorestartfile"
vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile)

for i in range(1,nsteps):
	fp.step_vpm(vpmdata, vinf, alpha_rad)
