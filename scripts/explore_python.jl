import FLOWPy as fp
using PyPlot
# prep restart file

rpm = 1800.0
nrevs = 1/5 # 100
duration = nrevs / (rpm/60)
dt = 0.001 # approximately
nsteps = Int(ceil(duration/dt))

vinf = 67.0 # m/s
alpha_rad = 10.0 * 3.1415926/180.0 # aoa in radians
restart_vpmfile = nothing
run_name="FLOWUnsteadyPy"
vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile, add_rotors=true, run_name)

# fp.update_ecrm_002(vpmdata, data_list)

for i in 0:nsteps
    save_file = i == nsteps ? true : false
    fp.step_vpm(vpmdata, vinf, alpha_rad; save_file)
end
