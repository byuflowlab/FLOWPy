import FLOWPy as fp

# prep restart file

rpm = 1800.0
nrevs = 1/5 * 100
duration = nrevs / (rpm/60)
dt = 0.001 # approximately
nsteps = Int(ceil(duration/dt))

vinf = 67.0 # m/s
alpha_rad = 4.0 * 3.1415926/180.0 # aoa in radians
restart_vpmfile = nothing
vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile)

for _ in 0:nsteps
    fp.step_vpm(vpmdata, vinf, alpha_rad)
end

#####
##### continue stepping from restart file
#####
#=
restart_vpmfile = "second_run/FLOWUnsteadyPy_pfield.7.h5"
vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile)
for _ in 0:nsteps
    fp.step_vpm(vpmdata, vinf, alpha_rad)
end
=#