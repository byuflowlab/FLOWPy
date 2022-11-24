import FLOWPy as fp
using PyPlot
# prep restart file

rpm = 1800.0
nrevs = 1/5 * 35 # 100
duration = nrevs / (rpm/60)
dt = 0.001 # approximately
nsteps = Int(ceil(duration/dt))

vinf = 67.0 # m/s
alpha_rad = 10.0 * 3.1415926/180.0 # aoa in radians
restart_vpmfile = nothing
run_name="FLOWUnsteadyPy"
vpmdata = fp.VPMData(nsteps, duration; restart_vpmfile, add_rotors=false, run_name)

# fp.update_ecrm_002(vpmdata, data_list)

for i in 0:nsteps
    save_file = i == nsteps ? true : false
    fp.step_vpm(vpmdata, vinf, alpha_rad; save_file)
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



directory = ".."
filename = joinpath(directory, run_name * "_history.h5")
include("plot.jl")