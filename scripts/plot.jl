this_data = fp.read_history_h5(filename)

fig = figure()
fig.clear()
fig.add_subplot(211, ylabel="lift coefficient")
fig.add_subplot(212, ylabel="drag coefficient", xlabel="timestep")
ax_cl = fig.get_axes()[1]
ax_cd = fig.get_axes()[2]
ax_cl.plot(1:length(this_data["CL"]), this_data["CL"])
ax_cd.plot(1:length(this_data["CD"]), this_data["CD"])

fig2 = figure()
fig2.clear()
fig.add_subplot(211, ylabel="lift coefficient")
fig.add_subplot(212, ylabel="drag coefficient", xlabel="spanwise location")
ax_cl = fig.get_axes()[1]
ax_cd = fig.get_axes()[2]
for i in 1:strider:size(this_data["cl"])[2]
    ax_cl.plot(this_data["y2b"], this_data["cl"][:,i], color=(0.0,0.1, 1-i/size(this_data["cl"])[2]))
    ax_cd.plot(this_data["y2b"], this_data["cd"][:,i], color=(0.0,0.1, 1-i/size(this_data["cd"])[2]))
end