using HDF5
using PyPlot

function read_history_h5(filename)
    keys = ["cl", "cd", "CL", "CD", "y2b"]
    this_data = Dict()
    h5open(filename, "r") do io
        for key in keys
            this_data[key] = read(io, key)
        end
    end
    return this_data
end

function plot_this!(fig_name, series_name, filename, dt, hold, make_legend; strider_cl=1, color_cl=nothing)
    this_data = read_history_h5(filename)

    fig = figure(fig_name*"_CL")
    if !hold
        fig.clear()
        fig.add_subplot(211, ylabel="lift coefficient")
        fig.add_subplot(212, ylabel="drag coefficient", xlabel="time, seconds")
    end
    ax_CL = fig.get_axes()[1]
    ax_CD = fig.get_axes()[2]
    this_CL = this_data["CL"]
    this_CD = this_data["CD"]
    t = range(0,stop=dt*(length(this_CL)-1), length=length(this_CL))
    ax_CL.plot(t, this_CL)#, label=series_name)
    ax_CD.plot(t, this_CD, label=series_name)
    if make_legend
        ax_CD.legend()
    end

    fig2 = figure(fig_name*"_cldist")
    if !hold
        fig2.clear()
        fig2.add_subplot(211, ylabel="lift coefficient")
        fig2.add_subplot(212, ylabel="drag coefficient", xlabel="spanwise location")
    end
    ax_cl = fig2.get_axes()[1]
    ax_cd = fig2.get_axes()[2]
    this_index = strider_cl == nothing ? [size(this_data["cl"])[2]] : 1:strider_cl:size(this_data["cl"])[2]
    for i in this_index
        color_i = color_cl == nothing ? (0.0,1-i/size(this_data["cl"])[2], i/size(this_data["cl"])[2]) : color_cl
        ax_cl.plot(this_data["y2b"], this_data["cl"][:,i], color=color_i)#, label=series_name)
        ax_cd.plot(this_data["y2b"], this_data["cd"][:,i], color=color_i, label=series_name)
    end
    if make_legend
        ax_cd.legend()
    end
end
