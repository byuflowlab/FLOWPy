include("plot.jl")

function history_tag(alpha, dt, n, nr, rotors)
    alpha_tag = "a$(Int(round(alpha; digits=0)))"
    get_digits = dt > 0.0009999 ? 3 : 4
    dt_tag = "_dt$(round(dt; digits=get_digits))"
    n_tag = "_n$n"
    nr_tag = "_nr$nr"
    rotors_tag = "_rotors$rotors"
    return alpha_tag*dt_tag*n_tag*nr_tag*rotors_tag*"_history.h5"
end

#####
##### no rotors
#####
# dt
fig_name = "norotors_dt"
alpha = [6.0,6,6]
dt = [0.002, 0.001, 0.0005]
n = [10,10,10]
nr = [10,10,10]
rotors = [false, false, false]
striders = [1,2,4]
series_names = ["dt = $(round(dt[i]; digits=dt[i] > 0.0009999 ? 3 : 4))s" for i in 1:length(dt)]
for i in 1:length(alpha)
    filename = history_tag(alpha[i], dt[i], n[i], nr[i], rotors[i])
    hold = i == 1 ? false : true
    make_legend = i == length(alpha) ? true : false
    plot_this!(fig_name, series_names[i], "convergence_results/"*filename, dt[i], hold, make_legend; strider_cl=nothing, strider_CL=striders[i], color_cl=(0.0,1-i/length(alpha), i/length(alpha)))
end

# n
fig_name = "norotors_n"
alpha = [6.0,6,6,6]
dt = [0.001, 0.001, 0.001, 0.001]
n = [10,20,40,80]
nr = [10,10,10,10]
rotors = [false, false, false, false]
series_names = ["n = $(n[i])" for i in 1:length(n)]
for i in 1:length(alpha)
    filename = history_tag(alpha[i], dt[i], n[i], nr[i], rotors[i])
    @show filename
    hold = i == 1 ? false : true
    make_legend = i == length(alpha) ? true : false
    plot_this!(fig_name, series_names[i], "convergence_results/"*filename, dt[i], hold, make_legend; strider_cl=nothing, color_cl=(0.0,1-i/length(alpha), i/length(alpha)))
end

#####
##### rotors
#####

# n
fig_name = "rotors_n"
alpha = [6.0,6,6,6]
dt = [0.001, 0.001, 0.001, 0.001]
n = [10,20,40,80]
nr = [10,10,10,10]
rotors = [true, true, true, true]
series_names = ["n = $(n[i])" for i in 1:length(n)]
for i in 1:length(alpha)
    filename = history_tag(alpha[i], dt[i], n[i], nr[i], rotors[i])
    @show filename
    hold = i == 1 ? false : true
    make_legend = i == length(alpha) ? true : false
    plot_this!(fig_name, series_names[i], "convergence_results/"*filename, dt[i], hold, make_legend; strider_cl=nothing, color_cl=(0.0,1-i/length(alpha), i/length(alpha)))
end

# nr
fig_name = "rotors_nr"
alpha = [6.0,6,6]
dt = [0.001, 0.001, 0.001]
n = [80,80,80]
nr = [10,20,40]
rotors = [true, true, true]
series_names = ["n=80, nr = $(nr[i])" for i in 1:length(n)]
for i in 1:length(alpha)
    filename = history_tag(alpha[i], dt[i], n[i], nr[i], rotors[i])
    @show filename
    hold = i == 1 ? false : true
    make_legend = i == length(alpha) ? true : false
    plot_this!(fig_name, series_names[i], "convergence_results/"*filename, dt[i], hold, make_legend; strider_cl=nothing, color_cl=(0.0,1-i/length(alpha), i/length(alpha)))
end

# n2
fig_name = "rotors_n2"
alpha = [6.0,6,6]
dt = [0.001, 0.001, 0.001]
n = [20,40,80]
nr = [20,20,20]
rotors = [true, true, true]
series_names = ["n=$(n[i]), nr = 20" for i in 1:length(n)]
for i in 1:length(alpha)
    filename = history_tag(alpha[i], dt[i], n[i], nr[i], rotors[i])
    @show filename
    hold = i == 1 ? false : true
    make_legend = i == length(alpha) ? true : false
    plot_this!(fig_name, series_names[i], "convergence_results/"*filename, dt[i], hold, make_legend; strider_cl=nothing, color_cl=(0.0,1-i/length(alpha), i/length(alpha)))
end
