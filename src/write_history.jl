function write_history_h5(filename, vpmdata)
    data_history = vpmdata.data_history
    h5open(filename, "w") do io
        for this_key in keys(data_history)
            @show this_key
            this_data = data_history[this_key]
            write(io, this_key, this_data)
        end
    end
end

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

# # test
# write_history_h5("test.h5", vpmdata)

# stuff = read_history_h5("test.h5")