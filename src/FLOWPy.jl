module FLOWPy

import FLOWUnsteady as uns
vlm = uns.vlm
vpm = uns.vpm
using HDF5

for filename in ["build_ecrm_002", "vpmdata", "step_vpm", "update_ecrm_002", "forces", "write_history"]
	include(filename*".jl")
end

end # module
