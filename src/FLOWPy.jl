module FLOWPy

import FLOWUnsteady as uns
vlm = uns.vlm
vpm = uns.vpm

for filename in ["build_ecrm_002", "vpmdata", "step_vpm"]
	include(filename*".jl")
end

end # module
