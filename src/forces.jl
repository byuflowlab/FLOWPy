function forces(vpmdata::VPMData)
    return vpmdata.simulation.vehicle.system.wings[1].sol["Ftot"]
end
