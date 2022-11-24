function update_ecrm_002(vpmdata, data_list)
    wing = vpmdata.simulation.vehicle.system.wings[1]
    _xlwingdcr, _xtwingdcr, _ywingdcr, _zlwingdcr, _ztwingdcr, 
        _xn, _yn, _zn, _xm, _ym, _zm = data_list
    wing._xlwingdcr .+= _xlwingdcr
    wing._xtwingdcr .+= _xtwingdcr
    wing._ywingdcr .+= _ywingdcr
    wing._zlwingdcr .+= _zlwingdcr
    wing._ztwingdcr .+= _ztwingdcr
    wing._xn .+= _xn
    wing._yn .+= _yn
    wing._zn .+= _zn
    wing._xm .+= _xm
    wing._ym .+= _ym
    wing._zm .+= _zm
end