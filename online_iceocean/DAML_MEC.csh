#!/bin/csh -fx                                                                                                                                                                                                            

set ensemble_size  = `grep ensemble_size ../input.nml | cut -d'=' -f2`
foreach member (1 ${ensemble_size})
   set TMP = `printf %02d ${member}`
   touch -r ice_model.res.ens_$TMP.nc ice_model$TMP.timestamp
end

python update_model_state_MEC.py
rm update_model_state_MEC.py

foreach member (1 ${ensemble_size})
   set TMP = `printf %02d ${member}`
   touch -r ice_model$TMP.timestamp ice_model.res.ens_$TMP.nc
   rm ice_model$TMP.timestamp
end
