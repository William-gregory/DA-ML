#!/bin/csh -fx                                                                                                                                                                         

cp /ncrc/home2/William.Gregory/DA-ML/CNN_SIS2.py .

set ensemble_size  = `grep ensemble_size ../input.nml | cut -d'=' -f2`
foreach member (1 ${ensemble_size})
   set TMP = `printf %02d ${member}`
   touch -r ice_model.res.ens_$TMP.nc ice_model$TMP.timestamp
end

python CNN_SIS2.py
rm CNN_SIS2.py

foreach member (1 ${ensemble_size})
   set TMP = `printf %02d ${member}`
   touch -r ice_model$TMP.timestamp ice_model.res.ens_$TMP.nc
   rm ice_model$TMP.timestamp
end
