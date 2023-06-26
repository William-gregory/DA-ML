#!/bin/csh -fx                                                                                                                                                                         

cd ../
cp compute_increments_OPTp1.py .
cp update_model_state_CNN.py RESTART/.

python compute_increments_OPTp1.py
rm compute_increments_OPTp1.py
mv dSICN_increment.npy RESTART/.

cd RESTART
set ensemble_size  = `grep ensemble_size ../input.nml | cut -d'=' -f2`
foreach member (1 ${ensemble_size})
   set TMP = `printf %02d ${member}`
   touch -r ice_model.res.ens_$TMP.nc ice_model$TMP.timestamp
end

python update_model_state_CNN.py
rm update_model_state_CNN.py
rm dSICN_increment.npy

foreach member (1 ${ensemble_size})
   set TMP = `printf %02d ${member}`
   touch -r ice_model$TMP.timestamp ice_model.res.ens_$TMP.nc
   rm ice_model$TMP.timestamp
end
