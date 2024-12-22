### Data Assimilation with Kalman Filtering

This python code uses Ensemble Kalman Filtering to assimilate satellite observations of sea ice concentration into the SIS2 sea ice model. This procedure is performed offline during a numerical simulation, updating the model's sea ice concentration state variables within the output restart files. This is effectively a post-processing step which can be added to the model's run script (using freInclude for GFDL users of FRE) as folows:

         cd $work/RESTART                                                                                                                                                       
         if (! -e coupler.res) then                                                                                                                                            
             echo model still running                                                                                                                                           
         else                                                                                                                                                                   
             cp EAKF_sequential_mpi.py .                                                                                                          
             scontrol show hostnames | head -n 2 > hostfile.txt                                                                                                                 
             mpiexec -n 2 --hostfile hostfile.txt python3 EAKF_sequential_mpi.py >& EAKF_log.out || exit 12                                
             wait                                                                                                                                                               
         endif

I originally wrote the code `ENKF_singlenode.py`
