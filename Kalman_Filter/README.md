### Data Assimilation with Kalman Filtering

This python code was written to 

         cd $work/RESTART                                                                                                                                                                    
         if (! -e coupler.res) then                                                                                                                                                          
             echo model still running                                                                                                                                                        
         else                                                                                                                                                                                
             cp EAKF_sequential_mpi.py .                                                                                                          
             scontrol show hostnames | head -n 2 > hostfile.txt                                                                                                                              
             mpiexec -n 2 --hostfile hostfile.txt python3 EAKF_sequential_mpi.py >& EAKF_log.out || exit 12                                
             wait                                                                                                                                                                            
         endif
