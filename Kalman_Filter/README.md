### Data Assimilation with Kalman Filtering

This python code uses Ensemble Kalman Filtering to assimilate satellite observations of sea ice concentration into the SIS2 sea ice model. This procedure is performed offline during a numerical simulation, updating the model's sea ice concentration state variables within the output restart files. This is effectively a post-processing step which can be added to the model's run script (using `freInclude` for GFDL users of FRE) as folows:

         cd $work/RESTART                                                                                                                                                       
         if (! -e coupler.res) then                                                                                                                                            
             echo model still running                                                                                                                                           
         else                                                                                                                                                                   
             cp EAKF_sequential_mpi.py .                                                                                                          
             scontrol show hostnames | head -n 2 > hostfile.txt                                                                                                                 
             mpiexec -n 2 --hostfile hostfile.txt python3 EAKF_sequential_mpi.py >& EAKF_log.out || exit 12                                
             wait                                                                                                                                                               
         endif

I originally wrote the code `ENKF_mpi.py` to get the code structure largely in place. This isn't super flexible as it assumes that the observations and model are on the same 2D structured grid. It is however quite fast to run. It divides up the global domain into a series of tiles (appropriately padding each tile depending on the chosen localization radius) and performs the Ensemble Kalman Filter on each tile using parallel (mpi) processing. This code has a one-time cost associated with computing the tiles and their halos. This information is then saved to disk and used in subsequent DA calls.

The `EAKF_sequential_mpi.py` code is then implementations of the Ensemble Adjustment Kalman Filter, as outlined in [Jeff Anderson's 2003 paper](https://doi.org/10.1175/1520-0493(2003)131<0634:ALLSFF>2.0.CO;2). This code assimilates the observations on their native grid, and is therefore more flexible. However because of the sequential nature of the implementation, I haven't figured out the best way to parallelize this code with mpi--eventually it would be good to take the same approach as used for the Data Assimilation Research Testbed ([DART](https://docs.dart.ucar.edu/en/latest/guide/forward_operator.html#parallelism-implementation-details)). The current parallel implementation here just splits the job into 2 tasks to process the Arctic and Antarctic in parallel.

If you would like to test this code, I would recommend using the `EAKF_sequential_testing.py` code, which should run on your local machine (may be slow depending on the horizontal grid resolution and number of obs being assimilated). 
