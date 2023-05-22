#!/bin/bash                                                                                                                                                                                                                             
#SBATCH --nodes=1                                                                                                                                                                                                                       
#SBATCH --ntasks=1                                                                                                                                                                                                                      
#SBATCH --cpus-per-task=1                                                                                                                                                                                                               
#SBATCH --mem=250G                                                                                                                                                                                                                      
#SBATCH --time=04:00:00                                                                                                                                                                                                                 
#SBATCH --gres=gpu:1                                                                                                                                                                                                                    
#SBATCH --job-name=CNN                                                                                                                                                                                                                 

module load intel/2021.1.2 openmpi/intel-2021.1/4.1.0 hdf5/intel-2021.1/1.10.6 netcdf/intel-19.1/hdf5-1.10.6/4.7.4

python final_network.py
