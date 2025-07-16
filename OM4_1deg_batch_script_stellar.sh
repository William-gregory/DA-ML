#!/bin/csh -fx
#SBATCH --chdir=/scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp/stdout/run
#SBATCH --output=/scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp/stdout/run/%x.o%j
#SBATCH --job-name=OM4_hybrid_seaiceML
#SBATCH --mail-user=wg4031@princeton.edu
#SBATCH --mail-type=fail
#SBATCH --ntasks=72
#SBATCH --exclusive
#SBATCH --qos=cimes-short
#SBATCH --mem-per-cpu=2G
#SBATCH --time=10:00:00

set -r echoOn = $?echo
set -r runtimeBeg = `date "+%s"`

if ( $echoOn ) unset echo
echo '<NOTE> : ====== FRE RUNSCRIPT ======'
echo "<NOTE> : Starting at $HOST on `date`"
if ( $echoOn ) set echo

unalias *

if ( $echoOn ) unset echo
set -r modulesHomeDir = $MODULESHOME
source $modulesHomeDir/init/tcsh
if ( $echoOn ) set echo

if ( $?SLURM_JOB_ID ) then
   tty -s
   if ($status) then
      set -r FRE_JOBID = $SLURM_JOB_NAME:t.o$SLURM_JOB_ID
      set -r batch
   else
      set -r FRE_JOBID = $0:t.o`date +%s`
   endif
else
   set -r FRE_JOBID = $0:t.o`date +%s`
endif

if ( $echoOn ) unset echo
echo "################################################################################"
echo "# $FRE_JOBID"
echo "################################################################################"
if ( $echoOn ) set echo

################################################################################
#---------------- global constants and variables, set by frerun ----------------
################################################################################

set -r platform = ncrc5.intel-classic
set -r target = repro-openmp
set -r name = OM4_hybrid_seaiceML
set -r stdoutDir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp/stdout/run
set -r stdoutTmpDir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp/stdout/run
set -r stateDir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp/state/run
set -r workDir = /scratch/cimes/wg4031/OM4_1deg/work/$FRE_JOBID
set -r ptmpDir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ptmp
set -r archiveDir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp
set -r scriptName = /home/wg4031/OM4_runs/OM4_hybrid_seaiceML
set -r executable = /scratch/cimes/wg4031/MOM6-examples/build/ice_ocean_SIS2/OM4_seaiceML/OM4_seaiceML
set -r segmentsPerSimulation = 1
set segmentsPerPPCall = 0
set -r segmentsPerJob = 1
set -r monthslist = ( 12 )
set -r dayslist = ( 0 )
set -r hourslist = ( 0 )
set -r timeStampOptions = ( -f digital )
set -r baseDate = '2018 1 1 0 0 0'
set -r mailMode = fail
set -r mailList = wg4031@princeton.edu
set -r includeDir = /scratch/cimes/wg4031/SPEAR_include
set -r ardiffTmpdir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ptmp

set -r flagRunTypeProduction
set -r flagCheckSumOff
set -r flagWorkDirCleanOn
set -r flagOutputTypeOverwrite
set -r flagOutputFormat64Bit
set -r flagOutputStagingTypeChained
set -r flagOutputCacheHistoryOff
set -r flagOutputCombineHistoryOn
set -r flagOutputCompressAsciiOff
set -r flagOutputCompressRestartOff
set -r flagOutputCompressHistoryOff
set -r flagOutputArchiveOn
set -r flagOutputPostProcessOn
set -r flagOutputXferOn
set -r flagOutputCheckOff
set -r flagVerbosityOff
set -r flagOutputFillGridOff

set outputDir = /scratch/cimes/wg4031/OM4_1deg/OM4_hybrid_seaiceML/ncrc5.intel-classic-repro-openmp
set gridSpec = /scratch/cimes/wg4031/inputs/mosaic.tar
set initCond = /scratch/cimes/wg4031/initCond_OM4/20180101.tar

  set -r npes = 72
  set -r atm_ranks = 
  set -r tot_atm_ranks = 0
  set -r atm_threads = 1
  set -r atm_layout =
  set -r atm_io_layout =
  set -r atm_mask_table = 
  set -r atm_hyperthread = .true.
  set -r scheduler_atm_threads =
  set -r atm_nxblocks = 1
  set -r atm_nyblocks = 1
  set -r ocn_ranks = 72
  set -r tot_ocn_ranks = 72
  set -r ocn_threads = 1
  set -r ocn_layout = 12,6
  set -r ocn_io_layout = 1,1
  set -r ocn_mask_table = MOM_mask_table
  set -r ocn_hyperthread = .false.
  set -r scheduler_ocn_threads = 1
  set -r lnd_ranks = 
  set -r tot_lnd_ranks = 
  set -r lnd_threads = 
  set -r lnd_layout = 
  set -r lnd_io_layout =
  set -r lnd_mask_table = 
  set -r lnd_hyperthread = .false.
  set -r scheduler_lnd_threads = 
  set -r ice_ranks = 
  set -r tot_ice_ranks = 
  set -r ice_threads = 
  set -r ice_layout = 12,6
  set -r ice_io_layout = 1,1
  set -r ice_mask_table = 
  set -r ice_hyperthread = .false.
  set -r scheduler_ice_threads = 

alias runCommand time `which srun` --ntasks=$tot_ocn_ranks --cpus-per-task=$scheduler_ocn_threads --cpu-bind=cores --export=ALL,OMP_NUM_THREADS=$ocn_threads ./$executable:t
 
set -r mppnccombineOptsRestart = '-h 65536 -m'
set -r mppnccombineOptsHistory = '-h 65536 -m'

set -r FreCommandsSrcDir = /scratch/cimes/wg4031/MOM6-examples/src/
set -r FreCommandsBldDir = /scratch/cimes/wg4031/MOM6-examples/build/

################################################################################
#------------------------ global constants and aliases -------------------------
################################################################################

if ( $echoOn ) unset echo

# Platform environment defaults from /ncrc/home2/fms/local/opt/fre-commands/bronx-21/site/ncrc5/env.defaults.intel-classic
module load intel-rt/2021.1.2 intel-tbb/2021.1.1 intel-mkl/2021.1.1 intel-debugger/10.0.0 intel-dpl/2021.1.2 /opt/intel/oneapi/compiler/2021.1.2/linux/lib/oclfpga/modulefiles/oclfpga
module load intel/2021.1.2 openmpi/intel-2021.1/4.1.0 hdf5/intel-2021.1/1.10.6 netcdf/intel-19.1/hdf5-1.10.6/4.7.4 nco/netcdf-4.7.4/hdf5-1.10.6/5.0.3

setenv LD_LIBRARY_PATH "/usr/lib:/usr/local/openmpi/4.1.0/intel20211/lib64:/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64:${LD_LIBRARY_PATH}"
      
module list

set -r architecture = 'x86_64'
set machine = `uname -m`

if ( $machine != $architecture ) then
   if ( $echoOn ) unset echo
   echo "*ERROR*: The script '$scriptName' is intended for the machine architecture '$architecture'"
   if ( $echoOn ) set echo
   exit 1
endif

unset machine

set -r tmpOutputDir = $workDir/output.stager
set -r envFile = /tmp/shell.variables.$FRE_JOBID
set -r envFileDelay = 2

set -r patternGrepAscii = '\<out\>|\<results\>|\<log\>|\<timestats\>|\<stats\>|\<velocity_truncations\>'
set -r patternGrepRestart = '\<res\>|\<nc\>|\.input.\tgz$|\.ww3$'
set -r patternGrepRestartNextDrop = '\<res\>'
set -r patternGrepRestartNextMove = '\<res\>|\<nc\>|\.ww3$'
set -r patternGrepHistory = '\<nc\>|\.ww3$'
set -r patternGrepRegion = '^rregion'

set -r patternSedHome = 's/^\/(autofs|ncrc)\/.+\/'$USER'\//\/home\/'$USER'\/'
set -r patternSedF5 = 's|^/scratch/cimes/'$USER'/|/home/'$USER'/|'

set -r patternSedSCRATCH = "$patternSedF5"
set -r patternSedDEV = "$patternSedF5"
set -r patternSedCTMP = "$patternSedSCRATCH"
set -r patternSedCPERM = "$patternSedDEV"

set -r archExt = 'tar'

set -r outputStagingType = `set -r | grep '^flagOutputStagingType' | sed 's/flagOutputStagingType//'`

alias expandVariables /scratch/cimes/wg4031/fre_scripts/expand_variables --verbose
alias findModuleInfo /scratch/cimes/wg4031/fre_scripts/find_modules_info --verbose
alias findXIncludes /scratch/cimes/wg4031/fre_scripts/find_xincludes --verbose
alias prepareDir /scratch/cimes/wg4031/fre_scripts/prepare_dir.csh
alias timeStamp /scratch/cimes/wg4031/fre_scripts/time_stamp.csh $timeStampOptions
alias workDirCleaner /scratch/cimes/wg4031/fre_scripts/batch_rmdir.csh
alias adjust_dry_mass_tool /scratch/cimes/wg4031/fre_scripts/adjust_dry_mass.csh

set -r workDirCleaner = `alias workDirCleaner`

alias outputStager /scratch/cimes/wg4031/fre_scripts/output.stager

set -r outputStager = `alias outputStager`
@ outputStagerErrors = 0

################################################################################
#--------------------------------- environment ---------------------------------
################################################################################

limit stacksize unlimited
limit coredumpsize unlimited
limit

if ( $#dayslist != $segmentsPerJob || $#monthslist != $segmentsPerJob || $#hourslist != $segmentsPerJob ) then
   if ( $echoOn ) unset echo
   echo "*ERROR*: dayslist, monthslist and hourslist lengths must be equal to a number of segments per job"
   if ( $echoOn ) set echo
   exit 1
endif

  # NiNaC not loaded when script created

################################################################################
#----------------------------- global variables --------------------------------
################################################################################

set continueFlag = 1

set combineList = ( )
set saveJobIds = ( )
set argFiles = ( )

@ currentSeg = 1
@ irun = 1

################################################################################
#--------------------------- before the main loop ------------------------------
################################################################################

if ( $?flagRunTypeProduction ) then
   mkdir -p $stateDir
   if ( $status != 0 ) then
      if ( $echoOn ) unset echo
      echo "*ERROR*: Unable to setup the state directory '$stateDir'"
      if ( $echoOn ) set echo
      exit 1
   endif
   set reload_file = $stateDir/reload_commands

   if ( -f $reload_file ) then
      if ( -r $reload_file ) then
         source $reload_file
         if ( $#argFiles > 0 ) then
            if ( `echo $argFiles | tr ' ' "\n" | grep --count "^$FRE_JOBID"` != $#argFiles ) then
               set saveJobIds = ( )
               set argFiles = ( )
            endif
         endif
      else
         if ( $echoOn ) unset echo
         echo "*ERROR*: The reload file '$reload_file' is not readable"
         if ( $echoOn ) set echo
         exit 1
      endif
   endif

   set queue_file = $stateDir/queue_commands

   if ( -f $queue_file ) then
      if ( -r $queue_file ) then
         source $queue_file
      else
         if ( $echoOn ) unset echo
         echo "*ERROR*: The queue file '$queue_file' is not readable"
         if ( $echoOn ) set echo
         exit 1
      endif
   else
      touch $queue_file
      if ( $status == 0 ) then
         if ( $echoOn ) unset echo
         echo "<NOTE> : Writing queue information to the queue file '$queue_file' at `date +%s`"
         if ( $echoOn ) set echo
         echo "set continueFlag = $continueFlag" >> $queue_file
         chmod 644 $queue_file
      else
         if ( $echoOn ) unset echo
         echo "*ERROR*: The queue file '$queue_file' can't be saved"
         if ( $echoOn ) set echo
         exit 1
      endif
   endif

   if ( ! $continueFlag ) then
      if ( $echoOn ) unset echo
      echo "<NOTE> : Stopping execution"
      if ( $echoOn ) set echo
      exit 0
   endif
endif

if ( $?fyear ) then
   #remove leading zeros, fyear as integer
   set fyearint = `echo $fyear | sed 's/^0*//'`
   if ( ${fyearint} > 0 ) then
      @ fyearm1 = ${fyearint} - 1
      set fyearm1 = `printf "%04d" $fyearm1`
   else
      set fyearm1 = '0000'
   endif
   @ fyearp1 = ${fyearint} + 1
   set fyearp1 = `printf "%04d" $fyearp1`
   @ fyearp2 = ${fyearint} + 2
   set fyearp2 = `printf "%04d" $fyearp2`
endif

if ( $?ireload ) then
   # Using old way to calculate currentSeg
   # Get best guess --- may not be correct if user changed number of segments
   # per job after job started --- frerun -e does not update state file
   @ currentSeg = ( $ireload - 1 ) * $segmentsPerJob + $irun
endif

if ( -e $workDir ) then
   if ( -d $workDir ) then
      if ( -r $workDir ) then
         if ( -w $workDir ) then
            ls -1 --directory --file-type $workDir/* | grep --fixed-strings --invert-match $tmpOutputDir | xargs rm --force --recursive
            mkdir -p $workDir/INPUT 'clean'     || exit 1
            mkdir -p $workDir/RESTART 'clean'   || exit 1
         else
            if ( $echoOn ) unset echo
            echo "*ERROR*: The directory '$workDir' exists, but is not writable"
            if ( $echoOn ) set echo
            exit 1
         endif
      else
         if ( $echoOn ) unset echo
         echo "*ERROR*: The directory '$workDir' exists, but is not readable"
         if ( $echoOn ) set echo
         exit 1
      endif
   else
      if ( $echoOn ) unset echo
      echo "*ERROR*: The pathname '$workDir' exists, but is not a directory"
      if ( $echoOn ) set echo
      exit 1
   endif
else
   mkdir -p $workDir         || exit 1
   mkdir -p $workDir/INPUT   || exit 1
   mkdir -p $workDir/RESTART || exit 1
endif

cd $workDir

set dataFilesNotOK = ( )

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/MOM_input_P39_init $workDir/INPUT/MOM_input_P39_init
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/MOM_input_P39_init )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/salt_restore_PHC2.1degOM4.v20180716.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/salt_restore_PHC2.1degOM4.v20180716.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/MOM_saltrestore $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/MOM_saltrestore )
  endif
  
  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/tidal_amplitude.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/tidal_amplitude.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/MOM_channels_SPEAR $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/MOM_channels_SPEAR )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/seawifs_1998-2006_smoothed_2X.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/seawifs_1998-2006_smoothed_2X.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/hycom1_75_800m.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/hycom1_75_800m.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/All_edits.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/All_edits.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/woa13_decav_ptemp_monthly_fulldepth_01.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/woa13_decav_ptemp_monthly_fulldepth_01.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/woa13_decav_s_monthly_fulldepth_01.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025.JRA/INPUT/woa13_decav_s_monthly_fulldepth_01.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/SPEAR_include/datasets/1949_2040/landuse_transitions/transitions.Hist_SSP585_1940_2040.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/SPEAR_include/datasets/1949_2040/landuse_transitions/transitions.Hist_SSP585_1940_2040.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/geothermal_davies2013_v1.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/geothermal_davies2013_v1.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/layer_coord.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/layer_coord.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025/INPUT/diag_rho2.nc -C $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/MOM6-examples/ice_ocean_SIS2/OM4_025/INPUT/diag_rho2.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/vgrid_75_2m.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/vgrid_75_2m.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/Kh_background_5e3_50to90ns.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/Kh_background_5e3_50to90ns.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/topo_edits_011818.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/topo_edits_011818.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/diag_vgrid.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/diag_vgrid.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/WOA05_pottemp_salt.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/WOA05_pottemp_salt.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/edits_012016.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/edits_012016.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/edits_013016a.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/edits_013016a.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/ocean_static.nc $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/ocean_static.nc )
  endif

  if (! -d $workDir/INPUT/) mkdir -p $workDir/INPUT/ && \
  ln -f /scratch/cimes/wg4031/inputs/SIS_input_C39_init $workDir/INPUT/.
  if ( $status != 0 ) then
    set dataFilesNotOK = ( $dataFilesNotOK /scratch/cimes/wg4031/inputs/SIS_input_C39_init )
  endif

if ( $#dataFilesNotOK > 0 ) then
   if ( $echoOn ) unset echo
   foreach dataFile ( $dataFilesNotOK )
      echo "*ERROR*: A problem with the data file: $dataFile"
   end
   echo "*ERROR*: Failed to copy data files"
   if ( $echoOn ) set echo
   exit 1
endif

@ gridSpecStatus = 0

tar -xf $gridSpec -C $workDir/INPUT
@ gridSpecStatus = $status

if ( $gridSpecStatus ) then
   if ( $echoOn ) unset echo
   echo "*ERROR*: Failed to copy grid specification"
   if ( $echoOn ) set echo
   exit 1
endif

unset gridSpecStatus

@ initCondStatus = 0

tar -xf $initCond -C $workDir/INPUT
@ initCondStatus = $status

if ( $initCondStatus ) then
   if ( $echoOn ) unset echo
   echo "*ERROR*: Failed to copy initial conditions"
   if ( $echoOn ) set echo
   exit 1
endif

unset initCondStatus

set domains_stack_size = "2097152"
      
set domains_stack_size = "3097152"
#set ocn_mask_table = "MOM_mask_table"

cat > $workDir/INPUT/MOM_layout << LAYOUT_EOF
#override LAYOUT = $ocn_layout
#override IO_LAYOUT = $ocn_io_layout
MASKTABLE = $ocn_mask_table
LAYOUT_EOF

cat > $workDir/INPUT/SIS_layout << LAYOUT_EOF
#override LAYOUT = $ice_layout
#override IO_LAYOUT = $ice_io_layout
MASKTABLE = $ocn_mask_table
LAYOUT_EOF

cd $workDir

if ( $echoOn ) unset echo
ls -l INPUT/*
if ( $echoOn ) set echo

ln -f $executable . || cp -pf $executable .

if ( $status == 0 ) then
   if ( $echoOn ) unset echo
   echo "<NOTE> : Using the executable '$executable'"
   if ( $echoOn ) set echo
else
   if ( $echoOn ) unset echo
   echo "*ERROR*: Failed to copy the executable"
   if ( $echoOn ) set echo
   exit 1
endif

cat >> diag_table <<EOF
OM4_hybrid_seaiceML
2018 1 1 0 0 0
     "grid_spec",            -1,  "months",  1, "days", "time"
"ice_static",           -1, "months",  1, "days", "time"
"ice_daily",             24, "hours",  1, "days", "time"
"ice_month",              1, "months", 1, "days", "time",
"ocean_static",          -1, "months", 1, "days", "time"
"ocean_month",            1, "months", 1, "days", "time"
#
#===============================================================================
# ICE DIAGNOSTICS:
#=============================================================================== 
 "ice_model", "CELL_AREA",  "CELL_AREA",    "ice_static", "all", "none", "none", 2
 "ice_model", "COSROT",     "COSROT",       "ice_static", "all", "none", "none", 2
 "ice_model", "GEOLAT",     "GEOLAT",       "ice_static", "all", "none", "none", 2
 "ice_model", "GEOLON",     "GEOLON",       "ice_static", "all", "none", "none", 2
 "ice_model", "SINROT",     "SINROT",       "ice_static", "all", "none", "none", 2
# Daily sea-ice
 "ice_model", "SST",        "SST",          "ice_daily", "all", "mean", "none", 2
 "ice_model", "SSH",        "SSH",          "ice_daily", "all", "mean", "none", 2
 "ice_model", "siconc",     "siconc",       "ice_daily", "all", "mean", "none", 2
 "ice_model", "sithick",    "sithick",      "ice_daily", "all", "mean", "none", 2
 "ice_model", "dCN",        "dCN",          "ice_daily", "all", "mean", "none", 2
# Monthly sea-ice
 "ice_model", "CELL_AREA",  "CELL_AREA",    "ice_month", "all", "none", "none", 2
 "ice_model", "COSROT",     "COSROT",       "ice_month", "all", "none", "none", 2
 "ice_model", "GEOLAT",     "GEOLAT",       "ice_month", "all", "none", "none", 2
 "ice_model", "GEOLON",     "GEOLON",       "ice_month", "all", "none", "none", 2
 "ice_model", "SINROT",     "SINROT",       "ice_month", "all", "none", "none", 2
 "ice_model", "ALB",        "ALB",          "ice_month", "all", "mean", "none", 2
 "ice_model", "BHEAT",      "BHEAT",        "ice_month", "all", "mean", "none", 2
 "ice_model", "BMELT",      "BMELT",        "ice_month", "all", "mean", "none", 2
 "ice_model", "BSNK",       "BSNK",         "ice_month", "all", "mean", "none", 2
 "ice_model", "CALVING",    "CALVING",      "ice_month", "all", "mean", "none", 2
 "ice_model", "CALVING_HFLX","CALVING_HFLX","ice_month", "all", "mean", "none", 2
 "ice_model", "CN",         "CN",           "ice_month", "all", "mean", "none", 2
 "ice_model", "E2MELT",     "E2MELT",       "ice_month", "all", "mean", "none", 2
 "ice_model", "EVAP",       "EVAP",         "ice_month", "all", "mean", "none", 2
 "ice_model", "EXT",        "EXT",          "ice_month", "all", "mean", "none", 2
 "ice_model", "EXT",        "EXT_MIN",      "ice_month", "all", "min",  "none", 2
 "ice_model", "EXT",        "EXT_MAX",      "ice_month", "all", "max",  "none", 2
 "ice_model", "FA_X",       "FA_X",         "ice_month", "all", "mean", "none", 2
 "ice_model", "FA_Y",       "FA_Y",         "ice_month", "all", "mean", "none", 2
 "ice_model", "FI_X",       "FI_X",         "ice_month", "all", "mean", "none", 2
 "ice_model", "FI_Y",       "FI_Y",         "ice_month", "all", "mean", "none", 2
 "ice_model", "FRAZIL",     "FRAZIL",       "ice_month", "all", "mean", "none", 2
 "ice_model", "IX_TRANS",   "IX_TRANS",     "ice_month", "all", "mean", "none", 2
 "ice_model", "IY_TRANS",   "IY_TRANS",     "ice_month", "all", "mean", "none", 2
 "ice_model", "LH",         "LH",           "ice_month", "all", "mean", "none", 2
 "ice_model", "LSNK",       "LSNK",         "ice_month", "all", "mean", "none", 2
 "ice_model", "LSRC",       "LSRC",         "ice_month", "all", "mean", "none", 2
 "ice_model", "LW",         "LW",           "ice_month", "all", "mean", "none", 2
 "ice_model", "RAIN",       "RAIN",         "ice_month", "all", "mean", "none", 2
 "ice_model", "RUNOFF",     "RUNOFF",       "ice_month", "all", "mean", "none", 2
 "ice_model", "SALTF",      "SALTF",        "ice_month", "all", "mean", "none", 2
 "ice_model", "SH",         "SH",           "ice_month", "all", "mean", "none", 2
 "ice_model", "SNOWFL",     "SNOWFL",       "ice_month", "all", "mean", "none", 2
 "ice_model", "SN2IC",      "SN2IC",        "ice_month", "all", "mean", "none", 2
 "ice_model", "SSH",        "SSH",          "ice_month", "all", "mean", "none", 2
 "ice_model", "SSS",        "SSS",          "ice_month", "all", "mean", "none", 2
 "ice_model", "SST",        "SST",          "ice_month", "all", "mean", "none", 2
 "ice_model", "SW",         "SW",           "ice_month", "all", "mean", "none", 2
 "ice_model", "TMELT",      "TMELT",        "ice_month", "all", "mean", "none", 2
 "ice_model", "TSN",        "TSN",          "ice_month", "all", "mean", "none", 2
 "ice_model", "TS",         "TS",           "ice_month", "all", "mean", "none", 2
 "ice_model", "T1",         "T1",           "ice_month", "all", "mean", "none", 2
 "ice_model", "T2",         "T2",           "ice_month", "all", "mean", "none", 2
 "ice_model", "T3",         "T3",           "ice_month", "all", "mean", "none", 2
 "ice_model", "T4",         "T4",           "ice_month", "all", "mean", "none", 2
 "ice_model", "UO",         "UO",           "ice_month", "all", "mean", "none", 2
 "ice_model", "VO",         "VO",           "ice_month", "all", "mean", "none", 2
 "ice_model", "XPRT",       "XPRT",         "ice_month", "all", "mean", "none", 2
 "ice_model", "siu",        "siu",          "ice_month", "all", "mean", "none", 2
 "ice_model", "siv",        "siv",          "ice_month", "all", "mean", "none", 2
 "ice_model", "sispeed",    "sispeed",      "ice_month", "all", "mean", "none", 2
 "ice_model", "STRENGTH_hf","STRENGTH_hf",  "ice_month", "all", "mean", "none", 2
 "ice_model", "sitimefrac", "sitimefrac",   "ice_month", "all", "mean", "none", 2
 "ice_model", "sitemptop",  "sitemptop",    "ice_month", "all", "mean", "none", 2
 "ice_model", "siconc",     "siconc",       "ice_month", "all", "mean", "none", 2
 "ice_model", "dCN",        "dCN",          "ice_month", "all", "mean", "none", 2
 "ice_model", "sisnconc",   "sisnconc",     "ice_month", "all", "mean", "none", 2
 "ice_model", "simass",     "simass",       "ice_month", "all", "mean", "none", 2
 "ice_model", "sisnmass",   "sisnmass",     "ice_month", "all", "mean", "none", 2
 "ice_model", "sisnthick",  "sisnthick",    "ice_month", "all", "mean", "none", 2
 "ice_model", "sithick",    "sithick",      "ice_month", "all", "mean", "none", 2
 "ice_model", "sivol",      "sivol",        "ice_month", "all", "mean", "none", 2
 "ice_model", "MIB",        "MIB",          "ice_month", "all", "mean", "none", 2
"icebergs", "BERG_MELT",   "BERG_MELT",  "iceberg_month", "all", .true., "none", 2
"icebergs", "MASS",        "MASS",       "iceberg_month", "all", .true., "none", 2
"icebergs", "MELT",        "MELT",       "iceberg_month", "all", .true., "none", 2
#
#===============================================================================
# OCEAN DIAGNOSTICS:
#=============================================================================== 
"ocean_model", "Coriolis",               "Coriolis",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "area_t",                 "area_t",                 "ocean_static",        "all", "none", "none", 2
"ocean_model", "areacello",              "areacello",              "ocean_static",        "all", "none", "none", 2
"ocean_model", "areacello_bu",           "areacello_bu",           "ocean_static",        "all", "none", "none", 2
"ocean_model", "areacello_cu",           "areacello_cu",           "ocean_static",        "all", "none", "none", 2
"ocean_model", "areacello_cv",           "areacello_cv",           "ocean_static",        "all", "none", "none", 2
"ocean_model", "depth_ocean",            "depth_ocean",            "ocean_static",        "all", "none", "none", 2
"ocean_model", "deptho",                 "deptho",                 "ocean_static",        "all", "none", "none", 2
"ocean_model", "dxCu",                   "dxCu",                   "ocean_static",        "all", "none", "none", 2
"ocean_model", "dxCv",                   "dxCv",                   "ocean_static",        "all", "none", "none", 2
"ocean_model", "dxt",                    "dxt",                    "ocean_static",        "all", "none", "none", 2
"ocean_model", "dyCu",                   "dyCu",                   "ocean_static",        "all", "none", "none", 2
"ocean_model", "dyCv",                   "dyCv",                   "ocean_static",        "all", "none", "none", 2
"ocean_model", "dyt",                    "dyt",                    "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolat",                 "geolat",                 "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolat_c",               "geolat_c",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolat_u",               "geolat_u",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolat_v",               "geolat_v",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolon",                 "geolon",                 "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolon_c",               "geolon_c",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolon_u",               "geolon_u",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "geolon_v",               "geolon_v",               "ocean_static",        "all", "none", "none", 2
"ocean_model", "sftof",                  "sftof",                  "ocean_static",        "all", "none", "none", 2
"ocean_model", "wet",                    "wet",                    "ocean_static",        "all", "none", "none", 2
"ocean_model", "wet_c",                  "wet_c",                  "ocean_static",        "all", "none", "none", 2
"ocean_model", "wet_u",                  "wet_u",                  "ocean_static",        "all", "none", "none", 2
"ocean_model", "wet_v",                  "wet_v",                  "ocean_static",        "all", "none", "none", 2
"ocean_model", "Heat_PmE",                 "Heat_PmE",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "LW",                       "LW",                         "ocean_month",     "all", "mean", "none",2
"ocean_model", "LwLatSens",                "LwLatSens",                  "ocean_month",     "all", "mean", "none",2
"ocean_model", "MLD_003",                  "MLD_003",                    "ocean_month",     "all", "mean", "none",2
"ocean_model", "MLD_003",                  "MLD_003_max",                "ocean_month",     "all", "max",  "none",2
"ocean_model", "MLD_003",                  "MLD_003_min",                "ocean_month",     "all", "min",  "none",2
"ocean_model", "ML_buoy_restrat",          "ML_buoy_restrat",            "ocean_month",     "all", "mean", "none",2
"ocean_model", "PRCmE",                    "PRCmE",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "SSH",                      "ssh",                        "ocean_month",     "all", "mean", "none",2
"ocean_model", "SSS",                      "sss",                        "ocean_month",     "all", "mean", "none",2
"ocean_model", "SST",                      "SST",                        "ocean_month",     "all", "mean", "none",2
"ocean_model", "SST",                      "sst_max",                    "ocean_month",     "all", "max",  "none",2
"ocean_model", "SST",                      "sst_min",                    "ocean_month",     "all", "min",  "none",2
"ocean_model", "SW",                       "SW",                         "ocean_month",     "all", "mean", "none",2
"ocean_model", "S_adx_2d",                 "S_adx_2d",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "S_ady_2d",                 "S_ady_2d",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "TKE_tidal",                "TKE_tidal",                  "ocean_month",     "all", "mean", "none",2
"ocean_model", "T_adx_2d",                 "T_adx_2d",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "T_ady_2d",                 "T_ady_2d",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "evap",                     "evap",                       "ocean_month",     "all", "mean", "none",2
"ocean_model", "evs",                      "evs",                        "ocean_month",     "all", "mean", "none",2
"ocean_model", "ficeberg",                 "ficeberg",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "fprec",                    "fprec",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "friver",                   "friver",                     "ocean_month",     "all", "mean", "none",2
"ocean_model", "frunoff",                  "frunoff",                    "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_cond",        "heat_content_cond",          "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_fprec",       "heat_content_fprec",         "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_frunoff",     "heat_content_frunoff",       "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_lprec",       "heat_content_lprec",         "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_lrunoff",     "heat_content_lrunoff",       "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_massin",      "heat_content_massin",        "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_massout",     "heat_content_massout",       "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_content_surfwater",   "heat_content_surfwater",     "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfds",                     "hfds",                       "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfevapds",                 "hfevapds",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfibthermds",              "hfibthermds",                "ocean_month",     "all", "mean", "none",2
"ocean_model", "hflso",                    "hflso",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfrainds",                 "hfrainds",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfrunoffds",               "hfrunoffds",                 "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfsifrazil",               "hfsifrazil",                 "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfsnthermds",              "hfsnthermds",                "ocean_month",     "all", "mean", "none",2
"ocean_model", "hfsso",                    "hfsso",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "latent",                   "latent",                     "ocean_month",     "all", "mean", "none",2
"ocean_model", "lprec",                    "lprec",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "lrunoff",                  "lrunoff",                    "ocean_month",     "all", "mean", "none",2
"ocean_model", "net_heat_coupler",         "net_heat_coupler",           "ocean_month",     "all", "mean", "none",2
"ocean_model", "net_heat_surface",         "net_heat_surface",           "ocean_month",     "all", "mean", "none",2
"ocean_model", "net_massin",               "net_massin",                 "ocean_month",     "all", "mean", "none",2
"ocean_model", "net_massout",              "net_massout",                "ocean_month",     "all", "mean", "none",2
"ocean_model", "p_surf",                   "p_surf",                     "ocean_month",     "all", "mean", "none",2
"ocean_model", "prlq",                     "prlq",                       "ocean_month",     "all", "mean", "none",2
"ocean_model", "prsn",                     "prsn",                       "ocean_month",     "all", "mean", "none",2
"ocean_model", "rlntds",                   "rlntds",                     "ocean_month",     "all", "mean", "none",2
"ocean_model", "rsntds",                   "rsntds",                     "ocean_month",     "all", "mean", "none",2
"ocean_model", "salt_flux",                "salt_flux",                  "ocean_month",     "all", "mean", "none",2
"ocean_model", "sensible",                 "sensible",                   "ocean_month",     "all", "mean", "none",2
"ocean_model", "sfdsi",                    "sfdsi",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "speed",                    "speed",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "tauuo",                    "tauuo",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "tauvo",                    "tauvo",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "taux",                     "taux",                       "ocean_month",     "all", "mean", "none",2
"ocean_model", "tauy",                     "tauy",                       "ocean_month",     "all", "mean", "none",2
"ocean_model", "ustar",                    "ustar",                      "ocean_month",     "all", "mean", "none",2
"ocean_model", "wfo",                      "wfo",                        "ocean_month",     "all", "mean", "none",2
"ocean_model", "heat_added",               "heat_added",                 "ocean_month",     "all", "mean", "none",2
"ocean_model", "salt_flux_added",          "salt_flux_added",            "ocean_month",     "all", "mean", "none",2
EOF

#FRE table(diag_table.yaml)                                                                                                                                                                                      
cat >> field_table <<EOF                                                                                                                                                                                         
# specific humidity for moist runs                                                                                                                                                                               
 "TRACER", "atmos_mod", "sphum"                                                                                                                                                                                  
           "longname",     "specific humidity"                                                                                                                                                                   
           "units",        "kg/kg" /                                                                                                                                                                             
##         "profile_type", "fixed", "surface_value=3.e-6" /                                                                                                                                                      
# prognostic cloud scheme tracers                                                                                                                                                                                
  "TRACER", "atmos_mod", "liq_wat"                                                                                                                                                                               
            "longname",     "cloud liquid specific humidity"                                                                                                                                                     
            "units",        "kg/kg" /                                                                                                                                                                            
  "TRACER", "atmos_mod", "ice_wat"                                                                                                                                                                               
            "longname",     "cloud ice water specific humidity"                                                                                                                                                  
            "units",        "kg/kg" /                                                                                                                                                                            
  "TRACER", "atmos_mod", "cld_amt"                                                                                                                                                                               
            "longname",     "cloud fraction"                                                                                                                                                                     
            "units",        "none" /
# sphum must be present on land as well                                                                                                                                                                          
 "TRACER", "land_mod",     "sphum"                                                                                                                                                                               
           "longname",     "specific humidity"                                                                                                                                                                   
           "units",        "kg/kg" /                                                                                                                                                                             
# test tracer for radon                                                                                                                                                                                          
#                                                                                                                                                                                                                
# "TRACER", "atmos_mod", "radon"                                                                                                                                                                                 
#           "longname",     "radon test tracer"                                                                                                                                                                  
#           "units",        "kg/kg" /                                                                                                                                                                            
###.................................................                                                                                                                                                             
EOF

#FRE table(field_table.yaml)                                                                                                                                                                                     
cat >> data_table <<EOF                                                                                                                                                                                          
"ATM", "p_surf",             "psl",      "./INPUT/JRA_psl.nc",     "bilinear",   1.0                                                                                                                             
 "ATM", "p_bot",              "psl",      "./INPUT/JRA_psl.nc",     "bilinear",   1.0                                                                                                                            
 "ATM", "t_bot",              "tas",      "./INPUT/JRA_tas.nc",     "bilinear",   1.0                                                                                                                            
 "ATM", "sphum_bot",          "huss",     "./INPUT/JRA_huss.nc",    "bilinear",   1.0                                                                                                                            
 "ATM", "u_bot",              "uas",      "./INPUT/JRA_uas.nc",     "bicubic",    1.0                                                                                                                            
 "ATM", "v_bot",              "vas",      "./INPUT/JRA_vas.nc",     "bicubic",    1.0                                                                                                                            
 "ATM", "z_bot",              "",         "",                       "bilinear",  10.0                                                                                                                            
 "ATM", "gust",               "",         "",                       "bilinear",   1.0e-4                                                                                                                         
 "ICE", "lw_flux_dn",         "rlds",     "./INPUT/JRA_rlds.nc",    "bilinear",   1.0                                                                                                                            
 "ICE", "sw_flux_vis_dir_dn", "rsds",     "./INPUT/JRA_rsds.nc",    "bilinear",   0.285                                                                                                                          
 "ICE", "sw_flux_vis_dif_dn", "rsds",     "./INPUT/JRA_rsds.nc",    "bilinear",   0.285                                                                                                                          
 "ICE", "sw_flux_nir_dir_dn", "rsds",     "./INPUT/JRA_rsds.nc",    "bilinear",   0.215                                                                                                                          
 "ICE", "sw_flux_nir_dif_dn", "rsds",     "./INPUT/JRA_rsds.nc",    "bilinear",   0.215                                                                                                                          
 "ICE", "lprec",              "prra",     "./INPUT/JRA_prra.nc",    "bilinear",   1.0                                                                                                                            
 "ICE", "fprec",              "prsn",     "./INPUT/JRA_prsn.nc",    "bilinear",   1.0                                                                                                                            
 "ICE", "dhdt",               "",         "",                       "none",      80.0                                                                                                                            
 "ICE", "dedt",               "",         "",                       "none",       2.0e-6                                                                                                                         
 "ICE", "drdt",               "",         "",                       "none",      10.0                                                                                                                            
#JRA runoff                                                                                                                                                                                                      
"ICE" , "runoff"        , "friver"      , "./INPUT/JRA_friver_360x320.nc", "none" ,  1.0                                                                                                                         
"ICE" , "calving"       , "licalvf"     , "./INPUT/JRA_licalvf_360x320.nc", "none" ,  1.0                                                                                                                        
EOF                                                                                                                                                                                                              

#FRE table(data_table.yaml)

cat > input.nml.unexpanded <<\EOF
 &MOM_input_nml
         output_directory = './',
         input_filename = "MOM.res.ens_%2E.nc"
         restart_input_dir = 'INPUT/',
         restart_output_dir = 'RESTART/',
         parameter_filename = 'INPUT/MOM_input','INPUT/MOM_layout','INPUT/MOM_override'
/

 &SIS_input_nml
         output_directory = './',
         input_filename = '$sis_restart_flag'
         restart_input_dir = 'INPUT/',
         restart_output_dir = 'RESTART/',
         parameter_filename = 'INPUT/SIS_input','INPUT/SIS_layout','INPUT/SIS_override'
/

 &atmos_model_nml                                                                                                                                                                                                
       layout= 0, 0
 
 &coupler_nml
        months = $months,
        days   = $days,
        current_date = 2018,1,1,0,0,0,
        calendar = 'julian',
        dt_cpld = 1800,
        dt_atmos = 1800,
        do_atmos = .false.,
        do_land = .false.,
        do_ice = .true.,
        do_ocean = .true.,
        atmos_npes = 0,
	ocean_npes = 0,
        atmos_nthreads = $atm_threads
        ocean_nthreads = $ocn_threads
        concurrent = .false.
        use_lag_fluxes=.false.
/

 &diag_manager_nml
        mix_snapshot_average_fields = .false.,
        max_input_fields = 1200,
        max_output_fields = 1800
        max_axes = 400
        max_num_axis_sets = 400
        max_files = 400
        max_axis_attributes=4
/

&flux_exchange_nml                                                                                                                                                                                              
        debug_stocks = .FALSE.                                                                                                                                                                               
        divert_stocks_report = .TRUE.                                                                                                                                                                        
        do_area_weighted_flux = .FALSE.                                                                                                                                                                      
/

&fms_io_nml
 	fms_netcdf_restart=.true.
        threading_read  = 'multi',
        max_files_r = 800,
        max_files_w = 800,
/

&fms_nml
        domains_stack_size = 8000000
	stack_size = 0
        clock_grain = 'ROUTINE'
	clock_flags = 'NONE'
/

&ice_albedo_nml                                                                                                                                                                                                 
        t_range = 10.                                                                                                                                                                                        
/  

&ice_model_nml

/

&icebergs_nml                                                                                                                                                                                                   
       	make_calving_reproduce=.TRUE.                                                                                                                                                                             
	really_debug=.FALSE.                                                                                                                                                                                     
        debug=.FALSE.                                                                                                                                                                                            
        verbose=.FALSE.                                                                                                                                                                                          
        verbose_hrs=7200                                                                                                                                                                                         
        use_operator_splitting=.TRUE.                                                                                                                                                                            
        bergy_bit_erosion_fraction=0.0                                                                                                                                                                           
        sicn_shift=0.1                                                                                                                                                                                           
        parallel_reprod=.TRUE.                                                                                                                                                                                   
        traj_sample_hrs=0                                                                                                                                                                                        
        add_weight_to_ocean=.false.                                                                                                                                                                              
        tidal_drift = 0.005                                                                                                                                                                                      
        coastal_drift = 0.001
/

&monin_obukhov_nml                                                                                                                                                                                              
            neutral = .true.                                                                                                                                                                                     
/

&ocean_albedo_nml
      	ocean_albedo_option = 2
/

&ocean_rough_nml
        rough_scheme = 'beljaars'
/

&sat_vapor_pres_nml                                                                                                                                                                                             
     	construct_table_wrt_liq = .true.,                                                                                                                                                                          
     	construct_table_wrt_liq_and_ice = .true.,                                                                                                                                                                  
/                                                                                                                                                                                                                
                                                                                                                                                                                                                 
&surface_flux_nml                                                                                                                                                                                               
        ncar_ocean_flux = .true.                                                                                                                                                                             
	raoult_sat_vap = .true.                                                                                                                                                                              
/                                                                                                                                                                                                                
                                                                                                                                                                                                                 
&topography_nml                                                                                                                                                                                                 
        topog_file = 'INPUT/navy_topography.data.nc'                                                                                                                                                         
/

&xgrid_nml
        make_exchange_reproduce = .true.
        interp_method = 'second_order'
/

\EOF

set months = $monthslist[1]
set days = $dayslist[1]
set hours = $hourslist[1]
set adjust_dry_mass = `adjust_dry_mass_tool`

set | sort > $envFile
sleep $envFileDelay
set -r | sort >> $envFile
sleep $envFileDelay
env --unset=TERMCAP | grep -e '^[a-zA-Z0-9_]*=' | sort >> $envFile

expandVariables $envFile < input.nml.unexpanded > input.nml || exit 1

rm -f $envFile

if ( $?FRE_STAGE ) then
   if ( $FRE_STAGE == 'CHAIN' ) then
      if ( -f $scriptName ) then
         if ( -r $scriptName ) then
            set result = `submit $scriptName`
            if ( $status == 0 ) then
               if ( $echoOn ) unset echo
               echo "<NOTE> : The job '$result' to run the '$scriptName' has been submitted successfully"
               if ( $echoOn ) set echo
               workDirCleaner $workDir
               if ( $echoOn ) unset echo
               echo "<NOTE> : Natural end-of-input-chaining for '$scriptName'"
               if ( $echoOn ) set echo
               exit 0
            else
               workDirCleaner $workDir
               if ( $echoOn ) unset echo
               echo "*ERROR*: Can't submit the '$scriptName'"
               if ( $echoOn ) set echo
               exit 1
            endif
            unset result
         else
            workDirCleaner $workDir
            if ( $echoOn ) unset echo
            echo "*ERROR*: The script '$scriptName' exists, but is not readable - it can't be submitted"
            if ( $echoOn ) set echo
            exit 1
         endif
      else
         workDirCleaner $workDir
         if ( $echoOn ) unset echo
         echo "*ERROR*: The script '$scriptName' does not exist - it can't be submitted"
         if ( $echoOn ) set echo
         exit 1
      endif
   else if ( $FRE_STAGE == 'DEBUG' ) then
      if ( $echoOn ) unset echo
      echo "<NOTE> : The working directory '$workDir' is ready for debugging"
      echo "<NOTE> : Natural end-of-debug-staging for '$scriptName'"
      if ( $echoOn ) set echo
      exit 0
   else
      workDirCleaner $workDir
      if ( $echoOn ) unset echo
      echo "<NOTE> : Natural end-of-input-staging for '$scriptName'"
      if ( $echoOn ) set echo
      exit 0
   endif
endif

################################################################################
#------------------------------- the main loop ---------------------------------
################################################################################

while ( $irun <= $segmentsPerJob && $currentSeg <= $segmentsPerSimulation )
   if ( $echoOn ) unset echo
   echo "################################################################################"
   echo "# $currentSeg/$segmentsPerSimulation"
   echo "################################################################################"
   if ( $echoOn ) set echo

   # ---------------- reload the queue file and exit if it has been requested

   if ( $?flagRunTypeProduction ) then
      if ( -f $queue_file ) then
         if ( -r $queue_file ) then
            source $queue_file
         else
            if ( $echoOn ) unset echo
            echo "*ERROR*: The queue file '$queue_file' is not readable"
            if ( $echoOn ) set echo
            exit 1
         endif
      endif

      if ( ! $continueFlag ) then
         if ( $echoOn ) unset echo
         echo "<NOTE> : Stopping execution"
         if ( $echoOn ) set echo
         exit 0
      endif
   endif

   # ---------------- commands

#------------------------------------------                                                                                                                                                                      
## Find out whether to restart.                                                                                                                                                                                  
# MOM6 restart switch                                                                                                                                                                                            
if ( $currentSeg == 1 ) then                                                                                                                                                                                     
   set restart_flag = 'n'                                                                                                                                                                                        
else                                                                                                                                                                                                             
   set restart_flag = 'r'                                                                                                                                                                                        
endif                                                                                                                                                                                                            
   set restart_flag = 'r'                                                                                                                                                                                        
   set sis_restart_flag = 'r'                                                                                                                                                                                    
                                                                                                                                                                                                                 
   #Because of a technical issue FRE does not delete old uncombined restarts before copying the new ones for the next segment                                                                                    
   #The INPUT/ ends up with inconsistent restart files which at best cause the model to crash right away.                                                                                                        
   if(-e $workDir/INPUT/MOM.res.nc.0000 ) \rm -rf $workDir/INPUT/MOM.res*.nc                                                                                                                                           
endif                                                                                                                                                                                                            
                                                                                                                                                                                                                 
ln -s $workDir/INPUT/ocean_topog.nc $workDir/INPUT/topog.nc                                                                                                                                                            
                                                                                                                                                                                                                 
touch $workDir/INPUT/MOM_override                                                                                                                                                                                   
touch $workDir/INPUT/MOM_layout                                                                                                                                                                                     
touch $workDir/INPUT/SIS_layout  

# Record the job stdout location for later use timings database                                                                                                                                                  
cat >> /ncrc/home2/William.Gregory/frejobs_stdout <<EOF_frejobs                                                                                                                                                  
$stdoutDir/$FRE_JOBID                                                                                                                                                                                            
EOF_frejobs                                                                                                                                                                                                      
                                                                                                                                                                                                                 
cat > $work/INPUT/MOM_layout << MOM_LAYOUT_EOF                                                                                                                                                                   
#override IO_LAYOUT = $ocn_io_layout                                                                                                                                                                             
#override LAYOUT    = $ocn_layout                                                                                                                                                                                
#override MASKTABLE = $ocn_mask_table                                                                                                                                                                            
#override OCEAN_OMP_THREADS = $ocn_threads                                                                                                                                                                       
MOM_LAYOUT_EOF                                                                                                                                                                                                   
                                                                                                                                                                                                                 
cat > $work/INPUT/SIS_layout << SIS_LAYOUT_EOF                                                                                                                                                                   
#override IO_LAYOUT = $ice_io_layout                                                                                                                                                                             
#override LAYOUT    = $ice_layout                                                                                                                                                                                
#override MASKTABLE = $ice_mask_table                                                                                                                                                                            
SIS_LAYOUT_EOF 

#Note: No SIS_diurnal hence                                                                                                                                                                                      
#ADD_DIURNAL_SW = False                                                                                                                                                                                          
#since JRA forcings has a diurnal cycle                                                                                                                                                                          
                                                                                                                                                                                                                 
#######IAF cycle mechanism                                                                                                                                                                                       
echo "Model year = $fyear"                                                                                                                                                                                       
#Current JRA dataset starts at 1958 and ends at 2018                                                                                                                                                             
#These numbers should be adjusted when the datasets start or length changes.                                                                                                                                     
set JRA_START_YEAR = 1958                                                                                                                                                                                        
#set JRA_LEN = 64                                                                                                                                                                                                
#@ modulyr = ( $fyear - $JRA_START_YEAR ) % $JRA_LEN                                                                                                                                                             
#@ forceyr = $JRA_START_YEAR + $modulyr                                                                                                                                                                          
#if ( $modulyr < 0 ) then                                                                                                                                                                                        
#   @ forceyr = $forceyr + $JRA_LEN                                                                                                                                                                              
#endif                                                                                                                                                                                                           
#The above logic is not needed for cycles that start at 1958 (set by xml property fyear and reset in coupler_nml:current_date)                                                                                   
@ forceyr = $fyear                                                                                                                                                                                               
echo "Forcings file year = $forceyr"                                                                                                                                                                             
@ forceyrp1 = $forceyr + 1

cd $workDir/INPUT/                                                                                                                                                                                                  
set fetch_cmd = 'ln -sf ' #This might be the cause of frequent crashes with HDF errors?!                                                                                                                         
$fetch_cmd /scratch/cimes/wg4031/JRA/huss_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010000-${forceyr}12312100.padded.nc JRA_huss.nc                                                                                                                                                                                                         
$fetch_cmd /scratch/cimes/wg4031/JRA/prra_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010130-${forceyr}12312230.padded.nc JRA_prra.nc                                                                                                                                                                                                         
$fetch_cmd /scratch/cimes/wg4031/JRA/prsn_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010130-${forceyr}12312230.padded.nc JRA_prsn.nc                                                                                                                                                                                                         
$fetch_cmd /scratch/cimes/wg4031/JRA/psl_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010000-${forceyr}12312100.padded.nc JRA_psl.nc                                                                                                                                                                                                           
$fetch_cmd /scratch/cimes/wg4031/JRA/rlds_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010130-${forceyr}12312230.padded.nc JRA_rlds.nc                                                                                                                                                                                                         
$fetch_cmd /scratch/cimes/wg4031/JRA/rsds_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010130-${forceyr}12312230.padded.nc JRA_rsds.nc                                                                                                                                                                                                         
$fetch_cmd /scratch/cimes/wg4031/JRA/tas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010000-${forceyr}12312100.padded.nc JRA_tas.nc                                                                                                                                                                                                           
$fetch_cmd /scratch/cimes/wg4031/JRA/uas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010000-${forceyr}12312100.padded.nc JRA_uas.nc                                                                                                                                                                                                           
$fetch_cmd /scratch/cimes/wg4031/JRA/vas_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}01010000-${forceyr}12312100.padded.nc JRA_vas.nc                                                                                                                                                                                                           
#                                                                                                                                                                                                                
$fetch_cmd /scratch/cimes/wg4031/JRA/friver_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}0101-${forceyr}1231.padded.regrid360x320.nc JRA_friver_360x320.nc    
$fetch_cmd /scratch/cimes/wg4031/JRA/licalvf_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_${forceyr}0101-${forceyr}1231.padded.regrid360x320.nc JRA_licalvf_360x320.nc  
#                                                                                                                                                                                                                
$fetch_cmd /scratch/cimes/wg4031/OISST/sst_oidaily_v2p1_icecorr_icec30_tripolar_${forceyr}.nc temp_restore.nc                                                                             
cd $workDir


cd INPUT

# adjustment of dry mass first time only
set adjust_dry_mass = `/scratch/cimes/wg4031/fre_scripts/adjust_dry_mass.csh`

# create dummy MOM6 parameter file
touch MOM_input
cd ..


      

cat > ocean_obs_table <<EOF
EOF

######BEGIN MOM6 and SIS2 layout setup################
##THIS HAS TO BE PART OF csh type=always###

cat > $workDir/INPUT/MOM_override_00 << EOF
#override COORD_DEF = "FILE:vgrid_oda.nc,dz"
#override COORD_FILE = "coord_oda.nc"
#override REGRIDDING_COORDINATE_MODE = "Z*"
#override ALE_COORDINATE_CONFIG = "FILE:vgrid_oda.nc,dz"
#override NK = 75
#override MAXIMUM_INT_DEPTH_CONFIG = "NONE"
#override MAX_LAYER_THICKNESS_CONFIG = "NONE"
NJHALO_ODA = 12
NIHALO_ODA = 8
#override REMAPPING_SCHEME = "PPM_H4"

! === module ODA ===
ASSIM_METHOD = 'NONE'
ASSIM_FREQUENCY = 24
BASIN_FILE = "basin.nc"         ! default = "basin.nc"

#override TOPO_EDITS_FILE = "topo_edits_011818.nc"
#override CHANNEL_LIST_FILE = MOM_channels_SPEAR
RESTART_CONTROL = 2
NUM_DIAG_COORDS = 2             ! default = 1
                                ! The number of diagnostic vertical coordinates to use.
                                ! For each coordinate, an entry in DIAG_COORDS must be provided.
DIAG_COORDS = "z Z ZSTAR", "rho2 RHO2 RHO" !

DIAG_COORD_DEF_RHO2 = "RFNC1:35,999.5,1028,1028.5,8.,1038.,0.0078125" ! default = "WOA09"
                                ! Determines how to specify the coordinate
                                ! resolution. Valid options are:
                                !  PARAM       - use the vector-parameter DIAG_COORD_RES_RHO2
                                !  UNIFORM[:N] - uniformly distributed
                                !  FILE:string - read from a file. The string specifies
                                !                the filename and variable name, separated
                                !                by a comma or space, e.g. FILE:lev.nc,dz
                                !                or FILE:lev.nc,interfaces=zw
                                !  WOA09[:N]   - the WOA09 vertical grid (approximately)
                                !  FNC1:string - FNC1:dz_min,H_total,power,precision
                                !  HYBRID:string - read from a file. The string specifies
                                !                the filename and two variable names, separated
                                !                by a comma or space, for sigma-2 and dz. e.g.
                                !                HYBRID:vgrid.nc,sigma2,dz

EOF

cat > $workDir/INPUT/MOM_override_01 << EOF
#override TOPO_EDITS_FILE = "topo_edits_011818.nc"
#override CHANNEL_LIST_FILE = MOM_channels_SPEAR
#override MIXEDLAYER_RESTRAT = False
ENSEMBLE_OCEAN = True
RESTART_CONTROL = 2
NUM_DIAG_COORDS = 3             ! default = 1
                                ! The number of diagnostic vertical coordinates to use.
                                ! For each coordinate, an entry in DIAG_COORDS must be provided.
DIAG_COORDS = "z Z ZSTAR, rho2 RHO2 RHO, ztop ZTOP300 ZSTAR"

DIAG_COORD_DEF_RHO2 = "RFNC1:35,999.5,1028,1028.5,8.,1038.,0.0078125" ! default = "WOA09"
DIAG_COORD_DEF_ZTOP300 =  "UNIFORM:2,600."
                                ! Determines how to specify the coordinate
                                ! resolution. Valid options are:
                                !  PARAM       - use the vector-parameter DIAG_COORD_RES_RHO2
                                !  UNIFORM[:N] - uniformly distributed
                                !  FILE:string - read from a file. The string specifies
                                !                the filename and variable name, separated
                                !                by a comma or space, e.g. FILE:lev.nc,dz
                                !                or FILE:lev.nc,interfaces=zw
                                !  WOA09[:N]   - the WOA09 vertical grid (approximately)
                                !  FNC1:string - FNC1:dz_min,H_total,power,precision
                                !  HYBRID:string - read from a file. The string specifies
                                !                the filename and two variable names, separated
                                !                by a comma or space, for sigma-2 and dz. e.g.
                                !                HYBRID:vgrid.nc,sigma2,dz

EOF

set ENS_I = 2
while ( $ENS_I <= 9 )
	cp $workDir/INPUT/MOM_override_01 $workDir/INPUT/MOM_override_0$ENS_I
	@ ENS_I++
end
set ENS_I = 10
while ( $ENS_I <= 30 )
	cp $workDir/INPUT/MOM_override_01 $workDir/INPUT/MOM_override_$ENS_I
	@ ENS_I++
end

touch $workDir/INPUT/SIS_override

cat > $workDir/INPUT/SIS_override << EOF
#override CP_SEAWATER = 3992.
#override CP_BRINE = 3992.
#override ICE_BULK_SALINITY = 0.0
#override ICE_RELATIVE_SALINITY = 0.17
#override SIS2_FILLING_FRAZIL = T
#override SIS_THICKNESS_ADVECTION_SCHEME = "PCM"
#override SIS_CONTINUITY_SCHEME = "PCM"
#override SIS_TRACER_ADVECTION_SCHEME = "PPM:H3"
EOF

#scan auxiliary_nml
#set OCN_MASK_TABLE =
#set OCN_MASK_TABLE = `awk '/ocean_mask_table/ {gsub(/=/," ");gsub(/,/," ");print $2}' $workDir/input.nml`
#set OCN_LAYOUT     = `awk '/ocean_layout/     {gsub(/=/," ");           print $2,$3}' $workDir/input.nml`
#set OCN_IO_LAYOUT  = `awk '/ocean_io_layout/  {gsub(/=/," ");           print $2,$3}' $workDir/input.nml`

cat > $workDir/INPUT/MOM_layout << MOM_LAYOUT_EOF
#override IO_LAYOUT = $ocn_io_layout
#override LAYOUT    = $ocn_layout
#override MASKTABLE = $ocn_mask_table
MOM_LAYOUT_EOF

#scan auxiliary_nml
#set ICE_MASK_TABLE = `awk '/ice_mask_table/ {gsub(/=/," ");gsub(/,/," ");print $2}' $workDir/input.nml`
#set ICE_LAYOUT     = `awk '/ice_layout/     {gsub(/=/," ");           print $2,$3}' $workDir/input.nml`
#set ICE_IO_LAYOUT  = `awk '/ice_io_layout/  {gsub(/=/," ");           print $2,$3}' $workDir/input.nml`

cat > $workDir/INPUT/SIS_layout << SIS_LAYOUT_EOF
#override IO_LAYOUT = $ice_io_layout
#override LAYOUT    = $ice_layout
#override MASKTABLE = $ocn_mask_table
SIS_LAYOUT_EOF

######END MOM6 and SIS2 layout setup################

cd $workDir/INPUT

if ( $currentSeg == 1 ) then
     rm -f coupler.res
endif
cd $workDir

#------------------------------------------
## Find out whether to restart.
# MOM6 restart switch
if ( $currentSeg == 1 ) then
   set restart_flag = 'n'
else
   set restart_flag = 'r'
endif
   set sis_restart_flag = 'r'

cd $workDir

## Run the model
#------------------------------------------
###The following vars are set for the timing database ingestion tool only and are not necessary for a successful run
set OCN_MASK_TABLE  =  $ocn_mask_table
set OCN_LAYOUT    = $ocn_layout
set OCN_IO_LAYOUT = $ocn_io_layout
set ICE_LAYOUT    = $ice_layout
set ICE_IO_LAYOUT = $ice_io_layout
set atmos_npes = $atm_ranks
set ocean_npes = $ocn_ranks
set fv_layout = $atm_layout
set fv_io_layout  = $atm_io_layout
set land_layout   = $lnd_layout
set land_io_layout = $lnd_io_layout
set dummy  = `awk '/nxblocks/ {gsub(/=/," ");gsub(/,/," ");print $2}' $workDir/input.nml`
set nxblocks = $dummy
set dummy  = `awk '/nyblocks/ {gsub(/=/," ");gsub(/,/," ");print $2}' $workDir/input.nml`
set nyblocks = $dummy

   cd $workDir

   # ---------------- expand namelists
    
   set months = $monthslist[$irun]
   set days = $dayslist[$irun]
   set hours = $hourslist[$irun]
   set adjust_dry_mass = `adjust_dry_mass_tool`

   set | sort > $envFile
   sleep $envFileDelay
   set -r | sort >> $envFile
   sleep $envFileDelay
   env --unset=TERMCAP | grep -e '^[a-zA-Z0-9_]*=' | sort >> $envFile

   expandVariables $envFile < input.nml.unexpanded > input.nml || exit 1

   rm -f $envFile

   # ---------------- prepare MPI call, execute it, analyze results

   unsetenv OMP_NUM_THREADS

   echo "Time before runCommand"
   date

   runCommand |& tee fms.out

   if ( $status == 0 ) then
      if ( $echoOn ) unset echo
      echo "<NOTE> : The MPI launcher (srun) exited normally"
      if ( $echoOn ) set echo
   else if ( $status == 1 ) then
      set msg =       "*ERROR*: Automatic message from the job '$FRE_JOBID'\n"
      set msg = "${msg}*ERROR*: -----------------------------------------------------------------------\n"
      set msg = "${msg}*ERROR*: The MPI launcher (srun) exited with error status\n"
      set msg = "${msg}*ERROR*: \n"
      set msg = "${msg}*ERROR*: Possible Reasons: incorrect srun options (for example more cores specified\n"
      set msg = "${msg}*ERROR*: than available), node failure or untrapped srun error.\n"
      set msg = "${msg}*ERROR*: \n"
      set msg = "${msg}*ERROR*: Please see the job stdout, located at:\n"
      set msg = "${msg}*ERROR*: \n"
      set msg = "${msg}*ERROR*: \t$stdoutDir/$FRE_JOBID\n"
      set msg = "${msg}*ERROR*: \n"

      set MPI_FAIL
   else
      set msg =       "*ERROR*: Automatic message from the job '$FRE_JOBID'\n"
      set msg = "${msg}*ERROR*: -----------------------------------------------------------------------\n"
      set msg = "${msg}*ERROR*: The MPI launcher (srun) exited abnormally\n"
      set msg = "${msg}*ERROR*: \n"
      set msg = "${msg}*ERROR*: Possible Reasons: job cancelled or job ended through MPI_Abort or segfault.\n"
      set msg = "${msg}*ERROR*: \n"
      set msg = "${msg}*ERROR*: Please see the job stdout, located at:\n"
      set msg = "${msg}*ERROR*: \n"
      set msg = "${msg}*ERROR*: \t$stdoutDir/$FRE_JOBID\n"
      set msg = "${msg}*ERROR*: \n"

      set MPI_FAIL
   endif

   if ( $?MPI_FAIL ) then
      set coreFiles = ( `ls core*` )

      if ( $#coreFiles > 0 ) then
         if ( ! $?MPI_COREDUMP_DEBUGGER ) setenv MPI_COREDUMP_DEBUGGER 'gdb -batch'
         echo 'where' > .gdbinit

         set coreFileFirst = $coreFiles[1]
         echo "Dump of the core file '$coreFileFirst'" > $coreFileFirst.out
         $MPI_COREDUMP_DEBUGGER $executable:t $coreFileFirst >> $coreFileFirst.out
         cat $coreFileFirst.out >> fms.out
         cat $coreFileFirst.out
         unset coreFileFirst

         set msg = "${msg}*ERROR*: Your job has produced $#coreFiles core files (segment $currentSeg)\n"
         set msg = "${msg}*ERROR*: Please go to the working directory '$workDir' and issue the following command for each core file there:\n"
         set msg = "${msg}*ERROR*: \n"

         @ count = 0
         @ countMax = 20

         foreach coreFile ( $coreFiles )
            set msg = "${msg}*ERROR*: \t$MPI_COREDUMP_DEBUGGER $executable:t $coreFile >> $coreFile.out\n"
            if ( $count < $countMax ) then
               @ count++
            else
               break
            endif
         end

         set msg = "${msg}*ERROR*: \n"
         set msg = "${msg}*ERROR*: FRE has executed the command above for one core file and echoed the result to the job stdout.\n"

         if ( $count == $countMax ) then
            set msg = "${msg}*ERROR*: In order to save space only the first $countMax core files are presented in this list.\n"
            set msg = "${msg}*ERROR*: \n"
         endif

         unset countMax
         unset count
      else
         set cdsize = `limit coredumpsize`
         set msg = "${msg}*ERROR*: No core files produced (segment $currentSeg)\n"
         set msg = "${msg}*ERROR*: You are using the $cdsize\n"
         set msg = "${msg}*ERROR*: \n"
         unset cdsize
      endif

      set msg = "${msg}*ERROR*: -----------------------------------------------------------------------\n"
      set msg = "${msg}*ERROR*: This message has been generated by FRE\n"
      set msg = "${msg}*ERROR*: `date`"

      if ( $?batch ) then
         if ( $echoOn ) unset echo
         printf "$msg" | mailx -s "The MPI launcher has failed" $mailList
         printf "$msg"
         if ( $echoOn ) set echo
      endif

      unset coreFiles
      unset msg

      set outputDir = ${outputDir}_crash
   endif

   echo "Time after runCommand"
   date

   # ---------------- generate date for file names

   set begindate = `timeStamp -b`
   if ( $begindate == 'no_time_stamp' ) set begindate = tmp`date '+%j%H%M%S'`
   set enddate = `timeStamp -e`
   if ( $enddate == 'no_time_stamp' ) set enddate = tmp`date '+%j%H%M%S'`
   set fyear = `echo $enddate | timeStamp -y`

   # ---------------- commands, copied from XML (experiment/postProcess/csh)

         cd $workDir/RESTART
         ncrcat icebergs.res.nc.* icebergs.res.nc
         if(-e icebergs.res.nc) rm -f icebergs.res.nc.*
         cd $workDir
         ncrcat iceberg_trajectories.nc.* iceberg_trajectories.nc
         rm -f iceberg_trajectories.nc.*
         #Save the last line of timestats without the first (step number) column as a signature of MOM6 answers
         tail -1 timestats | cut -c8- > RESTART/$enddate.timestats.res
         #Save the whole timestats under ascii/
         mv timestats $enddate.timestats
         mv timestats.nc $enddate.timestats.nc

         #Make a directory to trick FRE to pick up and archive in ascii
         mkdir -p extra.results

         mv *velocity_truncations CPU_stats SIS_parameter_doc* MOM_parameter_doc* extra.results/         

   cd $workDir

   # ---------------- remove time_stamp.out file

   if ( -f time_stamp.out ) rm -f time_stamp.out

   # ---------------- save ascii files

   set asciiFiles = ( `ls -1 | egrep "$patternGrepAscii"` )

   if ( $#asciiFiles > 0 ) then
      set asciiSuffix = ascii/$begindate.ascii_out
      set asciiBackup = $begindate.ascii_out.tar
      set asciiArchDir = $outputDir/$asciiSuffix
      set asciiWorkDir = $tmpOutputDir$asciiArchDir

      mkdir -p $asciiWorkDir 'clean' || exit 1

      # include batch job stdout in ascii tarfile
      cp $stdoutDir/$FRE_JOBID $asciiWorkDir

      if ( ! $?MPI_FAIL ) then
         ls -1 | egrep "$patternGrepAscii" | xargs -I'{}' mv --force '{}' $asciiWorkDir/$begindate.'{}'

	 set asciiOutputDirRemote = ""
         set actionSaveOn         = 1
         set actionXferOn         = 0
         set paramArchiveOn       = 1
         set paramPtmpOn          = 0
         set paramCheckSumOn      = 0
         set paramCompressOn      = 0
      else
         ls -1 | egrep "$patternGrepAscii" | xargs -I'{}' ln --force '{}' $asciiWorkDir/$begindate.'{}'

         set asciiOutputDirRemote = ""
         set actionSaveOn         = 1
         set actionXferOn         = 0
         set paramArchiveOn       = 1
         set paramPtmpOn          = 0
         set paramCheckSumOn      = 0
         set paramCompressOn      = 0
      endif

      set asciiJobName = $FRE_JOBID.output.stager.$begindate.A
      set asciiArgFile = $stateDir/$asciiJobName.args

      echo "set expName                   =   $name"                       > $asciiArgFile
      echo "set type                      =   ascii"                      >> $asciiArgFile
      echo "set stagingType               =   $outputStagingType"         >> $asciiArgFile
      echo "set actionCombineOn           =   0"                          >> $asciiArgFile
      echo "set actionCheckOn             =   0"                          >> $asciiArgFile
      echo "set actionSaveOn              =   $actionSaveOn"              >> $asciiArgFile
      echo "set actionXferOn              =   $actionXferOn"              >> $asciiArgFile
      echo "set actionPPStartOn           =   0"                          >> $asciiArgFile
      echo "set paramArchiveOn            =   $paramArchiveOn"            >> $asciiArgFile
      echo "set paramPtmpOn               =   $paramPtmpOn"               >> $asciiArgFile
      echo "set paramCheckSumOn           =   $paramCheckSumOn"           >> $asciiArgFile
      echo "set paramCompressOn           =   $paramCompressOn"           >> $asciiArgFile
      echo "set workDir                   =   $tmpOutputDir"              >> $asciiArgFile
      echo "set ptmpDir                   =   $ptmpDir"                   >> $asciiArgFile
      echo "set archDir                   =   $asciiArchDir"              >> $asciiArgFile
      echo "set outputDirRemote           =   $asciiOutputDirRemote"      >> $asciiArgFile
      echo "set saveRetry                 =   0"                          >> $asciiArgFile
      echo "set xferRetry                 =   0"                          >> $asciiArgFile


      outputStager $asciiArgFile

      mkdir -p $outputDir/ascii
      tar -cf $outputDir/ascii/$asciiBackup -C $asciiWorkDir .

      rm $asciiArgFile

      unset asciiResult
      unset asciiBackup
	 
      unset asciiXferOptions
      unset asciiSaveOptions

      unset asciiArgFile
      unset asciiJobName

      unset paramCompressOn
      unset paramCheckSumOn
      unset paramPtmpOn
      unset paramArchiveOn
      unset actionXferOn
      unset actionSaveOn

      unset asciiOutputDirRemote
      unset asciiWorkDir
      unset asciiArchDir
      unset asciiSuffix
   endif

   unset asciiFiles

# ---------------- save restart files, namelist, tables etc. and move them from RESTART to INPUT

   pushd $workDir/RESTART

   set restartFiles = ( `ls -1 | egrep "$patternGrepRestart"` )

   if ( $#restartFiles > 0 ) then
      set restartSuffix = restart/$enddate
      set restartBackup = $enddate.tar
      set restartArchDir = $outputDir/$restartSuffix
      set restartWorkDir = $tmpOutputDir$restartArchDir

      mkdir -p $restartWorkDir 'clean' || exit 1

      ls -1 | egrep "$patternGrepRestart" | xargs ln --force --target-directory=$restartWorkDir

      cp --force --preserve=mode,ownership,timestamps $workDir/input.nml $restartWorkDir
      cp --force --preserve=mode,ownership,timestamps $workDir/*_table   $restartWorkDir
      cp --force --preserve=mode,ownership,timestamps $workDir/*_table.yaml $restartWorkDir
      cp --force --preserve=mode,ownership,timestamps $scriptName        $restartWorkDir

      if ( ! $?MPI_FAIL ) then
         set restartOutputDirRemote = ""

         if ( $currentSeg < $segmentsPerSimulation && $irun < $segmentsPerJob ) then
            find $workDir/INPUT   -maxdepth 1 -type f | egrep "$patternGrepRestartNextDrop" | xargs --no-run-if-empty rm --force
            find $workDir/RESTART -maxdepth 1 -type f | egrep "$patternGrepRestartNextMove" | xargs --no-run-if-empty mv --force --target-directory ../INPUT

            set actionCombineOn = $?flagRunTypeRegression
            set actionCheckOn   = $?flagOutputCheckOn
            set actionSaveOn    = 1
            set actionXferOn    = 0
            set paramArchiveOn  = $?flagOutputArchiveOn
            set paramPtmpOn     = 1
            set paramCheckSumOn = $?flagCheckSumOn
            set paramCompressOn = $?flagOutputCompressRestartOn
         else
            if ( $currentSeg < $segmentsPerSimulation && ( $?flagOutputStagingTypeStaged || $?flagOutputStagingTypeChained ) ) then
               if ( $status == 0 ) then
                  if ( $echoOn ) unset echo
                  echo "<NOTE> : The restart directory '$restartArchDir' has been saved successfully"
                  if ( $echoOn ) set echo
               else
                  if ( $echoOn ) unset echo
                  set msg =       "*ERROR*: Can't save the restart directory '$restartArchDir'\n"
                  set msg = "${msg}*ERROR*: restart files have not been saved.  Files will remain in the work directory\n\n"
                  set msg = "${msg}*ERROR*: $workDir\n\n"
                  set msg = "${msg}*ERROR*: To continue the model, you will need to recover the restart files manually.\n"
                  if ( $?batch ) then
                     printf "$msg" | mailx -s "Can't save the restart directory '$restartArchDir'\n" $mailList
                  endif
                  printf "$msg"
                  set restartSaveFailure = 1
                  if ( $echoOn ) set echo
               endif
            endif

            set actionCombineOn = $?flagRunTypeRegression
            set actionCheckOn   = $?flagOutputCheckOn
            set actionSaveOn    = 1
            set actionXferOn    = $?flagOutputXferOn
            set paramArchiveOn  = $?flagOutputArchiveOn
            set paramPtmpOn     = 1
            set paramCheckSumOn = $?flagCheckSumOn
            set paramCompressOn = $?flagOutputCompressRestartOn
         endif
      else
         set restartOutputDirRemote = ""
         set actionCombineOn        = 0
         set actionCheckOn          = 0
         set actionSaveOn           = 1
         set actionXferOn           = 0
         set paramArchiveOn         = 1
         set paramPtmpOn            = 0
         set paramCheckSumOn        = 0
         set paramCompressOn        = 0
      endif

      set restartJobName = $FRE_JOBID.output.stager.$enddate.R
      set restartArgFile = $stateDir/$restartJobName.args

      echo "set expName                   =   $name"                       > $restartArgFile
      echo "set type                      =   restart"                    >> $restartArgFile
      echo "set stagingType               =   $outputStagingType"         >> $restartArgFile
      echo "set actionCombineOn           =   $actionCombineOn"           >> $restartArgFile
      echo "set actionCheckOn             =   $actionCheckOn"             >> $restartArgFile
      echo "set actionSaveOn              =   $actionSaveOn"              >> $restartArgFile
      echo "set actionXferOn              =   $actionXferOn"              >> $restartArgFile
      echo "set actionPPStartOn           =   0"                          >> $restartArgFile
      echo "set paramArchiveOn            =   $paramArchiveOn"            >> $restartArgFile
      echo "set paramPtmpOn               =   $paramPtmpOn"               >> $restartArgFile
      echo "set paramCheckSumOn           =   $paramCheckSumOn"           >> $restartArgFile
      echo "set paramCompressOn           =   $paramCompressOn"           >> $restartArgFile
      echo "set paramVerbosityOn          =   $?flagVerbosityOn"          >> $restartArgFile
      echo "set workDir                   =   $tmpOutputDir"              >> $restartArgFile
      echo "set ptmpDir                   =   $ptmpDir"                   >> $restartArgFile
      echo "set archDir                   =   $restartArchDir"            >> $restartArgFile
      echo "set outputDirRemote           =   $restartOutputDirRemote"    >> $restartArgFile
      echo "set saveRetry                 =   0"                          >> $restartArgFile
      echo "set xferRetry                 =   0"                          >> $restartArgFile
      echo "set mppnccombineOptString     =  '$mppnccombineOptsRestart'"  >> $restartArgFile
      echo "set ardiffTmpdir              =   $ardiffTmpdir"              >> $restartArgFile

      if ( $?flagOutputStagingTypeOnline ) then
         if ( $?MPICH_RANK_REORDER_METHOD ) then
            set mpiRankReorderMethod = $MPICH_RANK_REORDER_METHOD
            unsetenv MPICH_RANK_REORDER_METHOD
         endif

         outputStager $restartArgFile

         if ( $status == 0 ) then
            if ( $echoOn ) unset echo
            echo "<NOTE> : The restart directory '$restartArchDir' has been processed successfully"
            if ( $echoOn ) set echo
         else
            @ outputStagerErrors += 1
            if ( $echoOn ) unset echo
            set msg =       "*WARNING*: Can't save the restart directory '$restartArchDir'\n"
            set msg = "${msg}*WARNING*: restart files have not been saved, you may need to transfer them manually.\n\n"
            set msg = "${msg}*WARNING*: The restart ArgFile has been saved at $restartArgFile.  You may be able\n"
            set msg = "${msg}*WARNING*: use the following command:\n\n"
            set msg = "${msg}*WARNING*: $outputStager $restartArgFile\n"
            if ( $?batch ) then
               printf "$msg" | mailx -s "Can't save the restart directory '$restartArchDir'\n" $mailList
            endif
            printf "$msg"
            if ( $echoOn ) set echo
         endif

         if ( $?mpiRankReorderMethod ) then
            setenv MPICH_RANK_REORDER_METHOD $mpiRankReorderMethod
            unset $mpiRankReorderMethod
         endif
      endif
	 
      mkdir -p $outputDir/restart
      tar -cf $outputDir/restart/$restartBackup -C $restartWorkDir .

      rm $restartArgFile
   
      unset restartResult
      unset restartBackup

      unset restartArgFile
      unset restartJobName

      unset paramCompressOn
      unset paramCheckSumOn
      unset paramPtmpOn
      unset paramArchiveOn
      unset actionXferOn
      unset actionSaveOn
      unset actionCheckOn
      unset actionCombineOn

      unset restartOutputDirRemote
      unset restartWorkDir
      unset restartSuffix
   endif
   popd

   # ---------------- rename region history files

   set regionHistoryFiles = ( `ls -1 | egrep "$patternGrepRegion"` )

   if ( $#regionHistoryFiles > 0 ) then
      if ( ! $?MPI_FAIL ) then
         foreach file ( $regionHistoryFiles )
            mv -f $file `echo $file | sed -r "s/$patternGrepRegion//"`
         end
      endif
   endif

   unset regionHistoryFiles

   # ---------------- combine, save and post-process history files

   set historyFiles = ( `ls -1 | egrep "$patternGrepHistory"` )

   if ( $#historyFiles > 0 ) then
      if ( $?flagOutputCombineHistoryOn ) then
         set historySuffix = history/$begindate.nc
      else
         set historySuffix = history/$begindate.raw.nc
      endif

      set historyBackup = $begindate.nc.tar
      set historyArchDir = $outputDir/$historySuffix
      set historyWorkDir = $tmpOutputDir$historyArchDir

      mkdir -p $historyWorkDir 'clean' || exit 1

      if ( ! $?MPI_FAIL ) then
         ls -1 | egrep "^[0-9][0-9][0-9][0-9][0-9][0-9].+$patternGrepHistory" | xargs -I'{}' mv --force '{}' $historyWorkDir/'{}'
         ls -1 | egrep "$patternGrepHistory" | xargs -I'{}' mv --force '{}' $historyWorkDir/$begindate.'{}'

         set historyOutputDirRemote = ""

         set actionCombineOn = $?flagOutputCombineHistoryOn
         set actionCheckOn   = $?flagOutputCheckOn
         set actionSaveOn    = 1
         set actionXferOn    = $?flagOutputXferOn
         set actionFillGridOn = $?flagOutputFillGridOn

         set actionPPStartOn = 0
         set ppStarterOptions = ( )

         set actionRetryOn   =   1
         set paramArchiveOn  =   $?flagOutputArchiveOn
         @   paramPtmpOn     = ! $?flagOutputArchiveOn
         set paramCheckSumOn =   $?flagCheckSumOn
         set paramCompressOn =   $?flagOutputCompressHistoryOn
      else
         ls -1 | egrep "^[0-9][0-9][0-9][0-9][0-9][0-9].+$patternGrepHistory" | xargs -I'{}' mv --force '{}' $historyWorkDir/$begindate.'{}'
         ls -1 | egrep "$patternGrepHistory" | xargs -I'{}' ln --force '{}' $historyWorkDir/$begindate.'{}'

         set historyOutputDirRemote = ""
         set actionCombineOn        = 0
         set actionCheckOn          = 0
         set actionSaveOn           = 1
         set actionXferOn           = 0
         set actionPPStartOn        = 0
         set actionRetryOn          = 0
         set ppStarterOptions       = ( )
         set paramArchiveOn         = 1
         set paramPtmpOn            = 0
         set paramCheckSumOn        = 0
         set paramCompressOn        = 0
         set actionFillGridOn       = 0
      endif

      set historyJobName = $FRE_JOBID.output.stager.$begindate.H
      set historyArgFile = $stateDir/$historyJobName.args

      echo "set expName                   =   $name"                       > $historyArgFile
      echo "set type                      =   history"                    >> $historyArgFile
      echo "set stagingType               =   $outputStagingType"         >> $historyArgFile
      echo "set actionCombineOn           =   $actionCombineOn"           >> $historyArgFile
      echo "set actionCheckOn             =   $actionCheckOn"             >> $historyArgFile
      echo "set actionSaveOn              =   $actionSaveOn"              >> $historyArgFile
      echo "set actionXferOn              =   $actionXferOn"              >> $historyArgFile
      echo "set actionPPStartOn           =   $actionPPStartOn"           >> $historyArgFile
      echo "set actionRetryOn             =   $actionRetryOn"             >> $historyArgFile
      echo "set actionFillGridOn          =   $actionFillGridOn"          >> $historyArgFile
      echo "set paramArchiveOn            =   $paramArchiveOn"            >> $historyArgFile
      echo "set paramPtmpOn               =   $paramPtmpOn"               >> $historyArgFile
      echo "set paramCheckSumOn           =   $paramCheckSumOn"           >> $historyArgFile
      echo "set paramCompressOn           =   $paramCompressOn"           >> $historyArgFile
      echo "set paramVerbosityOn          =   $?flagVerbosityOn"          >> $historyArgFile
      echo "set workDir                   =   $tmpOutputDir"              >> $historyArgFile
      echo "set ptmpDir                   =   $ptmpDir"                   >> $historyArgFile
      echo "set archDir                   =   $historyArchDir"            >> $historyArgFile
      echo "set outputDirRemote           =   $historyOutputDirRemote"    >> $historyArgFile
      echo "set saveRetry                 =   0"                          >> $historyArgFile
      echo "set xferRetry                 =   0"                          >> $historyArgFile
      echo "set includeDir                =   $includeDir"                >> $historyArgFile
      echo "set mppnccombineOptString     =  '$mppnccombineOptsHistory'"  >> $historyArgFile
      echo "set gridSpec                  =   $gridSpec"                  >> $historyArgFile
      echo "set ardiffTmpdir              =   $ardiffTmpdir"              >> $historyArgFile

      if ( $?flagOutputStagingTypeOnline ) then
         if ( $?MPICH_RANK_REORDER_METHOD ) then
            set mpiRankReorderMethod = $MPICH_RANK_REORDER_METHOD
            unsetenv MPICH_RANK_REORDER_METHOD
         endif

         outputStager $historyArgFile

         if ( $status == 0 ) then
            if ( $echoOn ) unset echo
            echo "<NOTE> : The history directory '$historyArchDir' has been processed successfully"
            if ( $echoOn ) set echo
         else
            @ outputStagerErrors += 1
            if ( $echoOn ) unset echo
            set msg =       "*WARNING*: Can't process the history directory '$historyArchDir'\n"
            set msg = "${msg}*WARNING*: history files have not been saved, you may need to transfer them manually.\n\n"
            set msg = "${msg}*WARNING*: The history ArgFile has been saved at $historyArgFile.  You may be able\n"
            set msg = "${msg}*WARNING*: use the following command:\n\n"
            set msg = "${msg}*WARNING*: $outputStager $historyArgFile\n"
            if ( $?batch ) then
               printf "$msg" | mailx -s "Can't process the history directory '$historyArchDir'" $mailList
            endif
            printf "$msg"
            if ( $echoOn ) set echo
         endif

         if ( $?mpiRankReorderMethod ) then
            setenv MPICH_RANK_REORDER_METHOD $mpiRankReorderMethod
            unset $mpiRankReorderMethod
         endif
      endif

      mkdir -p $outputDir/history
      tar -cf $outputDir/history/$historyBackup -C $historyWorkDir .

      rm $historyArgFile
      
      unset historyResult
      unset historyXferOptions
      unset historySaveOptions

      unset historyArgFile
      unset historyJobName

      unset paramCompressOn
      unset paramCheckSumOn
      unset paramPtmpOn
      unset paramArchiveOn
      unset ppStarterOptions
      unset actionPPStartOn
      unset actionXferOn
      unset actionSaveOn
      unset actionCheckOn
      unset actionCombineOn

      unset historyOutputDirRemote
      unset historyWorkDir
      unset historyArchDir
      unset historySuffix
   endif

   unset historyFiles

   # ---------------- terminate script if MPI failed

   if ( $?MPI_FAIL ) then
      if ( $echoOn ) unset echo
      echo "*ERROR*: The MPI failed (segment $currentSeg)"
      echo "*ERROR*: Any output that may have been generated is in the '$outputDir'"
      echo "*ERROR*: The '$workDir' is being kept for possible debugging"
      if ( $echoOn ) set echo

      exit 1
   endif

   # ---------------- unset remaining restart variables

   unset restartArchDir
   unset restartFiles

   # ---------------- increment loop counters

   @ currentSeg++
   @ irun++
end

################################################################################
#--------------------------- after the main loop -------------------------------
################################################################################

if ( $echoOn ) unset echo
echo "################################################################################"
echo "# ending"
echo "################################################################################"
if ( $echoOn ) set echo

if ( $?flagRunTypeProduction ) then
   if ( $currentSeg <= $segmentsPerSimulation ) then
      if ( -f $queue_file ) then
         if ( -r $queue_file ) then
            source $queue_file
         else
            if ( $echoOn ) unset echo
            echo "*ERROR*: The queue file '$queue_file' is not readable"
            if ( $echoOn ) set echo
            exit 1
         endif
      endif

      if ( ! $continueFlag ) then
         if ( $echoOn ) unset echo
         echo "<NOTE> : Stopping execution"
         if ( $echoOn ) set echo
         exit 0
      endif

      if ( -f $scriptName ) then
         if ( -r $scriptName ) then
            set nextOptions = ( $submitOptionsProject )
            set result = `submit -O "$nextOptions" $scriptName`
            if ( $status == 0 ) then
               if ( $echoOn ) unset echo
               echo "<NOTE> : The job '$result' to run the '$scriptName' has been submitted successfully"
               if ( $echoOn ) set echo
            else
               if ( $echoOn ) unset echo
               echo "*ERROR*: Can't submit the '$scriptName'"
               if ( $echoOn ) set echo
               exit 1
            endif
            unset result
            unset nextOptions
         else
            if ( $echoOn ) unset echo
            @ lastSeg = $currentSeg - 1
            echo "WARNING: The script '$scriptName' exists, but is not readable (run $lastSeg) - it can't be resubmitted"
            unset lastSeg
            if ( $echoOn ) set echo
         endif
      else
         if ( $echoOn ) unset echo
         @ lastSeg = $currentSeg - 1
         echo "WARNING: The script '$scriptName' does not exist (run $lastSeg) - it can't be resubmitted"
         unset lastSeg
         if ( $echoOn ) set echo
      endif
   endif
endif

set -r runtimeEnd = `date "+%s"`
set -r runtime = `echo "$runtimeEnd - $runtimeBeg" | bc -l`

rm -rf $workDir

if ( $echoOn ) unset echo
echo "<NOTE> : Finishing on `date`"
echo "<NOTE> : Runtime = '$runtime' (seconds)"
echo "<NOTE> : Natural end-of-script for '$scriptName'"
if ( $echoOn ) set echo

exit 0
