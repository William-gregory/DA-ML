#!/bin/tcsh -fx                                                                                                                                                                                                                          
#FRE scheduler-options                                                                                                                                                                                                                   
#SBATCH --chdir=/gpfs/f5/gfdl_o/scratch/William.Gregory/SPEAR_reforecast/SPEAR_Q_nonsymMOM6_OTA_FTorch/ncrc5.intel-classic-repro-openmp/exec                                                                                             
#SBATCH --output=/gpfs/f5/gfdl_o/scratch/William.Gregory/SPEAR_reforecast/SPEAR_Q_nonsymMOM6_OTA_FTorch/ncrc5.intel-classic-repro-openmp/exec/%x.o%j                                                                                     
#SBATCH --job-name=compile_SPEAR_Q_nonsymMOM6_OTA_FTorch.csh                                                                                                                                                                             
#SBATCH --comment=fre/bronx-21                                                                                                                                                                                                           
#SBATCH --time=120                                                                                                                                                                                                                       
#SBATCH --nodes=1 --ntasks=8                                                                                                                                                                                                             
#SBATCH --partition=eslogin_c5                                                                                                                                                                                                           
#SBATCH --mail-user=William.Gregory@noaa.gov                                                                                                                                                                                             
#SBATCH --mail-type=fail                                                                                                                                                                                                                 
#SBATCH --export=NONE                                                                                                                                                                                                                    
#SBATCH --clusters=es                                                                                                                                                                                                                    
#SBATCH --account=gfdl_o                                                                                                                                                                                                                 

# Compile Script for Experiment 'SPEAR_Q_nonsymMOM6_OTA_FTorch'                                                                                                                                                                          
# ------------------------------------------------------------------------------                                                                                                                                                         
# The script created at 2024-11-18T14:03:31 via:                                                                                                                                                                                         
# /ncrc/home2/fms/local/opt/fre-commands/bronx-21/bin/fremake --link --ncores=8 --platform=ncrc5.intel-classic --target=repro,openmp --walltime=120 --xmlfile=/autofs/ncrc-svm1_home2/William.Gregory/F5_coupled/SPEAR_Q_nonsymMOM6_OTA.\
xml SPEAR_Q_nonsymMOM6_OTA_FTorch                                                                                                                                                                                                        
# ------------------------------------------------------------------------------                                                                                                                                                         

set -r echoOn = $?echo

if ( $echoOn ) unset echo
echo "<NOTE> : Starting at $HOST on `date`"
if ( $echoOn ) set echo

unalias *

# ---------------- Set build, src and stage directories                                                                                                                                                                                  

set src_dir = /gpfs/f5/gfdl_o/scratch/William.Gregory/SPEAR_reforecast/SPEAR_Q_nonsymMOM6_OTA_FTorch/src
set bld_dir = /gpfs/f5/gfdl_o/scratch/William.Gregory/SPEAR_reforecast/SPEAR_Q_nonsymMOM6_OTA_FTorch/ncrc5.intel-classic-repro-openmp/exec
set ptmp_dir = /gpfs/f5/gfdl_o/scratch/William.Gregory/ptmp/SPEAR_reforecast/SPEAR_Q_nonsymMOM6_OTA_FTorch

# ---------------- Make template                                                                                                                                                                                                         

set mkmf_template = /ncrc/home2/fms/local/opt/fre-commands/bronx-21/site/ncrc5/intel-classic.mk

# ---------------- set environment                                                                                                                                                                                                       

source $MODULESHOME/init/tcsh
module use /ncrc/home2/fms/local/modulefiles
if ( $echoOn ) unset echo
source $bld_dir/env.cshrc
if ( $echoOn ) set echo

# ---------------- write main Makefile                                                                                                                                                                                                   

sed -e 's/<TAB>/\t/' >$bld_dir/Makefile <<END                                                                                                                                                                                            
# Makefile for Experiment 'SPEAR_Q_nonsymMOM6_OTA_FTorch'                                                                                                                                                                                
                                                                                                                                                                                                                                         
SRCROOT = $src_dir/                                                                                                                                                                                                                      
BUILDROOT = $bld_dir/                                                                                                                                                                                                                    
STAGEDIR = $ptmp_dir/$bld_dir/                                                                                                                                                                                                           
                                                                                                                                                                                                                                         
MK_TEMPLATE = $mkmf_template                                                                                                                                                                                                             
include \$(MK_TEMPLATE)                                                                                                                                                                                                                  
                                                                                                                                                                                                                                         
fms_SPEAR_Q_nonsymMOM6_OTA_FTorch.x: coupler/libcoupler.a atmos_dyn/libatmos_dyn.a sis2/libsis2.a icebergs/libicebergs.a atmos_phys/libatmos_phys.a ice_param/libice_param.a mom6/libmom6.a land_lad2/libland_lad2.a FMS/libFMS.a        
<TAB>\$(LD) \$^ \$(LDFLAGS) -o \$@ \$(STATIC_LIBS)                                                                                                                                                                                       
                                                                                                                                                                                                                                         
FMS/libFMS.a:  FORCE                                                                                                                                                                                                                     
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=FMS \$(@F)                                                                                                                            
                                                                                                                                                                                                                                         
atmos_phys/libatmos_phys.a: FMS/libFMS.a FORCE                                                                                                                                                                                           
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=atmos_phys \$(@F)                                                                                                                     
                                                                                                                                                                                                                                         
atmos_dyn/libatmos_dyn.a: atmos_phys/libatmos_phys.a FMS/libFMS.a FORCE                                                                                                                                                                  
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=atmos_dyn \$(@F)                                                                                                                      
                                                                                                                                                                                                                                         
land_lad2/libland_lad2.a: FMS/libFMS.a FORCE                                                                                                                                                                                             
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=land_lad2 \$(@F)                                                                                                                      
                                                                                                                                                                                                                                         
mom6/libmom6.a: FMS/libFMS.a FORCE                                                                                                                                                                                                       
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE) OPENMP="" --directory=mom6 \$(@F)

ice_param/libice_param.a: FMS/libFMS.a FORCE                                                                                                                                                                                             
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=ice_param \$(@F)                                                                                                                      
                                                                                                                                                                                                                                         
sis2/libsis2.a: mom6/libmom6.a icebergs/libicebergs.a ice_param/libice_param.a FMS/libFMS.a FORCE                                                                                                                                        
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=sis2 \$(@F)                                                                                                                           
                                                                                                                                                                                                                                         
icebergs/libicebergs.a: FMS/libFMS.a FORCE                                                                                                                                                                                               
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=icebergs \$(@F)                                                                                                                       
                                                                                                                                                                                                                                         
coupler/libcoupler.a: sis2/libsis2.a atmos_dyn/libatmos_dyn.a mom6/libmom6.a atmos_phys/libatmos_phys.a icebergs/libicebergs.a land_lad2/libland_lad2.a FMS/libFMS.a FORCE                                                               
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=coupler \$(@F)                                                                                                                        
                                                                                                                                                                                                                                         
FORCE:                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                         
stage:                                                                                                                                                                                                                                   
<TAB>install -d \$(STAGEDIR)                                                                                                                                                                                                             
<TAB>install -m 555 fms_SPEAR_Q_nonsymMOM6_OTA_FTorch.x \$(STAGEDIR)                                                                                                                                                                     
                                                                                                                                                                                                                                         
clean:                                                                                                                                                                                                                                   
<TAB>\$(MAKE) --directory=FMS clean                                                                                                                                                                                                      
<TAB>\$(MAKE) --directory=atmos_phys clean                                                                                                                                                                                               
<TAB>\$(MAKE) --directory=atmos_dyn clean                                                                                                                                                                                                
<TAB>\$(MAKE) --directory=land_lad2 clean                                                                                                                                                                                                
<TAB>\$(MAKE) --directory=mom6 clean                                                                                                                                                                                                     
<TAB>\$(MAKE) --directory=ice_param clean                                                                                                                                                                                                
<TAB>\$(MAKE) --directory=sis2 clean                                                                                                                                                                                                     
<TAB>\$(MAKE) --directory=icebergs clean                                                                                                                                                                                                 
<TAB>\$(MAKE) --directory=coupler clean

localize:                                                                                                                                                                                                                                
<TAB>\$(MAKE) -f \$(BUILDROOT)FMS/Makefile localize                                                                                                                                                                                      
<TAB>\$(MAKE) -f \$(BUILDROOT)atmos_phys/Makefile localize                                                                                                                                                                               
<TAB>\$(MAKE) -f \$(BUILDROOT)atmos_dyn/Makefile localize                                                                                                                                                                                
<TAB>\$(MAKE) -f \$(BUILDROOT)land_lad2/Makefile localize                                                                                                                                                                                
<TAB>\$(MAKE) -f \$(BUILDROOT)mom6/Makefile localize                                                                                                                                                                                     
<TAB>\$(MAKE) -f \$(BUILDROOT)ice_param/Makefile localize                                                                                                                                                                                
<TAB>\$(MAKE) -f \$(BUILDROOT)sis2/Makefile localize                                                                                                                                                                                     
<TAB>\$(MAKE) -f \$(BUILDROOT)icebergs/Makefile localize                                                                                                                                                                                 
<TAB>\$(MAKE) -f \$(BUILDROOT)coupler/Makefile localize                                                                                                                                                                                  
                                                                                                                                                                                                                                         
distclean:                                                                                                                                                                                                                               
<TAB>\$(RM) -r FMS                                                                                                                                                                                                                       
<TAB>\$(RM) -r atmos_phys                                                                                                                                                                                                                
<TAB>\$(RM) -r atmos_dyn                                                                                                                                                                                                                 
<TAB>\$(RM) -r land_lad2                                                                                                                                                                                                                 
<TAB>\$(RM) -r mom6                                                                                                                                                                                                                      
<TAB>\$(RM) -r ice_param                                                                                                                                                                                                                 
<TAB>\$(RM) -r sis2                                                                                                                                                                                                                      
<TAB>\$(RM) -r icebergs                                                                                                                                                                                                                  
<TAB>\$(RM) -r coupler                                                                                                                                                                                                                   
<TAB>\$(RM) fms_SPEAR_Q_nonsymMOM6_OTA_FTorch.x                                                                                                                                                                                          
<TAB>\$(RM) Makefile                                                                                                                                                                                                                     
                                                                                                                                                                                                                                         
END

# ---------------- create component Makefiles                                                                                                                                                                                            

mkdir -p $bld_dir/FMS
list_paths -l -o $bld_dir/FMS/pathnames_FMS $src_dir/FMS
cd $bld_dir
pushd FMS
mkmf -m Makefile -a $src_dir -b $bld_dir -p libFMS.a -t $mkmf_template -c "-DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF -DMAXFIELDMETHODS_=500 -D__PGI" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framewo\
rk $bld_dir/FMS/pathnames_FMS
popd

mkdir -p $bld_dir/atmos_phys
list_paths -l -o $bld_dir/atmos_phys/pathnames_atmos_phys $src_dir/atmos_param $src_dir/atmos_shared
cd $bld_dir
pushd atmos_phys
mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_phys.a -t $mkmf_template -c "-DINTERNAL_FILE_NML -DAM3_CHEM" -o "-I$bld_dir/FMS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir/atmos_p\
hys/pathnames_atmos_phys
popd

mkdir -p $bld_dir/atmos_dyn
list_paths -l -o $bld_dir/atmos_dyn/pathnames_atmos_dyn $src_dir/atmos_drivers/coupled $src_dir/atmos_cubed_sphere/driver/coupled $src_dir/atmos_cubed_sphere/model $src_dir/atmos_cubed_sphere/model_nh_null $src_dir/atmos_cubed_spher\
e/tools $src_dir/atmos_cubed_sphere/GFDL_tools
cd $bld_dir
pushd atmos_dyn
mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_dyn.a -t $mkmf_template -c "-DINTERNAL_FILE_NML -DSPMD -D__PGI" -o "-I$bld_dir/atmos_phys -I$bld_dir/FMS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/f\
ramework $bld_dir/atmos_dyn/pathnames_atmos_dyn
popd

mkdir -p $bld_dir/land_lad2
list_paths -l -o $bld_dir/land_lad2/pathnames_land_lad2 $src_dir/land_lad2
cd $bld_dir
pushd land_lad2
mkmf -m Makefile -a $src_dir -b $bld_dir -p libland_lad2.a -t $mkmf_template --use-cpp -c "-DINTERNAL_FILE_NML -nostdinc" -o "-I$bld_dir/FMS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir\
/land_lad2/pathnames_land_lad2
popd

mkdir -p $bld_dir/mom6
list_paths -l -o $bld_dir/mom6/pathnames_mom6 $src_dir/MOM6/config_src/dynamic $src_dir/MOM6/config_src/coupled_driver $src_dir/MOM6/src/*/ $src_dir/MOM6/src/*/*/ $src_dir/MOM6/ECDA/
cd $bld_dir
pushd mom6
mkmf -m Makefile -a $src_dir -b $bld_dir -p libmom6.a -t $mkmf_template -c "-DINTERNAL_FILE_NML -DMAX_FIELDS_=100" -o "-I$bld_dir/FMS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir/mom6/p\
athnames_mom6
popd

mkdir -p $bld_dir/ice_param
list_paths -l -o $bld_dir/ice_param/pathnames_ice_param $src_dir/ice_param
cd $bld_dir
pushd ice_param
mkmf -m Makefile -a $src_dir -b $bld_dir -p libice_param.a -t $mkmf_template -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/FMS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir/ice_param/pathnames\
_ice_param
popd

mkdir -p $bld_dir/sis2
list_paths -l -o $bld_dir/sis2/pathnames_sis2 $src_dir/SIS2/config_src/dynamic $src_dir/SIS2/src
cd $bld_dir
pushd sis2
mkmf -m Makefile -a $src_dir -b $bld_dir -p libsis2.a -t $mkmf_template -c "-DINTERNAL_FILE_NML -DNONSYMMETRIC_MEMORY_" -o "-I$bld_dir/ice_param -I$bld_dir/icebergs -I$bld_dir/mom6 -I$bld_dir/FMS -I/gpfs/f5/gfdl_o/scratch/William.Gr\
egory/FTorch/src/build/torch_v2p1_Inteloneapi2023p2p0/include/ftorch" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir/sis2/pathnames_sis2
popd

mkdir -p $bld_dir/icebergs
list_paths -l -o $bld_dir/icebergs/pathnames_icebergs $src_dir/icebergs
cd $bld_dir
pushd icebergs
mkmf -m Makefile -a $src_dir -b $bld_dir -p libicebergs.a -t $mkmf_template -o "-I$bld_dir/FMS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir/icebergs/pathnames_icebergs
popd

mkdir -p $bld_dir/coupler
list_paths -l -o $bld_dir/coupler/pathnames_coupler $src_dir/FMScoupler/full $src_dir/FMScoupler/shared
cd $bld_dir
pushd coupler
mkmf -m Makefile -a $src_dir -b $bld_dir -p libcoupler.a -t $mkmf_template -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/atmos_dyn -I$bld_dir/sis2 -I$bld_dir/atmos_phys -I$bld_dir/icebergs -I$bld_dir/mom6 -I$bld_dir/land_lad2 -I$bld_dir/F\
MS" -IFMS/include -IFMS/mpp/include -IMOM6/pkg/CVMix-src/include -IMOM6/src/framework $bld_dir/coupler/pathnames_coupler
popd

# ---------------- call make on the main Makefile

make  REPRO=on OPENMP=on NETCDF=3 LDFLAGS="-L/gpfs/f5/gfdl_o/scratch/William.Gregory/FTorch/src/build/torch_v2p1_Inteloneapi2023p2p0/lib64 -lftorch" fms_SPEAR_Q_nonsymMOM6_OTA_FTorch.x

if ( $status == 0 ) then
  if ( $?NiNaC_LVL ) then
    if ( $NiNaC_LVL > 0 ) then
      # Run NiNaC                                                                                                                                                                                                                        
      $NiNaC_BldRx $src_dir $bld_dir
      if ( $status != 0 ) then
        if ( $echoOn ) unset echo
        echo "NiNaC Note: While NiNaC loaded attempt at NiNaC_BldRx failed with exit status $status : FRE continuing as normal."
        if ( $echoOn ) set echo
      endif
    endif
  endif

  if ( $echoOn ) unset echo
  echo "<NOTE> : make succeeded for SPEAR_Q_nonsymMOM6_OTA_FTorch."
  if ( $echoOn ) set echo
else
  if ( $echoOn ) unset echo
  echo "*ERROR*: make failed for SPEAR_Q_nonsymMOM6_OTA_FTorch."
  if ( $echoOn ) set echo
  exit 1
endif

exit 0
