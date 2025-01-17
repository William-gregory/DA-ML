#!/bin/bash

EXP=SPEAR_Q_nonsymMOM6_OTA_seaiceML

module load intel-rt/2021.1.2 intel-tbb/2021.1.1 intel-mkl/2021.1.1 intel-debugger/10.0.0 intel-dpl/2021.1.2 /opt/intel/oneapi/compiler/2021.1.2/linux/lib/oclfpga/modulefiles/oclfpga
module load intel/2021.1.2 openmpi/intel-2021.1/4.1.0 hdf5/intel-2021.1/1.10.6 netcdf/intel-19.1/hdf5-1.10.6/4.7.4

export LD_LIBRARY_PATH=/usr/lib:/usr/local/openmpi/4.1.0/intel20211/lib64:/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64:$LD_LIBRARY_PATH

if [ ! -d /scratch/cimes/wg4031/$EXP ]; then
    cd /scratch/cimes/wg4031/
    mkdir $EXP
    cd $EXP
    mkdir src
    mkdir build
    mkdir ptmp
fi

if [ ! -d /scratch/cimes/wg4031/$EXP/src ]; then
    mkdir /scratch/cimes/wg4031/$EXP/src
fi

if [ ! -d /scratch/cimes/wg4031/$EXP/build ]; then
    mkdir /scratch/cimes/wg4031/$EXP/build
fi

if [ ! -d /scratch/cimes/wg4031/$EXP/ptmp ]; then
    mkdir /scratch/cimes/wg4031/$EXP/ptmp
fi

src_dir=/scratch/cimes/wg4031/$EXP/src
bld_dir=/scratch/cimes/wg4031/$EXP/build
ptmp_dir=/scratch/cimes/wg4031/$EXP/ptmp
cd $src_dir

### CLONE REPOS ###
if [ ! -d FMS ]; then
    git clone -b main --recursive http://github.com/NOAA-GFDL/FMS.git FMS
    cd FMS
    git checkout cb85999e0f9c168f41a169548d8272e738f5eb1a
    cd ..
fi

if [ ! -d atmos_shared ]; then
    git clone -b xanadu http://egitlab.gfdl.noaa.gov/fms/atmos_shared.git atmos_shared
    cd atmos_shared
    git checkout 23054287c2edd2f4148e6ea8fd2db152ec5b470c
    cd ..
fi

if [ ! -d atmos_param ]; then
    git clone -b xanadu http://egitlab.gfdl.noaa.gov/fms/atmos_param.git atmos_param
    cd atmos_param
    git checkout 8814f6a412bd00a6b79b0b184e643e393f9c719a
    cd ..
fi

if [ ! -d atmos_drivers ]; then
    git clone -b xanadu http://egitlab.gfdl.noaa.gov/fms/atmos_drivers.git atmos_drivers
    cd atmos_drivers
    git checkout 3be6ed406de2db29766746a69115fd6a47048692
    cd ..
fi

if [ ! -d atmos_cubed_sphere ]; then
    git clone -b xanadu http://egitlab.gfdl.noaa.gov/fms/atmos_cubed_sphere.git atmos_cubed_sphere
    cd atmos_cubed_sphere
    git checkout 9a92a8fa0551ca27a8382a046ed72d33f13f40e3
    cd ..
fi

if [ ! -d land_lad2 ]; then
    git clone -b xanadu http://egitlab.gfdl.noaa.gov/fms/land_lad2.git land_lad2
    cd land_lad2
    git checkout 7a8a4fcaf6b49b02926cd739f225e20f42d01615
    cd ..
fi

if [ ! -d ice_param ]; then
    git clone -b dev/master http://egitlab.gfdl.noaa.gov/fms/ice_param.git ice_param
fi

if [ ! -d MOM6 ]; then
    git clone -b dev/SPEAR_ODA --recursive http://github.com/feiyulu/MOM6.git MOM6
    cd MOM6
    git clone -b GFDL https://egitlab.gfdl.noaa.gov/Feiyu.Lu/ECDA.git ECDA
    cd ..
fi

if [ ! -d SIS2 ]; then
    git clone -b ML_pure_fortran --recursive http://github.com/William-gregory/SIS2.git SIS2
fi

if [ ! -d icebergs ]; then
    git clone -b dev/gfdl --recursive http://github.com/NOAA-GFDL/icebergs.git icebergs
    cd icebergs
    git checkout b114a809187317909e19fdca7ff843f2a603f011
    cd ..
fi

if [ ! -d FMScoupler ]; then
    git clone -b main --recursive http://github.com/NOAA-GFDL/FMScoupler.git FMScoupler
    cd FMScoupler
    git checkout f48e67ba9c343152045182e80df5c14021130d47
    cd ..
fi

if [ ! -d mkmf ]; then
    git clone https://github.com/NOAA-GFDL/mkmf.git mkmf
    cd mkmf/templates
    sed -i 's/ftn/ifort/g' ncrc5-intel-classic.mk
    sed -i 's/cc/icc/g' ncrc5-intel-classic.mk
    cd ../../
fi
### FINISH CLONE REPOS ###


### COMPILE STEP ###

# ---------------- write main Makefile
sed -e 's/<TAB>/\t/' >$bld_dir/Makefile <<END
# Makefile for Experiment '$EXP'

SRCROOT = $src_dir/
BUILDROOT = $bld_dir/
STAGEDIR = $ptmp_dir/$bld_dir/

MK_TEMPLATE = $src_dir/mkmf/templates/ncrc5-intel-classic.mk
include \$(MK_TEMPLATE)

fms_$EXP.x: coupler/libcoupler.a atmos_dyn/libatmos_dyn.a sis2/libsis2.a icebergs/libicebergs.a atmos_phys/libatmos_phys.a mom6/libmom6.a ice_param/libice_param.a land_lad2/libland_lad2.a FMS/libFMS.a
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

sis2/libsis2.a: icebergs/libicebergs.a mom6/libmom6.a ice_param/libice_param.a FMS/libFMS.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=sis2 \$(@F)

icebergs/libicebergs.a: FMS/libFMS.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=icebergs \$(@F)

coupler/libcoupler.a: sis2/libsis2.a atmos_dyn/libatmos_dyn.a atmos_phys/libatmos_phys.a mom6/libmom6.a land_lad2/libland_lad2.a icebergs/libicebergs.a FMS/libFMS.a FORCE
<TAB>\$(MAKE) SRCROOT=\$(SRCROOT) BUILDROOT=\$(BUILDROOT) MK_TEMPLATE=\$(MK_TEMPLATE)  --directory=coupler \$(@F)

FORCE:

stage:
<TAB>install -d \$(STAGEDIR)
<TAB>install -m 555 fms_$EXP.x \$(STAGEDIR)

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
<TAB>\$(RM) fms_$EXP.x
<TAB>\$(RM) Makefile

END

# ---------------- create component Makefiles

cd $bld_dir

#FMS
mkdir -p $bld_dir/FMS
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/FMS/pathnames_FMS $src_dir/FMS
cd $bld_dir
pushd FMS
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libFMS.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF -DMAXFIELDMETHODS_=500 -D__PGI" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include $bld_dir/FMS/pathnames_FMS
popd

#ATMOS PHYSICS
mkdir -p $bld_dir/atmos_phys
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/atmos_phys/pathnames_atmos_phys $src_dir/atmos_param $src_dir/atmos_shared
cd $bld_dir
pushd atmos_phys
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_phys.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-DINTERNAL_FILE_NML -DAM3_CHEM" -o "-I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/atmos_phys/pathnames_atmos_phys
popd

#ATMOS DYNAMICS
mkdir -p $bld_dir/atmos_dyn
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/atmos_dyn/pathnames_atmos_dyn $src_dir/atmos_drivers/coupled $src_dir/atmos_cubed_sphere/driver/coupled $src_dir/atmos_cubed_sphere/model $src_dir/atmos_cubed_sphere/model_nh_null $src_dir/atmos_cubed_sphere/tools $src_dir/atmos_cubed_sphere/GFDL_tools
cd $bld_dir
pushd atmos_dyn
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libatmos_dyn.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-DINTERNAL_FILE_NML -DSPMD -D__PGI" -o "-I$bld_dir/atmos_phys -I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/atmos_dyn/pathnames_atmos_dyn
popd

#LAND
mkdir -p $bld_dir/land_lad2
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/land_lad2/pathnames_land_lad2 $src_dir/land_lad2
cd $bld_dir
pushd land_lad2
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libland_lad2.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk --use-cpp -c "-DINTERNAL_FILE_NML -nostdinc" -o "-I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/land_lad2/pathnames_land_lad2
popd

#OCEAN
mkdir -p $bld_dir/mom6
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/mom6/pathnames_mom6 $src_dir/MOM6/config_src/dynamic $src_dir/MOM6/config_src/coupled_driver $src_dir/MOM6/src/*/ $src_dir/MOM6/src/*/*/ $src_dir/MOM6/ECDA/
cd $bld_dir
pushd mom6
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libmom6.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-UHAVE_GETTID -DINTERNAL_FILE_NML -DMAX_FIELDS_=100" -o "-I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/mom6/pathnames_mom6
popd

#ICE PARAMETERS
mkdir -p $bld_dir/ice_param
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/ice_param/pathnames_ice_param $src_dir/ice_param
cd $bld_dir
pushd ice_param
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libice_param.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/ice_param/pathnames_ice_param
popd

#SEA ICE
mkdir -p $bld_dir/sis2
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/sis2/pathnames_sis2 $src_dir/SIS2/config_src/dynamic $src_dir/SIS2/src
cd $bld_dir
pushd sis2
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libsis2.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-DINTERNAL_FILE_NML -DNONSYMMETRIC_MEMORY_" -o "-I$bld_dir/mom6 -I$bld_dir/icebergs -I$bld_dir/ice_param -I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/sis2/pathnames_sis2
popd

#ICEBERGS
mkdir -p $bld_dir/icebergs
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/icebergs/pathnames_icebergs $src_dir/icebergs
cd $bld_dir
pushd icebergs
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libicebergs.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -o "-I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/icebergs/pathnames_icebergs
popd

#COUPLER
mkdir -p $bld_dir/coupler
$src_dir/mkmf/bin/list_paths -l -o $bld_dir/coupler/pathnames_coupler $src_dir/FMScoupler/full $src_dir/FMScoupler/shared
cd $bld_dir
pushd coupler
$src_dir/mkmf/bin/mkmf -m Makefile -a $src_dir -b $bld_dir -p libcoupler.a -t $src_dir/mkmf/templates/ncrc5-intel-classic.mk -c "-DINTERNAL_FILE_NML" -o "-I$bld_dir/atmos_dyn -I$bld_dir/sis2 -I$bld_dir/icebergs -I$bld_dir/mom6 -I$bld_dir/land_lad2 -I$bld_dir/atmos_phys -I$bld_dir/FMS" -I$src_dir/FMS/include -I$src_dir/FMS/mpp/include -I$src_dir/MOM6/pkg/CVMix-src/include -I$src_dir/MOM6/src/framework -I/usr/local/openmpi/4.1.0/intel20211/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/include -I/usr/local/hdf5/intel-2021.1/1.10.6/include -I/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/include $bld_dir/coupler/pathnames_coupler
popd

#MAKE EXECUTABLE
make REPRO=on OPENMP=on NETCDF=3 LDFLAGS="-L/usr/lib -lrt -lpthread -L/opt/intel/oneapi/compiler/2021.1.2/linux/compiler/lib/intel64 -liomp5 -L/usr/local/openmpi/4.1.0/intel20211/lib64 -lmpi -lmpi_mpifh -L/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/4.7.4/lib64 -lnetcdf -lnetcdff -L/usr/local/hdf5/intel-2021.1/1.10.6/lib64 -lhdf5 -lhdf5_fortran" fms_$EXP.x




    
