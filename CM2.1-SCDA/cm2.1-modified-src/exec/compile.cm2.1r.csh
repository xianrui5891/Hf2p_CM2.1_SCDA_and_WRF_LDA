#!/bin/tcsh -f

# Compile Script for Experiment 'CM2.1R_exec'
# ------------------------------------------------------------------------------
# The script created at 2010-04-30T13:57:41 via:
# /home/decp/fjz/CM2.1R/bin/fremake --npes=8 --platform=hpcs.hpcs --target=prod --walltime=120 --xmlfile=/home/snz/CM2.1R/CM2.1R.xml
#
# ------------------------------------------------------------------------------

# Change Log
# ------------------------------------------------------------------------------
# 20101027 - Seth Underwood <sunderwood@hpti.com>
#     - Modified for use with fre 4 intel makefile template on gaea
# ------------------------------------------------------------------------------

set -r echoOn = $?echo

if ( $echoOn ) unset echo
echo "<NOTE> : Starting at $HOST on `date`"
if ( $echoOn ) set echo

unalias *

#set root = /Users/lmkli/work/zhangsq
set root = /home/ouc/zhuxianrui/CM2-f2py/CM2.1-python
set exec_dir = $root/exec
set src_dir =  $root/src_v5.0

# ---------------- set environment

if ( $echoOn ) unset echo
source $exec_dir/env.cshrc
if ( $echoOn ) set echo

# ---------------- write main Makefile
sed -e 's/<TAB>/\t/' >${exec_dir}/Makefile <<END
# Makefile for Experiment 'CM2.1U_Control-1990_E1.M_3A'
include $exec_dir/intel.mk

fms_CM2.1R_exec_v5.0_ps.x: libcoupler.a libland.a libecda.a libfms.a
<TAB>\$(LD) \$^ \$(LDFLAGS) -o \$@ 

libfms.a:  FORCE
<TAB>make  VERBOSE=TRUE NETCDF=3 -f Makefile.fms \$@

libecda.a: libfms.a  FORCE
<TAB>make  VERBOSE=TRUE NETCDF=3 -f Makefile.ecda \$@

libland.a: libfms.a FORCE
<TAB>make VERBOSE=TRUE NETCDF=3 -f Makefile.land \$@

libcoupler.a: libfms.a libecda.a libland.a FORCE
<TAB>make VERBOSE=TRUE NETCDF=3 -f Makefile.coupler \$@

FORCE:

clean:
<TAB>make  NETCDF=3 -f Makefile.ecda clean
<TAB>make  NETCDF=3 -f Makefile.fms clean
<TAB>make  NETCDF=3 -f Makefile.land clean
<TAB>make  NETCDF=3 -f Makefile.coupler clean

localize:
<TAB>make  NETCDF=3 -f Makefile.ecda localize
<TAB>make  NETCDF=3 -f Makefile.fms localize
<TAB>make  NETCDF=3 -f Makefile.land localize
<TAB>make  NETCDF=3 -f Makefile.coupler localize

END

# ---------------- Local variables to for ECDA source, and src/exec directories
set cppDefs = ( "-DMAX_LINKS_=10 -DENABLE_ECDA -Duse_netCDF -Duse_libMPI -DSPMD -DLAND_BND_TRACERS -DUSE_OCEAN_BGC -Duse_shared_pointers -DNI_=360 -DNJ_=200 -DNK_=50 -DNI_LOCAL_=60 -DNJ_LOCAL_=40 -DENABLE_ODA -DENABLE_ADA -DENABLE_IDA -from-rtn1 -to-rtn999" )
set srcList = "( ida/{ida_driver.F90,ida_types.F90,ida_core.F90,eakf_mtv_ida_dn1.f90,eakf_mtv_ida_up.f90} ida/model/{model_ida_dn1.f90,model_ida_up.f90} ida/obs/{obs_eakf_ida_dn.F90,obs_eakf_ida_up.F90} )"

# ---------------- create component Makefile's

cd $src_dir
list_paths -o pathnames_ecda ecda/shared oda ice_sis ice_param mom4p1 ocean_shared atmos_param atmos_shared atmos_coupled atmos_fv_dynamics ada
cd $exec_dir
mkmf -m Makefile.ecda -a $src_dir -p libecda.a -t $exec_dir/intel.mk -c "$cppDefs" $srcList pathnames_ecda shared/mpp/include shared/include 

cd $src_dir
list_paths -o pathnames_coupler coupler
cd $exec_dir
mkmf -m Makefile.coupler -a $src_dir -p libcoupler.a -t $exec_dir/intel.mk -c "$cppDefs" pathnames_coupler shared/mpp/include shared/include

cd $src_dir
list_paths -o pathnames_land land_lad land_param
cd $exec_dir
mkmf -m Makefile.land -a $src_dir -p libland.a -t $exec_dir/intel.mk -c "$cppDefs" pathnames_land shared/mpp/include shared/include

cd $src_dir
list_paths -o pathnames_fms shared
cd $exec_dir
mkmf -m Makefile.fms -a $src_dir -p libfms.a -t $exec_dir/intel.mk -c "$cppDefs" pathnames_fms shared/mpp/include shared/include 

# ---------------- adjust the main Makefile

#cat Makefile | sed -e 's/<TAB>/\t/' > Makefile.$$ && mv -f Makefile.$$ Makefile

# ---------------- call make on the main Makefile

make  NETCDF=3 fms_CM2.1R_exec_v5.0_ps.x

if ( $status ) then
  unset echo
  echo ERROR: make failed for CM2.1R_exec_ECDA_v5.0_ps
  exit 1
else
  unset echo
  echo NOTE: make succeeded for CM2.1R_exec_ECDA_v5.0_ps
endif
