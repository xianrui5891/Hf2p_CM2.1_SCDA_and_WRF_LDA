# Makefile for Experiment 'CM2.1U_Control-1990_E1.M_3A'
include /home/ouc/zhuxianrui/CM2-f2py/CM2.1-python/exec/intel.mk

fms_CM2.1R_exec_v5.0_ps.x: libcoupler.a libland.a libecda.a libfms.a
	$(LD) $^ $(LDFLAGS) -o $@ 

libfms.a:  FORCE
	make  VERBOSE=TRUE NETCDF=3 -f Makefile.fms $@

libecda.a: libfms.a  FORCE
	make  VERBOSE=TRUE NETCDF=3 -f Makefile.ecda $@

libland.a: libfms.a FORCE
	make VERBOSE=TRUE NETCDF=3 -f Makefile.land $@

libcoupler.a: libfms.a libecda.a libland.a FORCE
	make VERBOSE=TRUE NETCDF=3 -f Makefile.coupler $@

FORCE:

clean:
	make  NETCDF=3 -f Makefile.ecda clean
	make  NETCDF=3 -f Makefile.fms clean
	make  NETCDF=3 -f Makefile.land clean
	make  NETCDF=3 -f Makefile.coupler clean

localize:
	make  NETCDF=3 -f Makefile.ecda localize
	make  NETCDF=3 -f Makefile.fms localize
	make  NETCDF=3 -f Makefile.land localize
	make  NETCDF=3 -f Makefile.coupler localize

