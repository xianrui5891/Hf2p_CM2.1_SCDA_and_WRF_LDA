module ocean_tracer_util_mod
!
! <CONTACT EMAIL="Stephen.Griffies@noaa.gov">
! S. M. Griffies 
! </CONTACT>
!
!<CONTACT EMAIL="Ronald.Pacanowski@noaa.gov"> Ron Pacanowski
!</CONTACT>
!
!<OVERVIEW>
! This module contains many routines of use for tracer diagnostics in mom4. 
!</OVERVIEW>
!
!<DESCRIPTION>
! Tracer utility module for mom4. Of use for tracer diagnostics. 
!</DESCRIPTION>
!
use constants_mod,        only: epsln
use mpp_mod,              only: stdout, stdlog, FATAL
use mpp_mod,              only: mpp_error, mpp_chksum, mpp_pe, mpp_min, mpp_max
use platform_mod,         only: i8_kind
  
use ocean_domains_mod,    only: get_local_indices
use ocean_parameters_mod, only: ADVECT_PSOM
use ocean_types_mod,      only: ocean_grid_type, ocean_domain_type
use ocean_types_mod,      only: ocean_prog_tracer_type, ocean_diag_tracer_type
use ocean_types_mod,      only: ocean_time_type, ocean_thickness_type 
use ocean_types_mod,      only: ocean_density_type
use ocean_util_mod,       only: write_timestamp
use ocean_workspace_mod,  only: wrk1 

implicit none

private

character(len=256) :: version='CVS $Id'
character(len=256) :: tagname='Tag $Name'

! for output
integer :: unit=6


#include <ocean_memory.h>

logical :: module_is_initialized = .false.

type(ocean_domain_type), pointer :: Dom =>NULL()
type(ocean_grid_type), pointer   :: Grd =>NULL()

public ocean_tracer_util_init
public tracer_min_max
public dzt_min_max
public tracer_prog_chksum 
public tracer_diag_chksum 
public tracer_psom_chksum
public sort_pick_array
public sort_shell_array
public remap_onto_nrho


contains


!#######################################################################
! <SUBROUTINE NAME="ocean_tracer_util_init">
!
! <DESCRIPTION>
! Initialize mom4 tracer utilities.
! </DESCRIPTION>
!
subroutine ocean_tracer_util_init (Grid, Domain)

  type(ocean_grid_type),   intent(in), target :: Grid
  type(ocean_domain_type), intent(in), target :: Domain
  integer :: stdlogunit

  if (module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error in ocean_tracer_util_mod (ocean_tracer_util_init): module already initialized')
  endif 

  module_is_initialized = .true.
  stdlogunit=stdlog()
  write( stdlogunit,'(/a/)') trim(version)

  Dom => Domain
  Grd => Grid

#ifndef MOM4_STATIC_ARRAYS
  call get_local_indices(Domain, isd, ied, jsd, jed, isc, iec, jsc, jec)
  nk=Grd%nk
#endif


end subroutine ocean_tracer_util_init
! </SUBROUTINE> NAME="ocean_tracer_util_init">



!#######################################################################
! <SUBROUTINE NAME="tracer_min_max">
!
! <DESCRIPTION>
! Compute the global min and max for tracers.  
!
! Vectorized using maxloc() and minloc() intrinsic functions by 
! Russell.Fiedler@csiro.au (May 2005).
!
! Modified by Zhi.Liang@noaa.gov (July 2005)
!          
! </DESCRIPTION>
!
subroutine tracer_min_max(Time, Thickness, Tracer)
  
  type(ocean_time_type),        intent(in) :: Time
  type(ocean_thickness_type),   intent(in) :: Thickness
  type(ocean_prog_tracer_type), intent(in) :: Tracer

  real    :: tmax, tmin, tmax0, tmin0
  integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: tau 
  real    :: fudge

  ! arrays to enable vectorization
  integer :: iminarr(3),imaxarr(3)

  if (.not. module_is_initialized) then
     call mpp_error(FATAL,&
     '==>Error from ocean_tracer_util_mod (tracer_min_max): module not initialized')
  endif 

  tmax=-1.e10;tmin=1.e10
  itmax=0;jtmax=0;ktmax=0
  itmin=0;jtmin=0;ktmin=0

  tau = Time%tau

  call write_timestamp(Time%model_time)
  
  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:,tau)
  
  if(ANY(Grd%tmask(isc:iec,jsc:jec,:) > 0.)) then
     iminarr=minloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     imaxarr=maxloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     itmin=iminarr(1)+isc-1
     jtmin=iminarr(2)+jsc-1
     ktmin=iminarr(3)
     itmax=imaxarr(1)+isc-1
     jtmax=imaxarr(2)+jsc-1 
     ktmax=imaxarr(3)
     tmin=wrk1(itmin,jtmin,ktmin)
     tmax=wrk1(itmax,jtmax,ktmax)
  end if

  ! use "fudge" to distinguish processors when tracer extreme is independent of processor
  fudge = 1.0 + 1.e-12*mpp_pe() 
  tmax = tmax*fudge
  tmin = tmin*fudge
  if(tmax == 0.0) then 
    tmax = tmax + 1.e-12*mpp_pe() 
  endif 
  if(tmin == 0.0) then 
    tmin = tmin + 1.e-12*mpp_pe() 
  endif 
  

  tmax0=tmax;tmin0=tmin

  call mpp_max(tmax)
  call mpp_min(tmin)

  if (tmax0 == tmax) then
      if (trim(Tracer%name) == 'temp') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The maximum T is ', tmax,&
               ' deg C at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
            Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
            Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod: The maximum temperature is outside allowable range')
          endif
      else if (trim(Tracer%name) == 'salt') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The maximum S is ', tmax,&
               ' psu at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
            Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
            Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL,&
              '==>Error from ocean_tracer_util_mod: The maximum salinity is outside allowable range')
          endif
      else
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The maximum '//trim(Tracer%name)// ' is ', tmax,&
                ' '// trim(Tracer%units)// ' at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,','&
               ,ktmax,'),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
               Grd%yt(itmax,jtmax),',', Grd%zt(ktmax),' m)'
          write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
            Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
             Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
          if (tmax > Tracer%max_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod: The maximum tracer is outside allowable range')
          endif 
      endif
  endif
  
  if (tmin0 == tmin) then
      if (trim(Tracer%name) == 'temp') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The minimum T is ', tmin,&
               ' deg C at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
            Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
             Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod (tracer_min_max): minimum temp outside allowable range')
          endif
      else if (trim(Tracer%name) == 'salt') then
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The minimum S is ', tmin,&
               ' psu at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
               '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
            Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
             Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod: The minimum salinity is outside allowable range')
          endif
      else
          write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
               ' The minimum '//trim(Tracer%name)// ' is ', tmin,&
                ' '// trim(Tracer%units)// ' at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,','&
               ,ktmin,'),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
               Grd%yt(itmin,jtmin),',', Grd%zt(ktmin),' m)'
          write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
            Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
          write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
             Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)
          write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
          if (tmin < Tracer%min_tracer) then
              call mpp_error(FATAL, &
               '==>Error from ocean_tracer_util_mod (tracer_min_max): minimum tracer outside allowable range')
          endif
      endif
  endif

  return


end subroutine tracer_min_max
! </SUBROUTINE>  NAME="tracer_min_max"


!#######################################################################
! <SUBROUTINE NAME="dzt_min_max">
!
! <DESCRIPTION>
! Compute the global min and max for dzt.  
!
! Modified by Stephen.Griffies@noaa.gov from subroutine tracer_min_max
!          
! </DESCRIPTION>
!
subroutine dzt_min_max(Time, Thickness, filecaller)
  
  type(ocean_time_type),      intent(in) :: Time
  type(ocean_thickness_type), intent(in) :: Thickness
  character(len=*),           intent(in) :: filecaller

  real    :: dztmax, dztmin, dztmax0, dztmin0
  integer :: itmax, jtmax, ktmax, itmin, jtmin, ktmin
  integer :: tau
  real    :: fudge

  ! arrays to enable vectorization
  integer :: iminarr(3),imaxarr(3)

  tau = Time%tau

  dztmax=-1.e10;dztmin=1.e10
  itmax=0;jtmax=0;ktmax=0
  itmin=0;jtmin=0;ktmin=0

  call write_timestamp(Time%model_time)
  
  wrk1(isc:iec,jsc:jec,:) = Thickness%dzt(isc:iec,jsc:jec,:)
  
  if(ANY(Grd%tmask(isc:iec,jsc:jec,:) > 0.)) then
     iminarr=minloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     imaxarr=maxloc(wrk1(isc:iec,jsc:jec,:),Grd%tmask(isc:iec,jsc:jec,:) > 0.)
     itmin=iminarr(1)+isc-1
     jtmin=iminarr(2)+jsc-1
     ktmin=iminarr(3)
     itmax=imaxarr(1)+isc-1
     jtmax=imaxarr(2)+jsc-1 
     ktmax=imaxarr(3)
     dztmin=wrk1(itmin,jtmin,ktmin)
     dztmax=wrk1(itmax,jtmax,ktmax)
  end if

  ! use "fudge" to distinguish processors when extreme is independent of processor
  fudge = 1.0 + 1.e-12*mpp_pe() 
  dztmax = dztmax*fudge
  dztmin = dztmin*fudge
  if(dztmax == 0.0) then 
    dztmax = dztmax + 1.e-12*mpp_pe() 
  endif 
  if(dztmin == 0.0) then 
    dztmin = dztmin + 1.e-12*mpp_pe() 
  endif 

  dztmax0=dztmax;dztmin0=dztmin

  call mpp_max(dztmax)
  call mpp_min(dztmin)

  if (dztmax0 == dztmax) then

      write(unit,'(/a)') trim(filecaller)
      call write_timestamp(Time%model_time)

      write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
           ' The maximum dzt is ', dztmax,&
           ' metre at (i,j,k) = (',itmax+Dom%ioff,',',jtmax+Dom%joff,',',ktmax,&
           '),  (lon,lat,dpt) = (',Grd%xt(itmax,jtmax),',', &
           Grd%yt(itmax,jtmax),',', Thickness%depth_zt(itmax,jtmax,ktmax),' m)'
      write(unit,'(a,es22.12,a,es22.12,a,es22.12)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
           Grd%dxt(itmax,jtmax),',', Grd%dyt(itmax,jtmax),',', Thickness%dzt(itmax,jtmax,ktmax),' m)'
      write(unit,'(a,i6)') ' The number of cells in the column are kmt = ', Grd%kmt(itmax,jtmax)
      write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
           Thickness%dst(itmax,jtmax,ktmax),', ',Thickness%rho_dzt(itmax,jtmax,ktmax,tau)
  endif

  if (dztmin0 == dztmin) then
      write(unit,'(/,a,es22.12,a,i4,a1,i4,a1,i3,a,f10.4,a,f10.4,a,f10.4,a)') &
           ' The minimum dzt is ', dztmin,&
           ' metre at (i,j,k) = (',itmin+Dom%ioff,',',jtmin+Dom%joff,',',ktmin,&
           '),  (lon,lat,dpt) = (',Grd%xt(itmin,jtmin),',', &
           Grd%yt(itmin,jtmin),',', Thickness%depth_zt(itmin,jtmin,ktmin),' m)'
      write(unit,'(a,es22.12,a,es22.12,a,es22.12,a)')      &
           ' The grid dimensions are (dxt, dyt, dzt) = (', &
           Grd%dxt(itmin,jtmin),',', Grd%dyt(itmin,jtmin),',', Thickness%dzt(itmin,jtmin,ktmin),' m)'
      write(unit,'(a,i6)') ' And the number of cells in the column are kmt = ', Grd%kmt(itmin,jtmin)
      write(unit,'(a,es22.12,a,es22.12)') ' The grid dimensions (dst,rho_dzt) = ', &
           Thickness%dst(itmin,jtmin,ktmin),', ',Thickness%rho_dzt(itmin,jtmin,ktmin,tau)

      if(dztmin < 0.0) then 
          write(unit,'(a)') '==>Error in ocean_tracer_util_mod (dzt_min_max): dzt < 0 detected.'
          call mpp_error(FATAL,&
               '==>Error in ocean_tracer_util_mod (dzt_min_max): dzt < 0 detected.')
      endif

  endif


  return


end subroutine dzt_min_max
! </SUBROUTINE>  NAME="dzt_min_max"


!#######################################################################
! <SUBROUTINE NAME="tracer_prog_chksum">
!
! <DESCRIPTION>
! Compute checksums for prognostic tracers 
! </DESCRIPTION>
subroutine tracer_prog_chksum(Time, Tracer, index, chksum)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: Tracer
  integer,                      intent(in)  :: index
  integer(i8_kind), optional, intent(inout) :: chksum
  integer(i8_kind)                          :: chk_sum
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (tracer_prog_chksum): module not initialized ')
  endif 

  write(stdoutunit,*) ' '
  write(stdoutunit,*) '=== Prognostic tracer checksum follows ==='
  write(stdoutunit,*) 'Tracer name = ', Tracer%name

  call write_timestamp(Time%model_time)

  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:,index)*Grd%tmask(isc:iec,jsc:jec,:)

  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  
  write(stdoutunit,*) 'Tracer chksum = ',  chk_sum

  if (PRESENT(chksum)) chksum = chk_sum

end subroutine tracer_prog_chksum
! </SUBROUTINE>  NAME="tracer_prog_chksum"


!#######################################################################
! <SUBROUTINE NAME="tracer_diag_chksum">
!
! <DESCRIPTION>
! Compute checksums for diagnostic tracers 
! </DESCRIPTION>
subroutine tracer_diag_chksum(Time, Tracer, chksum)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_diag_tracer_type), intent(in)    :: Tracer
  integer(i8_kind), optional,   intent(inout) :: chksum
  integer(i8_kind)                            :: chk_sum
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (tracer_diag_chksum): module not initialized ')
  endif 

  write(stdoutunit,*) '=== Diagnostic tracer checksum follows ==='
  write(stdoutunit,*) 'Tracer name = ', Tracer%name

  call write_timestamp(Time%model_time)

  wrk1(isc:iec,jsc:jec,:) = Tracer%field(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)

  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  
  write(stdoutunit,*) 'Tracer chksum = ',  chk_sum

  if (PRESENT(chksum)) chksum = chk_sum

end subroutine tracer_diag_chksum
! </SUBROUTINE>  NAME="tracer_diag_chksum"


!#######################################################################
! <SUBROUTINE NAME="tracer_psom_chksum">
!
! <DESCRIPTION>
! Compute checksums for PSOM advection second order moments. 
! </DESCRIPTION>
subroutine tracer_psom_chksum(Time, Tracer)

  type(ocean_time_type),        intent(in)  :: Time
  type(ocean_prog_tracer_type), intent(in)  :: Tracer
  integer(i8_kind)                          :: chk_sum
  
  integer :: stdoutunit 
  stdoutunit=stdout() 

  if (.not. module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (tracer_psom_chksum): module not initialized ')
  endif 

  write(stdoutunit,*) ' '
  write (stdoutunit,*) 'Writing psom moments for tracer ',trim(Tracer%name)
  call write_timestamp(Time%model_time)

  wrk1(isc:iec,jsc:jec,:) = Tracer%s0(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%s0 chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%sx(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%sx chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%sxx(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%sxx chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%sy(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%sy chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%syy(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%syy chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%sz(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%sz chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%szz(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%szz chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%sxy(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%sxy chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%sxy(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%sxz chksum = ',  chk_sum

  wrk1(isc:iec,jsc:jec,:) = Tracer%syz(isc:iec,jsc:jec,:)*Grd%tmask(isc:iec,jsc:jec,:)
  chk_sum = mpp_chksum(wrk1(isc:iec,jsc:jec,:))
  write(stdoutunit,*) 'Tracer%syz chksum = ',  chk_sum


end subroutine tracer_psom_chksum
! </SUBROUTINE>  NAME="tracer_psom_chksum"



!#######################################################################
! <SUBROUTINE NAME="sort_pick_array">
!
! <DESCRIPTION>
! Simplest, and slowest, sorting algorithm from Numerical Recipes.
! Called "sort_pick" in Numerical Recipes.  
!
! Input are two arrays, first array defines the ascending sort 
! and second is a slave to the sort. 
!
! Typical example is sorting a vector of water parcels lightest
! to densest, with slave being volume of the parcels.  
!
! More sophisticated sorting algorithms exist, and may need to 
! be coded should this method prove too slow. 
!
! This scheme has order N^2 operations, which is a lot. 
!
! output has array(1) smallest and a(nsortpts) largest 
! with corresponding slave array.  
!
! coded Stephen.Griffies@noaa.gov June 2005
!
! </DESCRIPTION>
!
subroutine sort_pick_array(array, slave)

  real, dimension(:), intent(inout) :: array
  real, dimension(:), intent(inout) :: slave

  real    :: tmp_a, tmp_s
  integer :: m, n, nsortpts

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (sort_pick_array): module not initialized')
  endif 

  nsortpts = size(array)

  do n=2,nsortpts
      tmp_a = array(n)
      tmp_s = slave(n)
     do m=n-1,1,-1
        if(array(m) <= tmp_a) exit
        array(m+1)=array(m)
        slave(m+1)=slave(m)
     enddo
     array(m+1) = tmp_a
     slave(m+1) = tmp_s    
  enddo


end subroutine sort_pick_array
! </SUBROUTINE>  NAME="sort_pick_array"



!#######################################################################
! <SUBROUTINE NAME="sort_shell_array">
!
! <DESCRIPTION>
! Shell (or diminishing increment) sort from Numerical Recipes.
! Called "sort_shell" in Numerical Recipes.  
!
! Input are two arrays, first array defines the ascending sort 
! and second is a slave to the sort array. 
!
! Typical example is sorting a vector of water parcels lightest  
! to densest, with slave being volume of the parcels.  
!
! More sophisticated sorting algorithms exist, and may need to 
! be coded should this method prove too slow. 
!
! This scheme has order N^(5/4) operations. 
!
! output has array(1) smallest and a(nsortpts) largest,
! with corresponding ordering for slave array.  
!
! coded Stephen.Griffies@noaa.gov June 2005
!
! </DESCRIPTION>
!
subroutine sort_shell_array(array, slave)

  real, dimension(:), intent(inout) :: array
  real, dimension(:), intent(inout) :: slave

  real    :: tmp_a, tmp_s
  integer :: m, n, inc, nsortpts

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL,&
    '==>Error from ocean_tracer_util_mod (sort_shell_array): module not initialized')
  endif 

  nsortpts = size(array)

  inc=1
  do 
     inc=3*inc+1
     if(inc > nsortpts) exit
  enddo

  do 
     inc=inc/3

     do m=inc+1,nsortpts
        tmp_a = array(m)
        tmp_s = slave(m)

        n=m
        do
           if(array(n-inc) <= tmp_a) exit
           array(n) = array(n-inc)
           slave(n) = slave(n-inc)
           n=n-inc
           if(n <= inc) exit
        enddo

        array(n)=tmp_a
        slave(n)=tmp_s

     enddo

     if(inc <=1) exit
  enddo


end subroutine sort_shell_array
! </SUBROUTINE>  NAME="sort_shell_array"


!#######################################################################
! <SUBROUTINE NAME="remap_onto_nrho">
!
! <DESCRIPTION>
! Bin a level input tendency field according to neutral density classes. 
! The binning is meant for tendencies and transports, as used in particular
! for the neutral_rho and wdian diagnostics.  
! </DESCRIPTION>
!
subroutine remap_onto_nrho (Time, Dens, infield_level, outfield_nrho)

  type(ocean_time_type),        intent(in)    :: Time
  type(ocean_density_type),     intent(in)    :: Dens
  real, dimension(isd:,jsd:,:), intent(in)    :: infield_level
  real, dimension(isd:,jsd:,:), intent(inout) :: outfield_nrho

  integer :: i, j, k, k_rho, neutralrho_nk
  real    :: W1, W2
  real    :: weight 

  if (.not.module_is_initialized) then 
    call mpp_error(FATAL, &
    '==>Error ocean_neutral_util_mod (remap_onto_nrho): module needs initialization ')
  endif 

  neutralrho_nk = size(Dens%neutralrho_ref(:))
  outfield_nrho(:,:,:) = 0.0


  ! interpolate infield_level to get outfield_nrho
!!$  do k_rho=1,neutralrho_nk
!!$     do k=1,nk-1
!!$        do j=jsc,jec
!!$           do i=isc,iec
!!$              if(     Dens%neutralrho_ref(k_rho) >  Dens%neutralrho(i,j,k)  ) then
!!$                  if( Dens%neutralrho_ref(k_rho) <= Dens%neutralrho(i,j,k+1)) then 
!!$                      W1= Dens%neutralrho_ref(k_rho)- Dens%neutralrho(i,j,k)
!!$                      W2= Dens%neutralrho(i,j,k+1)  - Dens%neutralrho_ref(k_rho)
!!$                      outfield_nrho(i,j,k_rho) = (infield_level(i,j,k+1)*W1 +infield_level(i,j,k)*W2) &
!!$                                                 /(W1 + W2 + epsln)
!!$                  endif
!!$              endif
!!$           enddo
!!$        enddo
!!$     enddo
!!$  enddo


! rebin infield_level to get outfield_nrho
      do k=1,nk
         do j=jsc,jec
            do i=isc,iec
               if(Dens%neutralrho(i,j,k) <= Dens%neutralrho_ref(1)) then 
                   outfield_nrho(i,j,1) = outfield_nrho(i,j,1) + infield_level(i,j,k)
               endif
               if(Dens%neutralrho_ref(neutralrho_nk) < Dens%neutralrho(i,j,k)) then 
                   outfield_nrho(i,j,neutralrho_nk) = outfield_nrho(i,j,neutralrho_nk) + infield_level(i,j,k)             
               endif
            enddo
         enddo
      enddo

      do k_rho=1,neutralrho_nk-1
         do k=1,nk
            do j=jsc,jec
               do i=isc,iec
                  if(     Dens%neutralrho_ref(k_rho) <  Dens%neutralrho(i,j,k)  ) then
                      if( Dens%neutralrho(i,j,k)     <= Dens%neutralrho_ref(k_rho+1)) then 
                          weight = (Dens%neutralrho(i,j,k)-Dens%neutralrho_ref(k_rho)) &
                                  /(Dens%neutralrho_ref(k_rho+1)-Dens%neutralrho_ref(k_rho))
                          outfield_nrho(i,j,k_rho)   = outfield_nrho(i,j,k_rho)   + (1.0-weight)*infield_level(i,j,k)
                          outfield_nrho(i,j,k_rho+1) = outfield_nrho(i,j,k_rho+1) + weight*infield_level(i,j,k)
                      endif
                  endif
               enddo
            enddo
         enddo
      enddo

!!$      do k_rho=2,neutralrho_nk-1
!!$         do k=1,nk
!!$            do j=jsc,jec
!!$               do i=isc,iec
!!$                  if(     Dens%neutralrho_ref(k_rho) <  Dens%neutralrho(i,j,k)  ) then
!!$                      if( Dens%neutralrho(i,j,k)     <= Dens%neutralrho_ref(k_rho+1)) then 
!!$                          weight = (Dens%neutralrho(i,j,k)-Dens%neutralrho_ref(k_rho)) &
!!$                                  /(Dens%neutralrho_ref(k_rho+1)-Dens%neutralrho_ref(k_rho))
!!$                          outfield_nrho(i,j,k_rho)   = outfield_nrho(i,j,k_rho)   + weight*infield_level(i,j,k)
!!$                          outfield_nrho(i,j,k_rho+1) = outfield_nrho(i,j,k_rho+1) + (1.0-weight)*infield_level(i,j,k)
!!$                      endif
!!$                  endif
!!$               enddo
!!$            enddo
!!$         enddo
!!$      enddo
!!$

  do k_rho=1,neutralrho_nk
     do j=jsc,jec
        do i=isc,iec
           outfield_nrho(i,j,k_rho) = outfield_nrho(i,j,k_rho)*Grd%tmask(i,j,1)
        enddo
     enddo
  enddo

end subroutine remap_onto_nrho
! </SUBROUTINE> NAME="remap_onto_nrho"



end module ocean_tracer_util_mod

