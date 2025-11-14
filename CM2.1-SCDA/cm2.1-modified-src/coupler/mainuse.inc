!
!  coupler_main couples component models and controls the time integration
!
!program coupler_main
! mainuse.inc start>
!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
! <CONTACT EMAIL="Bruce.Wyman@noaa.gov"> Bruce Wyman </CONTACT>
! <CONTACT EMAIL="V.Balaji@noaa.gov"> V. Balaji </CONTACT>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!  A main program that couples component models for atmosphere, ocean, land, 
!  and sea ice on independent grids. 
! </OVERVIEW>

! <DESCRIPTION>
!  This version couples model components representing atmosphere, ocean, land 
!  and sea ice on independent grids. Each model component is represented by a 
!  data type giving the instantaneous model state.
!
!  The component models are coupled to allow implicit vertical diffusion of 
!  heat and moisture at the interfaces of the atmosphere, land, and ice models. 
!  As a result, the atmosphere, land, and ice models all use the same time step. 
!  The atmospheric model has been separated into down and up calls that 
!  correspond to the down and up sweeps of the standard tridiagonal elimination.
!
!  The ocean interface uses explicit mixing. Fluxes to and from the ocean must
!  be passed through the ice model. This includes atmospheric fluxes as well as 
!  fluxes from the land to the ocean (runoff).
!
!  This program contains the model's main time loop. Each iteration of the 
!  main time loop is one coupled (slow) time step. Within this slow time step 
!  loop is a fast time step loop, using the atmospheric time step, where the
!  tridiagonal vertical diffusion equations are solved. Exchange between sea 
!  ice and ocean occurs once every slow timestep.
!
! <PRE>
!      MAIN PROGRAM EXAMPLE
!      --------------------
!
!         DO slow time steps (ocean)
!
!              call flux_ocean_to_ice
!
!              call ICE_SLOW_UP
!
!              DO fast time steps (atmos)
!
!                   call flux_calculation
!
!                   call ATMOS_DOWN
!
!                   call flux_down_from_atmos
!
!                   call LAND_FAST
!
!                   call ICE_FAST
!
!                   call flux_up_to_atmos
!
!                   call ATMOS_UP
!
!              END DO
!
!              call ICE_SLOW_DN
!
!              call flux_ice_to_ocean
!
!              call OCEAN
!
!         END DO

!  </PRE>

! </DESCRIPTION>
! <INFO>
!   <NOTE>
!     <PRE>
!   1.If no value is set for current_date, start_date, or calendar (or default value 
!     specified) then the value from restart file "INPUT/coupler.res" will be used. 
!     If neither a namelist value or restart file value exist the program will fail. 
!   2.The actual run length will be the sum of months, days, hours, minutes, and 
!     seconds. A run length of zero is not a valid option. 
!   3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>

!   <ERROR MSG="no namelist value for current_date " STATUS="FATAL">
!     A namelist value for current_date must be given if no restart file for
!     coupler_main (INPUT/coupler.res) is found. 
!   </ERROR>
!   <ERROR MSG="invalid namelist value for calendar" STATUS="FATAL">
!     The value of calendar must be 'julian', 'noleap', or 'thirty_day'. 
!     See the namelist documentation. 
!   </ERROR>
!   <ERROR MSG="no namelist value for calendar" STATUS="FATAL">
!     If no restart file is present, then a namelist value for calendar 
!     must be specified. 
!   </ERROR>
!   <ERROR MSG="initial time is greater than current time" STATUS="FATAL">
!     If a restart file is present, then the namelist value for either 
!     current_date or start_date was incorrectly set. 
!   </ERROR>
!   <ERROR MSG="run length must be multiple of ocean time step " STATUS="FATAL">
!     There must be an even number of ocean time steps for the requested run length. 
!   </ERROR>
!   <ERROR MSG="final time does not match expected ending time " STATUS="WARNING">
!     This error should probably not occur because of checks done at initialization time. 
!   </ERROR>

! </INFO>

  use constants_mod,           only: constants_init

  use time_manager_mod,        only: time_type, set_calendar_type, set_time
  use time_manager_mod,        only: set_date, get_date, days_in_month, month_name
  use time_manager_mod,        only: operator(+), operator(-), operator (<)
  use time_manager_mod,        only: operator (>), operator ( /= ), operator ( / )
  use time_manager_mod,        only: operator (*), THIRTY_DAY_MONTHS, JULIAN
  use time_manager_mod,        only: NOLEAP, NO_CALENDAR, INVALID_CALENDAR
  use time_manager_mod,        only: date_to_string, increment_date
  use time_manager_mod,        only: operator(>=), operator(<=), operator(==)

  use fms_mod,                 only: open_namelist_file, field_exist, file_exist, check_nml_error
  use fms_mod,                 only: uppercase, error_mesg, write_version_number
  use fms_mod,                 only: fms_init, fms_end, stdout
  use fms_mod,                 only: read_data, write_data

  use fms_io_mod,              only: fms_io_exit
  use fms_io_mod,              only: restart_file_type, register_restart_field, save_restart

  use diag_manager_mod,        only: diag_manager_init, diag_manager_end, diag_grid_end
  use diag_manager_mod,        only: DIAG_OCEAN, DIAG_OTHER, DIAG_ALL, get_base_date

  use field_manager_mod,       only: MODEL_ATMOS, MODEL_LAND, MODEL_ICE

  use tracer_manager_mod,      only: tracer_manager_init, get_tracer_index
  use tracer_manager_mod,      only: get_number_tracers, get_tracer_names, NO_TRACER

  use coupler_types_mod,       only: coupler_types_init

  use data_override_mod,       only: data_override_init

!
! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!

  use atmos_model_mod,         only: atmos_model_init, atmos_model_end
  use atmos_model_mod,         only: update_atmos_model_down
  use atmos_model_mod,         only: update_atmos_model_up
  use atmos_model_mod,         only: atmos_data_type
  use atmos_model_mod,         only: land_ice_atmos_boundary_type
  use atmos_model_mod,         only: atmos_model_restart
  use ada_types_mod,           only: scda_type ! snz 202250426

  use land_model_mod,          only: land_model_init, land_model_end
  use land_model_mod,          only: land_data_type, atmos_land_boundary_type
  use land_model_mod,          only: update_land_model_fast, update_land_model_slow
  use land_model_mod,          only: land_model_restart

  use ice_model_mod,           only: ice_model_init, ice_model_end
  use ice_model_mod,           only: update_ice_model_slow_up
  use ice_model_mod,           only: update_ice_model_fast
  use ice_model_mod,           only: update_ice_model_slow_dn
  use ice_model_mod,           only: ice_data_type, land_ice_boundary_type
  use ice_model_mod,           only: ocean_ice_boundary_type, atmos_ice_boundary_type
  use ice_model_mod,           only: ice_model_restart

  use ocean_model_mod,         only: update_ocean_model, ocean_model_init
  use ocean_model_mod,         only: ocean_model_end, ocean_public_type, ocean_state_type, ice_ocean_boundary_type
  use ocean_model_mod,         only: ocean_model_restart
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!

  use flux_exchange_mod,       only: flux_exchange_init
  use flux_exchange_mod,       only: sfc_boundary_layer
  use flux_exchange_mod,       only: generate_sfc_xgrid
  use flux_exchange_mod,       only: flux_down_from_atmos
  use flux_exchange_mod,       only: flux_up_to_atmos
  use flux_exchange_mod,       only: flux_land_to_ice
  use flux_exchange_mod,       only: flux_ice_to_ocean
  use flux_exchange_mod,       only: flux_ocean_to_ice
  use flux_exchange_mod,       only: flux_check_stocks, flux_init_stocks
  use flux_exchange_mod,       only: flux_ice_to_ocean_stocks, flux_ocean_from_ice_stocks

  use atmos_tracer_driver_mod, only: atmos_tracer_driver_gather_data

  use mpp_mod,                 only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_chksum
  use mpp_mod,                 only: mpp_init, mpp_pe, mpp_npes, mpp_root_pe, MAXPES
  use mpp_mod,                 only: stderr, stdlog, mpp_error, NOTE, FATAL, WARNING
  use mpp_mod,                 only: mpp_set_current_pelist, mpp_declare_pelist

  use mpp_io_mod,              only: mpp_open, mpp_close, mpp_io_clock_on
  use mpp_io_mod,              only: MPP_NATIVE, MPP_RDONLY, MPP_DELETE

  use mpp_domains_mod,         only: mpp_broadcast_domain

  use memutils_mod,            only: print_memuse_stats

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.F90,v 18.0 2010/03/02 23:35:50 fms Exp $'
  character(len=128) :: tag = '$Name: riga $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land
  type   (ice_data_type) :: Ice
  ! allow members of ocean type to be aliased (ap)
  type (ocean_public_type), target :: Ocean
  type (ocean_state_type),  pointer :: Ocean_state => NULL()

  type(atmos_land_boundary_type)     :: Atmos_land_boundary
  type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
  type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
  type(land_ice_boundary_type)       :: Land_ice_boundary
  type(ice_ocean_boundary_type)      :: Ice_ocean_boundary
  type(ocean_ice_boundary_type)      :: Ocean_ice_boundary

  type(scda_type)      :: CBVar ! snz 20250426

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type (time_type) :: Time, Time_init, Time_end, &
                      Time_step_atmos, Time_step_cpld
  type(time_type) :: Time_atmos, Time_ocean
  integer :: num_atmos_calls, na
  integer :: num_cpld_calls, nc

!------ for intermediate restart
  type(restart_file_type), allocatable :: Ice_bc_restart(:), Ocn_bc_restart(:)
  character(len=64),       allocatable :: ice_bc_restart_file(:), ocn_bc_restart_file(:) 
  integer                              :: num_ice_bc_restart=0, num_ocn_bc_restart=0
  type(time_type)                      :: Time_restart, Time_restart_current, Time_start
  character(len=32)                    :: timestamp

! ----- coupled model initial date -----

  integer :: date_init(6)
  integer :: calendar_type = INVALID_CALENDAR

!-----------------------------------------------------------------------
!------ namelist interface -------

! <NAMELIST NAME="coupler_nml">
!   <DATA NAME="current_date"  TYPE="integer, dimension(6)"  DEFAULT="0">
!     The date that the current integration starts with. 
!   </DATA>
!   <DATA NAME="force_date_from_namelist"  TYPE="logical"  DEFAULT=".false.">
!     Flag that determines whether the namelist variable current_date should 
!     override the date in the restart file INPUT/coupler.res. If the restart 
!     file does not exist then force_date_from_namelist has not effect, the value of current_date 
!     will be used.
!   </DATA>
!   <DATA NAME="calendar"  TYPE="character(maxlen=17)"  DEFAULT="''">
!     The calendar type used by the current integration. Valid values are consistent 
!     with the time_manager module: 'julian', 'noleap', or 'thirty_day'. The value 
!     'no_calendar' can not be used because the time_manager's date  function are used. 
!     All values must be lowercase.
!   </DATA>
!   <DATA NAME="months "  TYPE="integer"  DEFAULT="0">
!     The number of months that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="days "  TYPE="integer"  DEFAULT="0">
!     The number of days that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="hours"  TYPE="integer"  DEFAULT="0">
!     The number of hours that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="minutes "  TYPE="integer"  DEFAULT="0">
!     The number of minutes that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="seconds"  TYPE="integer"  DEFAULT="0">
!     The number of seconds that the current integration will be run for. 
!   </DATA>
!   <DATA NAME="dt_atmos"  TYPE="integer"  DEFAULT="0">
!     Atmospheric model time step in seconds, including the fast coupling with 
!     land and sea ice. 
!   </DATA>
!   <DATA NAME="dt_cpld"  TYPE="integer"  DEFAULT="0">
!     Time step in seconds for coupling between ocean and atmospheric models: 
!     must be an integral multiple of dt_atmos and dt_ocean. This is the "slow" timestep.
!   </DATA>
!  <DATA NAME="do_atmos, do_ocean, do_ice, do_land, do_flux" TYPE="logical">
!  If true (default), that particular model component (atmos, etc.) is run.
!  If false, the execution of that component is skipped. This is used when
!  ALL the output fields sent by that component to the coupler have been
!  overridden using the data_override feature. For advanced users only:
!  if you're not sure, you should leave these values at TRUE.
!  </DATA> 
!  <DATA NAME="concurrent" TYPE="logical">
!  If true, the ocean executes concurrently with the atmosphere-land-ocean
!   on a separate set of PEs.
!  If false (default), the execution is serial: call atmos... followed by
!  call ocean...
!  If using concurrent execution, you must set one of
!   atmos_npes and ocean_npes, see below.
!  </DATA> 
!  <DATA NAME="atmos_npes, ocean_npes" TYPE="integer">
!  If concurrent is set to true, we use these to set the list of PEs on which
!   each component runs.
!  At least one of them must be set to a number between 0 and NPES.
!  If exactly one of these two is set non-zero, the other is set to the
!   remainder from NPES.
!  If both are set non-zero they must add up to NPES.
!  </DATA> 
!  <DATA NAME="use_lag_fluxes" TYPE="logical">
!  If true, then mom4 is forced with SBCs from one coupling timestep ago
!  If false, then mom4 is forced with most recent SBCs.
!  For a leapfrog MOM coupling with dt_cpld=dt_ocean, lag fluxes
!  can be shown to be stable and current fluxes to be unconditionally unstable.
!  For dt_cpld>dt_ocean there is probably sufficient damping.
!  use_lag_fluxes is set to TRUE by default.
!  </DATA>
!  <DATA NAME="n_mask" TYPE="integer">
!    number of region to be masked out. Its value should be less than MAX_PES.
!  </DATA>
!  <DATA NAME="mask_list(2,MAXPES)" TYPE="integer, dimension(2,MAX_MASK_REGION)">
!    The position of the region to be masked out. mask_list(1,:) is the x-layout position
!    and mask_list(2,:) is y-layout position.  
!  </DATA>
!  <DATA NAME="layout_mask" TYPE="integer, dimension(2)">
!   Processor domain layout for all the component model. layout_mask need to be set when and only 
!   when n_mask is greater than 0 ( some domain region is masked out ). When this namelist is set,
!   it will overload the layout in each component model. The default value is (0,0).
!   Currently we require all the component model has the same layout and same grid size.
!  </DATA>
!  <DATA NAME="restart_interval" TYPE="integer, dimension(6)"  DEFAULT="0">
!     The time interval that write out intermediate restart file. The format is (yr,mo,day,hr,min,sec).
!     When restart_interval is all zero, no intermediate restart file will be written out.
!   </DATA>
!   <NOTE>
!     <PRE>
!     1.If no value is set for current_date, start_date, or calendar (or default value specified) then the value from restart
!       file "INPUT/coupler.res" will be used. If neither a namelist value or restart file value exist the program will fail. 
!     2.The actual run length will be the sum of months, days, hours, minutes, and seconds. A run length of zero is not a
!       valid option. 
!     3.The run length must be an intergal multiple of the coupling timestep dt_cpld. 
!     </PRE>
!   </NOTE>
! </NAMELIST>

  integer, dimension(6) :: restart_interval = (/ 0, 0, 0, 0, 0, 0/)
  integer, dimension(6) :: current_date     = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '
  logical :: force_date_from_namelist = .false.  ! override restart values for date
  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & ice/land
  integer :: dt_cpld  = 0  ! fluxes passed between ice & ocean
  integer :: atm_step =1 !lulv 20200608
  integer :: atm_end_flag =0 !lulv 20201111
  integer :: ocn_end_flag =0 !lulv 20201118

  integer ::atmos_npes=0, ocean_npes=0, ice_npes=0, land_npes=0
  logical :: do_atmos =.true., do_land =.true., do_ice =.true., do_ocean=.true.
  logical :: do_flux =.true.
  logical :: concurrent=.FALSE.
  logical :: use_lag_fluxes=.TRUE.
  logical :: do_chksum=.FALSE.
  integer :: layout_mask(2) = (/0 , 0/)
  integer :: n_mask = 0
  integer :: mask_list(2, MAXPES), n, m 
  integer :: check_stocks = 0 ! -1: never 0: at end of run only n>0: every n coupled steps
  integer, parameter :: mp = 2*MAXPES
  data ((mask_list(n,m),n=1, 2),m=1,MAXPES) /mp*0/

  namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, months, days, hours,      &
                         minutes, seconds, dt_cpld, dt_atmos, do_atmos,              &
                         do_land, do_ice, do_ocean, do_flux, atmos_npes, ocean_npes, &
                         ice_npes, land_npes, concurrent, use_lag_fluxes, do_chksum, &
                         n_mask, layout_mask, mask_list, check_stocks, restart_interval

  integer :: initClock, mainClock, termClock

  character(len=80) :: text
  character(len=48), parameter                    :: mod_name = 'coupler_main_mod'
 
  integer :: ensemble_id = 1 , outunit
  integer, allocatable :: ensemble_pelist(:, :) 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_main'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!#######################################################################

! mainuse.inc end>
