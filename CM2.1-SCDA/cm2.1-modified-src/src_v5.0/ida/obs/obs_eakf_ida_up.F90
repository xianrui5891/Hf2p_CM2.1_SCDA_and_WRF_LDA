module obs_eakf_sice_up_mod

use obs_tools_mod, only : conv_state_to_obs, obs_def_type, def_single_obs, &
   def_single_obs_end
use loc_and_dist_mod, only : loc_type
use fms_mod, only : file_exist, open_file, open_namelist_file, &
   check_nml_error, write_version_number, close_file
use mpp_mod, only : mpp_error, FATAL, mpp_pe, mpp_root_pe

private
public take_single_obs_sice, sice_obs_init_up, sice_obs_end_up, &
       get_close_grids_sice_up

integer :: no_obs = 0
integer :: no_atm_obs = 0

real, parameter :: pai = 3.1415926
real, parameter :: radian = pai/180.0

! Following is to allow initialization of obs_def_type
logical :: first_run_call = .true.

#ifdef STATIC_MEMORY

integer, parameter :: max_sice_obs = 50000 ! must be same as eakf_atm
type (obs_def_type) :: sice_obs_def(max_sice_obs)

#else

type (obs_def_type), allocatable :: sice_obs_def(:)

#endif

!=======================================================================

!---- namelist with default values
! Set a cut-off for lon and lat for close obs search
real :: close_lat_window = 10.0
real :: close_lon_window = 10.0

namelist /sice_obs_wdw_nml/ close_lat_window, close_lon_window

!--- module name and version number
character(len = 18), parameter :: module_name = 'obs:cm2.0_verssion'
character(len = 12), parameter :: vers_num = '12/18/2005'

!=======================================================================

contains

!=======================================================================

subroutine sice_obs_init_up(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, &
           sice_obs, num_obs, T_grid, list_loc_halo_obs, num_prfs_loc_halo)

use ida_types_mod, only : ida_grid_type, sice_obs_type
use ada_types_mod, only : ada_grid_type

! Initializes the description a linear observations operator. For each
! observation, a list of the state variables on which it depends and the
! coefficients for each state variable is passed to def_single_obs
! which establishes appropriate data structures.

implicit none

type(ada_grid_type), intent(in) :: T_grid
type(sice_obs_type), intent(in) :: sice_obs(:)
integer, intent(in) :: num_obs
integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy
integer, intent(inout) :: list_loc_halo_obs(:), num_prfs_loc_halo
integer :: i, j, k, k0, ii, jj, kk, kk0, blk, unit_num, i_o, lon_len
INTEGER :: state_index(8), unit, ierr, io, idx_obs
!::sdu integer :: no_obs
INTEGER :: ni, nj
real :: olon, olat, coef(6), temp_lon, frac_lon, frac_lat

ni = T_grid%ni
nj = T_grid%nj

lon_len = ied_ens-isd_ens+2*halox+1
blk = (jed_ens-jsd_ens+2*haloy+1)*lon_len

! Read namelist for run time control
if(file_exist('input.nml')) then
   unit = open_namelist_file()
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = sice_obs_wdw_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'sice_obs_wdw_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Write the namelist to a log file
unit = open_file(file = 'logfile.out', action = 'append')
call write_version_number(vers_num, module_name, unit)
write(unit, nml = sice_obs_wdw_nml)
call close_file(unit)


! Initialization for identity observations
if (mpp_pe() == mpp_root_pe() .and. first_run_call) then
  write(*, *) 'close_lat_window, close_lon_window ', close_lat_window, close_lon_window
end if

num_prfs_loc_halo = 0
list_loc_halo_obs(:) = 0
do i = 1, num_obs
   ii = sice_obs(i)%i_index_up
   jj = sice_obs(i)%j_index_up

   if (jsd_ens > 1 .and. jed_ens < nj) then ! for m-middle domains

     if ( jj >= jsd_ens-haloy .and. jj <= jed_ens+haloy ) then

       if (isd_ens > 1 .and. ied_ens < ni) then ! for z-middle domains
         if ( ii >= isd_ens-halox .and. ii <= ied_ens+halox ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       elseif (ied_ens == ni) then ! for z-rightmost domains
         if ( (ii >= isd_ens-halox .and. ii <= ni) .or. &
              (ii >= 1 .and. ii <= halox) ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       elseif (isd_ens == 1) then ! for z-leftmost domains
         if ( (ii >= 1 .and. ii <= ied_ens+halox) .or. &
              (ii >= ni-halox+1 .and. ii <= ni) ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       end if

     end if

   elseif (jed_ens == nj) then ! for m-northmost domains

     if (jj >= jsd_ens-haloy .and. jj <= nj) then

       if (ied_ens < ni .and. isd_ens > 1 ) then ! for z-middle domains
         if (ii >= isd_ens-halox .and. ii <= ied_ens+halox) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       elseif (ied_ens == ni) then ! for z-rightmost domains
         if ( (ii >= isd_ens-halox .and. ii <= ni) .or. &
              (ii >= 1 .and. ii<= halox) ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       elseif (isd_ens == 1) then ! for z-leftmost domains
         if ( (ii >= ni-halox+1 .and. ii <= ni) .or. &
              (ii >= 1 .and. ii <= ied_ens+halox) ) then
          
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       end if

     end if
       
   elseif (jsd_ens == 1) then ! for m-southmost domains

     if ( jj >= 1 .and. jj <= jed_ens+haloy ) then

       if (isd_ens > 1 .and. ied_ens < ni ) then ! for z-middle domains
         if ( ii > isd_ens-halox .and. ii < ied_ens+halox ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       elseif (ied_ens == ni) then ! for z-rightmost domains
         if ( (ii >= isd_ens-halox .and. ii <= ni ) .or. &
              (ii >= 1 .and. ii <= halox) ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       elseif (isd_ens == 1) then ! for z-leftmost domains
         if ( (ii >= 1 .and. ii <= ied_ens+halox) .or. &
              (ii >= ni-halox+1 .and. ii <= ni) ) then
           num_prfs_loc_halo = num_prfs_loc_halo + 1
           list_loc_halo_obs(num_prfs_loc_halo) = i
         end if
       end if

     end if

   end if

end do

!print*,'in pe ',mpp_pe(),'num_prfs_loc,num_prfs_loc_halo=',num_prfs_loc,num_prfs_loc_halo

no_obs = 0
do i = 1, num_prfs_loc_halo
  no_obs = no_obs+1
end do

#ifndef STATIC_MEMORY
allocate(sice_obs_def(no_obs))
#endif

!print*,'pe ',mpp_pe(),'running here 0'
idx_obs = 0
do i = 1, num_prfs_loc_halo

   i_o = list_loc_halo_obs(i)
   ii = sice_obs(i_o)%i_index_up
   jj = sice_obs(i_o)%j_index_up

   if (1 < isd_ens .and. ied_ens < ni) then ! 4 i-interior sub-domain
   if (ii < (isd_ens-halox) .or. ii > (ied_ens+halox) .or. &
       jj < (jsd_ens-haloy) .or. jj > (jed_ens+haloy)) then
     if (ii < (isd_ens-halox))print*,'ii=',ii,'less than isd_ens-halox=',isd_ens-halox
     if (ii > (ied_ens+halox))print*,'ii=',ii,'greater than ied_ens+halox=',ied_ens+halox
     stop
   end if
   end if
   if (jj < 1 .or. jj > nj) then
     print*,'jj=',jj,'less than 1 or greater than nj'
     stop
   end if

   frac_lat = sice_obs(i_o)%j_index_up - jj
   frac_lon = sice_obs(i_o)%i_index_up - ii

   coef(1) = (1.0 - frac_lon) * (1.0 - frac_lat)
   coef(2) = frac_lon * (1.0 - frac_lat)
   coef(3) = (1.0 - frac_lon) * frac_lat
   coef(4) = frac_lon * frac_lat

   state_index(1) = (jj-jsd_ens+haloy)*lon_len + ii-isd_ens+halox + 1
   state_index(2) = (jj-jsd_ens+haloy)*lon_len + ii-isd_ens+halox + 2
   state_index(3) = (jj-jsd_ens+haloy+1)*lon_len + ii-isd_ens+halox + 1
   state_index(4) = (jj-jsd_ens+haloy+1)*lon_len + ii-isd_ens+halox + 2

   idx_obs = idx_obs + 1
   call def_single_obs(6, state_index(1:4), coef(1:4), sice_obs_def(idx_obs))

end do

first_run_call = .false.

end subroutine sice_obs_init_up

!------------------------------------------------------------------------

subroutine sice_obs_end_up()

implicit none

integer :: i

#ifndef STATIC_MEMORY

do i = 1, no_obs
  call def_single_obs_end(sice_obs_def(i))
end do

deallocate(sice_obs_def)

#endif

end subroutine sice_obs_end_up
!=======================================================================

function take_single_obs_sice(x, index)

implicit none
      
real :: take_single_obs_sice
real, intent(in) :: x(:)
integer, intent(in) :: index
            
real :: take(1)
         
! Given a model state, x, returns the expection of observations for 
! assimilation. For perfect model, take_obs is just state_to_obs
take = conv_state_to_obs(x, sice_obs_def(index:index), 1, index)
take_single_obs_sice = take(1)
       
end function take_single_obs_sice

subroutine get_close_grids_sice_up(obs_loc, isd_ens, ied_ens, jsd_ens, jed_ens, &
                           halox, haloy, T_grid, list, num)

use ida_types_mod, only : ida_grid_type
use ada_types_mod, only : ada_grid_type

! Computes a list of grids 'close' to the index_obs'th' atm_obs and
! returns the list. Num is the number of close grids requested for each state
! variable; list is the returned array of close grid points that contain the
! grid idx in global grid index from 1:model_size_atm. SNZ -- 12/27/05

implicit none

type(loc_type), intent(in) :: obs_loc
integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy
type(ada_grid_type), intent(in) :: T_grid
integer, intent(inout) :: list(:)
integer, intent(inout) :: num

integer :: i, j, i0, i_m, j_m, ni, nj
real :: olon, olat, low_lat, hi_lat, low_lon, hi_lon, olat0
type(loc_type) :: model_loc

ni = T_grid%ni
nj = T_grid%nj

olon = obs_loc%lon
olat = obs_loc%lat
olat0 = olat
if (abs(olat0) > 80.0) olat0 = 80.0

! Get the latitudinal arrange first

low_lat = olat - close_lat_window
if(low_lat < -90.0)low_lat = -90.0
hi_lat = olat + close_lat_window
if(hi_lat > 90.0)hi_lat = 90.0

num = 0
do j = jsd_ens-haloy, jed_ens+haloy
do i = isd_ens-halox, ied_ens+halox

  i_m = i
  j_m = j
  if (i_m <= 0) i_m = i_m + ni
  if (i_m > ni) i_m = i_m - ni
  if (j_m <= 0) j_m = 1
  if (j_m > nj) j_m = nj
  model_loc%lon = T_grid%x(i_m)
  model_loc%lat = T_grid%y(j_m)
  if(model_loc%lon < 0.0) model_loc%lon = model_loc%lon + 360.0
  if(model_loc%lon > 360.0) model_loc%lon = 360.0 - model_loc%lon

  if(low_lat < model_loc%lat .and. model_loc%lat < hi_lat) then
    low_lon = olon - close_lon_window/cos(olat0*radian)
    if(low_lon < 0.0) low_lon = low_lon + 360.0
    hi_lon  = olon + close_lon_window/cos(olat0*radian)
    if(hi_lon >= 360.5) hi_lon = 360.0 - hi_lon
	        
    if (low_lon < hi_lon) then 
      if ((model_loc%lon > low_lon) .and. (model_loc%lon < hi_lon)) then
        num = num + 1
        list(num) = (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
      end if 
    elseif (hi_lon < low_lon) then
      if ((model_loc%lon > low_lon .and. model_loc%lon < 360.0) .or. &
          (model_loc%lon >= 0.0 .and. model_loc%lon < hi_lon )) then
        num = num + 1
        list(num) = (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
      end if 
    else
      num = num + 1
      list(num) = (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
    end if  
					               
  end if

end do
end do

end subroutine get_close_grids_sice_up

end module obs_eakf_sice_up_mod
