module eakf_ida_up_mod

! this module is produced by snz for parallelism of ensemble-based filtering
! data assimilation algorithm, starting from August 9, 2002. A linear 
! regression eakf is used as a primary algorithm to design the parallelism.
! Communications are broadcasting the information at the observation location.

!----------------------------------------------------------------------

!use eakf_tab_mod, only : rrt, rrr, rrh
use assim_tools_mod, only : assim_tools_init, &
                            obs_increment, update_from_obs_inc
use model_sice_up_mod, only : init_model_sice_up, ens_ics_sice_up, red_ens_sice_up, &
                           get_model_size_sice_up
use obs_eakf_sice_up_mod, only : sice_obs_init_up, sice_obs_end_up, & 
                              get_close_grids_sice_up
use fms_mod, only: open_namelist_file, file_exist, check_nml_error, &
                   write_version_number, close_file
use fms_io_mod, only: open_file
use loc_and_dist_mod, only : loc_type, get_dist
use cov_cutoff_mod, only : comp_cov_factor

use ida_types_mod, only : ida_grid_type, sice_obs_type
use ada_types_mod, only : ada_field_type, ada_grid_type
use time_manager_mod, only : time_type, get_time

use mpp_mod, only : mpp_error, FATAL, mpp_sync_self, mpp_set_stack_size, &
                    mpp_pe, mpp_npes, mpp_broadcast, mpp_clock_id, &
                    mpp_clock_begin, mpp_clock_end, mpp_root_pe

use rand_no_mod, only : gau0 ! for adding errors on idealized obs
!-----------------------------------------------------------------------

implicit none

integer, parameter, private :: CN_ID = 1, UI_ID = 2, VI_ID = 3, HI_ID = 4
integer, parameter, private :: max_sice_obs = 50000 ! must be same as obs

logical, private :: first_run_call = .true.

real, allocatable :: ens(:, :), ens_mean(:)
real, allocatable :: enso_ti(:), obs_inc_eakf_ti(:), &
        obs_inc_oi_ti(:), ens_inc(:)
real, allocatable :: enso_ui(:), obs_inc_eakf_ui(:), &
        obs_inc_oi_ui(:)
real, allocatable :: enso_vi(:), obs_inc_eakf_vi(:), &
        obs_inc_oi_vi(:)

public ensemble_filter_sice_up

contains

subroutine ensemble_filter_sice_up(u_ens_tau, v_ens_tau, t_ens_tau, &
             q_ens_tau, sice_obs, num_obs,isd, ied, jsd, jed, halox, haloy, nk, &
             T_grid, iass, m_time)

!---- namelist with default values

real :: std_cut_b = 0.031623
real :: cutoff_vd = 5.0 ! level numbers
real :: std_cut_t = 1.00
real :: tq_dist = 200.0e3
real :: uv_dist = 100.0e3
real :: sigma_o_ti = 1.0
real :: sigma_o_ui = 1.0
real :: sigma_o_vi = 1.0

integer :: assim_freq = 24
integer :: ass_method = 1 ! 0 for snz-oi, 1 for eakf inv, 2 for eakf multv

logical :: ass_ti = .false.
logical :: ti_impact_t = .false., ti_impact_uv = .false., ti_impact_q = .false.
logical :: ass_uv = .false.
logical :: ui_impact_u = .false., ui_impact_v = .false., ui_impact_t = .false., &
           ui_impact_q = .false.
logical :: vi_impact_u = .false., vi_impact_v = .false., vi_impact_t = .false., &
           vi_impact_q = .false.

real :: sice_ass_up_lat = 40.0
real :: std_oi_fraction = 0.1
real :: time_window = 1.0 ! (days)

namelist /eakf_ida_up_nml/sigma_o_ti, sigma_o_ui, sigma_o_vi, &
          cutoff_vd, assim_freq, ass_method, std_cut_b, std_cut_t, &
          tq_dist, uv_dist, sice_ass_up_lat, ass_ti, ti_impact_t, ti_impact_uv, &
	  ti_impact_q, ass_uv, ui_impact_u, ui_impact_v, ui_impact_t, ui_impact_q, &
	  vi_impact_u, vi_impact_v, vi_impact_t, vi_impact_q, &
	  std_oi_fraction, time_window

!--- module name and version number ----
character(len = 11), parameter :: module_name = 'eakf'
character(len = 5), parameter :: vers_num = 'x100.0'

!=======================================================================
! input variables

type(ada_field_type), intent(inout) :: u_ens_tau(:), v_ens_tau(:), &
                                       t_ens_tau(:), q_ens_tau(:)
type(sice_obs_type), intent(in) :: sice_obs(:)
type(ada_grid_type) :: T_grid
type(time_type), intent(in) :: m_time
integer, intent(inout) :: num_obs
integer, intent(in) :: isd, ied, jsd, jed, halox, haloy, nk
integer, intent(inout) :: iass
integer :: ass_variable = 1 ! 1 for temperature, 2 for salinity
integer :: ni, nj

!=======================================================================

integer :: num_prfs_loc_halo
integer :: list_loc_halo_obs(max_sice_obs)
integer :: list_close_grids(100*100), index_obs(max_sice_obs)

!=======================================================================

integer :: id_eakf_total

integer :: isd_ens, ied_ens, jsd_ens, jed_ens, num_obsb, ngrids, num_obst
real, parameter :: radius = 6370.e+3
real, parameter :: pai = 3.1415926
real, parameter :: radian = pai/180.0

real :: cor_oi, e_flder_aed, coslat, ramp_coef
real :: cov_factor, cov_factor_v, cov_factor_t, cov_factor_h

integer :: num_close, assim_flag
type(loc_type) :: model_loc, obs_loc, model_loc_u

integer :: ii_ens, jj_ens, kk_ens, nv
integer :: i0, i, j, k, k0, kk, num, blk, i_idx, k_as
integer :: i_t, i_u, i_v, i_q
integer :: idx_obs, idx_buf, idx_k, lji0, npe, npes, kk0, kk1, kk2, iii, jjj
integer :: ind_unit(20), lji, model_size_sice, ens_size
integer :: unit, ierr, io, pe_curr, j_ens, i_h
integer :: m_days, m_hours, m_seconds, o_days, o_hours, o_seconds, nk_adj

!---------------------------------------------------------------------------

real :: dist, dist0, obs_value, dist_uv
real :: obs_var_ti, obs_var_ui, obs_var_vi
real :: std_oi, std_oi_o, std_oi_g, std_c

!---------------------------------------------------------------------------

character*40 :: file_name
character*40 :: diag_file

!---------------------------------------------------------------------------

ni = T_grid%ni
nj = T_grid%nj

std_c = 0.5
!---------------------------------------------------------------------------

id_eakf_total = mpp_clock_id('(total eakf computation) ')

!---------------------------------------------------------------------------

call mpp_clock_begin(id_eakf_total)

isd_ens = isd
ied_ens = ied
jsd_ens = jsd
jed_ens = jed
ens_size = size(t_ens_tau)

if (iass <= 1) return
iass = iass + 1
 
npes = mpp_npes()
pe_curr = mpp_pe()

! Read namelist for run time control
if(file_exist('input.nml')) then
   unit = open_namelist_file()
   ierr = 1
   do while(ierr /= 0)
      read(unit, nml = eakf_ida_up_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'eakf_ida_up_nml')
   enddo
 11 continue
   call close_file(unit)
endif

! Write the namelist to a log file
unit = open_file(file = 'logfile.out', action = 'write')
call write_version_number(vers_num, module_name, unit)
write(unit, nml = eakf_ida_up_nml)
call close_file(unit)

if ( num_obs > max_sice_obs ) num_obs = max_sice_obs

if (mpp_pe() == mpp_root_pe())print*,'num_obs for atmosphere is ', num_obs

blk = (jed_ens-jsd_ens+2*haloy+1)*(ied_ens-isd_ens+2*halox+1)

!--------------------------------------------------------------------------
call init_model_sice_up(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, ass_method)
model_size_sice = get_model_size_sice_up()

! Begin by initializing the observations

call sice_obs_init_up(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, &
     sice_obs, num_obs, T_grid, list_loc_halo_obs, num_prfs_loc_halo)

if (first_run_call) then
allocate(ens(model_size_sice, ens_size), ens_mean(model_size_sice))
allocate(enso_ti(ens_size), obs_inc_eakf_ti(ens_size), &
         obs_inc_oi_ti(ens_size), ens_inc(ens_size))
allocate(enso_ui(ens_size), obs_inc_eakf_ui(ens_size), obs_inc_oi_ui(ens_size))
allocate(enso_vi(ens_size), obs_inc_eakf_vi(ens_size), obs_inc_oi_vi(ens_size))

! for special handling on corrections in vertical direction

end if

! Initialize assim tools module
call assim_tools_init()

! print namelist

if (pe_curr == mpp_root_pe() .and. first_run_call) then

write(*, *) 'no of sice obs up is ', num_obs
write(*, *) 'sice_up model size is ', model_size_sice, 'ensemble size is ', ens_size
write(*, *) 'assim_freq, ass_method is ', assim_freq, ass_method
write(*, *) 'std_cut_b, std_cut_t is', std_cut_b, std_cut_t
write(*, *) '{tq,uv}_dist are',tq_dist, uv_dist
write(*, *) 'ass_ti(uv) is', ass_ti, ass_uv
write(*, *) 'sice_ass_up_lat is', sice_ass_up_lat
write(*, *) 'ti_impact_t(uv,q) is', ti_impact_t,ti_impact_uv,ti_impact_q
write(*, *) 'ui_impact_t(uv,q) is', ui_impact_t,ui_impact_u,ui_impact_v,ui_impact_q
write(*, *) 'vi_impact_t(uv,q) is', vi_impact_t,vi_impact_u,vi_impact_v,vi_impact_q
write(*, *) 'std_oi_fraction is', std_oi_fraction
write(*, *) 'time_window is', time_window

end if

!write(*,*)'pe ',mpp_pe(),'finish eakf initialization'
! Form the ensemble state vector: ens(:, :)
call ens_ics_sice_up(t_ens_tau, u_ens_tau, v_ens_tau, q_ens_tau, &  
     isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, ens, ass_method)

!write(*,*)'pe ',mpp_pe(),'finish ens_ics_atm'
! ###########################################################
! The assimilation main part starts here

! Compute the ensemble mean of the initial ensemble before assimilation
  ens_mean = sum(ens, dim=2) / ens_size

call mpp_sync_self()

! Loop through each observation location available at this time
  index_obs = 0
  idx_obs = 0
  do lji = 1, num_prfs_loc_halo
    lji0 = list_loc_halo_obs(lji)
    idx_obs = idx_obs + 1
    index_obs(lji0) = idx_obs
  end do

!write(*,*)'pe ',mpp_pe(),'finish index_obs'

obs_var_ti = sigma_o_ti**2
obs_var_ui = sigma_o_ui**2
obs_var_vi = sigma_o_vi**2

! Section to do adjustment point by point
! coding for cov_factor

  call get_time(m_time, m_seconds, m_days)
  m_hours = m_seconds/3600 + m_days * 24

if(mpp_pe() == mpp_root_pe())write(*,*)'m_hours=',m_hours

  ngrids = 0

!===== Eakf assim start =====================================
! for special handling on corrections in vertical direction

call mpp_sync_self()

!write(*,*)'pe ',mpp_pe(),'start lji loop'

do lji = 1, num_prfs_loc_halo ! (1) loop for num_obs in the analysis domain

lji0 = list_loc_halo_obs(lji)

idx_obs = index_obs(lji0)
if(idx_obs == 0)print*,'lji0 4ps =',lji0

if (sice_obs(lji0)%flag) then ! general flag

call get_time(sice_obs(lji0)%time, o_seconds, o_days)
o_hours = o_seconds/3600 + o_days * 24

o_hours = abs(o_hours - m_hours)

cov_factor_t = comp_cov_factor(real(o_hours), 24.0*time_window)

if(cov_factor_t > 0.0) then ! control 1d time window (1+ and 1-) 

obs_loc%lon = sice_obs(lji0)%lon
obs_loc%lat = sice_obs(lji0)%lat

if (abs(obs_loc%lat) <= 80.0) then
  coslat = cos(obs_loc%lat*radian)
else
  coslat = cos(80.0*radian)
end if

std_oi_o = 0.0
obs_inc_oi_ti = 0.0
obs_value = sice_obs(lji0)%ti
enso_ti(:) = sice_obs(lji0)%enso_ti(:)
if (sice_obs(lji0)%flag_ti) then
  call obs_increment(enso_ti, ens_size, obs_value, obs_var_ti, &
                     obs_inc_eakf_ti, obs_inc_oi_ti, std_oi_o)
end if

obs_inc_oi_ui = 0.0
obs_inc_oi_vi = 0.0
if (sice_obs(lji0)%flag_uv) then
  obs_value = sice_obs(lji0)%ui
  enso_ui(:) = sice_obs(lji0)%enso_ui(:)
  call obs_increment(enso_ui, ens_size, obs_value, obs_var_ui, &
                     obs_inc_eakf_ui, obs_inc_oi_ui, std_oi_o)
  obs_value = sice_obs(lji0)%vi
  enso_vi(:) = sice_obs(lji0)%enso_vi(:)
  call obs_increment(enso_vi, ens_size, obs_value, obs_var_vi, &
                     obs_inc_eakf_vi, obs_inc_oi_vi, std_oi_o)
end if

call get_close_grids_sice_up(obs_loc, isd_ens, ied_ens, jsd_ens, jed_ens, &
                     halox, haloy, T_grid, list_close_grids, num_close)
ngrids = ngrids + num_close

do 10000 k = 1, num_close ! (2) loop for close gridpoints for obs lji0

  j = list_close_grids(k)

  jj_ens = (j-1)/(ied_ens-isd_ens+2*halox+1)+1 + (jsd_ens-1-haloy)

  ii_ens = mod(j, ied_ens-isd_ens+2*halox+1)
  if (ii_ens == 0 ) ii_ens = ied_ens-isd_ens+2*halox+1
  ii_ens = ii_ens + (isd_ens-1-halox)

  i_h = (jj_ens-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+ii_ens-isd_ens+halox+1

  if (ii_ens <= 0) ii_ens = ii_ens + ni
  if (ii_ens > ni) ii_ens = ii_ens - ni
  if (jj_ens <= 0) jj_ens = 1
  if (jj_ens > nj) jj_ens = ni

  model_loc%lon = T_grid%x(ii_ens)
  model_loc%lat = T_grid%y(jj_ens)
  model_loc%lon = model_loc%lon+360.0
  if(model_loc%lon > 360.0) model_loc%lon = model_loc%lon-360.0
  if(obs_loc%lon > 360.0) obs_loc%lon = obs_loc%lon-360.0

  if ( abs(model_loc%lat) > sice_ass_up_lat ) then ! (4) within sice assim lats

    if (abs(model_loc%lat) > sice_ass_up_lat+10.0) then
      ramp_coef = 1.0
    else
      ramp_coef = (abs(model_loc%lat)-sice_ass_up_lat)/10.0
    end if

    dist = get_dist(model_loc, obs_loc)
    dist = radius * sqrt(dist) 

! observed ti corrects t, u, v, q
	  
    if (ass_ti) then

! observed ti corrects t
    if (ti_impact_t) then

    ass_variable = 2
    dist0 = tq_dist*coslat
    cov_factor_h = comp_cov_factor(dist, dist0)
    do kk_ens = nk, 1, -1
    i_t = (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if (cov_factor /= 0.0 .and. sum(ens(i_t, :)) /= 0.0) then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ti, obs_inc_eakf_ti, &
             obs_inc_oi_ti, ens(i_t, :), ens_size, &
             ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
             ass_method, ass_variable)
      ens(i_t, :)   = ens(i_t, :) + ens_inc(:)
    end if
    end do

    end if

! observed ti corrects uv
    if (ti_impact_uv) then

    ass_variable = 2
    dist0 = uv_dist*coslat
    cov_factor_h = comp_cov_factor(dist, dist0)
    do kk_ens = nk, 1, -1
    i_u = blk * nk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_u, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ti, obs_inc_eakf_ti, &
           obs_inc_oi_ti, ens(i_u, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_u, :)   = ens(i_u, :) + ens_inc(:) 
    end if
    i_v = 2 * blk * nk + (kk_ens-1)*blk + i_h
    if(sum(ens(i_v, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ti, obs_inc_eakf_ti, &
           obs_inc_oi_ti, ens(i_v, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_v, :)   = ens(i_v, :) + ens_inc(:) 
    end if

    end do

    end if

! observed ti corrects q
    if (ti_impact_q) then

    ass_variable = 2
    dist0 = tq_dist*coslat
    cov_factor_h = comp_cov_factor(dist, dist0)
    do kk_ens = nk, 1, -1
    i_q = 3* blk * nk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_q, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ti, obs_inc_eakf_ti, &
           obs_inc_oi_ti, ens(i_q, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_q, :)   = ens(i_q, :) + ens_inc(:) 
    end if
    end do

    end if
    
! end of assimilating ti

    end if

! observed uv corrects uv
    if (ass_uv) then

! observed ui corrects t
    if (ui_impact_t) then

    ass_variable = 2
    dist0 = uv_dist*coslat
    cov_factor_h = comp_cov_factor(dist, dist0)
    do kk_ens = nk, 1, -1
    i_t  = (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_t, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ui, obs_inc_eakf_ui, &
           obs_inc_oi_ui, ens(i_t, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_t, :)   = ens(i_t, :) + ens_inc(:) 
    end if
    end do

    end if
    
! observed ui corrects u
    if (ui_impact_u) then

    ass_variable = 2
    do kk_ens = nk, 1, -1
    i_u  = nk*blk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_u, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ui, obs_inc_eakf_ui, &
           obs_inc_oi_ui, ens(i_u, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_u, :)   = ens(i_u, :) + ens_inc(:) 
    end if
    end do

    end if 

! observed ui corrects v
    if (ui_impact_v) then

    ass_variable = 2
    do kk_ens = nk, 1, -1
    i_v  = 2* nk*blk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_v, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ui, obs_inc_eakf_ui, &
           obs_inc_oi_ui, ens(i_v, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_v, :)   = ens(i_v, :) + ens_inc(:) 
    end if
    end do

    end if

! observed ui corrects q
    if (ui_impact_q) then

    ass_variable = 2
    do kk_ens = nk, 1, -1
    i_q  = 3* nk*blk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_q, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_ui, obs_inc_eakf_ui, &
           obs_inc_oi_ui, ens(i_q, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_q, :)   = ens(i_q, :) + ens_inc(:) 
    end if
    end do

    end if

! observed vi corrects t
    if (vi_impact_t) then

    ass_variable = 2
    dist0 = uv_dist*coslat
    cov_factor_h = comp_cov_factor(dist, dist0)
    do kk_ens = nk, 1, -1
    i_t  = (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_t, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_vi, obs_inc_eakf_vi, &
           obs_inc_oi_vi, ens(i_t, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_t, :)   = ens(i_t, :) + ens_inc(:) 
    end if
    end do

    end if
    
! observed vi corrects u
    if (vi_impact_u) then

    ass_variable = 2
    do kk_ens = nk, 1, -1
    i_u  = nk*blk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_u, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_vi, obs_inc_eakf_vi, &
           obs_inc_oi_vi, ens(i_u, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_u, :)   = ens(i_u, :) + ens_inc(:) 
    end if
    end do

    end if 

! observed ui corrects v
    if (vi_impact_v) then

    ass_variable = 2
    do kk_ens = nk, 1, -1
    i_v  = 2* nk*blk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_v, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_vi, obs_inc_eakf_vi, &
           obs_inc_oi_vi, ens(i_v, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_v, :)   = ens(i_v, :) + ens_inc(:) 
    end if
    end do

    end if

! observed vi corrects q
    if (vi_impact_q) then

    ass_variable = 2
    do kk_ens = nk, 1, -1
    i_q  = 3* nk*blk + (kk_ens-1)*blk + i_h
    cov_factor_v = comp_cov_factor(real(nk-kk_ens), cutoff_vd)
    cov_factor = cov_factor_h * cov_factor_t * cov_factor_v
    if(sum(ens(i_q, :)) /= 0.0 .and. cov_factor /= 0.0)then
      std_oi_g = 0.0
      std_oi_o = 0.0
      cor_oi = 0.0
      ens_inc(:) = 0.0
      call update_from_obs_inc(enso_vi, obs_inc_eakf_vi, &
           obs_inc_oi_vi, ens(i_q, :), ens_size, &
           ens_inc, cov_factor, cor_oi, std_oi_o, std_oi_g, &
           ass_method, ass_variable)
      ens(i_q, :)   = ens(i_q, :) + ens_inc(:) 
    end if
    end do

    end if

    end if ! end of ass_uv

  end if ! South or North sice assimilation domain

10000 end do ! finish all close gridpoints 

end if ! time_window

end if ! general flag

end do ! (1) finish all obs location

!print*,'ngrids=',ngrids,'in pe ',mpp_pe()

!===== Eakf assim finish =====================================

call mpp_sync_self()

! Redistribute the sub ensemble state vector ens(:, :) back to the model grids
! in the local-domain.
call red_ens_sice_up(t_ens_tau, u_ens_tau, v_ens_tau, q_ens_tau, &
             isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, ens, ass_method)

!write(*,*)'pe ',mpp_pe(),'finish red_ens'

!--------------------------------------------------------------------------

call sice_obs_end_up()

first_run_call = .false.

call mpp_clock_end(id_eakf_total)

return

end subroutine ensemble_filter_sice_up

end module eakf_ida_up_mod
