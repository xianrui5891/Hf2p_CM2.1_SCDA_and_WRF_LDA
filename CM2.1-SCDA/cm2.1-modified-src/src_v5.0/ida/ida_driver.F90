module ida_driver_mod

!
! This is the top-level module for ocean data assimilation.
!

  use ida_types_mod, only: ida_grid_type, ida_field_type, sice_obs_type
  use ida_core_mod, only : ida_core_init
  use mpp_domains_mod, only : domain2d, mpp_update_domains, &
              mpp_define_domains, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, &
              FOLD_NORTH_EDGE, mpp_get_compute_domain, mpp_get_data_domain, &
              null_domain2d, mpp_redistribute, mpp_broadcast_domain, &
              mpp_update_domains, mpp_global_field
  use mpp_mod, only : mpp_npes, stdout, stdlog, mpp_error, FATAL, &
              mpp_pe, mpp_declare_pelist, mpp_sync_self, mpp_sum, &
              mpp_set_current_pelist, mpp_root_pe, mpp_set_stack_size, &
	      mpp_broadcast
  use time_manager_mod, only : time_type, decrement_time, increment_time, &
              get_time
  use constants_mod, only : radius, epsln

  use axis_utils_mod, only : frac_index

  use ice_type_mod, only: ice_data_type, hlim
  use ice_grid_mod, only: Domain, isc, iec, jsc, jec, isd, ied, jsd, jed, im,&
                          jm, km, xb1d, yb1d, geo_lon, geo_lat

  use ensemble_manager_mod, only : get_ensemble_id, get_ensemble_size, &
                                   get_ensemble_pelist

  IMPLICIT NONE
  
  private

  integer, parameter :: NO_ASSIM = 0, OI=1, EAKF=2
  integer, parameter :: max_sice_obs = 20000

  integer :: filter_halo_x = 1, filter_halo_y = 1

  integer :: asm_method = NO_ASSIM
  character(len=8) :: assim_method = 'NO_ASSIM'

  type (ida_grid_type), private, save :: sice_grid

  integer :: ensemble_size = 1
  integer, allocatable, dimension(:,:) :: ensemble_pelist
  integer :: pe, npes, npes_pm, ensemble_id, ens_siz(4)
  logical :: sequential_filter = .false.

  integer :: assim_layout(2) ! snz

  real :: covar_cutoff = 0.0
  integer :: assim_frequency = 24

  type(sice_obs_type) :: sice_obs(max_sice_obs)

  public :: init_ida

contains  
  
  subroutine init_ida(sice_time)
    ! initialize First_guess and Analysis grid and domain information
    ! grids are identical to ocean_model grid.  Halo size may differ.

    use fms_mod, only : open_namelist_file,close_file,check_nml_error
    use mpp_mod, only : mpp_set_current_pelist
    use diag_manager_mod, only : register_diag_field, diag_axis_init
    use mpp_domains_mod, only : mpp_global_field

    type(time_type), intent(in), target :: sice_Time

    type(time_type), pointer :: Time

! snz tries to add the time period to limit prfs num.
    type(time_type)  :: time_s, time_e

    integer :: ioun, io_status, ierr
    integer :: m

!-----------------------------------------------------------------------

    namelist /ida_nml / assim_method, assim_layout, &
			covar_cutoff, assim_frequency, &
			sequential_filter, filter_halo_x, filter_halo_y

    ioun = open_namelist_file()
    read(ioun,nml=ida_nml,iostat = io_status)
    ierr = check_nml_error(io_status,'ida_nml')
    call close_file(ioun)

    if (assim_frequency < 1) then
        call mpp_error(FATAL,'invalid assim frequency')
    else
      ensemble_id = get_ensemble_id()
      ens_siz = get_ensemble_size()
      ensemble_size = ens_siz(1)
      npes_pm = ens_siz(2)
      allocate(ensemble_pelist(1:ensemble_size,1:npes_pm))

      call get_ensemble_pelist(ensemble_pelist)

      write(stdout(),*) 'Assimilating data every ', assim_frequency, ' timesteps'
    endif

    Time => sice_time

    call mpp_set_stack_size(500000)

! snz tries to add the time period to limit prfs num.
    time_s = decrement_time(Time,0,6)
    time_e = increment_time(Time,0,36)

! get global grid information from ocean_model

    sice_grid%ni = im
    sice_grid%nj = jm
    sice_grid%nk = km
    allocate(sice_grid%x_ocn(im, jm), sice_grid%y_ocn(im, jm))

    call mpp_global_field(Domain, geo_lon(:, :), sice_grid%x_ocn)
    call mpp_global_field(Domain, geo_lat(:, :), sice_grid%y_ocn)

    call ida_core_init(sice_grid, time_s, time_e, localize=.false.)       

    return

  end subroutine init_ida

end module ida_driver_mod
