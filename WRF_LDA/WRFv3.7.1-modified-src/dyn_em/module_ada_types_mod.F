module module_ada_types_mod
  ! This module contains type declarations and default values
  ! for oda modules.  
  !
  ! Contact: S. Zhang Shaoqing.Zhang@noaa.gov

!!  use time_manager_mod, only : time_type

  implicit none

  private

  ! codes for modeling error distributions
  type, public :: ada_grid_type
     real, pointer, dimension(:) :: x => NULL(), y => NULL()
     real, pointer, dimension(:) :: dx => NULL(), dy => NULL()
     real, pointer, dimension(:,:,:) :: z => NULL()
     real, pointer, dimension(:) :: dz => NULL()
     real, pointer, dimension(:) :: pfull => NULL()
     real, pointer, dimension(:,:,:) :: mask => NULL()
     real,pointer, dimension(:) :: xu => NULL(),yu=> NULL()
     real,pointer, dimension(:) :: xv => NULL(),yv=> NULL()
     integer :: ni, nj, nk
  end type ada_grid_type


  type, public :: ada_field_type
     type(ada_grid_type) :: grid
     real, pointer, dimension(:,:,:) :: data => NULL()
  end type ada_field_type

  type, public :: atmos_obs_type
     integer :: nlev
     logical :: uv_flag
     real :: lon, lon_uv
     real :: lat, lat_uv
     real :: ps
     real :: i_index, j_index
     real :: i_index_uv, j_index_uv
     real, pointer, dimension(:) :: k_index => NULL()
     real, pointer, dimension(:) :: z => NULL()
     real, pointer, dimension(:) :: temp => NULL()
     real, pointer, dimension(:) :: mois => NULL()
     real, pointer, dimension(:) :: u => NULL()
     real, pointer, dimension(:) :: v => NULL()
!!     type(time_type) :: time
  end type atmos_obs_type

!x1y add for ps data
 type, public :: ps_obs_type
!     integer :: nvar
!     character(len=19) :: reana_id
!     character(len=1) :: obs_tp
!     integer :: ncep_tp
!     real    :: database_id
!     real :: time_offset
     real :: elevation
!     integer :: qc_id
     real :: lat, lon
!     real :: depth
!     logical :: accepted
!     integer :: levels
     real :: ps_data
!     real, pointer :: ps_data
!     real :: model_ps, model_slp
!     real :: bias
!     integer  :: flag_ps ! flags for ps
!     integer  :: flag_assi ! flags for assimilation
!     integer  :: flag_reana, flag_bc, flag_buddy
     real    :: ps_error ! measurement error
!     real    :: guess_ensm, guess_std, omf_ps, ana_ensm, ana_std, oma_ps
!     integer :: melevation
!     character(len=13) :: obs_id
!     character(len=30) :: stn_name
!!     type(time_type) :: time
!     integer         :: yyyy
!     integer         :: mmdd
!     type(time_type), pointer :: Model_time ! each profile can be associated with a first-guess field with an associated time and grid
!     type(grid_type), pointer :: Model_grid
     real :: i_index, j_index ! model longitude and latitude indices respectively
!     type(forward_model_type) :: Forward_model  ! linear operation from model to observation
!!     type(time_type) :: tdiff      ! positive difference between model time and observation time
  end type ps_obs_type
!x1y end for ps data

  type, public :: interp_coef_type
     integer :: i_ls, i_rn
     real :: frac_ls, frac_rn
  end type interp_coef_type

  type, public :: assim_data_type
     real, pointer, dimension(:,:) :: x => NULL(), y => NULL()
     real, pointer, dimension(:,:,:) :: z => NULL()
     real, pointer, dimension(:,:,:) :: u_2 => NULL()
     real, pointer, dimension(:,:,:) :: t_2 => NULL()
     real, pointer, dimension(:,:,:) :: u => NULL()
     real, pointer, dimension(:,:,:) :: v => NULL()
     real, pointer, dimension(:,:,:) :: t => NULL()
     real, pointer, dimension(:,:,:) :: q => NULL()
     real, pointer, dimension(:,:,:) :: dp => NULL()
     real, pointer, dimension(:,:) :: ps => NULL()
     real, pointer, dimension(:,:) :: phis => NULL() !lulv 20200731
!     real, pointer, dimension(:,:,:) :: u_adj => NULL()
!     real, pointer, dimension(:,:,:) :: v_adj => NULL()
!     real, pointer, dimension(:,:,:) :: t_adj => NULL()
!     real, pointer, dimension(:,:,:) :: q_adj => NULL()
!     real, pointer, dimension(:,:,:) :: dp_adj => NULL()
!     real, pointer, dimension(:,:) :: ps_adj => NULL()
  end type assim_data_type

contains

  subroutine ada_types_init()
    ! Place holder routine
  end subroutine ada_types_init
    
end module module_ada_types_mod
