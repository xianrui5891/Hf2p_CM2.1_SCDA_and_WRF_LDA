module ida_types_mod

! This module contains type declarations and default values
! for oda modules.  
!
! Contact: S. Zhang Shaoqing.Zhang@noaa.gov

  use time_manager_mod, only : time_type

  implicit none

  private

  ! codes for modeling error distributions
  
  type, public :: ida_grid_type
     real, pointer, dimension(:) :: x_atm, y_atm
     real, pointer, dimension(:, :) :: x_ocn, y_ocn
     real, pointer, dimension(:,:,:) :: mask
     integer :: ni, nj, nk
  end type ida_grid_type

  type, public :: ida_field_type
     type(ida_grid_type) :: grid
     real, pointer, dimension(:,:,:) :: data
  end type ida_field_type

  type, public :: sice_obs_type
     logical :: flag = .false.
     logical :: flag_cn = .false.
     logical :: flag_uv = .false.
     logical :: flag_ti = .false.
     integer :: ens_size
     real :: lon
     real :: lat
     real :: i_index_up, j_index_up, i_index_dn, j_index_dn
     real :: i_index_up_uv, j_index_up_uv, i_index_dn_uv, j_index_dn_uv
     real :: cn, ui, vi, ti
     real, pointer, dimension(:) :: enso_ti, enso_ui, enso_vi
     type(time_type) :: time
  end type sice_obs_type

  contains

    subroutine ida_types_init()

      use fms_mod, only : open_namelist_file, check_nml_error, close_file
      

    end subroutine ida_types_init
    
end module ida_types_mod
