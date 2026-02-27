module obs_tools_mod
  ! Provides a set of tools that are used by observations modules to generate
  ! data structures and information needed to do ensemble assimilation.
  ! These should be common to a variety of models and observational operator
  ! types.

  use mpp_mod, only : mpp_pe
  use fms_mod, only : stdout

  private
  public conv_state_to_obs, obs_def_type, def_single_obs, def_single_obs_end

  ! Define a type for linear model operator
  type obs_def_type
     integer :: num
     integer, pointer :: state_var_index(:)
     real, pointer :: coef(:)
  end type obs_def_type

contains

  !=======================================================================
  ! Given a state vector x, returns corresponding observations given the
  ! observations data definition structure obs_def.
  function conv_state_to_obs(x, obs_def, num_obs, idx)
    real, intent(in) :: x(:)
    type (obs_def_type), dimension(num_obs), intent(in) :: obs_def
    integer, intent(in) :: num_obs
    integer, intent(in) :: idx

    real, dimension(num_obs) :: conv_state_to_obs

    real :: value_h, value_l, prdt
    integer :: i, j, idx_x, size_x, idx_x1, idx0, stdout_unit

    stdout_unit = stdout()
    
    idx0 = idx
    conv_state_to_obs = 0.0
    do i=1, num_obs
       prdt = 1.0
       do j=1, obs_def(i)%num
          idx_x = obs_def(i)%state_var_index(j)
          prdt = prdt * x(idx_x)
          size_x = size(x)
          if ( idx_x > size_x .or. idx_x <= 0 ) then
             write (UNIT=stdout_unit, FMT='("PE ",I5,": idx_x = ",I8,", size_x = ",I8,", idx0 = ",I8)') mpp_pe(), idx_x, size_x, idx0
          end if
       end do
       if ( prdt == 0.0 ) then
          conv_state_to_obs(i) = 0.0
       else
          value_h = 0.0
          value_l = 0.0
          do j=1, obs_def(i)%num/2
             idx_x = obs_def(i)%state_var_index(j)
             idx_x1 = obs_def(i)%state_var_index(j+obs_def(i)%num/2)
             value_h = value_h + obs_def(i)%coef(j)*x(idx_x)
             value_l = value_l + obs_def(i)%coef(j)*x(idx_x1)
          end do
          conv_state_to_obs(i) = obs_def(i)%coef(5)*value_h + obs_def(i)%coef(6)*value_l
       end if
    end do
  end function conv_state_to_obs


!!$  real function conv_state_to_obs_snz(x, coef)
!!$    real, dimension(:), intent(in) :: x, coef
!!$
!!$    real :: value_h, value_l, prdt
!!$    integer :: i, j
!!$
!!$    prdt = 1.0
!!$    do i=1, 8
!!$       prdt = prdt * x(i)
!!$    end do
!!$
!!$    if ( prdt == 0.0 ) then
!!$       conv_state_to_obs_snz = 0.0
!!$    else
!!$       value_h = 0.0
!!$       value_l = 0.0
!!$       do j=1, 4
!!$          value_h = value_h + x(j) * coef(j)
!!$          value_l = value_l + x(j+4)*coef(j)
!!$       end do
!!$       conv_state_to_obs_snz = value_h * coef(5) + value_l * coef(6)
!!$    end if
!!$  end function conv_state_to_obs_snz


  !=======================================================================
  ! Puts definition of a single observation into an obs_def data structure.
  subroutine def_single_obs(num_state, state_ind, coef, obs_def)
    integer, intent(in) :: num_state
    integer, dimension(num_state), intent(in) :: state_ind
    real, dimension(num_state-2), intent(in) :: coef
    type(obs_def_type), intent(inout) :: obs_def

    integer :: i

    ! Set aside storage for defining this ob
    obs_def%num = num_state
    allocate(obs_def%state_var_index(num_state), obs_def%coef(num_state-2))

    ! Load the state variable index and coefficient for each state variable
    do i = 1, num_state
       obs_def%state_var_index(i) = state_ind(i)
    end do

    do i=1, num_state-2
       obs_def%coef(i) = coef(i)
    end do
  end subroutine def_single_obs

  !=======================================================================
  subroutine def_single_obs_end(obs_def)
    type(obs_def_type), intent(inout) :: obs_def

    !::sdu:: Deallocate, and nullify pointers if allocated
    if ( associated(obs_def%state_var_index) ) then
       deallocate(obs_def%state_var_index)
       nullify(obs_def%state_var_index)
    end if

    if ( associated(obs_def%coef) ) then
       deallocate(obs_def%coef)
       nullify(obs_def%coef)
    end if
  end subroutine def_single_obs_end
end module obs_tools_mod
