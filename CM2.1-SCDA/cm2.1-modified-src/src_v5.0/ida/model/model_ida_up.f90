
module model_sice_up_mod

!-----------------------------------------------------------------------
!
!   model_mod that couples component models for the atmosphere,
!   ocean (amip), land, and sea-ice using the exchange module
!
!-----------------------------------------------------------------------

use ada_types_mod, only : ada_field_type

implicit none
public init_model_sice_up, ens_ics_sice_up, red_ens_sice_up, &
       get_model_size_sice_up

!-----------------------------------------------------------------------

integer, private :: model_size_sice

contains

!#######################################################################

subroutine init_model_sice_up(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, &
                      ass_method)

implicit none

integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk
integer, intent(in) :: ass_method

model_size_sice = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*nk*4

end subroutine init_model_sice_up

!-------------------------------------------------------------------------

function get_model_size_sice_up()

implicit none

integer :: get_model_size_sice_up

get_model_size_sice_up = model_size_sice

end function get_model_size_sice_up

!-----------------------------------------------------------------------

subroutine ens_ics_sice_up(t_ens_tau, u_ens_tau, v_ens_tau, q_ens_tau, &
                           isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, &
			   nk, as, ass_method)

implicit none

type(ada_field_type), intent(in) :: t_ens_tau(:), u_ens_tau(:), &
                                    v_ens_tau(:), q_ens_tau(:)
integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, &
                       ass_method
real, intent(inout) :: as(:, :)

integer :: m, i, j, k, it, idx, blk_h

!  Get initial state for ensemble members

   blk_h = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)

   as(:, :) = 0.0

   do m = 1, size(as, 2)

   do k = 1, nk
   do j = jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
   do i = isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
   idx =  (k-1)*blk_h + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   as(idx, m) = t_ens_tau(m)%data(i,j,k)
   idx = blk_h*nk + idx
   as(idx, m) = u_ens_tau(m)%data(i,j,k)
   idx = blk_h*nk + idx
   as(idx, m) = v_ens_tau(m)%data(i,j,k)
   idx = blk_h*nk + idx
   as(idx, m) = q_ens_tau(m)%data(i,j,k)
   end do
   end do
   end do

   end do

end subroutine ens_ics_sice_up

!=======================================================================

subroutine red_ens_sice_up(t_ens_tau, u_ens_tau, v_ens_tau, q_ens_tau, &
                           isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, &
			   nk, as, ass_method)

implicit none

type(ada_field_type), intent(inout) :: t_ens_tau(:), u_ens_tau(:), &
                                       v_ens_tau(:), q_ens_tau(:)
integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, &
                       ass_method
real, intent(in) :: as(:, :)

integer :: m, n, i, j, k, it, idx, blk_h

!  Get initial state for ensemble members

   blk_h = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)

   do m = 1, size(as, 2)

   do k = 1, nk
   do j = jsd_ens, jed_ens ! drop halo
   do i = isd_ens, ied_ens ! drop halo
   idx = (k-1)*blk_h + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   t_ens_tau(m)%data(i,j,k) = as(idx, m)
   idx = blk_h*nk + idx
   u_ens_tau(m)%data(i,j,k) = as(idx, m)
   idx = blk_h*nk + idx
   v_ens_tau(m)%data(i,j,k) = as(idx, m)
   idx = blk_h*nk + idx
   q_ens_tau(m)%data(i,j,k) = as(idx, m)
   end do
   end do
   end do

   end do

end subroutine red_ens_sice_up
!#######################################################################

end module model_sice_up_mod
