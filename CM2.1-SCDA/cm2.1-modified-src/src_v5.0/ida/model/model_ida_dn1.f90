
module model_sice_dn_mod

!-----------------------------------------------------------------------
!
!   model_mod that couples component models for the atmosphere,
!   ocean (amip), land, and sea-ice using the exchange module
!
!-----------------------------------------------------------------------

use ida_types_mod, only : ida_field_type

implicit none
public init_model_sice, ens_ics_sice, red_ens_sice, get_model_size_sice

!-----------------------------------------------------------------------

integer, private :: model_size_sice

contains

!#######################################################################

subroutine init_model_sice(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, &
                      nk_ice, ass_method)

implicit none

integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, nk_ice
integer, intent(in) :: ass_method

model_size_sice = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*nk*4 + &
                  (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*4 + &
                  (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*nk_ice*2+&
		  (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*2

end subroutine init_model_sice

!-------------------------------------------------------------------------

function get_model_size_sice()

implicit none

integer :: get_model_size_sice

get_model_size_sice = model_size_sice

end function get_model_size_sice

!-----------------------------------------------------------------------

subroutine ens_ics_sice(u_ens_tau, v_ens_tau, t_ens_tau, s_ens_tau, &
                   c_ens_ice, u_ens_ice, v_ens_ice, m_ens_ice, &
                   t1_ens_ice, t2_ens_ice, &
                   isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, nk_ice, &
                   as, ass_method)

implicit none

type(ida_field_type), intent(in) :: u_ens_tau(:), v_ens_tau(:), t_ens_tau(:), &
                                    s_ens_tau(:)
type(ida_field_type), intent(in) :: c_ens_ice(:), u_ens_ice(:), v_ens_ice(:), &
                                    m_ens_ice(:)
type(ida_field_type), intent(in) :: t1_ens_ice(:), t2_ens_ice(:)
integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, nk_ice, &
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
   idx = (k-1)*blk_h + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   as(idx, m) = t_ens_tau(m)%data(i,j,k)
   idx = blk_h*nk + idx
   as(idx, m) = u_ens_tau(m)%data(i,j,k)
   idx = blk_h*nk + idx
   as(idx, m) = v_ens_tau(m)%data(i,j,k)
   idx = blk_h*nk + idx
   as(idx, m) = s_ens_tau(m)%data(i,j,k)
   end do
   end do
   end do

   do j = jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
   do i = isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
   idx = 4*blk_h*nk + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   as(idx, m) = c_ens_ice(m)%data(i,j,1)
   idx =  blk_h + idx
   as(idx, m) = u_ens_ice(m)%data(i,j,1)
   idx =  blk_h + idx
   as(idx, m) = v_ens_ice(m)%data(i,j,1)
   idx =  blk_h + idx
   as(idx, m) = m_ens_ice(m)%data(i,j,1)
   idx = 4*blk_h*nk + 4*blk_h + 2*blk_h*nk_ice + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   as(idx, m) = u_ens_ice(m)%data(i,j,1)
   idx =  blk_h + idx
   as(idx, m) = v_ens_ice(m)%data(i,j,1)
   end do
   end do

   do k = 1, nk_ice
   do j = jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
   do i = isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
   idx = 4*blk_h*nk + 4*blk_h + (k-1)*blk_h + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   as(idx, m) = t1_ens_ice(m)%data(i,j,k)
   idx = blk_h*nk_ice + idx
   as(idx, m) = t2_ens_ice(m)%data(i,j,k)
   end do
   end do
   end do

   end do

end subroutine ens_ics_sice

!=======================================================================

subroutine red_ens_sice(u_ens_tau, v_ens_tau, t_ens_tau, s_ens_tau, &
                   u_ens_ice, v_ens_ice, t1_ens_ice, t2_ens_ice, &
                   isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, nk_ice,&
                   as, ass_method)

implicit none

type(ida_field_type), intent(inout) :: u_ens_tau(:), v_ens_tau(:), &
                                    t_ens_tau(:), s_ens_tau(:)
type(ida_field_type), intent(inout) :: u_ens_ice(:), v_ens_ice(:), &
                                    t1_ens_ice(:), t2_ens_ice(:)
integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, nk_ice,&
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
   s_ens_tau(m)%data(i,j,k) = as(idx, m)
   end do
   end do
   end do

   do k = 1, nk_ice
   do j = jsd_ens, jed_ens ! drop halo
   do i = isd_ens, ied_ens ! drop halo
   idx = 4*blk_h*nk + 4*blk_h + (k-1)*blk_h + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   t1_ens_ice(m)%data(i,j,k) = as(idx, m)
   idx = blk_h*nk_ice + idx
   t2_ens_ice(m)%data(i,j,k) = as(idx, m)
   end do
   end do
   end do

   do j = jsd_ens, jed_ens ! drop halo
   do i = isd_ens, ied_ens ! drop halo
   idx = 4*blk_h*nk + 4*blk_h + 2*blk_h*nk_ice + &
         (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1)+i-isd_ens+halox+1
   u_ens_ice(m)%data(i,j,1) = as(idx, m)
   idx = blk_h + idx
   v_ens_ice(m)%data(i,j,1) = as(idx, m)
   end do
   end do

   end do

end subroutine red_ens_sice
!#######################################################################

end module model_sice_dn_mod
