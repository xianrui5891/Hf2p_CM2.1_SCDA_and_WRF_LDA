module model_oda_mod
  !-----------------------------------------------------------------------
  !
  !   model_mod that couples component models for the atmosphere,
  !   ocean (amip), land, and sea-ice using the exchange module
  !
  !-----------------------------------------------------------------------

  use oda_types_mod, only : field_type
  use mpp_mod, only : mpp_pe !lulv
  implicit none
  public init_model, ens_ics, red_ens, get_model_size, ens_ics_stn !lulv add ens_ics_stn

  integer :: model_size

contains

  subroutine init_model(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, ass_method)
    integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk
    integer, intent(in) :: ass_method

    model_size = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*nk*2 +&
         & 2*(ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)
  end subroutine init_model

  integer function get_model_size()

    get_model_size = model_size
  end function get_model_size

  subroutine ens_ics(temp_ens_tau, salt_ens_tau, &
       & uflx_ens, vflx_ens, &
       & isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, as, ass_method)
    type(field_type), intent(in), dimension(:) :: temp_ens_tau, salt_ens_tau
    type(field_type), intent(in), dimension(:) :: uflx_ens, vflx_ens
    integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk
    real, intent(inout), dimension(:,:) :: as
    integer, intent(in) :: ass_method

    integer :: m, i, j, k, idx, blk_h

    !  Get initial state for ensemble members
    blk_h = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)

    as(:,:) = 0.0

    do m=1, size(as, 2)

       do k=1, nk
          do j=jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
             do i=isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
                idx = (k-1)*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
                as(idx, m) = temp_ens_tau(m)%data(i,j,k)
                idx = blk_h*nk + idx
                as(idx, m) = salt_ens_tau(m)%data(i,j,k)
             end do
          end do
       end do

       do j=jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
          do i=isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
             idx =  2*nk*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
             as(idx, m) = uflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = vflx_ens(m)%data(i,j,1)
          end do
       end do
    end do
  end subroutine ens_ics

  !----------------------------
  !lulv add subroutine ens_ics_stn below
  subroutine ens_ics_stn(T_ens_stn,temp_ens_tau, salt_ens_tau, &
       & uflx_ens, vflx_ens, &
       & isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, as, asb, ass_method)
    type(field_type), intent(in), dimension(:) :: temp_ens_tau, salt_ens_tau
    type(field_type), intent(in), dimension(:) :: uflx_ens, vflx_ens
    type(field_type), intent(in), dimension(:,:) :: T_ens_stn
    integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk
    real, intent(inout), dimension(:,:) :: as, asb
    integer, intent(in) :: ass_method

    integer :: m, i, j, k, idx, blk_h, blk, nuv

    !  Get initial state for ensemble members
    blk_h = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)
    blk=blk_h*nk

    as(:,:) = 0.0 
    asb(:,:) = 0.0

!convert temp_ens_tau and salt_ens_tau to as
    do m=1, size(as, 2)
       do k=1, nk
          do j=jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
             do i=isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
                idx = (k-1)*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
                as(idx, m) = temp_ens_tau(m)%data(i,j,k)
                idx = blk_h*nk + idx 
                as(idx, m) = salt_ens_tau(m)%data(i,j,k)
             end do
          end do
       end do

       do j=jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
          do i=isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
             idx =  2*nk*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
             as(idx, m) = uflx_ens(m)%data(i,j,1)
             idx = blk_h + idx 
             as(idx, m) = vflx_ens(m)%data(i,j,1)
          end do
       end do
    end do

!convert T_ens_stn to asb
    do nuv=1,2
      do m=1, size(as, 2)
         do k=1, nk
            do j=jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
               do i=isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
                  idx = (nuv-1)*blk + (k-1)*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
                  asb(idx, m) = T_ens_stn(nuv,m)%data(i,j,k)
               end do
            end do
         end do
      end do
    end do
  end subroutine ens_ics_stn
  !lulv add subroutine ens_ics_stn above
  !----------------------------
  subroutine red_ens(temp_ens_tau, salt_ens_tau, &
       & uflx_ens, vflx_ens, &
       & isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, as, ass_method)
    type(field_type), intent(inout), dimension(:) :: temp_ens_tau, salt_ens_tau
    type(field_type), intent(inout), dimension(:) :: uflx_ens, vflx_ens
    integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk
    real, intent(in), dimension(:,:) :: as
    integer, intent(in) :: ass_method

    integer :: m, i, j, k, idx, blk_h

    !  Get initial state for ensemble members
    blk_h = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)

    do m=1, size(as, 2)
       do k=1, nk
          do j=jsd_ens, jed_ens ! drop halo
             do i=isd_ens, ied_ens ! drop halo
                idx = (k-1)*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
                temp_ens_tau(m)%data(i,j,k) = as(idx, m)
                idx = blk_h*nk + idx
                salt_ens_tau(m)%data(i,j,k) = as(idx, m)
             end do
          end do
       end do

       do j=jsd_ens, jed_ens ! drop halo
          do i=isd_ens, ied_ens ! drop halo
             idx = 2*nk*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
             uflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             vflx_ens(m)%data(i,j,1) = as(idx, m)
          end do
       end do
    end do
  end subroutine red_ens
end module model_oda_mod
