module model_oda_mod
  !-----------------------------------------------------------------------
  !
  !   model_mod that couples component models for the atmosphere,
  !   ocean (amip), land, and sea-ice using the exchange module
  !
  !-----------------------------------------------------------------------

  use oda_types_mod, only : field_type

  implicit none
  public init_model, ens_ics, red_ens, get_model_size 

  integer :: model_size

contains

  subroutine init_model(isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, ass_method)
    integer, intent(in) :: isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk
    integer, intent(in) :: ass_method

    model_size = (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)*nk*4 +&
         & 6*(ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1) +&
         & (ied_ens-isd_ens+2*halox+1)*(jed_ens-jsd_ens+2*haloy+1)
  end subroutine init_model

  integer function get_model_size()

    get_model_size = model_size
  end function get_model_size

  subroutine ens_ics(temp_ens_tau, salt_ens_tau, u_ens_tau, v_ens_tau,&
       & uflx_ens, vflx_ens, tflx_ens, qflx_ens, lwflx_ens, swflx_ens,&
       & eta_ens,&
       & isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, as, ass_method)
    type(field_type), intent(in), dimension(:) :: temp_ens_tau, salt_ens_tau
    type(field_type), intent(in), dimension(:) :: u_ens_tau, v_ens_tau
    type(field_type), intent(in), dimension(:) :: uflx_ens, vflx_ens
    type(field_type), intent(in), dimension(:) :: tflx_ens, qflx_ens
    type(field_type), intent(in), dimension(:) :: lwflx_ens, swflx_ens, eta_ens
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
                idx = blk_h*nk + idx
                as(idx, m) = u_ens_tau(m)%data(i,j,k)
                idx = blk_h*nk + idx
                as(idx, m) = v_ens_tau(m)%data(i,j,k)
             end do
          end do
       end do

       do j=jsd_ens-haloy, jed_ens+haloy ! add halo grids 4 each pedomain in y
          do i=isd_ens-halox, ied_ens+halox ! add halo grids 4 each pedomain in x
             idx =  4*nk*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
             as(idx, m) = uflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = vflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = tflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = qflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = lwflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = swflx_ens(m)%data(i,j,1)
             idx = blk_h + idx
             as(idx, m) = eta_ens(m)%data(i,j,1)
          end do
       end do
    end do
  end subroutine ens_ics

  subroutine red_ens(temp_ens_tau, salt_ens_tau, u_ens_tau, v_ens_tau,&
       & uflx_ens, vflx_ens, tflx_ens, qflx_ens, lwflx_ens, swflx_ens,&
       & eta_ens,&
       & isd_ens, ied_ens, jsd_ens, jed_ens, halox, haloy, nk, as, ass_method)
    type(field_type), intent(inout), dimension(:) :: temp_ens_tau, salt_ens_tau
    type(field_type), intent(inout), dimension(:) :: u_ens_tau, v_ens_tau
    type(field_type), intent(inout), dimension(:) :: uflx_ens, vflx_ens
    type(field_type), intent(inout), dimension(:) :: tflx_ens, qflx_ens
    type(field_type), intent(inout), dimension(:) :: lwflx_ens, swflx_ens, eta_ens
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
                idx = blk_h*nk + idx
                u_ens_tau(m)%data(i,j,k) = as(idx, m)
                idx = blk_h*nk + idx
                v_ens_tau(m)%data(i,j,k) = as(idx, m)
             end do
          end do
       end do

       do j=jsd_ens, jed_ens ! drop halo
          do i=isd_ens, ied_ens ! drop halo
             idx = 4*nk*blk_h + (j-jsd_ens+haloy)*(ied_ens-isd_ens+2*halox+1) + i-isd_ens+halox+1
             uflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             vflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             tflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             qflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             lwflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             swflx_ens(m)%data(i,j,1) = as(idx, m)
             idx = blk_h + idx
             eta_ens(m)%data(i,j,1) = as(idx, m)
          end do
       end do
    end do
  end subroutine red_ens
end module model_oda_mod
