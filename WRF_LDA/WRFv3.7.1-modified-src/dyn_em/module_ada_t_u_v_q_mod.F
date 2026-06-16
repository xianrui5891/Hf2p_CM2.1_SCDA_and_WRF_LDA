module module_ada_t_u_v_q_mod
!USE module_utility
!USE module_ada_types_mod, only : ada_grid_type, ada_field_type, assim_data_type
!USE module_loc_and_dist_mod, only : loc_type,get_dist
!USE netcdf
USE module_machine  !wangj0710
USE module_domain, ONLY : domain
USE module_configure, ONLY : grid_config_rec_type
  implicit none
! INCLUDE 'mpif.h'
  private
  public :: mk_ens_ada_dayt2  !write subroutine name

contains

   subroutine  mk_ens_ada_dayt2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top,vari)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_day
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_day = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440
        character(len=1),intent(in) :: vari
  
   write(*,*) 'z_top',z_top

         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60 



     write(*,*) 'tensemble_min_day,time_min_ens',ensemble_min_day,time_min_ens,ps_ens_tmp_min(1)%data(isd+10,jsd+10,1)
     write(*,*) 'tit_step,i_a,ps_ens_tmp_min',it_step,i_a,ps_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

! 1 gather  t/u/v for min 1-60
  select case(vari)
     case('t')
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
   end select

! 2 gather t/u/v for hour 1-24
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = ps_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

     write(*,*) 'ti_ensb,ps_ens_tmp_hour',i_ensb,ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
                ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

! 3  loop to replace 2-20 with 1-19  then calculate new day20 
        if (i_day==0) then
!3.1
          write(*,*) 'tit_step,i_a,i_day',it_step,i_a,i_day
           do m = 1, heaada_enssize-1
               do kk = 1,z_top
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)= ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
          end do
!3.2
          ps_ens_tau(heaada_enssize)%data(:,:,:) =0.
!3.3
            do m=1,24
            do jj = jsd,jed
            do kk = 1.z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              ps_ens_tmp_hour(m)%data(ii,kk,jj)*1./24.
            enddo
            enddo
            enddo
            enddo

        write(*,*) 'tit_step,ps_ens_tau',it_step,ps_ens_tau(heaada_enssize)%data(isd+10,1,jsd+10),&
                ps_ens_tau(heaada_enssize)%data(isd+10,27,jsd+10)

        endif


      write(*,*) 'ps_ens_tau20',ps_ens_tau(heaada_enssize)%data(isd+10,1,jsd+10)
  

    end subroutine mk_ens_ada_dayt2


end module module_ada_t_u_v_q_mod
