module module_da_add_mod
  USE module_domain, ONLY : domain
  USE module_configure, ONLY : grid_config_rec_type
  USE module_utility

private
public da_ptuq,rest_ada,da_real_ptuq

contains

SUBROUTINE da_ptuq(grid,config_flags,intq,qintegal,psfc_new_dry_day,corr_all, &
        mu_new_w_day, mu_2_old, &
      corr_ps_dayo,corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo,&
      corr_t_dayo, corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo,&
      corr_u_dayo, corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo,&
      corr_v_dayo, corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo,&
      PARAM_FIRST_SCALAR,num_moist,iii_end,jjj_end,kkk_end,&
        ims,  ime,  jms,  jme,  kms,  kme,current_step)


   TYPE ( domain ), INTENT(INOUT) :: grid
   TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
   INTEGER , INTENT(IN) :: ims,  ime,  jms,  jme,  kms, kme,iii_end,jjj_end,kkk_end
   REAL , DIMENSION(ims:ime,jms:jme ) , INTENT(inout) :: intq,mu_new_w_day,mu_2_old
   real,INTENT(inout) ::  qintegal
   real, dimension(:,:), allocatable,intent(inout) :: psfc_new_dry_day
   real, dimension(:,:,:), allocatable,intent(inout) :: corr_all
   real, dimension(:,:), allocatable,intent(inout) :: corr_ps_dayo,corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo
   real, dimension(:,:,:), allocatable,intent(inout) :: corr_t_dayo, corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo
   real, dimension(:,:,:), allocatable,intent(inout) :: corr_u_dayo,corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo
   real, dimension(:,:,:), allocatable,intent(inout) :: corr_v_dayo,corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo
  INTEGER , INTENT(IN) :: PARAM_FIRST_SCALAR,num_moist,current_step
  integer :: iii,jjj,kkk
  integer :: dt_const
  real,save :: muw_2



    if (config_flags%hea_day .or. config_flags%hea_5day .or. config_flags%hea_6h &
                 .or.  config_flags%hea_3h  .or. config_flags%hea_1h)  then !for ps da
 
   write(*,*)   'corrps6h,corrps3h',corr_psfc6ho(grid%i_start(1)+10,grid%j_start(1)+10),corr_psfc3ho(grid%i_start(1)+10,grid%j_start(1)+10)
 
    do jjj = grid%j_start(1),jjj_end
    do iii = grid%i_start(1),iii_end

     psfc_new_dry_day(iii,jjj) = 0.
      corr_all(iii,1,jjj) = 0.

     corr_all(iii,1,jjj) = grid%dt_const_ps_day*corr_ps_dayo(iii,jjj) +  &
                        & grid%dt_const_ps_6h*corr_psfc6ho(iii,jjj) + &
                        & grid%dt_const_ps_3h*corr_psfc3ho(iii,jjj) + grid%dt_const_ps_1h*corr_psfc1ho(iii,jjj) &
                        & + grid%dt_const_ps_5day*corr_psfc5dayo(iii,jjj)
     enddo
     enddo

    do jjj = grid%j_start(1),jjj_end
    do iii = grid%i_start(1),iii_end

         do kkk = kkk_end-1,2,-1
         qintegal=0.0
         do n=PARAM_FIRST_SCALAR,num_moist
         qintegal=    qintegal +0.5*(grid%moist(iii,kkk,jjj,n)+grid%moist(iii,kkk+1,jjj,n))
         enddo !n

        intq(iii,jjj) = intq(iii,jjj) + qintegal*(grid%rho(iii,kkk,jjj)+grid%rho(iii,kkk+1,jjj))*0.5 /(1+qintegal)*&
                  &(grid%php(iii,kkk+1,jjj)-grid%php(iii,kkk,jjj))
         enddo !k

        qintegal=0.0
        do n=PARAM_FIRST_SCALAR,num_moist
            qintegal= qintegal+0.5*(grid%moist(iii,1,jjj,n)+grid%moist(iii,2,jjj,n))
        enddo !n

            intq(iii,jjj) = intq(iii,jjj) + qintegal*(grid%rho(iii,1,jjj)+grid%rho(iii,2,jjj))*0.5 /(1+qintegal)*&
                  &(grid%php(iii,2,jjj)-grid%php(iii,1,jjj))

    enddo !iii
    enddo !jjj

    do jjj = grid%j_start(1),jjj_end
    do iii = grid%i_start(1),iii_end
            mu_2_old(iii,jjj) = grid%mu_2(iii,jjj)
    enddo !iii
    enddo !jjj

         dt_const = 1.

        do jjj = grid%j_start(1),jjj_end
        do iii = grid%i_start(1),iii_end

           psfc_new_dry_day(iii,jjj) =   grid%psfc(iii,jjj) + corr_all(iii,1,jjj) -intq(iii,jjj)
            grid%psfc(iii,jjj) = grid%psfc(iii,jjj) + corr_all(iii,1,jjj)
            mu_new_w_day(iii,jjj) = psfc_new_dry_day(iii,jjj) - grid%mub(iii,jjj) -grid%p_top -mu_2_old(iii,jjj)
                muw_2= grid%mu_2(iii,jjj)
                grid%mu_2(iii,jjj)= mu_new_w_day(iii,jjj)*dt_const+muw_2
                do kkk = 1, kkk_end
                     pressure_new = 0.01*(grid%znu(kkk)*((grid%psfc(iii,jjj) + corr_all(iii,1,jjj))-grid%p_top)+grid%p_top)
                     pressure_old = 0.01*(grid%znu(kkk)*( grid%psfc(iii,jjj) -grid%p_top)+grid%p_top)
                     grid%t_2(iii,kkk,jjj) = grid%t_2(iii,kkk,jjj) * (pressure_new/pressure_old)**0.286 
                enddo !kkk

       enddo !iii
       enddo !jjj

    write(*,*) 'mu_day',grid%mu_2(grid%i_start(1)+10,grid%j_start(1)+10),corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10)

        grid%mu_change(:,:) = mu_new_w_day(:,:)*dt_const

         grid%corr_all(:,1,:) = corr_all(:,1,:)



     endif !da ps

! t
  if (config_flags%hea_t5day .or. config_flags%hea_tday .or. config_flags%hea_t6h .or. config_flags%hea_t3h ) then
!  write(*,*) 'bt_2',grid%t_2(grid%i_start(1)+10,1,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10),&
!&grid%t_2(grid%i_start(1)+10,15,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,15,grid%j_start(1)+10)

    do jjj = grid%j_start(1),jjj_end
    do kkk = 1,kkk_end
    do iii = grid%i_start(1),iii_end
         corr_all(iii,kkk,jjj) = 0.     
         corr_all(iii,kkk,jjj) = grid%dt_const_t_day*corr_t_dayo(iii,kkk,jjj) + grid%dt_const_t_6h*corr_t6ho(iii,kkk,jjj) + &
                        & grid%dt_const_t_3h*corr_t3ho(iii,kkk,jjj)+ grid%dt_const_t_5day*corr_t5dayo(iii,kkk,jjj)
                        
    grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj)
    grid%t_2(iii,kkk,jjj) = grid%t_2(iii,kkk,jjj) +corr_all(iii,kkk,jjj)
    
    enddo
    enddo
    enddo

  write(*,*) 'at_2',grid%t_2(grid%i_start(1)+10,1,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10),&
        grid%t_2(grid%i_start(1)+10,15,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,15,grid%j_start(1)+10)
  write(*,*) 'corr_tday,corr_t5day',corr_t_dayo(grid%i_start(1)+10,15,grid%j_start(1)+10),corr_t5dayo(grid%i_start(1)+10,15,grid%j_start(1)+10)
  write(*,*) 'corr_t3h',corr_t3ho(grid%i_start(1)+10,15,grid%j_start(1)+10)


  endif ! da t
!u
  if (config_flags%hea_u5day .or. config_flags%hea_uday .or.config_flags%hea_u6h .or. config_flags%hea_u3h ) then

!  write(*,*) 'bu_2',grid%u_2(grid%i_start(1)+10,1,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10),&
!  &grid%u_2(grid%i_start(1)+10,15,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,15,grid%j_start(1)+10)

    do jjj = grid%j_start(1),jjj_end
    do kkk = 1,kkk_end
    do iii = grid%i_start(1),iii_end
         corr_all(iii,kkk,jjj) = 0.
         corr_all(iii,kkk,jjj) = grid%dt_const_u_day*corr_u_dayo(iii,kkk,jjj) + grid%dt_const_u_6h*corr_u6ho(iii,kkk,jjj) + &
                        & grid%dt_const_u_3h*corr_u3ho(iii,kkk,jjj)+ grid%dt_const_u_5day*corr_u5dayo(iii,kkk,jjj)

    grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj)
    grid%u_2(iii,kkk,jjj) = grid%u_2(iii,kkk,jjj) +corr_all(iii,kkk,jjj)

    enddo
    enddo
    enddo
  write(*,*) 'au_2',grid%u_2(grid%i_start(1)+10,1,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10),&
        grid%u_2(grid%i_start(1)+10,15,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,15,grid%j_start(1)+10)
  write(*,*) 'corr_uday,corr_u5day',corr_u_dayo(grid%i_start(1)+10,15,grid%j_start(1)+10),corr_u5dayo(grid%i_start(1)+10,15,grid%j_start(1)+10)
  write(*,*) 'corr_u6h,corr_u3h',corr_u6ho(grid%i_start(1)+10,15,grid%j_start(1)+10),corr_u3ho(grid%i_start(1)+10,15,grid%j_start(1)+10)

  endif
!v

  if (config_flags%hea_v5day .or. config_flags%hea_vday .or.config_flags%hea_v6h .or. config_flags%hea_v3h ) then

  write(*,*) 'v_2',grid%v_2(grid%i_start(1)+10,1,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10),&
  &grid%v_2(grid%i_start(1)+10,15,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,15,grid%j_start(1)+10)

    do jjj = grid%j_start(1),jjj_end
    do kkk = 1,kkk_end
    do iii = grid%i_start(1),iii_end
         corr_all(iii,kkk,jjj) = 0.
         corr_all(iii,kkk,jjj) = grid%dt_const_v_day*corr_v_dayo(iii,kkk,jjj) + grid%dt_const_v_6h*corr_v6ho(iii,kkk,jjj) + &
                        & grid%dt_const_v_3h*corr_v3ho(iii,kkk,jjj)+ grid%dt_const_v_5day*corr_v5dayo(iii,kkk,jjj)

    grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj)
    grid%v_2(iii,kkk,jjj) = grid%v_2(iii,kkk,jjj) +corr_all(iii,kkk,jjj)

    enddo
    enddo
    enddo
  write(*,*) 'v_2',grid%v_2(grid%i_start(1)+10,1,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,1,grid%j_start(1)+10),&
        grid%v_2(grid%i_start(1)+10,15,grid%j_start(1)+10),grid%corr_all(grid%i_start(1)+10,15,grid%j_start(1)+10)
  write(*,*) 'corr_vday,corr_v5day',corr_v_dayo(grid%i_start(1)+10,15,grid%j_start(1)+10),corr_v5dayo(grid%i_start(1)+10,15,grid%j_start(1)+10)
  write(*,*) 'corr_v6h,corr_v3h',corr_v6ho(grid%i_start(1)+10,15,grid%j_start(1)+10),corr_v3ho(grid%i_start(1)+10,15,grid%j_start(1)+10)

  endif


END SUBROUTINE da_ptuq


subroutine da_real_ptuq(grid,config_flags,corr_all,&
          corr_ps_dayo_r,corr_psfc6ho_r,corr_psfc3ho_r,corr_psfc1ho_r,corr_psfc5dayo_r,&
          corr_t_dayo_r, corr_t6ho_r,corr_t3ho_r,corr_t1ho_r,corr_t5dayo_r,&
          corr_u_dayo_r, corr_u6ho_r,corr_u3ho_r,corr_u1ho_r,corr_u5dayo_r,&
          corr_v_dayo_r, corr_v6ho_r,corr_v3ho_r,corr_v1ho_r,corr_v5dayo_r,&
          corr_q_dayo_r, corr_q6ho_r,corr_q3ho_r,corr_q1ho_r,corr_q5dayo_r,&
            iii_end,jjj_end,kkk_end, ims,  ime,  jms,  jme,  kms,kme,current_step)

TYPE ( domain ), INTENT(INOUT) :: grid
TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
INTEGER , INTENT(IN) :: ims,  ime,  jms,  jme,kms,kme,iii_end,jjj_end,kkk_end
real, dimension(:,:,:), allocatable,intent(inout) :: corr_all
real, dimension(:,:), allocatable,intent(inout) :: corr_ps_dayo_r,corr_psfc6ho_r,corr_psfc3ho_r,corr_psfc1ho_r,corr_psfc5dayo_r
real, dimension(:,:,:), allocatable,intent(inout) :: corr_t_dayo_r, corr_t6ho_r,corr_t3ho_r,corr_t1ho_r,corr_t5dayo_r
real, dimension(:,:,:), allocatable,intent(inout) :: corr_u_dayo_r, corr_u6ho_r,corr_u3ho_r,corr_u1ho_r,corr_u5dayo_r
real, dimension(:,:,:), allocatable,intent(inout) :: corr_v_dayo_r, corr_v6ho_r,corr_v3ho_r,corr_v1ho_r,corr_v5dayo_r
real, dimension(:,:,:), allocatable,intent(inout) :: corr_q_dayo_r, corr_q6ho_r,corr_q3ho_r,corr_q1ho_r,corr_q5dayo_r
INTEGER , INTENT(IN) :: current_step
integer :: iii,jjj,kkk
integer :: dt_const


! da t
if ( config_flags%hea_ot6h .or. config_flags%hea_ot3h .or. &
  config_flags%hea_otday .or. config_flags%hea_ot5day ) then
  do jjj = grid%j_start(1),jjj_end
  do kkk = 1,kkk_end
  do iii = grid%i_start(1),iii_end
      corr_all(iii,kkk,jjj) = 0.
      corr_all(iii,kkk,jjj) = grid%dt_const_ot_day *corr_t_dayo_r(iii,kkk,jjj) + grid%dt_const_ot_3h*corr_t3ho_r(iii,kkk,jjj) + &
                    & grid%dt_const_ot_6h*corr_t6ho_r(iii,kkk,jjj)+ grid%dt_const_ot_5day*corr_t5dayo_r(iii,kkk,jjj)

      grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj) 
      grid%t_2(iii,kkk,jjj) = grid%t_2(iii,kkk,jjj) +corr_all(iii,kkk,jjj)

      if (abs(corr_all(iii,kkk,jjj)) >0) then
      write(*,*) "update T values (corr_all) = ", corr_all(iii,kkk,jjj) , "ii = " , iii,"jj = ",jjj,"kk = ", kkk
      endif
  enddo
  enddo
  enddo
endif ! da t


! da q
if ( config_flags%hea_oq6h .or. config_flags%hea_oq3h .or. &
  config_flags%hea_oqday .or. config_flags%hea_oq5day ) then
  do jjj = grid%j_start(1),jjj_end
  do kkk = 1,kkk_end
  do iii = grid%i_start(1),iii_end

      corr_all(iii,kkk,jjj) = 0.
      corr_all(iii,kkk,jjj) = grid%dt_const_oq_day *corr_q_dayo_r(iii,kkk,jjj) + grid%dt_const_oq_3h*corr_q3ho_r(iii,kkk,jjj) + &
                    & grid%dt_const_oq_6h*corr_q6ho_r(iii,kkk,jjj)+ grid%dt_const_oq_5day*corr_q5dayo_r(iii,kkk,jjj)

      grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj)
      ! grid%q(iii,kkk,jjj) = grid%q(iii,kkk,jjj) +corr_all(iii,kkk,jjj)  ! need modify

      if (abs(corr_all(iii,kkk,jjj)) >0) then
      write(*,*) "update q values (corr_all) = ", corr_all(iii,kkk,jjj) , "ii = " , iii,"jj = ",jjj,"kk = ", kkk
      endif
  enddo
  enddo
  enddo
endif ! da q

! da u_2
if ( config_flags%hea_ou6h .or. config_flags%hea_ou3h .or. &
  config_flags%hea_ouday .or. config_flags%hea_ou5day ) then
  do jjj = grid%j_start(1),jjj_end
  do kkk = 1,kkk_end
  do iii = grid%i_start(1),iii_end

      corr_all(iii,kkk,jjj) = 0.
      corr_all(iii,kkk,jjj) = grid%dt_const_ou_day *corr_u_dayo_r(iii,kkk,jjj) + grid%dt_const_ou_3h*corr_u3ho_r(iii,kkk,jjj) + &
                    & grid%dt_const_ou_6h*corr_u6ho_r(iii,kkk,jjj)+ grid%dt_const_ou_5day*corr_u5dayo_r(iii,kkk,jjj)

      grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj) 
      grid%u_2(iii,kkk,jjj) = grid%u_2(iii,kkk,jjj) +corr_all(iii,kkk,jjj)

      if (abs(corr_all(iii,kkk,jjj)) >0) then
      write(*,*) "update U values (corr_all) = ", corr_all(iii,kkk,jjj) , "ii = " , iii,"jj = ",jjj,"kk = ", kkk
      endif
  enddo
  enddo
  enddo
endif ! da u_2


! da v_2
if ( config_flags%hea_ov6h .or. config_flags%hea_ov3h .or. &
  config_flags%hea_ovday .or. config_flags%hea_ov5day ) then
  do jjj = grid%j_start(1),jjj_end
  do kkk = 1,kkk_end
  do iii = grid%i_start(1),iii_end

      corr_all(iii,kkk,jjj) = 0.
      corr_all(iii,kkk,jjj) = grid%dt_const_ov_day *corr_v_dayo_r(iii,kkk,jjj) + grid%dt_const_ov_3h*corr_v3ho_r(iii,kkk,jjj) + &
                    & grid%dt_const_ov_6h*corr_v6ho_r(iii,kkk,jjj)+ grid%dt_const_ov_5day*corr_v5dayo_r(iii,kkk,jjj)

      grid%corr_all(iii,kkk,jjj) = corr_all(iii,kkk,jjj) 
      grid%v_2(iii,kkk,jjj) = grid%v_2(iii,kkk,jjj) +corr_all(iii,kkk,jjj)

      if (abs(corr_all(iii,kkk,jjj)) >0) then
      write(*,*) "update V values (corr_all) = ", corr_all(iii,kkk,jjj) , "ii = " , iii,"jj = ",jjj,"kk = ", kkk
      endif
  enddo
  enddo
  enddo
endif ! da v_2


end subroutine da_real_ptuq

SUBROUTINE rest_ada(grid,config_flags,currentTime,atm_step,&
                ims, ime, jms, jme, kms, kme,its, ite, jts, jte,kkk_end)

   TYPE ( domain ), INTENT(INOUT) :: grid
   TYPE (grid_config_rec_type) , INTENT(IN)          :: config_flags
   TYPE(WRFU_Time)                :: currentTime
   INTEGER,intent(in) :: atm_step
    inTEGER                         :: yr, day, hr, minute, sec, rc
   

   CALL WRFU_TimeGet( currentTime, YY=yr, dayOfYear=day,H=hr,M=minute,S=sec,rc=rc)
   
   write(*,*) 'call rest_ada',atm_step,yr,day,hr,minute

END SUBROUTINE rest_ada

end module module_da_add_mod
