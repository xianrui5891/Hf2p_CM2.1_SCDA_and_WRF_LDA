!WRF:MEDIATION_LAYER:SOLVER

MODULE module_after_all_rk_steps

USE module_da_add_mod,only : rest_ada

CONTAINS

   !  This subroutine is called once per domain per time step.  It is outside
   !  of and after the end of the Runge-Kutta time steps, after the calls to 
   !  the explicit moisture driver, and after the polar filtering calls.  The
   !  variables in here are all up-to-date with the end of this current time 
   !  step.


   SUBROUTINE after_all_rk_steps ( grid, config_flags,                  &
                                   moist, chem, tracer, scalar,         &
                                   th_phy, pi_phy, p_phy, rho_phy,      & 
                                   p8w, t8w, dz8w,                      &
                                   curr_secs2,                          &
                                   diag_flag,                           &
                                   ids,  ide,  jds,  jde,  kds,  kde,   &
                                   ims,  ime,  jms,  jme,  kms,  kme,   &
                                   ips,  ipe,  jps,  jpe,  kps,  kpe,   &
                                   imsx, imex, jmsx, jmex, kmsx, kmex,  &
                                   ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                                   imsy, imey, jmsy, jmey, kmsy, kmey,  &
                                   ipsy, ipey, jpsy, jpey, kpsy, kpey,  &
                                        its, ite, jts, jte,kts,kte, &
                       currentTime,Varl,Varl0, corr_ps_hfqo,corr_ps_mino,corr_ps_dayo,&
                         corr_psfc6ho,corr_psfc3ho,corr_psfc1ho, corr_psfc5dayo,&
                        corr_t_dayo, corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo,&
                        corr_u_dayo, corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo,&
                        corr_v_dayo, corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo,&
                        corr_q_dayo, corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo,&
                        corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r,&
                        atm_step,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,ens_hfq_pso04,&
                       ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,&
                        ens_hfq_pso17,ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,kkk_end,var_t,var_u,&
                        corr_all)


      !=============================================================
      !  USE Association for Generic WRF Infrastructure
      !=============================================================

      !  Pick up the number of members for each of the 4d arrays - for declaration purposes.

      USE module_state_description, ONLY: num_moist, num_chem, num_tracer, num_scalar

      !  This gives us the type definition for grid (domain)

      USE module_domain, ONLY : domain

      !  All of the information from the namelist is in config_flags.  The
      !  type declaration for this puppy must be available.

      USE module_configure, ONLY : grid_config_rec_type

#ifdef DM_PARALLEL
      !  Ensure some of the fancy diagnostics variables that need to
      !  talk to other patches can do so.

      USE module_dm, ONLY : &
                  local_communicator, mytask, ntasks, ntasks_x, ntasks_y                   &
                 ,local_communicator_periodic, wrf_dm_maxval

      USE module_comm_dm, ONLY : &
                  halo_em_phys_w_sub
#endif

      !=============================================================
      !  USE Association for the Diagnostic Packages
      !=============================================================
   USE module_machine  !wangj0710
   USE module_driver_constants  !wangj0710      
    USE module_dm
   USE module_utility
      USE module_diagnostics_driver, ONLY : diagnostics_driver
      USE module_ada_types_mod, only : assim_data_type
      USE module_ada_driver_fv_mod, only : init_fv_ada,fv_ada,mda_o

      IMPLICIT NONE


      !=============================================================
      !  Subroutine Arguments
      !=============================================================

      !  Arguments passed in.  All of the diagnostics are part of the grid structure, so
      !  even though we are not changing any of the fundamental variables, we are computing
      !  the diagnostics.  Therefore grid is INOUT.
    INCLUDE 'mpif.h'

      TYPE ( domain ), INTENT(INOUT) :: grid

      !  We are not changing any of the namelist settings.

      TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags

      !  The 4d arrays are input only, no mods to them.

      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_chem  ) , INTENT(IN) :: chem
      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer) , INTENT(IN) :: tracer
      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar) , INTENT(IN) :: scalar

      !  A few handy 3d arrays computed for the physics scheme: pressure (Pa) and
      !  temperature (K), on both half (_phy) and full levels.

      REAL , DIMENSION(ims:ime,kms:kme,jms:jme)            , INTENT(IN) :: th_phy  , &
                                                                           p_phy   , &
                                                                           pi_phy  , &
                                                                           rho_phy , &
                                                                           dz8w    , &
                                                                           p8w     , &
                                                                           t8w

      !  Time (s) since the beginning of the restart.

      REAL :: curr_secs2

      !  Is this to be a history output time?  If so, compute the diagnostics.

      LOGICAL :: diag_flag

      !  The sundry dimensions required to keep a model running smoothly:
      !     The first letter:
      !        i: refers to the nominally west east direction, the inner-most (fastest)
      !           incrementing index
      !        j: refers to the nominally south north direction, the outer-most (slowest)
      !           incrementing index
      !        k: refers to the vertical direction form bottom to top, the second dimension
      !           in all 3d arrays
      !     The second letter: 
      !        d: refers to the domain size, the geophysical extent of the entire domain,
      !           not used in dimensions or looping, used to determine when we are close to
      !           the edge of the boundary
      !        m: refers to the memory size size, all 2d and 3d arrays from the Registry
      !           (passed into here via the grid structure or the I1 variables [such as
      !           p_phy, for example]) use these values for dimensioning
      !        p: refers to the patch size, the extent over which computational loops run

      INTEGER , INTENT(IN) :: ids, ide, jds, jde, kds, kde,     &
                              ims, ime, jms, jme, kms, kme,     &
                              ips, ipe, jps, jpe, kps, kpe

      !  Hopefully unnecessary, these are the filtered dimensions.

      INTEGER , INTENT(IN) :: imsx,imex,jmsx,jmex,kmsx,kmex,    &
                              ipsx,ipex,jpsx,jpex,kpsx,kpex,    &
                              imsy,imey,jmsy,jmey,kmsy,kmey,    &
                              ipsy,ipey,jpsy,jpey,kpsy,kpey

        integer :: i,j,k,i_start, i_end, j_start,j_end,iii,jjj,kkk,iii_end,jjj_end
        integer, intent(in) :: its, ite, jts, jte,kts,kte
        TYPE(WRFU_Time)                :: currentTime
        INTEGER,intent(in) :: atm_step 
 inTEGER                         :: yr, day, hr, minute, sec, rc
        type(assim_data_type), INTENT(INOUT) :: Varl,Varl0
   real, dimension(:,:), allocatable,intent(inout) :: corr_ps_hfqo
   real, dimension(:,:), allocatable,intent(inout) :: corr_ps_mino
   real, dimension(:,:), allocatable,intent(inout) :: corr_ps_dayo,corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_t_dayo,corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_u_dayo,corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_v_dayo,corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_q_dayo,corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo
   real,dimension(:,:,:), allocatable,intent(inout) :: corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso07,ens_hfq_pso08,ens_hfq_pso09
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso16,ens_hfq_pso17,ens_hfq_pso18
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso19,ens_hfq_pso20
!   real,DIMENSION(grid%sm31:grid%em31,grid%sm32:grid%em32,grid%sm33:grid%em33), intent(inout) :: var_t
   REAL , DIMENSION(ims:ime,kms:kme,jms:jme), intent(inout) :: var_t,var_u
   INTEGER , INTENT(IN) :: kkk_end
      REAL , DIMENSION(ims:ime,kms:kme,jms:jme,num_moist ) , INTENT(IN) :: moist
     real, dimension(:,:,:), allocatable, intent(inout) ::  corr_all
#ifdef DM_PARALLEL
      !=============================================================
      !  Include patch communications
      !=============================================================
#     include "HALO_EM_PHYS_W.inc"
#endif

      !=============================================================
      !  Start of executable code
      !=============================================================


      CALL wrf_debug ( 100 , '--> TOP OF AFTER ALL RK STEPS' ) 
      CALL wrf_debug ( 100 , '--> CALLING DIAGNOSTICS DRIVER' )

      CALL diagnostics_driver ( grid, config_flags,               &
                                moist, chem, tracer, scalar,         &
                                th_phy, pi_phy, p_phy, rho_phy,      & 
                                p8w, t8w, dz8w,                      &
                                curr_secs2,                          &
                                diag_flag,                           &
                                ids,  ide,  jds,  jde,  kds,  kde,   &
                                ims,  ime,  jms,  jme,  kms,  kme,   &
                                ips,  ipe,  jps,  jpe,  kps,  kpe,   &
                                imsx, imex, jmsx, jmex, kmsx, kmex,  &
                                ipsx, ipex, jpsx, jpex, kpsx, kpex,  &
                                imsy, imey, jmsy, jmey, kmsy, kmey,  &
                                ipsy, ipey, jpsy, jpey, kpsy, kpey   )

    i_start = its
    i_end   = min( ite,ide-1 )
    j_start = jts
    j_end   = min( jte,jde-1 )

!   do j = jms,jme
!   do i = ims,ime
!      write(*,*) 'U_2_after',i,j,grid%u_2(i,1,j)
!       write(*,*) 'psfc_after',i,j,grid%psfc(i,j)
!   enddo
!   enddo


CALL WRFU_TimeGet( currentTime, YY=yr, dayOfYear=day, H=hr,M=minute,S=sec,rc=rc)
!    if (grid%hea_hfq_1 .or. grid%hea_hfq_min) then
!        do j = j_start, j_end
!        do i = i_start, i_end

! if ( atm_step >= 28800 )  then
        do j = jms,jme
        do i = ims,ime

        Varl%ps(i,j)=grid%psfc(i,j)
!        write(*,*) 'varps',atm_step,i,j,Varl%ps(i,j)
        enddo
        enddo


   do j = jms,jme
   do k = 1,kkk_end
   do i = ims,ime
!        var_t(i,k,j) = grid%t_2(i,k,j)
!        var_u(i,k,j) = grid%u_2(i,k,j)
         Varl%t(i,k,j)=grid%t_2(i,k,j)
 
         Varl%u(i,k,j)=grid%u_2(i,k,j)
         Varl%v(i,k,j)=grid%v_2(i,k,j)

   enddo
   enddo
   enddo



  if (config_flags%restart)  call rest_ada(grid,config_flags,currentTime,atm_step,&
                                ims, ime, jms, jme, kms, kme,its, ite, jts,jte,kkk_end)

  if (.not. config_flags%restart  )    call fv_ada(grid,config_flags,currentTime,Varl,corr_ps_hfqo,&
                         corr_ps_mino,corr_ps_dayo,corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo, &
                        corr_t_dayo, corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo,&
                        corr_u_dayo, corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo,&
                        corr_v_dayo, corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo,&
                        corr_q_dayo, corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo,&
                        atm_step,ims, ime, jms, jme, kms, kme,its, ite, jts, jte,kkk_end,&
                        grid%hea_hfq_1,grid%hea_hfq_min,grid%hea_day,grid%hea_5day,grid%hea_6h,grid%hea_3h,grid%hea_1h, &
                        grid%dt_const_ps_hfq,grid%dt_const_ps_day,grid%dt_const_ps_5day,&
                       ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,ens_hfq_pso04,&
                       ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,&
                        ens_hfq_pso17,ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,var_t,var_u )

  if (.not. config_flags%restart  )    call mda_o(grid,config_flags,currentTime,Varl,&
                                corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r,&
                        atm_step,ims, ime, jms, jme, kms, kme,&
                                 its, ite, jts, jte, kts, kte, kkk_end,&
                        ids,  ide,  jds,  jde,  kds,  kde, &
                       ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,ens_hfq_pso04,&
                       ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,&
                        ens_hfq_pso17,ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20)


write(*,*) 'after_mda'

  if (.not. config_flags%restart  )  then

    if (config_flags%debug_ps) then

    do jjj = grid%j_start(1),j_end
    do iii = grid%i_start(1),i_end
!    grid%corr_d(iii,1,jjj)     = corr_ps_dayo(iii,jjj)
!    grid%corr_6h(iii,1,jjj)    = corr_psfc6ho(iii,jjj)
!    grid%corr_3h(iii,1,jjj)    = corr_psfc3ho(iii,jjj)
!    grid%corr_1h(iii,1,jjj)    = corr_psfc1ho(iii,jjj)
!    grid%corr_5day(iii,1,jjj)  = corr_psfc5dayo(iii,jjj)

    grid%ens1(iii,1,jjj) = ens_hfq_pso01(iii,1,jjj)
    grid%ens2(iii,1,jjj) = ens_hfq_pso02(iii,1,jjj)
    grid%ens3(iii,1,jjj) = ens_hfq_pso03(iii,1,jjj)
    grid%ens4(iii,1,jjj) = ens_hfq_pso04(iii,1,jjj)
    grid%ens5(iii,1,jjj) = ens_hfq_pso05(iii,1,jjj)
    grid%ens6(iii,1,jjj) = ens_hfq_pso06(iii,1,jjj)
    grid%ens7(iii,1,jjj) = ens_hfq_pso07(iii,1,jjj)
    grid%ens8(iii,1,jjj) = ens_hfq_pso08(iii,1,jjj)
    grid%ens9(iii,1,jjj) = ens_hfq_pso09(iii,1,jjj)
    grid%ens10(iii,1,jjj) = ens_hfq_pso10(iii,1,jjj)
    grid%ens11(iii,1,jjj) = ens_hfq_pso11(iii,1,jjj)
    grid%ens12(iii,1,jjj) = ens_hfq_pso12(iii,1,jjj)
    grid%ens13(iii,1,jjj) = ens_hfq_pso13(iii,1,jjj)
    grid%ens14(iii,1,jjj) = ens_hfq_pso14(iii,1,jjj)
    grid%ens15(iii,1,jjj) = ens_hfq_pso15(iii,1,jjj)
    grid%ens16(iii,1,jjj) = ens_hfq_pso16(iii,1,jjj)
    grid%ens17(iii,1,jjj) = ens_hfq_pso17(iii,1,jjj)
    grid%ens18(iii,1,jjj) = ens_hfq_pso18(iii,1,jjj)
    grid%ens19(iii,1,jjj) = ens_hfq_pso19(iii,1,jjj)
    grid%ens20(iii,1,jjj) = ens_hfq_pso20(iii,1,jjj)
    enddo
    enddo

  else
write(*,*) 'testaf'

    do jjj = grid%j_start(1),j_end
    do kkk = 1,kkk_end
    do iii = grid%i_start(1),i_end

    grid%ens1(iii,kkk,jjj)=ens_hfq_pso01(iii,kkk,jjj)
    grid%ens2(iii,kkk,jjj)=ens_hfq_pso02(iii,kkk,jjj)
    grid%ens3(iii,kkk,jjj)=ens_hfq_pso03(iii,kkk,jjj)
    grid%ens4(iii,kkk,jjj)=ens_hfq_pso04(iii,kkk,jjj)
    grid%ens5(iii,kkk,jjj)=ens_hfq_pso05(iii,kkk,jjj)
    grid%ens6(iii,kkk,jjj)=ens_hfq_pso06(iii,kkk,jjj)
    grid%ens7(iii,kkk,jjj)=ens_hfq_pso07(iii,kkk,jjj)
    grid%ens8(iii,kkk,jjj)=ens_hfq_pso08(iii,kkk,jjj)
    grid%ens9(iii,kkk,jjj)=ens_hfq_pso09(iii,kkk,jjj)
    grid%ens10(iii,kkk,jjj)=ens_hfq_pso10(iii,kkk,jjj)
    grid%ens11(iii,kkk,jjj)=ens_hfq_pso11(iii,kkk,jjj)
    grid%ens12(iii,kkk,jjj)=ens_hfq_pso12(iii,kkk,jjj)
    grid%ens13(iii,kkk,jjj)=ens_hfq_pso13(iii,kkk,jjj)
    grid%ens14(iii,kkk,jjj)=ens_hfq_pso14(iii,kkk,jjj)
    grid%ens15(iii,kkk,jjj)=ens_hfq_pso15(iii,kkk,jjj)
    grid%ens16(iii,kkk,jjj)=ens_hfq_pso16(iii,kkk,jjj)
    grid%ens17(iii,kkk,jjj)=ens_hfq_pso17(iii,kkk,jjj)
    grid%ens18(iii,kkk,jjj)=ens_hfq_pso18(iii,kkk,jjj)
    grid%ens19(iii,kkk,jjj)=ens_hfq_pso19(iii,kkk,jjj)
    grid%ens20(iii,kkk,jjj)=ens_hfq_pso20(iii,kkk,jjj)
   enddo
   enddo
   enddo

write(*,*) 'ens',ens_hfq_pso17(100,2,130)

   endif

! if (config_flags%debug_t) then

!   do jjj = grid%j_start(1),j_end
!    do kkk = 1,kkk_end
!    do iii = grid%i_start(1),i_end

!    grid%corr_d(iii,kkk,jjj)     = corr_t_dayo(iii,kkk,jjj)
!    grid%corr_6h(iii,kkk,jjj)    = corr_t6ho(iii,kkk,jjj)
!    grid%corr_3h(iii,kkk,jjj)    = corr_t3ho(iii,kkk,jjj)
!    grid%corr_1h(iii,kkk,jjj)    = corr_t1ho(iii,kkk,jjj)
!    grid%corr_5day(iii,kkk,jjj)  = corr_t5dayo(iii,kkk,jjj)

!    enddo
!    enddo
!    enddo
! endif


! if (config_flags%debug_u) then
!    do jjj = grid%j_start(1),j_end
!    do kkk = 1,kkk_end
!    do iii = grid%i_start(1),i_end

!    grid%corr_d(iii,kkk,jjj)     = corr_u_dayo(iii,kkk,jjj)
!    grid%corr_6h(iii,kkk,jjj)    = corr_u6ho(iii,kkk,jjj)
!    grid%corr_3h(iii,kkk,jjj)    = corr_u3ho(iii,kkk,jjj)
!    grid%corr_1h(iii,kkk,jjj)    = corr_u1ho(iii,kkk,jjj)
!    grid%corr_5day(iii,kkk,jjj)  = corr_u5dayo(iii,kkk,jjj)

!    enddo
!    enddo
!    enddo
! endif


! if (config_flags%debug_v) then
!    do jjj = grid%j_start(1),j_end
!    do kkk = 1,kkk_end
!    do iii = grid%i_start(1),i_end

!    grid%corr_d(iii,kkk,jjj)     = corr_v_dayo(iii,kkk,jjj)
!    grid%corr_6h(iii,kkk,jjj)    = corr_v6ho(iii,kkk,jjj)
!    grid%corr_3h(iii,kkk,jjj)    = corr_v3ho(iii,kkk,jjj)
!    grid%corr_1h(iii,kkk,jjj)    = corr_v1ho(iii,kkk,jjj)
!    grid%corr_5day(iii,kkk,jjj)  = corr_v5dayo(iii,kkk,jjj)

!    enddo
!    enddo
!    enddo
! endif
 endif !restart

! if (config_flags%debug_q) then
!write(*,*) 'test23'

!    do jjj = grid%j_start(1),j_end
!    do kkk = 1,kkk_end
!    do iii = grid%i_start(1),i_end

!    grid%corr_d(iii,kkk,jjj)     = corr_q_dayo_r(iii,kkk,jjj)
!    grid%corr_6h(iii,kkk,jjj)    = corr_q6ho_r(iii,kkk,jjj)
!    grid%corr_3h(iii,kkk,jjj)    = corr_q3ho_r(iii,kkk,jjj)
!    grid%corr_5day(iii,kkk,jjj)  = corr_q5dayo_r(iii,kkk,jjj)

!    enddo
!    enddo
!    enddo
!write(*,*) 'test24'

! endif



   END SUBROUTINE after_all_rk_steps





END MODULE module_after_all_rk_steps
