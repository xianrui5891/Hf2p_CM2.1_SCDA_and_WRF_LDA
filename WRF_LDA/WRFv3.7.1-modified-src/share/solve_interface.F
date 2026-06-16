!WRF:MEDIATION_LAYER:ADT_BARRIER
!

SUBROUTINE solve_interface ( grid ) 

   USE module_domain
   USE module_timing
   USE module_driver_constants
   USE module_configure
   USE module_wrf_error
   USE python_interface
   USE module_dm, ONLY : wrf_dm_max_real, wrf_dm_sum_real

   IMPLICIT NONE

   INTERFACE
#if (EM_CORE == 1 && DA_CORE != 1)
# include  <solve_em.int>
#endif
#if (NMM_CORE == 1)
# include  <solve_nmm.int>
#endif
#if (COAMPS_CORE == 1)
# include  <solve_coamps.int>
#endif
#if (EXP_CORE == 1)
# include  <solve_exp.int>
#endif
   END INTERFACE

   TYPE(domain) , INTENT(INOUT)  :: grid
   TYPE (grid_config_rec_type)   :: config_flags

   INTEGER     :: idum1, idum2
   CHARACTER(len=256) :: dbg_msg

   REAL :: minutesSinceSimulationStart_
   CHARACTER(len=256) :: current_timestr_
   INTEGER :: x_size, y_size, z_size

   REAL, DIMENSION(:,:,:), ALLOCATABLE :: errf_u,errf_v,errf_t,temp_t
   REAL, DIMENSION(:,:), ALLOCATABLE :: errf_q , xlat, xlon
   REAL, parameter :: r_d = 287.04
   REAL, parameter :: cp = 1004.67
   REAL, parameter :: p00 = 100000.0

   INTEGER :: rank, ierr
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   REAL :: t_lda_start, t_python_start
   REAL :: lda_total_time, lda_f2py_time
   REAL :: nproc_real
   REAL :: lda_total_max, lda_f2py_max
   REAL :: lda_total_avg, lda_f2py_avg
   INTEGER, SAVE :: lda_timing_count = 0
   REAL, SAVE :: lda_total_max_sum = 0.0
   REAL, SAVE :: lda_total_avg_sum = 0.0
   REAL, SAVE :: lda_f2py_max_sum = 0.0
   REAL, SAVE :: lda_f2py_avg_sum = 0.0
   INTEGER                         :: ids , ide , jds , jde , kds , kde , &
                                      ims , ime , jms , jme , kms , kme , &
                                      ips , ipe , jps , jpe , kps , kpe
   INTEGER ::              imsx, imex, jmsx, jmex, kmsx, kmex,    &
                           ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                           imsy, imey, jmsy, jmey, kmsy, kmey,    &
                           ipsy, ipey, jpsy, jpey, kpsy, kpey

    CHARACTER*255              :: mess
    

   CALL model_to_grid_config_rec ( grid%id , model_config_rec , config_flags )
   CALL set_scalar_indices_from_config ( grid%id , idum1 , idum2 )

#if (EM_CORE == 1 && DA_CORE != 1)
   CALL solve_em  ( grid , config_flags  &
!
# include <actual_new_args.inc>
!
               )
    write(*,*), "grid%id==================" , grid%id

!following is intended to call the python functions to do data assimilation
    IF( grid%id == 1) THEN
        
        CALL domain_clock_get(grid, minutesSinceSimulationStart = minutesSinceSimulationStart_)

        IF( MOD( minutesSinceSimulationStart_ , REAL(config_flags%history_interval) ) == 0. ) THEN
            WRITE(mess,'("WRF_LDA_START DOMAIN ID = ",I3)') grid%id
            CALL WRF_MESSAGE ( mess )
            t_lda_start = now_time()
            lda_f2py_time = 0.0

            CALL get_ijk_from_grid (  grid ,                   &
                             ids, ide, jds, jde, kds, kde,    &
                             ims, ime, jms, jme, kms, kme,    &
                             ips, ipe, jps, jpe, kps, kpe,    &
                             imsx, imex, jmsx, jmex, kmsx, kmex,    &
                             ipsx, ipex, jpsx, jpex, kpsx, kpex,    &
                             imsy, imey, jmsy, jmey, kmsy, kmey,    &
                             ipsy, ipey, jpsy, jpey, kpsy, kpey )

            allocate( errf_u( ids:ide-1, kds:kde-1, jds:jde-1 ) )
            allocate( errf_v( ids:ide-1, kds:kde-1, jds:jde-1 ) )
            allocate( errf_t( ids:ide-1, kds:kde-1, jds:jde-1 ) )
            allocate( temp_t( ims:ime, kms:kme, jms:jme) )
            allocate( errf_q( ids:ide-1, jds:jde-1 ) )
            allocate( xlat( ids:ide-1, jds:jde-1 ) ) 
            allocate( xlon( ids:ide-1, jds:jde-1 ) ) 

            errf_u = -9999.999
            errf_v = -9999.999
            errf_t = -9999.999
            errf_q = -9999.999
            xlat = -9999.999
            xlon = -9999.999

            WRITE(mess,*) "start_da_data_merge!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_patch_to_global_real(grid%u_2, errf_u, grid%domdesc, 'xyz', 'xzy', &
                                            ids, ide-1, kds, kde-1, jds, jde-1, &
                                            ims, ime, kms, kme, jms, jme, &
                                            ips, min(ipe,ide-1), kps, kpe-1, jps, min(jpe,jde-1) )
            WRITE(mess,*) "da_data_merge done 1!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_patch_to_global_real(grid%v_2, errf_v, grid%domdesc, 'xyz', 'xzy', &
                                            ids, ide-1, kds, kde-1, jds, jde-1, &
                                            ims, ime, kms, kme, jms, jme, &
                                            ips, min(ipe,ide-1), kps, kpe-1, jps, min(jpe,jde-1) )
            WRITE(mess,*) "da_data_merge done 2!"
            CALL WRF_MESSAGE ( mess )

            temp_t = (grid%t_2 + 300.0) * (((grid%p + grid%pb) / p00) ** (r_d / cp))

            CALL wrf_patch_to_global_real(temp_t, errf_t, grid%domdesc, 'xyz', 'xzy', &
                                            ids, ide-1, kds, kde-1, jds, jde-1, &
                                            ims, ime, kms, kme, jms, jme, &
                                            ips, min(ipe,ide-1), kps, kpe-1, jps, min(jpe,jde-1) )
            WRITE(mess,*) "da_data_merge done 3!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_patch_to_global_real(grid%q2, errf_q, grid%domdesc, 'xy', 'xy', &
                                            ids, ide-1, jds, jde-1, 1, 1, &
                                            ims, ime, jms, jme, 1, 1, &
                                            ips, min(ipe,ide-1), jps, min(jpe,jde-1), 1, 1 )
            WRITE(mess,*) "da_data_merge done 4!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_patch_to_global_real(grid%xlat, xlat, grid%domdesc, 'xy', 'xy', &
                                            ids, ide-1 , jds , jde-1 , 1 , 1 , &
                                            ims, ime   , jms , jme   , 1 , 1 , &
                                            ips, min(ipe,ide-1), jps , min(jpe,jde-1), 1 , 1 )
            WRITE(mess,*) "da_data_merge done 5!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_patch_to_global_real(grid%xlong, xlon, grid%domdesc, 'xy', 'xy', &
                                            ids, ide-1 , jds , jde-1 , 1 , 1 , &
                                            ims, ime   , jms , jme   , 1 , 1 , &
                                            ips, min(ipe,ide-1), jps , min(jpe,jde-1), 1 , 1 )
            WRITE(mess,*) "da_data_merge done 6!"
            CALL WRF_MESSAGE ( mess )
            WRITE(mess,*) "end_da_data_merge!"
            CALL WRF_MESSAGE ( mess )
            CALL domain_clock_get(grid, current_timestr = current_timestr_)

            x_size = size(errf_u, 1)
            y_size = size(errf_u, 2)
            z_size = size(errf_u, 3)
            
            WRITE(mess,*) "start_lda!"
            CALL WRF_MESSAGE ( mess )
            t_python_start = now_time()
            CALL python_da(current_timestr_, x_size, y_size, z_size, errf_u, errf_v, errf_t, errf_q, xlat, xlon)
            lda_f2py_time = now_time() - t_python_start
            WRITE(mess,*) "end_lda!"
            CALL WRF_MESSAGE ( mess )
            
            WRITE(mess,*) "start_lda_data_distribute!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_global_to_patch_real(errf_u, grid%u_2, grid%domdesc, 'xyz', 'xzy', &
                                            ids, ide-1, kds, kde-1, jds, jde-1, &
                                            ims, ime, kms, kme, jms, jme, &
                                            ips, min(ipe,ide-1), kps, kpe-1, jps, min(jpe,jde-1) )
            WRITE(mess,*) "lda_data_distribute done 1!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_global_to_patch_real(errf_v, grid%v_2, grid%domdesc, 'xyz', 'xzy', &
                                            ids, ide-1, kds, kde-1, jds, jde-1, &
                                            ims, ime, kms, kme, jms, jme, &
                                            ips, min(ipe,ide-1), kps, kpe-1, jps, min(jpe,jde-1) )
            WRITE(mess,*) "lda_data_distribute done 2!"
            CALL WRF_MESSAGE ( mess )
            CALL wrf_global_to_patch_real(errf_t, temp_t, grid%domdesc, 'xyz', 'xzy', &
                                            ids, ide-1, kds, kde-1, jds, jde-1, &
                                            ims, ime, kms, kme, jms, jme, &
                                            ips, min(ipe,ide-1), kps, kpe-1, jps, min(jpe,jde-1) )

            grid%t_2 = temp_t * (p00 / (grid%p + grid%pb)) ** (r_d / cp) - 300.0

            WRITE(mess,*) "lda_data_distribute done 3!"
            CALL WRF_MESSAGE ( mess )
            WRITE(mess,*) "end_lda_data_distribute!"
            CALL WRF_MESSAGE ( mess )

            DEALLOCATE(errf_u)
            DEALLOCATE(errf_v)
            DEALLOCATE(errf_t)
            DEALLOCATE(temp_t)
            DEALLOCATE(errf_q)
            DEALLOCATE(xlon)
            DEALLOCATE(xlat)
            lda_total_time = now_time() - t_lda_start

            nproc_real = wrf_dm_sum_real(1.0)
            lda_total_max = wrf_dm_max_real(lda_total_time)
            lda_f2py_max = wrf_dm_max_real(lda_f2py_time)
            lda_total_avg = wrf_dm_sum_real(lda_total_time) / MAX(nproc_real, 1.0)
            lda_f2py_avg = wrf_dm_sum_real(lda_f2py_time) / MAX(nproc_real, 1.0)
            lda_timing_count = lda_timing_count + 1
            lda_total_max_sum = lda_total_max_sum + lda_total_max
            lda_total_avg_sum = lda_total_avg_sum + lda_total_avg
            lda_f2py_max_sum = lda_f2py_max_sum + lda_f2py_max
            lda_f2py_avg_sum = lda_f2py_avg_sum + lda_f2py_avg

            IF ( wrf_dm_on_monitor() ) THEN
                WRITE(mess,'("WRF_LDA_TIMING call=",I6," da_total max/avg = ",F10.4," / ",F10.4," s")') &
                    lda_timing_count, lda_total_max, lda_total_avg
                CALL WRF_MESSAGE ( mess )
                WRITE(mess,'("WRF_LDA_TIMING call=",I6," f2py_python max/avg = ",F10.4," / ",F10.4," s")') &
                    lda_timing_count, lda_f2py_max, lda_f2py_avg
                CALL WRF_MESSAGE ( mess )
                WRITE(mess,'("WRF_LDA_TIMING_SUM calls=",I6," da_total max_sum/avg_sum = ",F12.4," / ",F12.4," s")') &
                    lda_timing_count, lda_total_max_sum, lda_total_avg_sum
                CALL WRF_MESSAGE ( mess )
                WRITE(mess,'("WRF_LDA_TIMING_SUM calls=",I6," f2py_python max_sum/avg_sum = ",F12.4," / ",F12.4," s")') &
                    lda_timing_count, lda_f2py_max_sum, lda_f2py_avg_sum
                CALL WRF_MESSAGE ( mess )
                WRITE(mess,'("WRF_LDA_TIMING_NOTE da_total includes gather/python/distribution.")')
                CALL WRF_MESSAGE ( mess )
                WRITE(mess,'("WRF_LDA_TIMING_NOTE f2py_python is only CALL python_da.")')
                CALL WRF_MESSAGE ( mess )
            END IF
        END IF
    END IF
# if ( WRF_CHEM == 1 )
   IF ( config_flags%chem_opt > 0 .or. config_flags%tracer_opt > 0 ) THEN

     CALL chem_driver  ( grid , config_flags  &
!
# include <actual_new_args.inc>
!
               )
   ENDIF
# endif
#endif
#if (NMM_CORE == 1)
   CALL solve_nmm  ( grid , config_flags   &
!
# include <actual_new_args.inc>
!
               )
# if ( WRF_CHEM == 1 )
   IF ( config_flags%chem_opt > 0  ) THEN

     CALL chem_driver  ( grid , config_flags  &
!
# include <actual_new_args.inc>
!  
               )
   ENDIF
# endif
#endif
#if (COAMPS_CORE == 1)
   CALL solve_coamps  ( grid , config_flags   &
!
# include <actual_new_args.inc>
!
               )
#endif

! ###### 4. Edit share/solve_interface.F to add call to experimental core

#if (EXP_CORE == 1)
   CALL solve_exp  ( grid              &
!
# include <exp_actual_args.inc>
!
               )
#endif

END SUBROUTINE solve_interface
