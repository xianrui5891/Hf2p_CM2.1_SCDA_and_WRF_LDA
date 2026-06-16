module module_eakf_ada_mod
USE module_utility
USE module_ada_types_mod, only : ada_grid_type, ada_field_type, assim_data_type
USE module_loc_and_dist_mod, only : loc_type,get_dist
USE netcdf
USE module_machine  !wangj0710
USE module_driver_constants  !wangj0710
USE module_dm
USE module_assim_tools_mod, only : obs_increment_hfq,update_from_obs_inc,obs_increment,update_from_obs_inc_k
USE module_cov_cutoff_mod, only  : comp_cov_factor
      USE module_domain, ONLY : domain
              USE module_configure, ONLY : grid_config_rec_type
  implicit none
 INCLUDE 'mpif.h'

private
!public ensemble_filter_atm_init,ensemble_filter_atm_hfq,ensemble_filter_atm_day,ensemble_filter_atm_min,ensemble_filter_atm_6h
public ensemble_filter_atm_init,ensemble_filter_atm_hfq,ensemble_filter_atm_day,ensemble_filter_atm_6h
public ensemble_filter_atm_6h_r,ensemble_filter_atm_5day_r,ensemble_filter_atm_day_r,ensemble_filter_atm_3h_r
public ensemble_filter_atm_3h,ensemble_filter_atm_5day
real, public, parameter :: RADIUS = 6371.0e3
real,public,parameter ::   ps_dist = 10.0e3
!real,public,parameter ::   ps_dist = 30.0e3
!  real :: t_dist = 1000.0e3
!  real :: q_dist = 500.0e3
!  real :: uv_dist = 200.0e3
  real :: t_dist = 10.0e3
  real :: q_dist = 10.0e3
  real :: uv_dist = 10.0e3
real,public ::   cov_factor = 0.0

logical :: obs_flag = .false.
logical :: ass_ps = .true.
logical :: assim_flag = .false.

real, public, parameter :: PI      = 3.14159265358979323846
real, public, parameter :: DEG_TO_RAD=PI/180.
real, dimension(:), allocatable :: enso_ps_tau,enso_t_now,enso_u_now,enso_v_now,enso_q_now
real, dimension(:), allocatable :: enso_ps_hfq,enso_t,enso_u,enso_v
real, dimension(:), allocatable :: enso_t_r_o,enso_u_r_o,enso_v_r_o,enso_ps_r_o,enso_q_r_o
real, dimension(:), allocatable :: obs_inc_eakf_ps_hfq,obs_inc_eakf_t,obs_inc_eakf_u,obs_inc_eakf_v,obs_inc_eakf_q
real, dimension(:), allocatable :: obs_inc_oi_ps_hfq,obs_inc_oi_t,obs_inc_oi_u,obs_inc_oi_v,obs_inc_oi_q
real, dimension(:), allocatable :: ens_in_hfq
real, dimension(:), allocatable :: ens_inc_hfq
real, dimension(:), allocatable :: obs_inc_eakf_tr,obs_inc_eakf_ur,obs_inc_eakf_vr
real :: obs_var_t, obs_var_q, obs_var_u, obs_var_v, obs_var_ps
real :: sigma_o_ps = 200.0
real :: sigma_o_t = 1.0
  real :: sigma_o_q = 1.0
  real :: sigma_o_u = 1.0
  real :: sigma_o_v = 1.0

real :: ass_ps_lat_start = -24.90977 ! ps=t=q
real :: ass_ps_lat_end   = 66.09582
integer :: ii_do = 10
integer :: jj_do = 10
integer :: ass_method = 2 ! 0 for snz-oi, 1 for eakf inv, 2 for eakf multv

contains


subroutine ensemble_filter_atm_init(e_size, hea_e_size)
 integer, intent(in) :: e_size, hea_e_size 
 integer :: unit, ierr, istat


allocate(enso_ps_hfq(hea_e_size))
allocate(enso_t(hea_e_size)) 
allocate(enso_u(hea_e_size))
allocate(enso_v(hea_e_size))
allocate(obs_inc_eakf_ps_hfq(hea_e_size))
allocate(obs_inc_eakf_t(hea_e_size))
allocate(obs_inc_eakf_u(hea_e_size))
allocate(obs_inc_eakf_v(hea_e_size))
allocate(obs_inc_eakf_q(hea_e_size))
allocate(obs_inc_oi_ps_hfq(hea_e_size))
allocate(obs_inc_oi_u(hea_e_size))
allocate(obs_inc_oi_v(hea_e_size))
allocate(obs_inc_oi_t(hea_e_size))
allocate(obs_inc_oi_q(hea_e_size))
allocate(enso_ps_tau(hea_e_size))
allocate(enso_t_now(hea_e_size))
allocate(enso_u_now(hea_e_size))
allocate(enso_v_now(hea_e_size))
allocate(enso_q_now(hea_e_size))
allocate(ens_in_hfq(hea_e_size), ens_inc_hfq(hea_e_size))
allocate(enso_t_r_o(hea_e_size))
allocate(enso_u_r_o(hea_e_size))
allocate(enso_v_r_o(hea_e_size))
allocate(enso_q_r_o(hea_e_size))
allocate(enso_ps_r_o(hea_e_size))
allocate(obs_inc_eakf_tr(hea_e_size))

enso_ps_hfq(:) = 0.0
enso_t(:) = 0.
enso_u(:) = 0. 
enso_v(:) = 0. 
enso_ps_tau(:) = 0.0
enso_t_now(:)  = 0.
enso_u_now(:)  = 0.
enso_v_now(:)  = 0.
enso_q_now(:)  = 0.

ens_in_hfq(:) = 0.0
ens_inc_hfq(:) = 0.0
obs_inc_eakf_ps_hfq(:) = 0.0
obs_inc_oi_ps_hfq(:) = 0.0
obs_inc_eakf_t(:) = 0.
obs_inc_eakf_u(:) = 0.
obs_inc_eakf_v(:) = 0.
obs_inc_eakf_q(:) = 0.

!obs_inc_eakf_t(:) = 0.
obs_inc_oi_u(:) = 0.
obs_inc_oi_v(:) = 0.
obs_inc_oi_t(:) = 0.
obs_inc_oi_q(:) = 0.

!
enso_t_r_o(:) = 0.
enso_u_r_o(:) = 0.
enso_v_r_o(:) = 0.
enso_q_r_o(:) = 0.
enso_ps_r_o(:) = 0.
enso_t_r_o(:) = 0.

obs_var_ps = (sigma_o_ps)**2 !Pa obs fangcha
    obs_var_t = (sigma_o_t)**2
    obs_var_q = (sigma_o_q)**2
    obs_var_u = (sigma_o_u)**2
    obs_var_v = (sigma_o_v)**2


return


end subroutine ensemble_filter_atm_init




 subroutine ensemble_filter_atm_3h(grid,config_flags,t_ens_tau,t_now,q_ens_tau,u_ens_tau,u_now,v_ens_tau,v_now,&
      & ps_ens, ps_tau, atm_var,atmos_obs, isc, iec, jsc, jec, halo,npx,npy,nk,&
        T_grid,m_time,atm_step,isd, ied, jsd, jed)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: q_ens_tau, t_ens_tau
    type(ada_field_type), dimension(:), intent(inout) :: u_ens_tau, v_ens_tau
    type(ada_field_type), dimension(:), intent(in) :: ps_ens
    type(ada_field_type), intent(inout) :: ps_tau
    type(ada_field_type), intent(inout) :: t_now,u_now,v_now !add
    type(assim_data_type), intent(in) :: atm_var
    type(assim_data_type), intent(in) :: atmos_obs
    type(ada_grid_type), intent(in) :: T_grid
    type(WRFU_Time), intent(in) :: m_time
    integer, intent(in) :: isc, iec, jsc, jec, halo, npx, npy, nk
    integer, intent(in) :: isd, ied, jsd, jed
   type(loc_type) :: model_loc, obs_loc_ptq,obs_loc_u,obs_loc_v  !add
    INTEGER             :: yr, month, day, hr, minute, sec, rc
    integer, intent(in)  :: atm_step
    integer :: i_o,j_o,k_o,i_m,j_m,k_m, j_ens, i_oo, i_mm,iii
    integer :: ens_size
    real :: cor_oi, coslat, rad2deg,coslatu,coslatv !add
    real :: dist, dist0, obs_value, obs_elev, mod_elev, obs_correct
    real :: std_oi_o, std_oi_g
    INTEGER myproc
    integer :: is_halo,js_halo,ie_halo,je_halo
    integer :: ass_variable = 1
    real :: cov_factor_v, cov_factor_t, cov_factor_h,sum_obsinc
    integer:: maxi,minj
        real :: sum_ens_in_hfq1,sum_ens_in_hfq2
        character(len=1) :: vari

        sum_ens_in_hfq1=0.
        sum_ens_in_hfq2=0.
      CALL WRFU_TimeGet( m_time, YY=yr, dayOfYear=day,H=hr,M=minute,S=sec,rc=rc)
      CALL wrf_get_myproc( myproc )
  write(*,*) 'running ensemble_filter_atm_3h'
        ens_size = size(ps_ens)
  is_halo=isd
  ie_halo=ied
  js_halo=jsd
  je_halo=jed
         do j_o =  js_halo,je_halo
         do i_o =  is_halo,ie_halo
                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
!                 if (i_o <= 0 .and. i_o > npx) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx

                 if (i_o <= 0) i_oo  = 1 
                 if (i_o > npx) i_oo = npx
                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        obs_value = atmos_obs%ps(i_oo,j_o)
!                          obs_elev = atmos_obs%phis(i_oo,j_o)
                           obs_elev = grid%ht(i_oo,j_o)
!          if (  obs_elev<=100.0 .and. obs_value <=96000.0) then
!                        else
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_ps_hfq(j_ens) = ps_ens(j_ens)%data(i_o,j_o,1)
                            enso_ps_tau(j_ens) = ps_tau%data(i_o,j_o,1)
                        end do
                        obs_correct = obs_value

            call  obs_increment_hfq(config_flags,enso_ps_hfq,enso_ps_tau,ens_size,obs_correct, obs_var_ps,&
              &obs_inc_eakf_ps_hfq,obs_inc_oi_ps_hfq,std_oi_o,i_o,j_o,atm_step)
         if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat>ass_ps_lat_end) obs_inc_eakf_ps_hfq(:) = 0.0
                        sum_obsinc=sum(obs_inc_eakf_ps_hfq(:))
              if (config_flags%debug_ps)      grid%obs_incre3h(i_o,1,j_o)=obs_inc_eakf_ps_hfq(20)
                 if (sum(obs_inc_eakf_ps_hfq(:)) /= 0.0) then
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.
                        if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0) i_mm = 1
                        if (i_m > npx) i_mm = npx

                           model_loc%lon = T_grid%x(i_mm)
                           if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                model_loc%lat = T_grid%y(j_m)
                                ass_variable = 1
                                ens_inc_hfq(:) = 0.0
                                dist = get_dist(model_loc, obs_loc_ptq, atm_step)

                                dist = RADIUS * sqrt(dist)
                                dist0 = ps_dist*coslat
                                cov_factor_h = comp_cov_factor(dist, dist0)
                                cor_oi = cov_factor_h
                                cov_factor = cov_factor_h
                                  do j_ens = 1, ens_size
                                    ens_in_hfq(j_ens) = ps_ens(j_ens)%data(i_m,j_m,1)
                                 end do
                sum_ens_in_hfq1=sum(ens_in_hfq)
                sum_ens_in_hfq2=sum(ens_in_hfq(:))
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                vari='p'
                                   call update_from_obs_inc(config_flags,enso_ps_hfq,obs_inc_eakf_ps_hfq,&
                                        & obs_inc_oi_ps_hfq,ens_in_hfq(:),ens_size,&
                                        & ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        & ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              end if

                              ps_tau%data(i_m,j_m,1) = ps_tau%data(i_m,j_m,1)+ens_inc_hfq(ens_size)
                       end if

                    enddo
                    enddo

                endif
                endif
!                endif

        end do
        end do

! process 3D obs now
! 1- 1.1   

         do j_o =  js_halo,je_halo
         do k_o = 1, nk
         do i_o =  is_halo,ie_halo

                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
!                 if (i_o <= 0 .and. i_o > npx) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx
                 if (i_o <= 0) i_oo  = 1
                 if (i_o > npx) i_oo = npx

                        obs_loc_u%lon   = T_grid%xu(i_oo)
                        obs_loc_v%lon   = T_grid%xv(i_oo)
                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_u%lon > 360.0 ) obs_loc_u%lon = obs_loc_u%lon-360.0
                        if ( obs_loc_v%lon > 360.0 ) obs_loc_v%lon = obs_loc_v%lon-360.0
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        obs_loc_u%lat = T_grid%yu(j_o)
                        obs_loc_v%lat = T_grid%yv(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        if ( abs(obs_loc_u%lat) <= 80.0 ) then
                         coslatu = cos(obs_loc_u%lat*DEG_TO_RAD)
                        else
                         coslatu = cos(80.0*DEG_TO_RAD)
                        end if

                        if ( abs(obs_loc_v%lat) <= 80.0 ) then
                         coslatv = cos(obs_loc_v%lat*DEG_TO_RAD)
                        else
                         coslatv = cos(80.0*DEG_TO_RAD)
                        endif
! 1.2 t
                        obs_value = atmos_obs%t(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_t(j_ens) = t_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_t_now(j_ens) = t_now%data(i_o,k_o,j_o)
                        end do
                        obs_correct = obs_value

                call  obs_increment_hfq(config_flags,enso_t,enso_t_now,ens_size,obs_correct,obs_var_t,&
                       &obs_inc_eakf_t,obs_inc_oi_t,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat> ass_ps_lat_end) obs_inc_eakf_t(:) = 0.0
                if (config_flags%debug_t) grid%obs_incre3h(i_o,k_o,j_o)=obs_inc_eakf_t(20)
! 1.3 q
! 1.4 u
                        obs_value = atmos_obs%u(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_u(j_ens) = u_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_u_now(j_ens) = u_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value

                call obs_increment_hfq(config_flags,enso_u,enso_u_now,ens_size,obs_correct,obs_var_u,&
                       &obs_inc_eakf_u,obs_inc_oi_u,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_u%lat < ass_ps_lat_start .or. obs_loc_u%lat> ass_ps_lat_end) obs_inc_eakf_u(:) = 0.0
                if (config_flags%debug_u) grid%obs_incre3h(i_o,k_o,j_o)=obs_inc_eakf_u(20)

!1.5 v

                        obs_value = atmos_obs%v(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_v(j_ens) = v_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_v_now(j_ens) = v_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value

                call  obs_increment_hfq(config_flags,enso_v,enso_v_now,ens_size,obs_correct,obs_var_v,&
                       &obs_inc_eakf_v,obs_inc_oi_v,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_v%lat < ass_ps_lat_start .or. obs_loc_v%lat> ass_ps_lat_end) obs_inc_eakf_v(:) = 0.0
                if (config_flags%debug_v) grid%obs_incre3h(i_o,k_o,j_o)=obs_inc_eakf_v(20)

! 2 
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.
                        if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0) i_mm = 1
                       if (i_m > npx) i_mm = npx

                           if (config_flags%hea_t3h) then
                            model_loc%lon = T_grid%x(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h

                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = t_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do
                        if (config_flags%debug_mda)  WRITE(*,*) 'ens_in_hfq,cov',sum(ens_in_hfq(:)),cov_factor

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                        vari='t'
                 if (config_flags%debug_mda) then 
                        WRITE(*,*) 'CALL update_from_obs_inc'
                 ENDIF


                                  call update_from_obs_inc(config_flags,enso_t,obs_inc_eakf_t,&
                                        &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              endif
                              t_now%data(i_m,k_o,j_m) = t_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                        endif !  t_tday 
!u
                         if (config_flags%hea_u3h) then !attention

                            model_loc%lon = T_grid%xu(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yu(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_u, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = u_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                          vari='u'
                                   call update_from_obs_inc(config_flags,enso_u,obs_inc_eakf_u,&
                                        &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_u,vari)
                              endif
                                u_now%data(i_m,k_o,j_m) =u_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_uday

!v
                         if (config_flags%hea_v3h) then !attention

                            model_loc%lon = T_grid%xv(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yv(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_v, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = v_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                          vari='v'
                                   call update_from_obs_inc(config_flags,enso_v,obs_inc_eakf_v,&
                                        &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_v,vari)
                              endif
                                v_now%data(i_m,k_o,j_m) =v_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_vday

                        endif !  (assim_flag)

                    enddo
                    enddo
                endif
        enddo !i_o
        enddo !k_o
        enddo !j_o


end  subroutine ensemble_filter_atm_3h

 subroutine ensemble_filter_atm_6h(grid,config_flags,t_ens_tau,t_now, q_ens_tau,u_ens_tau,u_now,v_ens_tau,v_now,&
      & ps_ens, ps_tau, atm_var,atmos_obs, isc, iec, jsc, jec, halo,npx,npy,nk,&
        T_grid,m_time,atm_step,isd, ied, jsd, jed)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: q_ens_tau, t_ens_tau
    type(ada_field_type), dimension(:), intent(inout) :: u_ens_tau, v_ens_tau
    type(ada_field_type), dimension(:), intent(in) :: ps_ens
    type(ada_field_type), intent(inout) :: ps_tau
    type(ada_field_type), intent(inout) :: t_now,u_now,v_now
    type(assim_data_type), intent(in) :: atm_var
    type(assim_data_type), intent(in) :: atmos_obs
    type(ada_grid_type), intent(in) :: T_grid
    type(WRFU_Time), intent(in) :: m_time
    integer, intent(in) :: isc, iec, jsc, jec, halo, npx, npy, nk
    integer, intent(in) :: isd, ied, jsd, jed
    type(loc_type) :: model_loc, obs_loc_ptq,obs_loc_u,obs_loc_v
    INTEGER             :: yr, month, day, hr, minute, sec, rc
    integer, intent(in)  :: atm_step
    integer :: i_o,j_o,k_o,i_m,j_m,k_m, j_ens, i_oo, i_mm,iii
    integer :: ens_size
    real :: cor_oi, coslat, rad2deg,coslatu,coslatv
    real :: dist, dist0, obs_value, obs_elev, mod_elev, obs_correct
    real :: std_oi_o, std_oi_g
    INTEGER myproc
    integer :: is_halo,js_halo,ie_halo,je_halo
    integer :: ass_variable = 1
    real :: cov_factor_v, cov_factor_t, cov_factor_h,sum_obsinc
    integer:: maxi,minj
    real :: sum_ens_in_hfq1,sum_ens_in_hfq2
    character(len=1) :: vari

        sum_ens_in_hfq1=0.
        sum_ens_in_hfq2=0.
      CALL WRFU_TimeGet( m_time, YY=yr, dayOfYear=day,H=hr,M=minute,S=sec,rc=rc)
      CALL wrf_get_myproc( myproc )
write(*,*) 'npx,  npy,  nk ',npx,npy,nk 
write(*,*) 'running ensemble_filter_atm_6h'
        ens_size = size(ps_ens)
  is_halo=isd
  ie_halo=ied
  js_halo=jsd
  je_halo=jed
!!!!!!!!!!!!!!!!!!!!!!!!!! 表层气压 (ps) 的同化  !!!!!!!!!!!!!!!!!!!!!!!!!
         do j_o =  js_halo,je_halo
         do i_o =  is_halo,ie_halo
                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx
                 if (i_o <= 0) i_oo  = 1
                 if (i_o > npx) i_oo = npx

                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        obs_value = atmos_obs%ps(i_oo,j_o)
!                          obs_elev = atmos_obs%phis(i_oo,j_o)
                           obs_elev = grid%ht(i_oo,j_o)
!          if (  obs_elev<=100.0 .and. obs_value <=96000.0) then
!                        else
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_ps_hfq(j_ens) = ps_ens(j_ens)%data(i_o,j_o,1)
                            enso_ps_tau(j_ens) = ps_tau%data(i_o,j_o,1)
                        end do
                        obs_correct = obs_value

                  call obs_increment_hfq(config_flags,enso_ps_hfq,enso_ps_tau,ens_size,obs_correct, obs_var_ps,&
                   &obs_inc_eakf_ps_hfq,obs_inc_oi_ps_hfq,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat>ass_ps_lat_end) obs_inc_eakf_ps_hfq(:) = 0.0
                        sum_obsinc=sum(obs_inc_eakf_ps_hfq(:))
                if (config_flags%debug_ps)    grid%obs_incre6h(i_o,1,j_o)=obs_inc_eakf_ps_hfq(20)

                 if (sum(obs_inc_eakf_ps_hfq(:)) /= 0.0) then
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.
                        if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0) i_mm = 1
                       if (i_m > npx) i_mm = npx

                           model_loc%lon = T_grid%x(i_mm)
                           if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                model_loc%lat = T_grid%y(j_m)
                                ass_variable = 1
                                ens_inc_hfq(:) = 0.0
                                dist = get_dist(model_loc, obs_loc_ptq, atm_step)

                                dist = RADIUS * sqrt(dist)
                                dist0 = ps_dist*coslat

                                cov_factor_h = comp_cov_factor(dist, dist0)
                                cor_oi = cov_factor_h
                                cov_factor = cov_factor_h
                                  do j_ens = 1, ens_size
                                    ens_in_hfq(j_ens) =  ps_ens(j_ens)%data(i_m,j_m,1)
                                 end do
                sum_ens_in_hfq1=sum(ens_in_hfq)
                sum_ens_in_hfq2=sum(ens_in_hfq(:))
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                            vari='p'
                                   call update_from_obs_inc(config_flags,enso_ps_hfq,obs_inc_eakf_ps_hfq,&
                                        & obs_inc_oi_ps_hfq,ens_in_hfq(:),ens_size,&
                                        & ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              end if

                              ps_tau%data(i_m,j_m,1) = ps_tau%data(i_m,j_m,1)+ens_inc_hfq(ens_size)
                       end if

                    enddo
                    enddo

                endif
                endif
!                endif

        end do
        end do

!!!!!!!!!!!!!!!!!!!!!! 三维场 (T, Q, U, V) 的同化 !!!!!!!!!!!!!!!!!!!!!!!!!
! process 3D obs now
! 1- 1.1  1. 三重循环遍历观测点        
         do j_o =  js_halo,je_halo   ! 水平方向 j
         do k_o = 1, nk              ! 垂直层 k
         do i_o =  is_halo,ie_halo   ! 水平方向 i
                !—————— 以下开始单点同化 ——————
                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.   ! j_o 越界处理
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx
                 if (i_o <= 0) i_oo = 1     ! i_o 越界处理
                 if (i_o > npx) i_oo = npx  ! i_o 越界处理
                        !在同化每个要素之前，需要先计算观测点在模型格点的地理位置与距离权重（经纬度余弦）
                        ! 观测点经度、纬度
                        obs_loc_u%lon   = T_grid%xu(i_oo)
                        obs_loc_v%lon   = T_grid%xv(i_oo)
                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_u%lon > 360.0 ) obs_loc_u%lon = obs_loc_u%lon-360.0
                        if ( obs_loc_v%lon > 360.0 ) obs_loc_v%lon = obs_loc_v%lon-360.0
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_u%lat = T_grid%yu(j_o)
                        obs_loc_v%lat = T_grid%yv(j_o)
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        ! 对纬度取绝对值后不超过 80°，再算余弦
                        ! 纬度余弦 用于后面协方差局部化时的“影响半径”缩放。
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        if ( abs(obs_loc_u%lat) <= 80.0 ) then
                         coslatu = cos(obs_loc_u%lat*DEG_TO_RAD)
                        else
                         coslatu = cos(80.0*DEG_TO_RAD)
                        end if

                        if ( abs(obs_loc_v%lat) <= 80.0 ) then
                         coslatv = cos(obs_loc_v%lat*DEG_TO_RAD)
                        else
                         coslatv = cos(80.0*DEG_TO_RAD)
                        endif

                       !!!!!!!!!! 1.2 t  温度 (T) 增量计算
                       ! ———— 提取观测值 ————
                        obs_value = atmos_obs%t(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        ! ———— 提取20个背景集合预报和背景场 ————
                        do j_ens = 1, ens_size
                            enso_t(j_ens) = t_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_t_now(j_ens) = t_now%data(i_o,k_o,j_o)
                        end do
                        obs_correct = obs_value
                ! ———— 计算观测增量（EAKF-HFQ 方法） ————
                call  obs_increment_hfq(config_flags,enso_t,enso_t_now,ens_size, obs_correct,obs_var_t,&
                       &obs_inc_eakf_t,obs_inc_oi_t,std_oi_o,i_o,j_o,atm_step)
                ! ———— 纬度范围外清零 & Debug 存储 ————
                if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat> ass_ps_lat_end) obs_inc_eakf_t(:) = 0.0
                if (config_flags%debug_t) grid%obs_incre6h(i_o,k_o,j_o)=obs_inc_eakf_t(20)
                ! enso_t：集合预报在观测点的阵列
                ! enso_t_now：上一次分析/背景场
                ! obs_increment_hfq：核心子例程，输出
                ! obs_inc_eakf_t(:)：每个集合成员的观测增量
                ! obs_inc_oi_t(:)：正交子空间增量

                        !!!!!!!!!!! 1.3 q
                        !!!!!!!!!!! 1.4 u  水平风 U 和 V 
                        ! ———— U 风观测增量 ————
                        obs_value = atmos_obs%u(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_u(j_ens) = u_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_u_now(j_ens) = u_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value

                call  obs_increment_hfq(config_flags,enso_u,enso_u_now,ens_size,obs_correct,obs_var_u,&
                       &obs_inc_eakf_u,obs_inc_oi_u,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_u%lat < ass_ps_lat_start .or. obs_loc_u%lat> ass_ps_lat_end) obs_inc_eakf_u(:) = 0.0
                if (config_flags%debug_u) grid%obs_incre6h(i_o,k_o,j_o)=obs_inc_eakf_u(20)

                        !!!!!!!!!!!! 1.5 v  ! ———— V 风观测增量 ————
                        obs_value = atmos_obs%v(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_v(j_ens) = v_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_v_now(j_ens) = v_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value
                call obs_increment_hfq(config_flags,enso_v,enso_v_now,ens_size,obs_correct,obs_var_v,&
                       &obs_inc_eakf_v,obs_inc_oi_v,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_v%lat < ass_ps_lat_start .or. obs_loc_v%lat> ass_ps_lat_end) obs_inc_eakf_v(:) = 0.0
                if (config_flags%debug_v) grid%obs_incre6h(i_o,k_o,j_o)=obs_inc_eakf_v(20)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2  邻域局部更新 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.
                        if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0)  i_mm = 1
                        if (i_m > npx) i_mm = npx
                             ! 本地化协方差
                           if (config_flags%hea_t6h) then
                            model_loc%lon = T_grid%x(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0

                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)  ! 水平给出距离衰减函数
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h  ! 距离系数，这里实际上只给了水平的系数
                            ! 提取邻点20个背景集合
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = t_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                        vari='t'

                                   call update_from_obs_inc(config_flags,enso_t,obs_inc_eakf_t,&
                                        &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              endif 
                               ! 将最后一个成员的增量累加到分析场
                              t_now%data(i_m,k_o,j_m) = t_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                        endif !  t_tday 

!u
                         if (config_flags%hea_u6h) then !attention

                            model_loc%lon = T_grid%xu(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yu(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_u, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = u_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0

                          vari='u'
                                   call update_from_obs_inc(config_flags,enso_u,obs_inc_eakf_u,&
                                        &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_u,vari)
                              endif
                                u_now%data(i_m,k_o,j_m) =u_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_uday
!v
                         if (config_flags%hea_v6h) then !attention

                            model_loc%lon = T_grid%xv(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yv(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_v, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h

                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = v_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                          vari='v'
                                   call update_from_obs_inc(config_flags,enso_v,obs_inc_eakf_v,&
                                        &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_v,vari)
                              endif

                                v_now%data(i_m,k_o,j_m) =v_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_vday


                        endif !  (assim_flag)

                    enddo
                    enddo

                endif
        enddo !i_o
        enddo !k_o
        enddo !j_o


end subroutine ensemble_filter_atm_6h






 subroutine ensemble_filter_atm_day(grid,config_flags,t_ens_tau, t_now,q_ens_tau, u_ens_tau,u_now,&
       & v_ens_tau, v_now,ps_ens, ps_tau, atm_var,atmos_obs, isc, iec, jsc, jec, halo, npx,npy,nk,&
        T_grid,m_time,atm_step,isd, ied, jsd, jed)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: q_ens_tau, t_ens_tau
    type(ada_field_type), dimension(:), intent(inout) :: u_ens_tau, v_ens_tau
    type(ada_field_type), dimension(:), intent(in) :: ps_ens
    type(ada_field_type), intent(inout) :: ps_tau
    type(ada_field_type), intent(inout) :: t_now,u_now,v_now
    type(assim_data_type), intent(in) :: atm_var
    type(assim_data_type), intent(in) :: atmos_obs
    type(ada_grid_type), intent(in) :: T_grid
    type(WRFU_Time), intent(in) :: m_time
    integer, intent(in) :: isc, iec, jsc, jec, halo, npx, npy, nk
    integer, intent(in) :: isd, ied, jsd, jed
    type(loc_type) :: model_loc, obs_loc_ptq,obs_loc_u,obs_loc_v
    INTEGER             :: yr, month, day, hr, minute, sec, rc
    integer, intent(in)  :: atm_step
    integer :: i_o,j_o,k_o,i_m,j_m,k_m, j_ens, i_oo, i_mm,iii
    integer :: ens_size
    real :: cor_oi, coslat, rad2deg,coslatu,coslatv
    real :: dist, dist0, obs_value, obs_elev, mod_elev, obs_correct
    real :: std_oi_o, std_oi_g
    INTEGER myproc
    integer :: is_halo,js_halo,ie_halo,je_halo
        integer :: ass_variable = 1
    real :: cov_factor_v, cov_factor_t, cov_factor_h,sum_obsinc
    integer:: maxi,minj
        real :: sum_ens_in_hfq1,sum_ens_in_hfq2
    character(len=1) :: vari


        sum_ens_in_hfq1=0.
        sum_ens_in_hfq2=0.

      CALL WRFU_TimeGet( m_time, YY=yr, dayOfYear=day,H=hr,M=minute,S=sec,rc=rc)
      CALL wrf_get_myproc( myproc )

!   write(*,*) 'z_level',nk
write(*,*) 'running ensemble_filter_atm_day'
        ens_size = size(ps_ens)

  is_halo=isd
  ie_halo=ied
  js_halo=jsd
  je_halo=jed

!  is_halo=max(isd,2)
!  ie_halo=min(ied,558)
!  js_halo=max(jsd,2)
!  je_halo=min(jed,408)


         do j_o =  js_halo,je_halo
         do i_o =  is_halo,ie_halo

                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx


                 if (i_o <= 0)  i_oo = 1
                 if (i_o > npx) i_oo = npx

                        obs_loc_u%lon   = T_grid%xu(i_oo)
                        obs_loc_v%lon   = T_grid%xv(i_oo)
                        obs_loc_ptq%lon = T_grid%x(i_oo)

                        if ( obs_loc_u%lon > 360.0 ) obs_loc_u%lon = obs_loc_u%lon-360.0
                        if ( obs_loc_v%lon > 360.0 ) obs_loc_v%lon = obs_loc_v%lon-360.0
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0

                        obs_loc_u%lat = T_grid%yu(j_o)
                        obs_loc_v%lat = T_grid%yv(j_o)
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if

                        if ( abs(obs_loc_u%lat) <= 80.0 ) then
                         coslatu = cos(obs_loc_u%lat*DEG_TO_RAD)
                        else
                         coslatu = cos(80.0*DEG_TO_RAD)
                        end if

                        if ( abs(obs_loc_v%lat) <= 80.0 ) then
                         coslatv = cos(obs_loc_v%lat*DEG_TO_RAD)
                        else
                         coslatv = cos(80.0*DEG_TO_RAD)
                        end if

                        obs_value = atmos_obs%ps(i_oo,j_o)
!                          obs_elev = atmos_obs%phis(i_oo,j_o)
                           obs_elev = grid%ht(i_oo,j_o) !HGT 
!                 if (  obs_elev<=100.0 .and. obs_value <=96000.0) then
!                write(*,*) 'obs_elevless100',j_o,i_o,obs_elev,obs_value

!                        else
!write(*,*) 'j_o, i_o,obs_loc_ptq%lat,obs_loc_ptq%lon, XLAT_U,XLONG_U,XLAT_V,XLONG_V', &
!          & j_o,i_o,obs_loc_ptq%lat,obs_loc_ptq%lon,grid%XLAT_U(i_o,j_o),grid%XLONG_U(i_o,j_o),grid%XLAT_V(i_o,j_o),grid%XLONG_V(i_o,j_o)

!write(*,*) '-----------------------------------------------------------------'
!write(*,*) ' i_o,j_o      ',i_o,j_o
!write(*,*) 'halo',is_halo,ie_halo,js_halo,je_halo,'calc',isc, iec, jsc, jec
!write(*,*) '-----------------------------------------------------------------'
                        std_oi_o = 0.0

                        do j_ens = 1, ens_size
                            enso_ps_hfq(j_ens) = ps_ens(j_ens)%data(i_o,j_o,1)
                            enso_ps_tau(j_ens) = ps_tau%data(i_o,j_o,1) 
                        end do

                        obs_correct = obs_value

!               write(*,'("ps_ens,ps_tau",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xf16.2,1xf16.2,1xf16.2)')&
!           atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,ps_ens(1)%data(i_o,j_o,1),&
!                ps_ens(20)%data(i_o,j_o,1), ps_tau%data(i_o,j_o,1)

!               write(*,'("enso_ps",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xf16.2,1xf16.2,1xf16.2)')&
!           atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,enso_ps_hfq(1),enso_ps_hfq(20),enso_ps_tau(20)

                        call obs_increment_hfq(config_flags,enso_ps_hfq,enso_ps_tau,ens_size, obs_correct, obs_var_ps,&
                       & obs_inc_eakf_ps_hfq,obs_inc_oi_ps_hfq,std_oi_o,i_o,j_o,atm_step)
                 if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat> ass_ps_lat_end) obs_inc_eakf_ps_hfq(:) = 0.0
                        sum_obsinc=sum(obs_inc_eakf_ps_hfq(:))

                if (config_flags%debug_ps)     grid%obs_incre1day(i_o,1,j_o)=obs_inc_eakf_ps_hfq(20)

                 if (sum(obs_inc_eakf_ps_hfq(:)) /= 0.0) then
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.

                        if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0)  i_mm = 1
                        if (i_m > npx) i_mm = npx

                           model_loc%lon = T_grid%x(i_mm)
                           if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                model_loc%lat = T_grid%y(j_m)
                                ass_variable = 1
                                ens_inc_hfq(:) = 0.0

!write(*,*) '-----------------------------------------------------------------'
!write(*,*) ' i_m,j_m      ',i_m,j_m,' i_mm,j_m     ',i_mm,j_m
!write(*,*) 'ii_do,jj_do',i_o-ii_do, i_o+ii_do,j_o-jj_do,j_o+jj_do
!write(*,*) '-----------------------------------------------------------------'

                                dist = get_dist(model_loc, obs_loc_ptq, atm_step)

!             write(*,'("loc",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf6.2,1xf6.2,1xf16.2)')&
!            atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,i_mm,j_m,model_loc%lon,model_loc%lat,dist

                                dist = RADIUS * sqrt(dist)
                                dist0 = ps_dist*coslat

!              write(*,'("dist0,coslat",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2,1xf16.2)')&
!                         atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,i_mm,j_m,dist0,coslat

                                cov_factor_h = comp_cov_factor(dist, dist0)
                                cor_oi = cov_factor_h
                                cov_factor = cov_factor_h
                                  do j_ens = 1, ens_size
                                    ens_in_hfq(j_ens) = ps_ens(j_ens)%data(i_m,j_m,1)
                                 end do
!         write(*,*) 'im,jm,ens_in_hfq',atm_step,j_m,i_m,ens_in_hfq(20)
                sum_ens_in_hfq1=sum(ens_in_hfq)
                sum_ens_in_hfq2=sum(ens_in_hfq(:))

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor /=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
        if (config_flags%debug_ps) then
!       write(*,'("cov_factor,ens_in_hfqin",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2,1xf16.2)')&
!                     atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,i_m,j_m,cov_factor,ens_in_hfq(20)
        endif
                vari='p'
                                   call update_from_obs_inc(config_flags,enso_ps_hfq,obs_inc_eakf_ps_hfq,&
                                        & obs_inc_oi_ps_hfq, ens_in_hfq(:),ens_size,&
                                        & ens_inc_hfq, cov_factor, cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step, &
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              end if
!                   write(*,'("ens_inc_hfq2",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2)')&
!                         atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,i_m,j_m,ens_inc_hfq(20)

!                    write(*,'("ps_tau1",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2)')&
!                         atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,i_m,j_m,ps_tau%data(i_m,j_m,1)
                              ps_tau%data(i_m,j_m,1) = ps_tau%data(i_m,j_m,1)+ens_inc_hfq(ens_size)

!                    write(*,'("ps_tau2",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2)')&
!                         atm_step,myproc,i_o,j_o,obs_loc_ptq%lon,obs_loc_ptq%lat,i_m,j_m,ps_tau%data(i_m,j_m,1)

                       end if

                    enddo
                    enddo


                endif
                endif
!                endif


        end do
        end do

! process 3D obs now
! 1- 1.1 
         do j_o =  js_halo,je_halo
         do k_o = 1, nk
         do i_o =  is_halo,ie_halo

                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx

                 if (i_o <= 0) i_oo  = 1
                 if (i_o > npx) i_oo = npx

                        obs_loc_u%lon   = T_grid%xu(i_oo)
                        obs_loc_v%lon   = T_grid%xv(i_oo)
                        obs_loc_ptq%lon = T_grid%x(i_oo)

                        if ( obs_loc_u%lon > 360.0 ) obs_loc_u%lon = obs_loc_u%lon-360.0
                        if ( obs_loc_v%lon > 360.0 ) obs_loc_v%lon = obs_loc_v%lon-360.0
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_u%lat = T_grid%yu(j_o)
                        obs_loc_v%lat = T_grid%yv(j_o)
                        obs_loc_ptq%lat = T_grid%y(j_o)

                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        if ( abs(obs_loc_u%lat) <= 80.0 ) then
                         coslatu = cos(obs_loc_u%lat*DEG_TO_RAD)
                        else
                         coslatu = cos(80.0*DEG_TO_RAD)
                        end if

                        if ( abs(obs_loc_v%lat) <= 80.0 ) then
                         coslatv = cos(obs_loc_v%lat*DEG_TO_RAD)
                        else
                         coslatv = cos(80.0*DEG_TO_RAD)
                        endif
! 1.2 t
                        obs_value = atmos_obs%t(i_oo,k_o,j_o)
                        std_oi_o = 0.0 
                        do j_ens = 1, ens_size
                            enso_t(j_ens) = t_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_t_now(j_ens) = t_now%data(i_o,k_o,j_o)
                        end do
                        obs_correct = obs_value
!                        call obs_increment(enso_t, ens_size, t_now,obs_var_t,&
!                        & obs_inc_eakf_t, obs_inc_oi_t, std_oi_o)
                call  obs_increment_hfq(config_flags,enso_t,enso_t_now, ens_size, obs_correct,obs_var_t,&
                       &obs_inc_eakf_t,obs_inc_oi_t,std_oi_o,i_o,j_o,atm_step)
!                if (obs_loc%lat < ass_t_lat_start .and. k_o > 16) obs_inc_eakf_t(:) = 0.0
                if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat> ass_ps_lat_end) obs_inc_eakf_t(:) = 0.0
                if (config_flags%debug_t)     grid%obs_incre1day(i_o,k_o,j_o)=obs_inc_eakf_t(20)

! 1.3 q
! 1.4 u
                        obs_value = atmos_obs%u(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_u(j_ens) = u_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_u_now(j_ens) = u_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value
                call  obs_increment_hfq(config_flags,enso_u,enso_u_now,ens_size, obs_correct,obs_var_u,&
                       &obs_inc_eakf_u,obs_inc_oi_u,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_u%lat < ass_ps_lat_start .or. obs_loc_u%lat> ass_ps_lat_end) obs_inc_eakf_u(:) = 0.0
                if (config_flags%debug_u) grid%obs_incre1day(i_o,k_o,j_o)=obs_inc_eakf_u(20)

! 1.5 v

                        obs_value = atmos_obs%v(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_v(j_ens) = v_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_v_now(j_ens) = v_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value
                call  obs_increment_hfq(config_flags,enso_v,enso_v_now,ens_size,obs_correct,obs_var_v,&
                       & obs_inc_eakf_v,obs_inc_oi_v,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_v%lat < ass_ps_lat_start .or. obs_loc_v%lat>ass_ps_lat_end) obs_inc_eakf_v(:) = 0.0
                if (config_flags%debug_v) grid%obs_incre1day(i_o,k_o,j_o)=obs_inc_eakf_v(20)


! 2 
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                          if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                         if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                          if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                          if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.        

                      if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0) i_mm = 1
                        if (i_m > npx) i_mm = npx

                         if (config_flags%hea_tday) then

                            model_loc%lon = T_grid%x(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j_m) 
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0

                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
!                write(*,*) 'dist,dist0,cov_factor',i_m,k_o,j_m,dist,dist0,cov_factor
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = t_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                        vari='t'
                                   call update_from_obs_inc(config_flags,enso_t,obs_inc_eakf_t,&
                                        &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              endif 
                              t_now%data(i_m,k_o,j_m) = t_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                        endif !  config_flags%hea_tday 
! u
                         if (config_flags%hea_uday) then

                            model_loc%lon = T_grid%xu(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yu(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_u, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = u_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0

                          vari='u'
                                   call update_from_obs_inc(config_flags,enso_u,obs_inc_eakf_u,&
                                        &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_u,vari)
                              endif
                                u_now%data(i_m,k_o,j_m) =u_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size) 
                            endif ! config_flags%hea_uday
! v
                         if (config_flags%hea_vday) then

                            model_loc%lon = T_grid%xu(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yv(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_v, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = v_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0

                          vari='v'
                                   call update_from_obs_inc(config_flags,enso_v,obs_inc_eakf_v,&
                                        &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_v,vari)
                              endif

                                v_now%data(i_m,k_o,j_m) =v_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_vday



                        endif !  (assim_flag)
                    enddo
                    enddo

!t_ens_tau, t_now,                
                endif !if (obs_flag)
        enddo !i_o
        enddo !k_o
        enddo !j_o

end subroutine ensemble_filter_atm_day


 subroutine ensemble_filter_atm_hfq(grid,config_flags,t_ens_tau, q_ens_tau, u_ens_tau, v_ens_tau,&
      & ps_ens, ps_tau, atm_var,&
      & atmos_obs, isc, iec, jsc, jec, halo, npx, npy, nk,T_grid,m_time,atm_step, &
      &  isd, ied, jsd, jed)

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: q_ens_tau, t_ens_tau
    type(ada_field_type), dimension(:), intent(inout) :: u_ens_tau, v_ens_tau
    type(ada_field_type), dimension(:), intent(in) :: ps_ens 
    type(ada_field_type), intent(inout) :: ps_tau 
    type(assim_data_type), intent(in) :: atm_var 
    type(assim_data_type), intent(in) :: atmos_obs
    type(ada_grid_type), intent(in) :: T_grid
    type(WRFU_Time), intent(in) :: m_time
    integer, intent(in) :: isc, iec, jsc, jec, halo, npx, npy, nk
    integer, intent(in) :: isd, ied, jsd, jed
    type(loc_type) :: model_loc, obs_loc_ptq
    INTEGER             :: yr, month, day, hr, minute, sec, rc
    integer, intent(in)  :: atm_step
    integer :: i_o,j_o,k_o,i_m,j_m,k_m, j_ens, i_oo, i_mm,iii
    integer :: ens_size
    real :: cor_oi, coslat, rad2deg
    real :: dist, dist0, obs_value, obs_elev, mod_elev, obs_correct 
    real :: std_oi_o, std_oi_g
!    real, dimension(:), intent(in) :: ak, bk !lulv 20200731
    INTEGER myproc
    integer :: is_halo,js_halo,ie_halo,je_halo
        integer :: ass_variable = 1
    real :: cov_factor_v, cov_factor_t, cov_factor_h,sum_obsinc
    integer:: maxi,minj
        real :: sum_ens_in_hfq1,sum_ens_in_hfq2
        character(len=1) :: vari
        sum_ens_in_hfq1=0.
        sum_ens_in_hfq2=0.

      CALL WRFU_TimeGet( m_time, YY=yr, dayOfYear=day, H=hr, M=minute,S=sec,rc=rc)
      CALL wrf_get_myproc( myproc ) 

        ens_size = size(ps_ens)
write(*,*) 'running ensemble_filter_atm_hfq'




!  is_halo=isc-halo
!  ie_halo=iec+halo
!  js_halo=jsc-halo
!  je_halo=jec+halo

! 1214 
  is_halo=isd
  ie_halo=ied
  js_halo=jsd
  je_halo=jed

!  is_halo=isc
!  ie_halo=iec
!  js_halo=jsc
!  je_halo=jec

         do j_o =  js_halo,je_halo
         do i_o =  is_halo,ie_halo
 
 
                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points

                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0)  i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx

                 if (i_o <= 0)  i_oo = i_o + 6
                 if (i_o > npx) i_oo = i_o - 6
                        
                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if

                
                           obs_value = atmos_obs%ps(i_oo,j_o)
!                          obs_elev = atmos_obs%phis(i_oo,j_o)
                           obs_elev = grid%ht(i_oo,j_o)
!                          mod_elev = atm_var%phis(i_oo,j_o)
!                          tm_lev = atm_var%t(i_oo,j_o,:)
        
  write(*,*) 'obs_elev',atm_step,i_oo,j_o,obs_elev
! if (obs_elev<=100.0 .and. obs_value <=96000.0) then 
!  if (  obs_elev<=100.0 .and. obs_value <=96000.0) then 

!                        else


                        std_oi_o = 0.0

                        do j_ens = 1, ens_size
                            enso_ps_hfq(j_ens) = ps_ens(j_ens)%data(i_o,j_o,1)
                            enso_ps_tau(j_ens) = ps_tau%data(i_o,j_o,1) !lulv 20200825
                        end do
                      
                        obs_correct = obs_value
!write(*,*) 'obs_var_ps',obs_var_ps



                        call obs_increment_hfq(config_flags,enso_ps_hfq,enso_ps_tau, ens_size, obs_correct, obs_var_ps,&
                       & obs_inc_eakf_ps_hfq, obs_inc_oi_ps_hfq,std_oi_o,i_o,j_o,atm_step) 

                       
!!                          if (obs_loc_ptq%lat < ass_ps_lat_start) obs_inc_eakf_ps_hfq(:) = 0.0
                 if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat> ass_ps_lat_end) obs_inc_eakf_ps_hfq(:) = 0.0

!isc, iec, jsc, jec
                        sum_obsinc=sum(obs_inc_eakf_ps_hfq(:))        

                 if (sum(obs_inc_eakf_ps_hfq(:)) /= 0.0) then
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do 
                        assim_flag = .true.

!                        maxi = max (is_halo, 0)
!                        minj = min (ie_halo,559)

                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false. 

                        if (assim_flag) then

                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
                        if (i_m <= 0) i_mm = i_m + npx
                        if (i_m > npx) i_mm = i_m - npx
                           model_loc%lon = T_grid%x(i_mm)
                           if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                model_loc%lat = T_grid%y(j_m)
                                ass_variable = 1
                                ens_inc_hfq(:) = 0.0
                                dist = get_dist(model_loc, obs_loc_ptq, atm_step)

                                dist = RADIUS * sqrt(dist)
                                dist0 = ps_dist*coslat


                                cov_factor_h = comp_cov_factor(dist, dist0)
                                cor_oi = cov_factor_h
                                cov_factor = cov_factor_h

!                                if ( ass_ps ) then ! take from nml todecide to do or not for ps

                                  do j_ens = 1, ens_size
                                    ens_in_hfq(j_ens) = ps_ens(j_ens)%data(i_m,j_m,1)
                                 end do


                sum_ens_in_hfq1=sum(ens_in_hfq)
                sum_ens_in_hfq2=sum(ens_in_hfq(:))



                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor /= 0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0

                                vari='p'
                                   call update_from_obs_inc(config_flags,enso_ps_hfq, obs_inc_eakf_ps_hfq,&
                                        & obs_inc_oi_ps_hfq, ens_in_hfq(:), ens_size,&
                                        & ens_inc_hfq, cov_factor, cor_oi, std_oi_o,std_oi_g,&
                                        & ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step, &
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              end if



                              ps_tau%data(i_m,j_m,1) = ps_tau%data(i_m,j_m,1)+ens_inc_hfq(ens_size)


                       end if

 !                       endif !assim_flag


                    enddo
                    enddo



                endif
                endif
!                endif


        end do
        end do



!deallocate(enso_ps_tau)
 end subroutine ensemble_filter_atm_hfq


 subroutine ensemble_filter_atm_5day(grid,config_flags,t_ens_tau,t_now,q_ens_tau,u_ens_tau,u_now,v_ens_tau,v_now,&
      & ps_ens, ps_tau, atm_var,atmos_obs, isc, iec, jsc, jec, halo,npx,npy,nk,&
        T_grid,m_time,atm_step,isd, ied, jsd, jed)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: q_ens_tau, t_ens_tau
    type(ada_field_type), dimension(:), intent(inout) :: u_ens_tau, v_ens_tau
    type(ada_field_type), dimension(:), intent(in) :: ps_ens
    type(ada_field_type), intent(inout) :: ps_tau
    type(ada_field_type), intent(inout) :: t_now,u_now,v_now
    type(assim_data_type), intent(in) :: atm_var
    type(assim_data_type), intent(in) :: atmos_obs
    type(ada_grid_type), intent(in) :: T_grid
    type(WRFU_Time), intent(in) :: m_time
    integer, intent(in) :: isc, iec, jsc, jec, halo, npx, npy, nk
    integer, intent(in) :: isd, ied, jsd, jed
    type(loc_type) :: model_loc, obs_loc_ptq,obs_loc_u,obs_loc_v
    INTEGER             :: yr, month, day, hr, minute, sec, rc
    integer, intent(in)  :: atm_step
    integer :: i_o,j_o,k_o,i_m,j_m,k_m, j_ens, i_oo, i_mm,iii
    integer :: ens_size
    real :: cor_oi, coslat, rad2deg,coslatu,coslatv
    real :: dist, dist0, obs_value, obs_elev, mod_elev, obs_correct
    real :: std_oi_o, std_oi_g
    INTEGER myproc
    integer :: is_halo,js_halo,ie_halo,je_halo
    integer :: ass_variable = 1
    real :: cov_factor_v, cov_factor_t, cov_factor_h,sum_obsinc
    integer:: maxi,minj
        real :: sum_ens_in_hfq1,sum_ens_in_hfq2
    character(len=1) :: vari

        sum_ens_in_hfq1=0.
        sum_ens_in_hfq2=0.
      CALL WRFU_TimeGet( m_time, YY=yr, dayOfYear=day,H=hr,M=minute,S=sec,rc=rc)
      CALL wrf_get_myproc( myproc )
  write(*,*) 'running ensemble_filter_atm_5day'
        ens_size = size(ps_ens)
  is_halo=isd
  ie_halo=ied
  js_halo=jsd
  je_halo=jed

         do j_o =  js_halo,je_halo
         do i_o =  is_halo,ie_halo
                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx


                 if (i_o <= 0)  i_oo = 1
                 if (i_o > npx) i_oo = npx

                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        obs_value = atmos_obs%ps(i_oo,j_o)
!                          obs_elev = atmos_obs%phis(i_oo,j_o)
                           obs_elev = grid%ht(i_oo,j_o)
!          if (  obs_elev<=100.0 .and. obs_value <=96000.0) then  ! not da sea low pressure
!                        else
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_ps_hfq(j_ens) = ps_ens(j_ens)%data(i_o,j_o,1)
                            enso_ps_tau(j_ens) = ps_tau%data(i_o,j_o,1)
                        end do
                        obs_correct = obs_value

            call  obs_increment_hfq(config_flags,enso_ps_hfq,enso_ps_tau,ens_size,obs_correct, obs_var_ps,&
              &obs_inc_eakf_ps_hfq,obs_inc_oi_ps_hfq,std_oi_o,i_o,j_o,atm_step)
         if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat>ass_ps_lat_end) obs_inc_eakf_ps_hfq(:) = 0.0
                        sum_obsinc=sum(obs_inc_eakf_ps_hfq(:))

                if (config_flags%debug_ps)       grid%obs_incre5day(i_o,1,j_o)=obs_inc_eakf_ps_hfq(20)
                 if (sum(obs_inc_eakf_ps_hfq(:)) /= 0.0) then
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
!                        if (i_m < maxi .or. i_m > minj ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.
                        if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                        if (i_m <= 0)  i_mm =1
                        if (i_m > npx) i_mm = npx

                           model_loc%lon = T_grid%x(i_mm)
                           if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                model_loc%lat = T_grid%y(j_m)
                                ass_variable = 1
                                ens_inc_hfq(:) = 0.0
                                dist = get_dist(model_loc, obs_loc_ptq, atm_step)

                                dist = RADIUS * sqrt(dist)
                                dist0 = ps_dist*coslat
                                cov_factor_h = comp_cov_factor(dist, dist0)
                                cor_oi = cov_factor_h
                                cov_factor = cov_factor_h
                                  do j_ens = 1, ens_size
                                    ens_in_hfq(j_ens) =  ps_ens(j_ens)%data(i_m,j_m,1)
                                 end do
                sum_ens_in_hfq1=sum(ens_in_hfq)
                sum_ens_in_hfq2=sum(ens_in_hfq(:))
                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                                vari='p'
                                   call update_from_obs_inc(config_flags,enso_ps_hfq,obs_inc_eakf_ps_hfq,&
                                        & obs_inc_oi_ps_hfq,ens_in_hfq(:),ens_size,&
                                        & ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        & ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              end if

                              ps_tau%data(i_m,j_m,1) = ps_tau%data(i_m,j_m,1)+ens_inc_hfq(ens_size)
                       end if

                    enddo
                    enddo

                endif
                endif
!                endif ! if (  obs_elev<=100.0 .and. obs_value <=96000.0) 

        end do
        end do

! process 3D obs now
! 1- 1.1

         do j_o =  js_halo,je_halo
         do k_o = 1, nk
         do i_o =  is_halo,ie_halo
                obs_flag = .true.
                if (j_o <= 0 .or. j_o > npy) obs_flag = .false.
                if (obs_flag) then ! only use for non-corner points
                 if (i_o >= 1 .and. i_o <= npx) i_oo = i_o
!                 if (i_o <= 0) i_oo = i_o + npx
!                 if (i_o > npx) i_oo = i_o - npx

                 if (i_o <= 0)  i_oo = 1
                 if (i_o > npx) i_oo = npx

                        obs_loc_u%lon   = T_grid%xu(i_oo)
                        obs_loc_v%lon   = T_grid%xv(i_oo)
                        obs_loc_ptq%lon = T_grid%x(i_oo)
                        if ( obs_loc_u%lon > 360.0 ) obs_loc_u%lon = obs_loc_u%lon-360.0
                        if ( obs_loc_v%lon > 360.0 ) obs_loc_v%lon = obs_loc_v%lon-360.0
                        if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon = obs_loc_ptq%lon-360.0
                        obs_loc_ptq%lat = T_grid%y(j_o)
                        obs_loc_u%lat = T_grid%yu(j_o)
                        obs_loc_v%lat = T_grid%yv(j_o)
                        if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
                         coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
                        else
                         coslat = cos(80.0*DEG_TO_RAD)
                        end if
                        if ( abs(obs_loc_u%lat) <= 80.0 ) then
                         coslatu = cos(obs_loc_u%lat*DEG_TO_RAD)
                        else
                         coslatu = cos(80.0*DEG_TO_RAD)
                        end if

                        if ( abs(obs_loc_v%lat) <= 80.0 ) then
                         coslatv = cos(obs_loc_v%lat*DEG_TO_RAD)
                        else
                         coslatv = cos(80.0*DEG_TO_RAD)
                        endif
! 1.2 t
                        obs_value = atmos_obs%t(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_t(j_ens) = t_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_t_now(j_ens) = t_now%data(i_o,k_o,j_o)
                        end do
                        obs_correct = obs_value
                call  obs_increment_hfq(config_flags,enso_t,enso_t_now,ens_size, obs_correct,obs_var_t,&
                       &obs_inc_eakf_t,obs_inc_oi_t,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_ptq%lat < ass_ps_lat_start .or. obs_loc_ptq%lat> ass_ps_lat_end) obs_inc_eakf_t(:) = 0.0
                if (config_flags%debug_t) grid%obs_incre5day(i_o,k_o,j_o)=obs_inc_eakf_t(20)
! 1.3 q
! 1.4 u
                        obs_value = atmos_obs%u(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_u(j_ens) = u_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_u_now(j_ens) = u_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value
                call obs_increment_hfq(config_flags,enso_u,enso_u_now,ens_size,obs_correct,obs_var_u,&
                       &obs_inc_eakf_u,obs_inc_oi_u,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_u%lat < ass_ps_lat_start .or. obs_loc_u%lat>ass_ps_lat_end) obs_inc_eakf_u(:) = 0.0
                if (config_flags%debug_u) grid%obs_incre5day(i_o,k_o,j_o)=obs_inc_eakf_u(20)

! 1.4 v
                        obs_value = atmos_obs%v(i_oo,k_o,j_o)
                        std_oi_o = 0.0
                        do j_ens = 1, ens_size
                            enso_v(j_ens) = v_ens_tau(j_ens)%data(i_o,k_o,j_o)
                            enso_v_now(j_ens) = v_now%data(i_o,k_o,j_o)
                        end do
                       obs_correct = obs_value
                call obs_increment_hfq(config_flags,enso_v,enso_v_now,ens_size,obs_correct,obs_var_v,&
                       &obs_inc_eakf_v,obs_inc_oi_v,std_oi_o,i_o,j_o,atm_step)
                if (obs_loc_v%lat < ass_ps_lat_start .or. obs_loc_v%lat>ass_ps_lat_end) obs_inc_eakf_v(:) = 0.0
                if (config_flags%debug_v) grid%obs_incre5day(i_o,k_o,j_o)=obs_inc_eakf_v(20)



! 2 
                    do j_m = j_o - jj_do, j_o + jj_do
                    do i_m = i_o - ii_do, i_o + ii_do
                        assim_flag = .true.
                        if (i_m < is_halo .or. i_m > ie_halo ) assim_flag = .false.
                        if (j_m < js_halo .or. j_m > je_halo ) assim_flag = .false.
                        if (j_m <= 0 .or. j_m > npy-1) assim_flag = .false.

                  if (assim_flag) then
                        if (i_m > 1 .and. i_m <= npx) i_mm = i_m
!                        if (i_m <= 0) i_mm = i_m + npx
!                        if (i_m > npx) i_mm = i_m - npx

                       if (i_m <= 0)  i_mm = 1
                       if (i_m > npx) i_mm = npx

                           if (config_flags%hea_t5day) then

                            model_loc%lon = T_grid%x(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
!                write(*,*)
!                'dist,dist0,cov_factor',i_m,k_o,j_m,dist,dist0,cov_factor
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = t_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and.cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                        vari='t'
                                   call update_from_obs_inc(config_flags,enso_t,obs_inc_eakf_t,&
                                        &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_ptq,vari)
                              endif

                              t_now%data(i_m,k_o,j_m) = t_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                        endif !  t_tday 

!u
                         if (config_flags%hea_u5day) then !attention

                            model_loc%lon = T_grid%xu(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yu(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_u, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = u_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                          vari='u'
                                   call update_from_obs_inc(config_flags,enso_u,obs_inc_eakf_u,&
                                        &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_u,vari)
                              endif

                                u_now%data(i_m,k_o,j_m) =u_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_u5day

                         if (config_flags%hea_v5day) then !attention

                            model_loc%lon = T_grid%xv(i_mm)
                            if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                            model_loc%lat = T_grid%yv(j_m)
                            ass_variable = 1
                            ens_inc_hfq(:) = 0.0
                            dist = get_dist(model_loc, obs_loc_v, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = ps_dist*coslat
                            cov_factor_h = comp_cov_factor(dist, dist0)
                            cor_oi = cov_factor_h
                            cov_factor = cov_factor_h
                              do j_ens = 1, ens_size
                                 ens_in_hfq(j_ens) = v_ens_tau(j_ens)%data(i_m,k_o,j_m)
                              end do

                              if ( sum(ens_in_hfq(:)) /= 0.0 .and. cov_factor/=0.0 ) then
                                  std_oi_g = 0.0
                                  std_oi_o = 0.0
                          vari='v'
                                   call update_from_obs_inc(config_flags,enso_v,obs_inc_eakf_v,&
                                        &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                        &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                        &ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
                                        & isc,iec,jsc,jec,obs_loc_v,vari)
                              endif
                                v_now%data(i_m,k_o,j_m) = v_now%data(i_m,k_o,j_m) + ens_inc_hfq(ens_size)
                            endif ! config_flags%hea_u5day

                        endif !  (assim_flag)

                    enddo
                    enddo
!t_ens_tau, t_now,                
                endif
        enddo !i_o
        enddo !k_o
        enddo !j_o

 end  subroutine ensemble_filter_atm_5day


 subroutine  ensemble_filter_atm_6h_r(grid,config_flags,ens_6h_q0,ens_6h_q,&
                   &ens_6h_t0,ens_6h_t,ens_6h_u0,ens_6h_u,&
                   &ens_6h_v0,ens_6h_v,ens_6h_ps0,ens_6h_ps,&
                  !  &atm_var, rio,rjo,rko,errf_u,errf_v,errf_t,errf_q,lat1,lon1,&
                  &atm_var,obs_z_number,lat_o,lon_o, rio,rjo,rko,errf_u,errf_v,errf_t,errf_q,&
                   &  T_grid, m_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts, jte,kts,kte,&
                     ids,  ide,  jds,  jde,  kds,  kde,dx)


  TYPE ( domain ), INTENT(INOUT) :: grid
  TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
  type(ada_field_type), dimension(:), intent(inout) :: ens_6h_q0,ens_6h_t0,ens_6h_u0,ens_6h_v0,ens_6h_ps0
  type(ada_field_type), intent(inout) :: ens_6h_q,ens_6h_t,ens_6h_u,ens_6h_v,ens_6h_ps
  type(assim_data_type), intent(in) :: atm_var
  integer, intent(in) :: isd, ied, jsd, jed,kms, kme,its, ite, jts,jte, kts, kte
  integer, intent(in) :: obs_z_number
  real,intent(in) :: rio(obs_z_number),rjo(obs_z_number)
  real,intent(in) :: rko(obs_z_number),errf_u(obs_z_number),errf_v(obs_z_number),errf_t(obs_z_number),errf_q(obs_z_number)
  real,intent(in) :: lat_o(obs_z_number),lon_o(obs_z_number)
  integer, intent(in) :: ids,  ide,  jds,  jde,  kds,  kde
  
  type(ada_grid_type), intent(in) :: T_grid
  integer, intent(in)  :: atm_step
  type(WRFU_Time), intent(in) :: m_time
  real, INTENT(IN)  :: dx
  integer :: n,nnn,k
  integer :: nsndlev,njcsnd
  integer nsta,I,J,maxj,minj,maxi,mini,igrid
  real  ::  ensem_stn_t(20,obs_z_number),ensem_stn_u(20,obs_z_number),ensem_stn_v(20,obs_z_number)
  real :: gfactor,rfactor,gridx,gridy,rindx,ris
  real :: rj,rx,ri,ry,cov_factor_h,cov_factor,cor_oi,cov_factor_v,cov_factor_t
  integer :: kob,kobp,kmin,kmax,ens_size,j_ens

  type(loc_type) :: model_loc,obs_loc_ptq,obs_loc_u,obs_loc_v
  real :: dist, dist0,coslat
  character(len=1) :: VARI
  real :: std_oi_o, std_oi_g
  integer :: ass_variable = 1
  real :: ra(obs_z_number),rb(obs_z_number),rc(obs_z_number)

  logical, allocatable :: assimilated_t(:,:,:)
  logical, allocatable :: assimilated_u(:,:,:)
  logical, allocatable :: assimilated_v(:,:,:)

ALLOCATE(assimilated_t(its:ite, kms:kme, jts:jte))
ALLOCATE(assimilated_u(its:ite, kms:kme, jts:jte))
ALLOCATE(assimilated_v(its:ite, kms:kme, jts:jte))


assimilated_t = .false.
assimilated_u = .true.
assimilated_v = .true.


  nsta = obs_z_number

  std_oi_o= 0.0
  std_oi_g=0.0

  ens_size = 20
  ensem_stn_t(:,:) = 0.
  ensem_stn_u(:,:) = 0.
  ensem_stn_v(:,:) = 0.
  write(*,*) 'running ensemble_filter_atm_6h_r'

  call ensem_stnt(grid,config_flags,ens_6h_t0,ensem_stn_t,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_6h_u0,ensem_stn_u,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_6h_v0,ensem_stn_v,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )

  GRIDX=0.5
  GRIDY=0.5
  IGRID=1

  write(*,*) "dx = " , dx
  RINDX=config_flags%e_r_6h *1000./dx   ! 影响半径对应的格点数

  DO J=jts,jte   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  水平方向 j 循环
        ! return if too close to j boundary
        if((j.le.2.or.j.ge.jde-1)) then
          ! return
          CYCLE
        endif
        ! MAXJ=J+IFIX(RINDX)  !J方向的上界（影响半径内）
        ! MINJ=J-IFIX(RINDX)  !J方向的下界（影响半径内）
        do n=1,nsta
           RA(n)=rio(n)-GRIDX  !i方向小半个网格
           RB(n)=rjo(n)-GRIDY  !j方向小半个网格
        enddo

        do n=1,nsta  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  观测点 n循环
              MAXJ=IFIX(RB(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia!qu zheng
              MAXJ=MIN0(jde,MAXJ)
              MINJ=IFIX(RB(N)-RINDX) !I方向的下界（影响半径内）
              MINJ=MAX0(1,MINJ) 
            cov_factor_t = 1.0
            IF( ABS(errf_u(n)) >= 99.0  .or. ABS(errf_v(n)) >= 99.0 .or. RKO(n) < 0.0   .or. RKO(n) >= 99.0 ) THEN  
            ! errf_u(n) 是第 n 个观测点的u/v误差估计值,如果是 无效值（大于99）或者垂直层k小于0，则认为该观测无效，跳过处理
            ! ELSEIF( RB(N).LT.FLOAT(MINJ) .OR. RB(N).GT.FLOAT(MAXJ) ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
            ELSEIF( J .LT. MINJ .OR. J .GT. MAXJ ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
            ELSE
              KOB=nint(RKO(N)) !k方向取整
              KOB=MIN0(kte,KOB) ! k方向取不大于kte
              KOB=MAX0(1,KOB)   !k方向取不小于1

              if(   RA(N).GE. ids   .and. RA(N).LE. float(ide) ) then         !      &
              ! .and. RB(N).GE. jds   .and. RB(N).LE.float(jde) ) then  !!! 检查该观测点是否落在当前子域（子网格）范围内，如果落在该进程负责范围，才继续处理
              ! RA(n)、RB(n) 是观测点的网格坐标（A方向 = X，B方向 = Y）
              ! ids ~ ide 是当前进程负责的网格 X 范围
              ! jds ~ jde 是当前进程负责的网格 Y 范围

              MAXI=IFIX(RA(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia
      !qu zheng
              MAXI=MIN0(ide,MAXI)
              MINI=IFIX(RA(N)-RINDX) !I方向的下界（影响半径内）
              MINI=MAX0(1,MINI)

                obs_loc_ptq%lon = lon_o(n) !观测点的经度  这里观测点经纬度为实际点的经纬度，不区分uvt的网格类型
              if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon  = obs_loc_ptq%lon-360.0
                obs_loc_ptq%lat = lat_o(n)
              if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
              coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)       
              else
              coslat = cos(80.0*DEG_TO_RAD)
              end if
              write(*,*) "coslat = " , coslat

               !!!!!!!!!!!!!!!! 1. 温度 (T) 增量计算
              if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                        DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                    model_loc%lat = T_grid%y(j)
                                  dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                  dist = RADIUS * sqrt(dist)
                                  dist0 = 1000*config_flags%e_r_6h*coslat
                                    cov_factor = 0
                                  cov_factor_h = comp_cov_factor(dist, dist0)
                                  cor_oi = cov_factor_h
                                  write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                  cov_factor = cov_factor_h* cov_factor_t * 0.1  !第一层不需要垂向衰减系数
                            do j_ens = 1, ens_size
                                enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)                              
                                ens_in_hfq(j_ens) =  ens_6h_t0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_t(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                             
                            if( assimilated_t(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then
                              ens_6h_t%data(i,kob,j) = ens_6h_t%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                              assimilated_t(i,kob,j) = .false.
                              write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                              write(*,*)  "at k=1 ens_6h_t%data(i,kob,j) = " , ens_6h_t%data(i,kob,j)
                            endif
                        ENDDO
              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                    kmin = IFIX(rko(n))
                    kmax = ceiling(rko(n)) !xiang shang qu zheng
                        DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j)
                            write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = 1000*config_flags%e_r_6h*coslat
                            write(*,*)  "dist =    dist0 = ",dist ,dist0
                            cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                            write(*,*)  "cov_factor_h" , cov_factor_h 
                            cov_factor = 0.0
                            if( kmax<=48 ) then
                            do k=kmin,kmax  !这里是
                                cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                write(*,*)  "cov_factor_v" , cov_factor_v
                                cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                write(*,*)  "cov_factor" , cov_factor                                
                                do j_ens = 1, ens_size
                                    enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_6h_t0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_t(n)
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_t(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then
                                    ens_6h_t%data(i,k,j) = ens_6h_t%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                    assimilated_t(i,k,j) = .false.
                                    ! write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_6h_t%data(i,k,j) = " , ens_6h_t%data(i,k,j)
                                endif 
                            enddo !k=kmin,kmax
                            endif ! if( kmax<=27 ) 
                      ENDDO
                endif !  if(kob .EQ.1 )
                !!!!!!!!!!!!!!!! 2. 风速 (u/v) 增量计算
                 if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                        DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                    model_loc%lat = T_grid%y(j)
                                  dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                  dist = RADIUS * sqrt(dist)
                                  dist0 = 1000*config_flags%e_r_6h*coslat
                                  cov_factor = 0
                                  cov_factor_h = comp_cov_factor(dist, dist0)
                                  cor_oi = cov_factor_h
                                  write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                  cov_factor = cov_factor_h* cov_factor_t* 0.1  !第一层不需要垂向衰减系数
                            !!!!!  U 
                            do j_ens = 1, ens_size
                                enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_6h_u0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_u(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_u(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0) then
                                ens_6h_u%data(i,kob,j) = ens_6h_u%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_u(i,kob,j) = .false.
                                write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                write(*,*)  "at k=1 ens_6h_u%data(i,kob,j) = " , ens_6h_u%data(i,kob,j)
                            endif
                            !!!!!  V 
                            do j_ens = 1, ens_size
                                enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_6h_v0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_v(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_v(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                ens_6h_v%data(i,kob,j) = ens_6h_v%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_v(i,kob,j) = .false.
                                write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                write(*,*)  "at k=1 ens_6h_v%data(i,kob,j) = " , ens_6h_v%data(i,kob,j)
                            endif
                        ENDDO
              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                    kmin = IFIX(rko(n))
                    kmax = ceiling(rko(n)) !xiang shang qu zheng
                        DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j)
                            write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = 1000*config_flags%e_r_6h*coslat
                            write(*,*)  "dist =    dist0 = ",dist ,dist0
                            cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                            write(*,*)  "cov_factor_h" , cov_factor_h 
                            cov_factor = 0.0

                            if( kmax<=48 ) then
                            do k=kmin,kmax  !这里是
                                cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                write(*,*)  "cov_factor_v" , cov_factor_v
                                cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                write(*,*)  "cov_factor" , cov_factor   
                                !!!!!  U                              
                                do j_ens = 1, ens_size
                                    enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_6h_u0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_u(n)  !!! 非常重要
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_u(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_6h_u%data(i,k,j) = ens_6h_u%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                    assimilated_u(i,k,j) = .false.
                                    write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_6h_u%data(i,k,j) = " , ens_6h_u%data(i,k,j)
                                endif
                                !!!!!  V                             
                                do j_ens = 1, ens_size
                                    enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_6h_v0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_v(n)  !!! 非常重要
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                 if( assimilated_v(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_6h_v%data(i,k,j) = ens_6h_v%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                    assimilated_v(i,k,j) = .false.
                                    write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_6h_v%data(i,k,j) = " , ens_6h_v%data(i,k,j)
                                endif
                            enddo !k=kmin,kmax
                            endif ! if( kmax<=27 ) 
                      ENDDO
                endif !  if(kob .EQ.1 )

              endif !ra
            ENDIF !IF
        enddo !n=1,6

  ENDDO  !J=jts,jte



 end subroutine ensemble_filter_atm_6h_r

 subroutine ensemble_filter_atm_3h_r(grid,config_flags,ens_3h_q0,ens_3h_q,&
                   &ens_3h_t0,ens_3h_t,ens_3h_u0,ens_3h_u,&
                   &ens_3h_v0,ens_3h_v,ens_3h_ps0,ens_3h_ps,&
                   &atm_var,obs_z_number,lat_o,lon_o, rio,rjo,rko,errf_u,errf_v,errf_t,errf_q,&
                  !  &Atm_var_g,obs_z_number,lat_o,lon_o,gridi,gridj,z_o2,errf_u,errf_v,errf_t,errf_q,&
                   &  T_grid, m_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts, jte,kts,kte,&
                     ids,  ide,  jds,  jde,  kds,  kde,dx)


    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: ens_3h_q0,ens_3h_t0,ens_3h_u0,ens_3h_v0,ens_3h_ps0
    type(ada_field_type), intent(inout) :: ens_3h_q,ens_3h_t,ens_3h_u,ens_3h_v,ens_3h_ps
    type(assim_data_type), intent(in) :: atm_var
    integer, intent(in) :: isd, ied, jsd, jed,kms, kme,its, ite, jts,jte, kts,kte
    integer, intent(in) :: obs_z_number
    real,intent(in) :: rio(obs_z_number),rjo(obs_z_number)
    real,intent(in) :: rko(obs_z_number),errf_u(obs_z_number),errf_v(obs_z_number),errf_t(obs_z_number),errf_q(obs_z_number)
    real,intent(in) :: lat_o(obs_z_number),lon_o(obs_z_number)
    integer, intent(in) :: ids,  ide,  jds,  jde,  kds,  kde

    type(ada_grid_type), intent(in) :: T_grid
    integer, intent(in)  :: atm_step
    type(WRFU_Time), intent(in) :: m_time
    real, INTENT(IN)  :: dx
    integer :: n,nnn,k
    integer :: nsndlev,njcsnd
    integer nsta,I,J,maxj,minj,maxi,mini,igrid
    real  ::  ensem_stn_t(20,obs_z_number),ensem_stn_u(20,obs_z_number),ensem_stn_v(20,obs_z_number)
    real :: gfactor,rfactor,gridx,gridy,rindx,ris
    real :: rj,rx,ri,ry,cov_factor_h,cov_factor,cor_oi,cov_factor_v,cov_factor_t
    integer :: kob,kobp,kmin,kmax,ens_size,j_ens


  type(loc_type) :: model_loc,obs_loc_ptq
  real :: dist, dist0,coslat
  character(len=1) :: VARI
  real :: std_oi_o, std_oi_g
  integer :: ass_variable = 1
  real :: ra(obs_z_number),rb(obs_z_number),rc(obs_z_number)

  logical, allocatable :: assimilated_t(:,:,:)
  logical, allocatable :: assimilated_u(:,:,:)
  logical, allocatable :: assimilated_v(:,:,:)

  ALLOCATE(assimilated_t(its:ite, kms:kme, jts:jte))
  ALLOCATE(assimilated_u(its:ite, kms:kme, jts:jte))
  ALLOCATE(assimilated_v(its:ite, kms:kme, jts:jte))

  assimilated_t = .false.
  assimilated_u = .true.
  assimilated_v = .true.

  nsta =  obs_z_number

  std_oi_o= 0.0
  std_oi_g=0.0

  ens_size = 20
  ensem_stn_t(:,:) = 0.
  ensem_stn_u(:,:) = 0.
  ensem_stn_v(:,:) = 0.
  write(*,*) 'running ensemble_filter_atm_3h_r'

  call ensem_stnt(grid,config_flags,ens_3h_t0,ensem_stn_t,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_3h_u0,ensem_stn_u,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_3h_v0,ensem_stn_v,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )

  GRIDX=0.5
  GRIDY=0.5
  IGRID=1

  RINDX=config_flags%e_r_3h *1000./dx   ! 37.037


  DO J=jts,jte  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  水平方向 j 循环
! return if too close to j boundary
        if((j.le.2.or.j.ge.jde-1)) then
          ! return
          CYCLE
        endif
     ! MAXJ=J+IFIX(RINDX)  !J方向的上界（影响半径内）
        ! MINJ=J-IFIX(RINDX)  !J方向的下界（影响半径内）
        do n=1,nsta
           RA(n)=rio(n)-GRIDX !i方向小半个网格
           RB(n)=rjo(n)-GRIDY !j方向小半个网格
        enddo

        do n=1,nsta    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  观测点 n循环
          MAXJ=IFIX(RB(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia!qu zheng
          MAXJ=MIN0(jde,MAXJ)
          MINJ=IFIX(RB(N)-RINDX) !I方向的下界（影响半径内）
          MINJ=MAX0(1,MINJ)
          cov_factor_t = 1.0
          IF( ABS(errf_u(n)) >= 99.0  .or. ABS(errf_v(n)) >= 99.0 .or. RKO(n) < 0.0  .or. RKO(n) >= 99.0  ) THEN  
              ! errf_u(n) 是第 n 个观测点的u/v误差估计值,如果是 无效值（大于99）或者垂直层k小于0，则认为该观测无效，跳过处理
          ! ELSEIF( RB(N).LT.FLOAT(MINJ) .OR. RB(N).GT.FLOAT(MAXJ) ) THEN ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
          ELSEIF( J .LT. MINJ .OR. J .GT. MAXJ ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
          ELSE
            KOB=nint(RKO(N))  !k方向取整
            KOB=MIN0(kte,KOB) ! k方向取不大于kte
            KOB=MAX0(1,KOB)   !k方向取不小于1
  
            if(   RA(N).GE. ids   .or. RA(N).LE. float(ide) )  then            ! &
            ! .and. RB(N).GE. jds   .and. RB(N).LE.float(jde) ) then!!! 检查该观测点是否落在当前子域（子网格）范围内，如果落在该进程负责范围，才继续处理
                ! RA(n)、RB(n) 是观测点的网格坐标（A方向 = X，B方向 = Y）
                ! ids ~ ide 是当前进程负责的网格 X 范围
                ! jds ~ jde 是当前进程负责的网格 Y 范围

            MAXI=IFIX(RA(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia
            MAXI=MIN0(ide,MAXI)
            MINI=IFIX(RA(N)-RINDX) !I方向的下界（影响半径内）
            MINI=MAX0(1,MINI)

              obs_loc_ptq%lon = lon_o(n)  !观测点的经度  这里观测点经纬度为实际点的经纬度，不区分uvt的网格类型
            if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon  = obs_loc_ptq%lon-360.0
              obs_loc_ptq%lat = lat_o(n)
            if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
            coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)
            else
            coslat = cos(80.0*DEG_TO_RAD)
            end if

            !!! 1. 温度 (T) 增量计算
              if(kob .EQ.1 ) then  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                      DO I=max0(its,MINI),min0(ite,MAXI)!!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                          model_loc%lon = T_grid%x(i)
                          if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                  model_loc%lat = T_grid%y(j)
                                dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                dist = RADIUS * sqrt(dist)
                                dist0 = 1000*config_flags%e_r_3h*coslat
                                  cov_factor = 0
                                cov_factor_h = comp_cov_factor(dist, dist0)
                                cor_oi = cov_factor_h
                                cov_factor = cov_factor_h* cov_factor_t* 0.1 !第一层不需要垂向衰减系数
                          do j_ens = 1, ens_size
                              enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)
                              ens_in_hfq(j_ens) =  ens_3h_t0(j_ens)%data(i,kob,j)
                              obs_inc_eakf_tr(j_ens) = errf_t(n)
                          enddo
                          ens_inc_hfq(:) = 0.0
                          call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                          &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                          &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                          &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                          & its, ite,jts, jte,obs_loc_ptq)
                          if( assimilated_t(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                              ens_3h_t%data(i,kob,j) = ens_3h_t%data(i,kob,j) + ens_inc_hfq(20)   !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_t(i,kob,j) = .false.
                          endif
                      ENDDO
                  else  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                  kmin = IFIX(rko(n))
                  kmax = ceiling(rko(n)) !xiang shang qu zheng
                      DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                          model_loc%lon = T_grid%x(i)
                          if ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                          model_loc%lat = T_grid%y(j)
                          dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                          dist = RADIUS * sqrt(dist)
                          dist0 = 1000*config_flags%e_r_3h*coslat
                          cov_factor_h = comp_cov_factor(dist, dist0)
                          cov_factor = 0
                          if( kmax<=48 ) then
                          do k=kmin,kmax
                              cov_factor_v =   1 - abs( (rko(n) - k) )  !垂直方向协方差系数
                              cov_factor   =  cov_factor_v * cov_factor_h * cov_factor_t
                              do j_ens = 1, ens_size
                                enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_3h_t0(j_ens)%data(i,k,j)
                                obs_inc_eakf_tr(j_ens) = errf_t(n)
                              enddo
                              ens_inc_hfq(:) = 0.0
                              call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                          &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                          &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                          &ass_method,ass_variable,i,j,n,k,atm_step,&
                                          & its, ite,jts, jte,obs_loc_ptq)
                              if( assimilated_t(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_3h_t%data(i,k,j) = ens_3h_t%data(i,k,j) + ens_inc_hfq(20)
                                    assimilated_t(i,k,j) = .false.
                              endif
                          enddo !k=kmin,kmax
                          endif ! if( kmax<=27 ) 
                    ENDDO
                  endif !  if(kob .EQ.1 )
                !!!!!!!!!!!!!!!!! 2. 风速 (u/v) 增量计算
                  if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                          DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                              model_loc%lon = T_grid%x(i)
                              if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                      model_loc%lat = T_grid%y(j)
                                    dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                    dist = RADIUS * sqrt(dist)
                                    dist0 = 1000*config_flags%e_r_3h*coslat
                                    cov_factor = 0
                                    cov_factor_h = comp_cov_factor(dist, dist0)
                                    cor_oi = cov_factor_h
                                    write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                    cov_factor = cov_factor_h* cov_factor_t* 0.1  !第一层不需要垂向衰减系数
                              !!!!!  U 
                              do j_ens = 1, ens_size
                                  enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                  ens_in_hfq(j_ens) =  ens_3h_u0(j_ens)%data(i,kob,j)
                                  obs_inc_eakf_tr(j_ens) = errf_u(n)
                              enddo
                              ens_inc_hfq(:) = 0.0
                              call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                              &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                              &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                              &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                              & its, ite,jts, jte,obs_loc_ptq)
                              if( assimilated_u(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                  ens_3h_u%data(i,kob,j) = ens_3h_u%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                  assimilated_u(i,kob,j) = .false.
                                  write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                  write(*,*)  "at k=1 ens_3h_u%data(i,kob,j) = " , ens_3h_u%data(i,kob,j)
                              endif
                              !!!!!  V 
                              do j_ens = 1, ens_size
                                  enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                  ens_in_hfq(j_ens) =  ens_3h_v0(j_ens)%data(i,kob,j)
                                  obs_inc_eakf_tr(j_ens) = errf_v(n)
                              enddo
                              ens_inc_hfq(:) = 0.0
                              call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                              &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                              &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                              &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                              & its, ite,jts, jte,obs_loc_ptq)
                              if( assimilated_v(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                  ens_3h_v%data(i,kob,j) = ens_3h_v%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                  assimilated_v(i,kob,j) = .false.
                                  write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                  write(*,*)  "at k=1 ens_3h_v%data(i,kob,j) = " , ens_3h_v%data(i,kob,j)
                              endif
                          ENDDO
                else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                      kmin = IFIX(rko(n))
                      kmax = ceiling(rko(n)) !xiang shang qu zheng
                          DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                              model_loc%lon = T_grid%x(i)
                              if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                              model_loc%lat = T_grid%y(j)
                              write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                              dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                              dist = RADIUS * sqrt(dist)
                              dist0 = 1000*config_flags%e_r_3h*coslat
                              write(*,*)  "dist =    dist0 = ",dist ,dist0
                              cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                              write(*,*)  "cov_factor_h" , cov_factor_h 
                              cov_factor = 0.0

                              if( kmax<=48 ) then
                              do k=kmin,kmax  !这里是
                                  cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                  write(*,*)  "cov_factor_v" , cov_factor_v
                                  cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                  ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                  write(*,*)  "cov_factor" , cov_factor   
                                  !!!!!  U                              
                                  do j_ens = 1, ens_size
                                      enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                      ens_in_hfq(j_ens) =  ens_3h_u0(j_ens)%data(i,k,j)
                                      obs_inc_eakf_tr(j_ens) = errf_u(n)  !!! 非常重要
                                  enddo
                                  ens_inc_hfq(:) = 0.0
                                  call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                                  &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                                  &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                  &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                  & its, ite,jts, jte,obs_loc_ptq)
                                  if( assimilated_u(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                      ens_3h_u%data(i,k,j) = ens_3h_u%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                      assimilated_u(i,k,j) = .false.
                                      write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                      write(*,*)  "ens_3h_u%data(i,k,j) = " , ens_3h_u%data(i,k,j)
                                  endif
                                  !!!!!  V                             
                                  do j_ens = 1, ens_size
                                      enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                      ens_in_hfq(j_ens) =  ens_3h_v0(j_ens)%data(i,k,j)
                                      obs_inc_eakf_tr(j_ens) = errf_v(n)  !!! 非常重要
                                  enddo
                                  ens_inc_hfq(:) = 0.0
                                  call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                                  &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                                  &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                  &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                  & its, ite,jts, jte,obs_loc_ptq)
                                  if( assimilated_v(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                      ens_3h_v%data(i,k,j) = ens_3h_v%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                      assimilated_v(i,k,j) = .false.
                                      write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                      write(*,*)  "ens_3h_v%data(i,k,j) = " , ens_3h_v%data(i,k,j)
                                  endif
                              enddo !k=kmin,kmax
                              endif ! if( kmax<=27 ) 
                        ENDDO
                  endif !  if(kob .EQ.1 )


            endif !ra
          ENDIF !IF
        enddo !n=1,6


    ENDDO  !J=jts,jte

 end subroutine ensemble_filter_atm_3h_r

 subroutine ensemble_filter_atm_day_r(grid,config_flags,ens_day_q0,ens_day_q,&
                   &ens_day_t0,ens_day_t,ens_day_u0,ens_day_u,&
                   &ens_day_v0,ens_day_v,ens_day_ps0,ens_day_ps,&
                  !  &atm_var, rio,rjo,rko,errf,lat1,lon1,&
                  &atm_var,obs_z_number,lat_o,lon_o, rio,rjo,rko,errf_u,errf_v,errf_t,errf_q,&
                   &  T_grid, m_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts, jte,kts,kte,&
                     ids,  ide,  jds,  jde,  kds,  kde,dx)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: ens_day_q0,ens_day_t0,ens_day_u0,ens_day_v0,ens_day_ps0
    type(ada_field_type), intent(inout) :: ens_day_q,ens_day_t,ens_day_u,ens_day_v,ens_day_ps
    type(assim_data_type), intent(in) :: atm_var
    integer, intent(in) :: isd, ied, jsd, jed,kms, kme,its, ite, jts,jte,kts,kte
    integer, intent(in) :: obs_z_number
    real,intent(in) :: rio(obs_z_number),rjo(obs_z_number)
    real,intent(in) :: rko(obs_z_number),errf_u(obs_z_number),errf_v(obs_z_number),errf_t(obs_z_number),errf_q(obs_z_number)
    real,intent(in) :: lat_o(obs_z_number),lon_o(obs_z_number)
    integer, intent(in) :: ids,  ide,  jds,  jde,  kds,  kde

    type(ada_grid_type), intent(in) :: T_grid
    integer, intent(in)  :: atm_step
    type(WRFU_Time), intent(in) :: m_time
    real, INTENT(IN)  :: dx
    integer :: n,nnn,k
    integer :: nsndlev,njcsnd
    integer nsta,I,J,maxj,minj,maxi,mini,igrid
    real  ::  ensem_stn_t(20,obs_z_number),ensem_stn_u(20,obs_z_number),ensem_stn_v(20,obs_z_number)
    real :: gfactor,rfactor,gridx,gridy,rindx,ris
    real :: rj,rx,ri,ry,cov_factor_h,cov_factor,cor_oi,cov_factor_v,cov_factor_t
    integer :: kob,kobp,kmin,kmax,ens_size,j_ens

  type(loc_type) :: model_loc,obs_loc_ptq
  real :: dist, dist0,coslat
  character(len=1) :: VARI
  real :: std_oi_o, std_oi_g
  integer :: ass_variable = 1
  real :: ra(obs_z_number),rb(obs_z_number),rc(obs_z_number)

  logical, allocatable :: assimilated_t(:,:,:)
  logical, allocatable :: assimilated_u(:,:,:)
  logical, allocatable :: assimilated_v(:,:,:)

  ALLOCATE(assimilated_t(its:ite, kms:kme, jts:jte))
  ALLOCATE(assimilated_u(its:ite, kms:kme, jts:jte))
  ALLOCATE(assimilated_v(its:ite, kms:kme, jts:jte))

  assimilated_t = .false.
  assimilated_u = .true.
  assimilated_v = .true.

  nsta =  obs_z_number

  std_oi_o= 0.0
  std_oi_g=0.0

  ens_size = 20
  ensem_stn_t(:,:) = 0.
  ensem_stn_u(:,:) = 0.
  ensem_stn_v(:,:) = 0.
  write(*,*) 'running ensemble_filter_atm_day_r'
  call ensem_stnt(grid,config_flags,ens_day_t0,ensem_stn_t,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_day_u0,ensem_stn_u,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_day_v0,ensem_stn_v,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )

  GRIDX=0.5
  GRIDY=0.5
  IGRID=1

  RINDX=config_flags%e_r_1day *1000./dx   ! 37.037

  DO J=jts,jte   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  水平方向 j 循环
        ! return if too close to j boundary
        if((j.le.2.or.j.ge.jde-1)) then
          ! return
          CYCLE
        endif
        ! MAXJ=J+IFIX(RINDX)  !J方向的上界（影响半径内）
        ! MINJ=J-IFIX(RINDX)  !J方向的下界（影响半径内）

        do n=1,nsta
           RA(n)=rio(n)-GRIDX  !i方向小半个网格
           RB(n)=rjo(n)-GRIDY  !j方向小半个网格
        enddo

        do n=1,nsta  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  观测点 n循环
              MAXJ=IFIX(RB(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia!qu zheng
              MAXJ=MIN0(jde,MAXJ)
              MINJ=IFIX(RB(N)-RINDX) !I方向的下界（影响半径内）
              MINJ=MAX0(1,MINJ)
            cov_factor_t = 1.0
            IF( ABS(errf_u(n)) >= 99.0  .or. ABS(errf_v(n)) >= 99.0 .or. RKO(n) < 0.0  .or. RKO(n) >= 99.0  ) THEN  
            ! errf_u(n) 是第 n 个观测点的u/v误差估计值,如果是 无效值（大于99）或者垂直层k小于0，则认为该观测无效，跳过处理
            ! ELSEIF( RB(N).LT.FLOAT(MINJ) .OR. RB(N).GT.FLOAT(MAXJ) ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
            ELSEIF( J .LT. MINJ .OR. J .GT. MAXJ ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
            ELSE
              KOB=nint(RKO(N))  !k方向取整
              KOB=MIN0(kte,KOB) ! k方向取不大于kte
              KOB=MAX0(1,KOB)   !k方向取不小于1

              if(   RA(N).GE. ids   .and. RA(N).LE. float(ide)  ) then  !           &
              ! .and. RB(N).GE. jds   .and. RB(N).LE.float(jde) ) then  !!! 检查该观测点是否落在当前子域（子网格）范围内，如果落在该进程负责范围，才继续处理
              ! RA(n)、RB(n) 是观测点的网格坐标（A方向 = X，B方向 = Y）
              ! ids ~ ide 是当前进程负责的网格 X 范围
              ! jds ~ jde 是当前进程负责的网格 Y 范围

              MAXI=IFIX(RA(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia
              !qu zheng
              MAXI=MIN0(ide,MAXI)
              MINI=IFIX(RA(N)-RINDX) !I方向的下界（影响半径内）
              MINI=MAX0(1,MINI)

                obs_loc_ptq%lon = lon_o(n) !观测点的经度  这里观测点经纬度为实际点的经纬度，不区分uvt的网格类型
              if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon  = obs_loc_ptq%lon-360.0
                obs_loc_ptq%lat = lat_o(n)
              if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
              coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)       
              else
              coslat = cos(80.0*DEG_TO_RAD)
              end if
              write(*,*) "coslat = " , coslat

               !!!!!!!!!!!!!!!! 1. 温度 (T) 增量计算
              if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                        DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                    model_loc%lat = T_grid%y(j)
                                  dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                  dist = RADIUS * sqrt(dist)
                                  dist0 = 1000*config_flags%e_r_1day*coslat
                                    cov_factor = 0
                                  cov_factor_h = comp_cov_factor(dist, dist0)
                                  cor_oi = cov_factor_h
                                  write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                  cov_factor = cov_factor_h* cov_factor_t* 0.1  !第一层不需要垂向衰减系数
                            do j_ens = 1, ens_size
                                enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)                              
                                ens_in_hfq(j_ens) =  ens_day_t0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_t(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_t(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                ens_day_t%data(i,kob,j) = ens_day_t%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_t(i,kob,j)  = .false.
                                write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                write(*,*)  "at k=1 ens_day_t%data(i,kob,j) = " , ens_day_t%data(i,kob,j)
                            endif
                        ENDDO
              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                    kmin = IFIX(rko(n))
                    kmax = ceiling(rko(n)) !xiang shang qu zheng
                        DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j)
                            write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = 1000*config_flags%e_r_1day*coslat
                            write(*,*)  "dist =    dist0 = ",dist ,dist0
                            cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                            write(*,*)  "cov_factor_h" , cov_factor_h 
                            cov_factor = 0.0
                            if( kmax<=48 ) then
                            do k=kmin,kmax  !这里是
                                cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                write(*,*)  "cov_factor_v" , cov_factor_v
                                cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                write(*,*)  "cov_factor" , cov_factor                                
                                do j_ens = 1, ens_size
                                    enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_day_t0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_t(n)
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_t(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_day_t%data(i,k,j) = ens_day_t%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                    assimilated_t(i,k,j) = .false.
                                    ! write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_day_t%data(i,k,j) = " , ens_day_t%data(i,k,j)
                                endif
                            enddo !k=kmin,kmax
                            endif ! if( kmax<=27 ) 
                      ENDDO
                endif !  if(kob .EQ.1 )
                 !!!!!!!!!!!!!!!! 2. 风速 (u/v) 增量计算
                 if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                        DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                    model_loc%lat = T_grid%y(j)
                                  dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                  dist = RADIUS * sqrt(dist)
                                  dist0 = 1000*config_flags%e_r_1day*coslat
                                  cov_factor = 0
                                  cov_factor_h = comp_cov_factor(dist, dist0)
                                  cor_oi = cov_factor_h
                                  write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                  cov_factor = cov_factor_h* cov_factor_t* 0.1  !第一层不需要垂向衰减系数
                            !!!!!  U 
                            do j_ens = 1, ens_size
                                enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_day_u0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_u(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                             if( assimilated_u(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                ens_day_u%data(i,kob,j) = ens_day_u%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_u(i,kob,j) = .false.
                                write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                write(*,*)  "at k=1 ens_day_u%data(i,kob,j) = " , ens_day_u%data(i,kob,j)
                            endif 
                            !!!!!  V 
                            do j_ens = 1, ens_size
                                enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_day_v0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_v(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_v(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                ens_day_v%data(i,kob,j) = ens_day_v%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_v(i,kob,j)  = .false.
                                write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                write(*,*)  "at k=1 ens_day_v%data(i,kob,j) = " , ens_day_v%data(i,kob,j)
                            endif
                        ENDDO
              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                    kmin = IFIX(rko(n))
                    kmax = ceiling(rko(n)) !xiang shang qu zheng
                        DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j)
                            write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = 1000*config_flags%e_r_1day*coslat
                            write(*,*)  "dist =    dist0 = ",dist ,dist0
                            cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                            write(*,*)  "cov_factor_h" , cov_factor_h 
                            cov_factor = 0.0

                            if( kmax<=48 ) then
                            do k=kmin,kmax  !这里是
                                cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                write(*,*)  "cov_factor_v" , cov_factor_v
                                cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                write(*,*)  "cov_factor" , cov_factor   
                                !!!!!  U                              
                                do j_ens = 1, ens_size
                                    enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_day_u0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_u(n)  !!! 非常重要
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_u(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_day_u%data(i,k,j) = ens_day_u%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                     assimilated_u(i,k,j) = .false.
                                    write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_day_u%data(i,k,j) = " , ens_day_u%data(i,k,j)
                                endif
                                !!!!!  V                             
                                do j_ens = 1, ens_size
                                    enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_day_v0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_v(n)  !!! 非常重要
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_v(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_day_v%data(i,k,j) = ens_day_v%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                    assimilated_v(i,k,j) = .false.
                                    write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_day_v%data(i,k,j) = " , ens_day_v%data(i,k,j)
                                endif 
                            enddo !k=kmin,kmax
                            endif ! if( kmax<=27 ) 
                      ENDDO
                endif !  if(kob .EQ.1 )

              endif !ra
            ENDIF !IF
        enddo !n=1,6

  ENDDO  !J=jts,jte


 end subroutine ensemble_filter_atm_day_r


 subroutine ensemble_filter_atm_5day_r(grid,config_flags,ens_5day_q0,ens_5day_q,&
                   &ens_5day_t0,ens_5day_t,ens_5day_u0,ens_5day_u,&
                   &ens_5day_v0,ens_5day_v,ens_5day_ps0,ens_5day_ps,&
                  !  &atm_var, rio,rjo,rko,errf,lat1,lon1,&
                  &atm_var,obs_z_number,lat_o,lon_o, rio,rjo,rko,errf_u,errf_v,errf_t,errf_q,&
                   &  T_grid, m_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts, jte,kts,kte,&
                     ids,  ide,  jds,  jde,  kds,  kde,dx)

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(ada_field_type), dimension(:), intent(inout) :: ens_5day_q0,ens_5day_t0,ens_5day_u0,ens_5day_v0,ens_5day_ps0
    type(ada_field_type), intent(inout) :: ens_5day_q,ens_5day_t,ens_5day_u,ens_5day_v,ens_5day_ps
    type(assim_data_type), intent(in) :: atm_var
    integer, intent(in) :: isd, ied, jsd, jed,kms, kme,its, ite, jts,jte,kts,kte
    integer, intent(in) :: obs_z_number
    real,intent(in) :: rio(obs_z_number),rjo(obs_z_number)
    real,intent(in) :: rko(obs_z_number),errf_u(obs_z_number),errf_v(obs_z_number),errf_t(obs_z_number),errf_q(obs_z_number)
    real,intent(in) :: lat_o(obs_z_number),lon_o(obs_z_number)
    integer, intent(in) :: ids,  ide,  jds,  jde,  kds,  kde

    type(ada_grid_type), intent(in) :: T_grid
    integer, intent(in)  :: atm_step
    type(WRFU_Time), intent(in) :: m_time
    real, INTENT(IN)  :: dx
    integer :: n,nnn,k
    integer :: nsndlev,njcsnd
    integer nsta,I,J,maxj,minj,maxi,mini,igrid
    real  ::  ensem_stn_t(20,obs_z_number),ensem_stn_u(20,obs_z_number),ensem_stn_v(20,obs_z_number)
    real :: gfactor,rfactor,gridx,gridy,rindx,ris
    real :: rj,rx,ri,ry,cov_factor_h,cov_factor,cor_oi,cov_factor_v,cov_factor_t
    integer :: kob,kobp,kmin,kmax,ens_size,j_ens

  type(loc_type) :: model_loc,obs_loc_ptq
  real :: dist, dist0,coslat
  character(len=1) :: VARI
  real :: std_oi_o, std_oi_g
  integer :: ass_variable = 1
  real :: ra(obs_z_number),rb(obs_z_number),rc(obs_z_number)

  logical, allocatable :: assimilated_t(:,:,:)
  logical, allocatable :: assimilated_u(:,:,:)
  logical, allocatable :: assimilated_v(:,:,:)

  ALLOCATE(assimilated_t(its:ite, kms:kme, jts:jte))
  ALLOCATE(assimilated_u(its:ite, kms:kme, jts:jte))
  ALLOCATE(assimilated_v(its:ite, kms:kme, jts:jte))

  assimilated_t = .false.
  assimilated_u = .true.
  assimilated_v = .true.

  nsta = obs_z_number

  std_oi_o= 0.0
  std_oi_g=0.0

  ens_size = 20
  ensem_stn_t(:,:) = 0.
  ensem_stn_u(:,:) = 0.
  ensem_stn_v(:,:) = 0.
  write(*,*) 'running ensemble_filter_atm_5day_r'

  call ensem_stnt(grid,config_flags,ens_5day_t0,ensem_stn_t,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_5day_u0,ensem_stn_u,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )
  call ensem_stnt(grid,config_flags,ens_5day_v0,ensem_stn_v,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
                ids,  ide,  jds,  jde,  kds,  kde )

  GRIDX=0.5
  GRIDY=0.5
  IGRID=1

  RINDX=config_flags%e_r_5day *1000./dx   ! 37.037


  DO J=jts,jte   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  水平方向 j 循环
        ! return if too close to j boundary
        if((j.le.2.or.j.ge.jde-1)) then
          ! return
          CYCLE
        endif
        ! MAXJ=J+IFIX(RINDX)  !J方向的上界（影响半径内）
        ! MINJ=J-IFIX(RINDX)  !J方向的下界（影响半径内）
        do n=1,nsta
           RA(n)=rio(n)-GRIDX  !i方向小半个网格
           RB(n)=rjo(n)-GRIDY  !j方向小半个网格
        enddo

        do n=1,nsta  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  观测点 n循环
              MAXJ=IFIX(RB(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia!qu zheng
              MAXJ=MIN0(jde,MAXJ)
              MINJ=IFIX(RB(N)-RINDX) !I方向的下界（影响半径内）
              MINJ=MAX0(1,MINJ)
            cov_factor_t = 1.0
            IF( ABS(errf_u(n)) >= 99.0  .or. ABS(errf_v(n)) >= 99.0 .or. RKO(n) < 0.0  .or. RKO(n) >= 99.0 ) THEN  
            ! errf_u(n) 是第 n 个观测点的u/v误差估计值,如果是 无效值（大于99）或者垂直层k小于0，则认为该观测无效，跳过处理
            ! ELSEIF( RB(N).LT.FLOAT(MINJ) .OR. RB(N).GT.FLOAT(MAXJ) ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
            ELSEIF( J .LT. MINJ .OR. J .GT. MAXJ ) THEN  ! 观测点不在当前处理范围内（MINJ ~ MAXJ），则跳过
            write(*,*)  "J .LT. MINJ .OR. J .GT. MAXJ",J , MINJ , MAXJ
            ELSE
              write(*,*) "RKO(N) = " , RKO(N)
              KOB=nint(RKO(N))  !k方向取整
              KOB=MIN0(kte,KOB) ! k方向取不大于kte
              KOB=MAX0(1,KOB)   !k方向取不小于1

              if(   RA(N).GE. ids   .and. RA(N).LE. float(ide) ) then    !           &
              ! .and. RB(N).GE. jds   .and. RB(N).LE.float(jde) ) then  !!! 检查该观测点是否落在当前子域（子网格）范围内，如果落在该进程负责范围，才继续处理
              ! RA(n)、RB(n) 是观测点的网格坐标（A方向 = X，B方向 = Y）
              ! ids ~ ide 是当前进程负责的网格 X 范围
              ! jds ~ jde 是当前进程负责的网格 Y 范围

              MAXI=IFIX(RA(N)+RINDX)  !I方向的上界（影响半径内） ! xiang 0 qu zheng,dui zheng shu ,ji xiang xia
              !qu zheng
              MAXI=MIN0(ide,MAXI)
              MINI=IFIX(RA(N)-RINDX) !I方向的下界（影响半径内）
              MINI=MAX0(1,MINI)

                obs_loc_ptq%lon = lon_o(n) !观测点的经度  这里观测点经纬度为实际点的经纬度，不区分uvt的网格类型
              if ( obs_loc_ptq%lon > 360.0 ) obs_loc_ptq%lon  = obs_loc_ptq%lon-360.0
                obs_loc_ptq%lat = lat_o(n)
              if ( abs(obs_loc_ptq%lat) <= 80.0 ) then
              coslat = cos(obs_loc_ptq%lat*DEG_TO_RAD)       
              else
              coslat = cos(80.0*DEG_TO_RAD)
              end if
              write(*,*) "coslat = " , coslat

              !!!!!!!!!!!!!!!! 1. 温度 (T) 增量计算
              if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                        DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                    model_loc%lat = T_grid%y(j)
                                  dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                  dist = RADIUS * sqrt(dist)
                                  dist0 = 1000*config_flags%e_r_5day*coslat
                                    cov_factor = 0
                                  cov_factor_h = comp_cov_factor(dist, dist0)
                                  cor_oi = cov_factor_h
                                  write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                  cov_factor = cov_factor_h* cov_factor_t* 0.1  !第一层不需要垂向衰减系数
                            do j_ens = 1, ens_size
                                enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)                              
                                ens_in_hfq(j_ens) =  ens_5day_t0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_t(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_t(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                              ens_5day_t%data(i,kob,j) = ens_5day_t%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_t(i,kob,j) = .false.
                              write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                              write(*,*)  "at k=1 ens_5day_t%data(i,kob,j) = " , ens_5day_t%data(i,kob,j)
                            endif
                        ENDDO
              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                    kmin = IFIX(rko(n))
                    kmax = ceiling(rko(n)) !xiang shang qu zheng
                        DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j)
                            write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = 1000*config_flags%e_r_5day*coslat
                            write(*,*)  "dist =    dist0 = ",dist ,dist0
                            cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                            write(*,*)  "cov_factor_h" , cov_factor_h 
                            cov_factor = 0.0
                            if( kmax<=48 ) then
                            do k=kmin,kmax  !这里是
                                cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                write(*,*)  "cov_factor_v" , cov_factor_v
                                cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                write(*,*)  "cov_factor" , cov_factor                                
                                do j_ens = 1, ens_size
                                    enso_t_r_o(j_ens) =  ensem_stn_t(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_5day_t0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_t(n)
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_t_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_t,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_t(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                  ens_5day_t%data(i,k,j) = ens_5day_t%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                  assimilated_t(i,k,j) = .false.
                                  ! write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                  write(*,*)  "ens_5day_t%data(i,k,j) = " , ens_5day_t%data(i,k,j)
                                endif 
                            enddo !k=kmin,kmax
                            endif ! if( kmax<=27 ) 
                      ENDDO
                endif !  if(kob .EQ.1 )
                !!!!!!!!!!!!!!!! 2. 风速 (u/v) 增量计算
                 if(kob .EQ.1 ) then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K=1时
                        DO I=max0(its,MINI),min0(ite,MAXI) !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if  ( model_loc%lon > 360.0 ) model_loc%lon = model_loc%lon-360.0
                                    model_loc%lat = T_grid%y(j)
                                  dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                                  dist = RADIUS * sqrt(dist)
                                  dist0 = 1000*config_flags%e_r_5day*coslat
                                  cov_factor = 0
                                  cov_factor_h = comp_cov_factor(dist, dist0)
                                  cor_oi = cov_factor_h
                                  write(*,*)  "at k=1 cov_factor_h" , cov_factor_h
                                  cov_factor = cov_factor_h* cov_factor_t* 0.1  !第一层不需要垂向衰减系数
                            !!!!!  U 
                            do j_ens = 1, ens_size
                                enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                write(*,*) "ensem_stn_u = " , ensem_stn_u(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_5day_u0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_u(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_u(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                              ens_5day_u%data(i,kob,j) = ens_5day_u%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                              assimilated_u(i,kob,j) = .false.
                              write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                              write(*,*)  "at k=1 ens_5day_u%data(i,kob,j) = " , ens_5day_u%data(i,kob,j)
                            endif
                            !!!!!  V 
                            do j_ens = 1, ens_size
                                enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                ens_in_hfq(j_ens) =  ens_5day_v0(j_ens)%data(i,kob,j)
                                obs_inc_eakf_tr(j_ens) = errf_v(n)
                            enddo
                            ens_inc_hfq(:) = 0.0
                            call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                            &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                            &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                            &ass_method,ass_variable,i,j,n,kob,atm_step,&
                                            & its, ite,jts, jte,obs_loc_ptq)
                            if( assimilated_v(i,kob,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                ens_5day_v%data(i,kob,j) = ens_5day_v%data(i,kob,j) + ens_inc_hfq(20)    !!! 为什么是20????  因为这里是20个ENS，第20个是最靠近同时次的ens
                                assimilated_v(i,kob,j) = .false. 
                                write(*,*)  "i=,kob=,j=,ens_inc_hfq(20)=" ,i,kob,j, ens_inc_hfq(20)
                                write(*,*)  "at k=1 ens_5day_v%data(i,kob,j) = " , ens_5day_v%data(i,kob,j)
                            endif
                        ENDDO
              else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 垂向 K/=1时
                    kmin = IFIX(rko(n))
                    kmax = ceiling(rko(n)) !xiang shang qu zheng
                        DO I=max0(its,MINI),min0(ite,MAXI)  !!!!!!!!!!!!!!!!!!!!!! 水平方向 i 循环
                            model_loc%lon = T_grid%x(i)
                            if ( model_loc%lon > 360.0 ) model_loc%lon =model_loc%lon-360.0
                            model_loc%lat = T_grid%y(j)
                            write(*,*)  "model_loc=, obs_loc_ptq=" , model_loc, obs_loc_ptq
                            dist = get_dist(model_loc, obs_loc_ptq, atm_step)
                            dist = RADIUS * sqrt(dist)
                            dist0 = 1000*config_flags%e_r_5day*coslat
                            write(*,*)  "dist =    dist0 = ",dist ,dist0
                            cov_factor_h = comp_cov_factor(dist, dist0)  !! ？ ==0?   !水平方向协方差系数
                            write(*,*)  "cov_factor_h" , cov_factor_h 
                            cov_factor = 0.0

                            if( kmax<=48 ) then
                            do k=kmin,kmax  !这里是
                                cov_factor_v =   1 - abs( (rko(n) - k) )    !垂直方向协方差系数
                                write(*,*)  "cov_factor_v" , cov_factor_v
                                cov_factor   =  cov_factor_v *cov_factor_h * cov_factor_t
                                ! cov_factor   =  (0.5*cov_factor_v + 0.5*cov_factor_h) * cov_factor_t  !! 水平和垂直的增益系数各自为0.5  
                                write(*,*)  "cov_factor" , cov_factor   
                                !!!!!  U                              
                                do j_ens = 1, ens_size
                                    enso_u_r_o(j_ens) =  ensem_stn_u(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_5day_u0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_u(n)  !!! 非常重要
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_u_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_u,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_u(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                  ens_5day_u%data(i,k,j) = ens_5day_u%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                  assimilated_u(i,k,j) = .false.
                                  write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                  write(*,*)  "ens_5day_u%data(i,k,j) = " , ens_5day_u%data(i,k,j)
                                endif 
                                !!!!!  V                             
                                do j_ens = 1, ens_size
                                    enso_v_r_o(j_ens) =  ensem_stn_v(j_ens,n)
                                    ens_in_hfq(j_ens) =  ens_5day_v0(j_ens)%data(i,k,j)
                                    obs_inc_eakf_tr(j_ens) = errf_v(n)  !!! 非常重要
                                enddo
                                ens_inc_hfq(:) = 0.0
                                call update_from_obs_inc_k(config_flags,enso_v_r_o,obs_inc_eakf_tr,&
                                                &obs_inc_oi_v,ens_in_hfq(:),ens_size,&
                                                &ens_inc_hfq,cov_factor,cor_oi,std_oi_o,std_oi_g,&
                                                &ass_method,ass_variable,i,j,n,k,atm_step,&
                                                & its, ite,jts, jte,obs_loc_ptq)
                                if( assimilated_v(i,k,j) .and. abs(ens_inc_hfq(20)) < 20.0 ) then 
                                    ens_5day_v%data(i,k,j) = ens_5day_v%data(i,k,j) + ens_inc_hfq(20) !!??==0
                                    assimilated_v(i,k,j) = .false.
                                    write(*,*)  "i=,k=,j=,ens_inc_hfq(20)=" ,i,k,j, ens_inc_hfq(20)
                                    write(*,*)  "ens_5day_v%data(i,k,j) = " , ens_5day_v%data(i,k,j)
                                endif 
                            enddo !k=kmin,kmax
                            endif ! if( kmax<=27 ) 
                      ENDDO
                endif !  if(kob .EQ.1 )

              endif !ra
            ENDIF !IF
        enddo !n=1,6

  ENDDO  !J=jts,jte

 end subroutine ensemble_filter_atm_5day_r


!   subroutine ensem_stnt(grid,config_flags,ens_6h_any,ensem_stn,atm_step,rio,rjo,rko, obs_z_number,&
!         & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
!         ids,  ide,  jds,  jde,  kds,  kde )
!     TYPE ( domain ), INTENT(INOUT) :: grid
!     TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
!     type(ada_field_type), dimension(:), intent(inout) :: ens_6h_any
!     integer, intent(in)  :: atm_step, obs_z_number
!     REAL, INTENT(IN)     :: rio(obs_z_number)                 ! Obs west-east            
!     REAL, INTENT(IN)     :: rjo(obs_z_number)                 ! Obs south-northcoordinate (non-stag grid).
!     REAL, INTENT(IN)     :: rko(obs_z_number)                 ! Obs vertical coordinate.
!     integer, intent(in)  :: isd, ied, jsd, jed,kms, kme,its, ite,jts,jte, kts, kte
!     integer, intent(in)  :: ids,  ide,  jds,  jde,  kds,  kde
!     real, intent(inout)  ::  ensem_stn(20,obs_z_number)
!   INTEGER :: IA(obs_z_number)
!   INTEGER :: IB(obs_z_number)
!   INTEGER :: IC(obs_z_number)
! ! new
!   integer nsta,n_obs,iob,job,kob,KOBP,iob_ms,job_ms,n
!   integer i_start,i_end,j_start,j_end    ! loop ranges for uratio,vratio calc.
!   integer k_start,k_end,n_en

!   real, parameter :: gridx_t = 0.5     ! Mass-point x stagger
!   real, parameter :: gridy_t = 0.5     ! Mass-point y stagger
!   real :: dummy = 999.
!   real :: ra(obs_z_number),rb(obs_z_number),rc(obs_z_number)
!    real gridx,gridy,dxob,dyob,dzob
!   real:: gridi(obs_z_number),gridj(obs_z_number)

!   gridi(:) = rio(:)
!   gridj(:) = rjo(:)
!   nsta = obs_z_number



!    do  n_obs = 1,nsta
!         ra(n_obs) = 0.
!         rb(n_obs) = 0.
!         rc(n_obs) = 0.
!    enddo
!    ensem_stn(:,1:nsta) = 0.0
! ! set boundary
!   i_start = max( its-1,ids )
!   i_end   = min( ite+1,ide-1 )
!   j_start = max( jts-1,jds )
!   j_end   = min( jte+1,jde-1 )
!   k_start = kts
!   k_end = min( kte, kde-1 )

!        GRIDX = gridx_t
!        GRIDY = gridy_t




!   do n = 1,nsta
!         RA(N)=rio(n) -GRIDX
!         RB(N)=rjo(n) -GRIDY
!         RC(N)=rko(n)

!         IA(N)=RA(N)
!         IB(N)=RB(N)
!         IC(N)=RC(N)  ! ic is integer ,rc is real 

!         RA(N)=RA(N)-FLOAT(IA(N))
!         RB(N)=RB(N)-FLOAT(IB(N))
!         RC(N)=RC(N)-FLOAT(IC(N))

!         IOB=MAX0(1,IA(N))
!         IOB=MIN0(IOB,ide-1)
!         JOB=MAX0(1,IB(N))
!         JOB=MIN0(JOB,jde-1)
!         DXOB=RA(N)-FLOAT(IA(N))
!         DYOB=RB(N)-FLOAT(IB(N))

! ! outside obs is big value
!     do n_en=1,20
!         IF((RA(N)-1.).LT.0)THEN
!            ensem_stn(n_en,N)=ensem_stn(n_en,N)+DUMMY
!         ENDIF
!         IF((RB(N)-1.).LT.0)THEN
!            ensem_stn(n_en,N)=ensem_stn(n_en,N)+DUMMY
!         ENDIF
!         IF((FLOAT(ide)-2.0*GRIDX-RA(N)-1.E-10).LT.0)THEN
!            ensem_stn(n_en,N)=ensem_stn(n_en,N)+DUMMY
!         ENDIF
!         IF((FLOAT(jde)-2.0*GRIDY-RB(N)-1.E-10).LT.0)THEN
!            ensem_stn(n_en,N)=ensem_stn(n_en,N)+DUMMY
!         ENDIF
!         if(rc(n).lt.1.) ensem_stn(n_en,n)=ensem_stn(n_en,n)+dummy
!      enddo


!   enddo !nsta

 
!   do n = 1,nsta

!        if ((isd <= gridi(n) .AND. gridi(n) <= ied) .AND. (jsd <= gridj(n) .AND. gridj(n) <= jed)) then

!             KOB=MIN0(k_end,IC(N))
!             KOBP=MIN0(KOB+1,K_END)
!             DZOB = RC(N)

!       do n_en=1,20
!          ensem_stn(n_en,n) =  ensem_stn(n_en,n) + ((1.-DZOB)*((1.-DyOB)*((1.-&
!                       DxOB)*ens_6h_any(n_en)%data(IOB,KOB,JOB)+DxOB* &
!                       ens_6h_any(n_en)%data(IOB+1,KOB,JOB))+DyOB*((1.-DxOB)* &
!                       ens_6h_any(n_en)%data(IOB,KOB,JOB+1)+DxOB* &
!                       ens_6h_any(n_en)%data(IOB+1,KOB,JOB+1)))+DZOB*((1.- &
!                       DyOB)*((1.-DxOB)*ens_6h_any(n_en)%data(IOB,KOBP,JOB)+DxOB &
!                       *ens_6h_any(n_en)%data(IOB+1,KOBP,JOB))+DyOB*((1.-DxOB &
!                       )*ens_6h_any(n_en)%data(IOB,KOBP,JOB+1)+DxOB* &
!                       ens_6h_any(n_en)%data(IOB+1,KOBP,JOB+1))))

!       write(*,*) 'n_en,n,ensem_stn',n_en,n,ensem_stn(n_en,n)

!       enddo
!         endif ! check obs location

!   enddo !nsta

!   end subroutine ensem_stnt




subroutine ensem_stnt(grid,config_flags,ens_6h_any,ensem_stn,atm_step,rio,rjo,rko, obs_z_number,&
        & isd, ied, jsd, jed,kms, kme,its, ite,jts, jte, kts,kte,&
        ids,  ide,  jds,  jde,  kds,  kde )

  use module_domain       ! 确保 grid 类型定义正确
  implicit none

  TYPE ( domain ), INTENT(INOUT) :: grid
  TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
  type(ada_field_type), dimension(:), intent(inout) :: ens_6h_any
  integer, intent(in)  :: atm_step, obs_z_number
  REAL, INTENT(IN)     :: rio(obs_z_number)
  REAL, INTENT(IN)     :: rjo(obs_z_number)
  REAL, INTENT(IN)     :: rko(obs_z_number)
  integer, intent(in)  :: isd, ied, jsd, jed,kms, kme,its, ite,jts,jte, kts, kte
  integer, intent(in)  :: ids,  ide,  jds,  jde,  kds,  kde
  real, intent(inout)  :: ensem_stn(20,obs_z_number)

  INTEGER :: IA(obs_z_number), IB(obs_z_number), IC(obs_z_number)
  INTEGER :: nsta, n_obs, n_en
  INTEGER :: IOB, JOB, KOB, KOBP
  REAL :: RA(obs_z_number), RB(obs_z_number), RC(obs_z_number)
  REAL :: DXOB, DYOB, DZOB
  REAL :: GRIDX, GRIDY
  REAL :: gridi(obs_z_number), gridj(obs_z_number)
  REAL :: dummy
  LOGICAL :: valid_obs(obs_z_number)

  ! 初始化
  GRIDX = 0.5
  GRIDY = 0.5
  dummy = 999.0
  nsta = obs_z_number
  valid_obs(:) = .true.
  gridi(:) = rio(:)
  gridj(:) = rjo(:)

  ensem_stn(:, :) = 0.0

  ! 第一遍：计算插值位置和合法性判断
  do n_obs = 1, nsta
     RA(n_obs) = rio(n_obs) - GRIDX
     RB(n_obs) = rjo(n_obs) - GRIDY
     RC(n_obs) = rko(n_obs)

     IA(n_obs) = int(RA(n_obs))
     IB(n_obs) = int(RB(n_obs))
     IC(n_obs) = int(RC(n_obs))

     RA(n_obs) = RA(n_obs) - float(IA(n_obs))
     RB(n_obs) = RB(n_obs) - float(IB(n_obs))
     RC(n_obs) = RC(n_obs) - float(IC(n_obs))

     ! 对应插值立方体起点格点
     IOB = IA(n_obs)
     JOB = IB(n_obs)

     ! 越界检查：必须确保 IOB~IOB+1, JOB~JOB+1, KOB~KOB+1 都在域内
     if (IOB < isd .or. IOB+1 > ied .or. JOB < jsd .or. JOB+1 > jed .or. &
         IC(n_obs) < kds .or. IC(n_obs)+1 > kde) then
         valid_obs(n_obs) = .false.
         ensem_stn(:, n_obs) = dummy
     endif
  enddo

  ! 第二遍：进行三维插值
  do n_obs = 1, nsta
     if (.not. valid_obs(n_obs)) cycle

     IOB = IA(n_obs)
     JOB = IB(n_obs)
     KOB = IC(n_obs)
     KOBP = min(KOB + 1, kde)

     DXOB = RA(n_obs)
     DYOB = RB(n_obs)
     DZOB = RC(n_obs)

     do n_en = 1, 20
        ensem_stn(n_en,n_obs) = (1.0 - DZOB) * &
            ((1.0 - DYOB) * ((1.0 - DXOB) * ens_6h_any(n_en)%data(IOB,KOB,JOB) + DXOB * ens_6h_any(n_en)%data(IOB+1,KOB,JOB)) + &
             DYOB * ((1.0 - DXOB) * ens_6h_any(n_en)%data(IOB,KOB,JOB+1) + DXOB * ens_6h_any(n_en)%data(IOB+1,KOB,JOB+1))) + &
          DZOB * &
            ((1.0 - DYOB) * ((1.0 - DXOB) * ens_6h_any(n_en)%data(IOB,KOBP,JOB) + DXOB * ens_6h_any(n_en)%data(IOB+1,KOBP,JOB)) + &
             DYOB * ((1.0 - DXOB) * ens_6h_any(n_en)%data(IOB,KOBP,JOB+1) + DXOB * ens_6h_any(n_en)%data(IOB+1,KOBP,JOB+1)))
          write(*,*) 'n_en,n_obs,ensem_stn',n_en,n_obs,ensem_stn(n_en,n_obs)
     enddo
  enddo

end subroutine ensem_stnt


end module module_eakf_ada_mod
