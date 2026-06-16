module module_ada_driver_fv_mod
  !
  ! This is the top-level module for ocean data assimilation.
  !
    USE module_utility  
! it seems like that this module will get type(wrfu_time)
    USE module_ada_types_mod, only : assim_data_type,ada_field_type,ada_grid_type
!    USE module_ada_core_mod, only : ada_core_init, get_ncep_field
    USE module_ada_core_mod, only :  get_ncep_field,get_obs_real
    USE module_eakf_ada_mod, only : ensemble_filter_atm_init,ensemble_filter_atm_6h_r,ensemble_filter_atm_3h_r,&
                                        & ensemble_filter_atm_5day_r,ensemble_filter_atm_day_r
    USE module_eakf_ada_mod, only : ensemble_filter_atm_hfq,ensemble_filter_atm_day,ensemble_filter_atm_6h,&
                                  & ensemble_filter_atm_3h,ensemble_filter_atm_5day
    USE netcdf
    USE module_machine  !wangj0710
    USE module_driver_constants  !wangj0710
    USE module_dm
    USE module_domain, ONLY : domain
    USE module_configure, ONLY : grid_config_rec_type
!  USE module_ada_t_u_v_q_mod
  implicit none
   INCLUDE 'mpif.h' 
  private
  public :: init_fv_ada, fv_ada,mda_o,get_z_o,get_obs_errf,isc_flth, iec_flth, jsc_flth,jec_flth,isd,&
            ied, jsd, jed, kms, kme,assim_frequency, heaada_enssize,&
                mpi_comm_here, ierr, myproc,nproc
  real, public, parameter :: PI      = 3.14159265358979323846
  real, public, parameter :: RAD_TO_DEG=PI/180.
  INTEGER mpi_comm_here, ierr, myproc,nproc
  integer, parameter :: NO_ASSIM = 0, OI=1, EAKF=2
  integer, parameter :: MAX_ATM_OBS = 20000
  integer, parameter :: MAX_SICE_OBS = 50000
  integer, parameter :: MAX_ENSEMBLE_SIZE = 50

  integer :: heaada_enssize = 20
  integer :: heahfq_minintv = 20
  integer :: isd, ied, jsd, jed, kms, kme !halo domain
  integer :: isc_flth, iec_flth, jsc_flth, jec_flth   !calculate domain
  integer :: ni, nj, nk
  integer, dimension(2) :: atm_layout, assim_layout
  integer :: nk_asm, nk_sice
  integer :: nk_inf = 1  ! x1y add for inflation levels
  integer :: filter_halo = 10
  integer :: assim_frequency = 6
  real    :: dt_const_ps = 1.0 !lulv 2020-05-24
  real    :: cov_inflate = 1.0
  real    :: corr_coef = 1.0 !lulv add for adjust dt_const in fv_physics.F90
  logical :: ida_up = .false.
  logical :: sequential_filter = .false.
!  logical :: hea_hfq = .true. !lulv 20200724
!  logical :: hea_day = .false. !lulv 20200724
!  logical :: hea_mon = .false. !lulv 20200724
  logical :: tra_eakf = .false. !lulv 20200606
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_hfq_ps0,ens_min_ps0,ens_day_ps0,ens_mon_ps0
! 20 ensb info
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_day_t0,ens_5day_t0,ens_6h_t0,ens_3h_t0
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_6h_ps0,ens_3h_ps0,ens_1h_ps0,ens_5day_ps0
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_day_u0,ens_5day_u0,ens_6h_u0,ens_3h_u0
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_day_v0,ens_5day_v0,ens_6h_v0,ens_3h_v0
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_day_q0,ens_5day_q0,ens_6h_q0,ens_3h_q0
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ps_ens_tmp,ps_ens_tmp_hour,ps_ens_tmp_day
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: t_ens_tmp,t_ens_tmp_hour,t_ens_tmp_day
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: u_ens_tmp,u_ens_tmp_hour,u_ens_tmp_day
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: v_ens_tmp,v_ens_tmp_hour,v_ens_tmp_day
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: q_ens_tmp,q_ens_tmp_hour,q_ens_tmp_day
  type(ada_field_type), dimension(60) :: ps_ens_tmp_min,t_ens_tmp_min,u_ens_tmp_min,v_ens_tmp_min,q_ens_tmp_min
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_hfq_t, ens_hfq_q !0621
  type(ada_field_type), dimension(MAX_ENSEMBLE_SIZE) :: ens_hfq_u, ens_hfq_v !0621
  type(ada_field_type) :: ens_hfq_ps,ens_day_ps,ens_min_ps !current step info
  type(ada_field_type) :: ens_6h_ps,ens_3h_ps,ens_1h_ps,ens_5day_ps
  type(ada_field_type) :: ens_day_t,ens_5day_t,ens_6h_t,ens_3h_t
  type(ada_field_type) :: ens_day_u,ens_5day_u,ens_6h_u,ens_3h_u
  type(ada_field_type) :: ens_day_v,ens_5day_v,ens_6h_v,ens_3h_v
  type(ada_field_type) :: ens_day_q,ens_5day_q,ens_6h_q,ens_3h_q
  type(ada_field_type) :: ens_6h_to,ens_3h_to,ens_day_to,ens_5day_to
  type(assim_data_type) :: Atmos_obs
!  real, dimension(:,:,:), allocatable, save :: corr_t, corr_q, corr_u, corr_v
!  real, dimension(:,:), allocatable, save :: corr_ps
!  real, dimension(:,:,:), allocatable, save :: corr_t_hfq, corr_q_hfq,corr_u_hfq, corr_v_hfq 
  real, dimension(:,:), allocatable, save :: corr_ps_hfq,corr_ps_day,corr_ps_min !corr_ps_mon
  real, dimension(:,:), allocatable, save :: corr_ps_6h,corr_ps_3h,corr_ps_1h,corr_ps_5day
  real, dimension(:,:,:), allocatable, save :: corr_t_day,corr_t_5day,corr_t_6h,corr_t_3h
  real, dimension(:,:,:), allocatable, save :: corr_u_day,corr_u_5day,corr_u_6h,corr_u_3h
  real, dimension(:,:,:), allocatable, save :: corr_v_day,corr_v_5day,corr_v_6h,corr_v_3h
  real, dimension(:,:,:), allocatable, save :: corr_q_day,corr_q_5day,corr_q_6h,corr_q_3h
  real, dimension(:,:,:), allocatable, save :: corr_to_6h,corr_to_3h,corr_to_day,corr_to_5day

  type(ada_grid_type) :: T_grid
  integer :: ensemble_size = 12

contains  
  
  subroutine init_fv_ada(grid,config_flags,Atm_time, nlon,mlat, & 
                        isd, ied, jsd, jed, kms, kme, &
                        isc_flth, iec_flth, jsc_flth, jec_flth, &
                        k_end,kde, atm_step, Var_xy,corr_ps_hfqo,corr_ps_mino,corr_ps_dayo, &
                        corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo, &
                        corr_t_dayo, corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo,&
                        corr_u_dayo, corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo,&
                        corr_v_dayo, corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo,&
                        corr_q_dayo, corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo,&
                       corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r,&
                       ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,ens_hfq_pso04,&
                       ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,& 
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,& 
                        ens_hfq_pso17,ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,&
                        hea_hfq_1,hea_hfq_min,hea_day,hea_5day,hea_6h,hea_3h,hea_1h)


    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    TYPE(WRFU_Time), intent(in), target        :: Atm_time
    integer, intent(in) :: k_end,kde,nlon,mlat
    integer, intent(in) :: isd, ied, jsd, jed, kms, kme !halo domain
    integer, intent(in) :: isc_flth, iec_flth, jsc_flth, jec_flth   !calculate domain
    integer :: m,i,j,nlev,atm_step,ii,jj
    type(assim_data_type), intent(in) :: Var_xy
    real,dimension(:,:), allocatable,intent(in) :: corr_ps_hfqo
    real,dimension(:,:), allocatable,intent(in) :: corr_ps_mino
    real,dimension(:,:), allocatable,intent(in) :: corr_ps_dayo,corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo
    real,dimension(:,:,:), allocatable,intent(in) :: corr_t_dayo,corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo
    real,dimension(:,:,:), allocatable,intent(in) :: corr_u_dayo,corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo
    real,dimension(:,:,:), allocatable,intent(in) :: corr_v_dayo,corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo
    real,dimension(:,:,:), allocatable,intent(in) :: corr_q_dayo,corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo
    real,dimension(:,:,:), allocatable,intent(in) :: corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso01,ens_hfq_pso02
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso03,ens_hfq_pso04
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso05,ens_hfq_pso06
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso07,ens_hfq_pso08
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso09,ens_hfq_pso10
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso11,ens_hfq_pso12
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso13,ens_hfq_pso14
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso15,ens_hfq_pso16
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso17,ens_hfq_pso18
   real, dimension(:,:,:), allocatable, intent(in) :: ens_hfq_pso19,ens_hfq_pso20
   logical,intent(in) :: hea_hfq_1,hea_hfq_min,hea_day,hea_5day,hea_6h,hea_3h,hea_1h
   real, pointer, dimension(:) :: lon => NULL()
   real, pointer, dimension(:) :: lat => NULL()
   real, pointer, dimension(:) :: lonu => NULL(),latu => NULL()
   real, pointer, dimension(:) :: lonv => NULL(),latv => NULL()
   character(len=200) :: file_obs1,file_obs2,file_obs3
   integer :: ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7,ierr8,varid,ierr9,ierr10,ierr11,ierr12,ncid2,ncid3
   integer :: ierr13,ierr14,ierr15,ierr16,ierr17,ierr18
   character(len=1) :: vari

       nlev = min (k_end, kde-1)
write(*,*) 'abb,nlev',nlev,k_end,kde-1
       ni = nlon
       nj = mlat
       nk = nlev

!write(*,*) 'nk',nk

           allocate(Atmos_obs%ps(1:nlon,1:mlat))
           allocate(Atmos_obs%phis(1:nlon,1:mlat))
           allocate(Atmos_obs%t(1:nlon,1:nk,1:mlat))
           allocate(Atmos_obs%u(1:nlon,1:nk,1:mlat))
           allocate(Atmos_obs%v(1:nlon,1:nk,1:mlat))
           allocate(Atmos_obs%u_2(1:nlon,1:nk,1:mlat))
           allocate(Atmos_obs%t_2(1:nlon,1:nk,1:mlat))

           allocate(T_grid%x(1:nlon),T_grid%y(1:mlat))
           allocate(T_grid%xu(1:nlon),T_grid%yu(1:mlat))
           allocate(T_grid%xv(1:nlon),T_grid%yv(1:mlat))
           allocate(lon(1:nlon))
           allocate(lat(1:mlat))
           allocate(lonu(1:nlon))
           allocate(latu(1:mlat))
           allocate(lonv(1:nlon))
           allocate(latv(1:mlat))
           
           Atmos_obs%ps(:,:) = 0.0
           Atmos_obs%phis(:,:) = 0.0
           Atmos_obs%t(:,:,:) = 0.0
           Atmos_obs%u(:,:,:) = 0.0
           Atmos_obs%v(:,:,:) = 0.0
           Atmos_obs%u_2(:,:,:) = 0.0
           Atmos_obs%t_2(:,:,:) = 0.0

           T_grid%x(:) = 0.0
           T_grid%y(:) = 0.0
           T_grid%xu(:) = 0.0
           T_grid%yu(:) = 0.0
           T_grid%xv(:) = 0.0
           T_grid%yv(:) = 0.0

           lon(:) = 0.0
           lat(:) = 0.0
        lonu(:) = 0.0
        latu(:) = 0.0
        lonv(:) = 0.0
        latu(:) = 0.0


     do m=1, heaada_enssize
             allocate(ens_hfq_t(m)%data(isd:ied,jsd:jed,nlev))
             allocate(ens_hfq_q(m)%data(isd:ied,jsd:jed,nlev))
             allocate(ens_hfq_u(m)%data(isd:ied,jsd:jed,nlev))
             allocate(ens_hfq_v(m)%data(isd:ied,jsd:jed,nlev))
             ens_hfq_t(m)%data(isd:ied,jsd:jed,nlev)= 0.0
             ens_hfq_q(m)%data(isd:ied,jsd:jed,nlev)= 0.0
             ens_hfq_u(m)%data(isd:ied,jsd:jed,nlev)= 0.0
             ens_hfq_v(m)%data(isd:ied,jsd:jed,nlev)= 0.0

        allocate(ens_hfq_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_hfq_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_min_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_min_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_day_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_day_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_6h_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_6h_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_5day_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_5day_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_3h_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_3h_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_1h_ps0(m)%data(isd:ied,jsd:jed,1))
        ens_1h_ps0(m)%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_day_t0(m)%data(isd:ied,1:27,jsd:jed))
        ens_day_t0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_5day_t0(m)%data(isd:ied,1:27,jsd:jed))
        ens_5day_t0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_6h_t0(m)%data(isd:ied,1:27,jsd:jed))
        ens_6h_t0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_3h_t0(m)%data(isd:ied,1:27,jsd:jed))
        ens_3h_t0(m)%data(isd:ied,1:27,jsd:jed) = 0.0




        allocate(ens_day_u0(m)%data(isd:ied,1:27,jsd:jed))
        ens_day_u0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_5day_u0(m)%data(isd:ied,1:27,jsd:jed))
        ens_5day_u0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_6h_u0(m)%data(isd:ied,1:27,jsd:jed))
        ens_6h_u0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_3h_u0(m)%data(isd:ied,1:27,jsd:jed))
        ens_3h_u0(m)%data(isd:ied,1:27,jsd:jed) = 0.0


        allocate(ens_day_v0(m)%data(isd:ied,1:27,jsd:jed))
        ens_day_v0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_5day_v0(m)%data(isd:ied,1:27,jsd:jed))
        ens_5day_v0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_6h_v0(m)%data(isd:ied,1:27,jsd:jed))
        ens_6h_v0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_3h_v0(m)%data(isd:ied,1:27,jsd:jed))
        ens_3h_v0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_day_q0(m)%data(isd:ied,1:27,jsd:jed))
        ens_day_q0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_5day_q0(m)%data(isd:ied,1:27,jsd:jed))
        ens_5day_q0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_6h_q0(m)%data(isd:ied,1:27,jsd:jed))
        ens_6h_q0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

        allocate(ens_3h_q0(m)%data(isd:ied,1:27,jsd:jed))
        ens_3h_q0(m)%data(isd:ied,1:27,jsd:jed) = 0.0

     enddo

        allocate(ens_hfq_ps%data(isd:ied,jsd:jed,1))
        ens_hfq_ps%data(isd:ied,jsd:jed,1) = 0.0 

        allocate(ens_day_ps%data(isd:ied,jsd:jed,1))
        ens_day_ps%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_min_ps%data(isd:ied,jsd:jed,1))
        ens_min_ps%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_6h_ps%data(isd:ied,jsd:jed,1))
        ens_6h_ps%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_3h_ps%data(isd:ied,jsd:jed,1))
        ens_3h_ps%data(isd:ied,jsd:jed,1) = 0.0
        allocate(ens_1h_ps%data(isd:ied,jsd:jed,1))
        ens_1h_ps%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_5day_ps%data(isd:ied,jsd:jed,1))
        ens_5day_ps%data(isd:ied,jsd:jed,1) = 0.0

        allocate(ens_day_t%data(isd:ied,1:nlev,jsd:jed))
        ens_day_t%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_5day_t%data(isd:ied,1:nlev,jsd:jed))
        ens_5day_t%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_6h_t%data(isd:ied,1:nlev,jsd:jed))
        ens_6h_t%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_3h_t%data(isd:ied,1:nlev,jsd:jed))
        ens_3h_t%data(isd:ied,1:nlev,jsd:jed) = 0.0


!write(*,*) 'ens_6h_t,ini',ens_6h_t%data(100,2,200)

        allocate(ens_day_to%data(isd:ied,1:nlev,jsd:jed))
        ens_day_to%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_5day_to%data(isd:ied,1:nlev,jsd:jed))
        ens_5day_to%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_6h_to%data(isd:ied,1:nlev,jsd:jed))
        ens_6h_to%data(isd:ied,1:nlev,jsd:jed) = 0.0

!write(*,*) 'ens_6h_to,ini',ens_6h_to%data(100,2,200)

        allocate(ens_3h_to%data(isd:ied,1:nlev,jsd:jed))
        ens_3h_to%data(isd:ied,1:nlev,jsd:jed) = 0.0


        allocate(ens_day_u%data(isd:ied,1:nlev,jsd:jed))
        ens_day_u%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_5day_u%data(isd:ied,1:nlev,jsd:jed))
        ens_5day_u%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_6h_u%data(isd:ied,1:nlev,jsd:jed))
        ens_6h_u%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_3h_u%data(isd:ied,1:nlev,jsd:jed))
        ens_3h_u%data(isd:ied,1:nlev,jsd:jed) = 0.0


        allocate(ens_day_v%data(isd:ied,1:nlev,jsd:jed))
        ens_day_v%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_5day_v%data(isd:ied,1:nlev,jsd:jed))
        ens_5day_v%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_6h_v%data(isd:ied,1:nlev,jsd:jed))
        ens_6h_v%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_3h_v%data(isd:ied,1:nlev,jsd:jed))
        ens_3h_v%data(isd:ied,1:nlev,jsd:jed) = 0.0


        allocate(ens_day_q%data(isd:ied,1:nlev,jsd:jed))
        ens_day_q%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_5day_q%data(isd:ied,1:nlev,jsd:jed))
        ens_5day_q%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_6h_q%data(isd:ied,1:nlev,jsd:jed))
        ens_6h_q%data(isd:ied,1:nlev,jsd:jed) = 0.0

        allocate(ens_3h_q%data(isd:ied,1:nlev,jsd:jed))
        ens_3h_q%data(isd:ied,1:nlev,jsd:jed) = 0.0


        do m=1,heahfq_minintv !wangj 202109
           allocate(ps_ens_tmp(m)%data(isd:ied,jsd:jed,1:1))
           ps_ens_tmp(m)%data(isd:ied,jsd:jed,1:1)=0.0
        enddo

          do m=1,60 !for daily mean, wangjj 20211130
             allocate(ps_ens_tmp_min(m)%data(isd:ied,jsd:jed,1))
             ps_ens_tmp_min(m)%data(isd:ied,jsd:jed,1) = 0.0
             allocate(t_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev))
             t_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev) = 0.0
             allocate(u_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev))
             u_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev) = 0.0
             allocate(v_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev))
             v_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev) = 0.0
             allocate(q_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev))
             q_ens_tmp_min(m)%data(isd:ied,jsd:jed,1:nlev) = 0.0
          enddo

          do m=1,24 !for daily mean, wangjj 20211130
             allocate(ps_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1))
             ps_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1) =0.0
             allocate(t_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev))
             t_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
             allocate(u_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev))
             u_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
             allocate(v_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev))
             v_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
             allocate(q_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev))
             q_ens_tmp_hour(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
          enddo


          do m=1,24 !for daily mean, wangjj 20211130
             allocate(ps_ens_tmp_day(m)%data(isd:ied,jsd:jed,1))
             ps_ens_tmp_day(m)%data(isd:ied,jsd:jed,1) =0.0
             allocate(t_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev))
             t_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
             allocate(u_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev))
             u_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
             allocate(v_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev))
             v_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
             allocate(q_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev))
             q_ens_tmp_day(m)%data(isd:ied,jsd:jed,1:nlev) =0.0
          enddo

!          allocate(corr_t_hfq(isd:ied,jsd:jed,1:nlev),corr_q_hfq(isd:ied,jsd:jed,1:nlev))
!          allocate(corr_u_hfq(isd:ied,jsd:jed,1:nlev),corr_v_hfq(isd:ied,jsd:jed,1:nlev))
          allocate(corr_ps_hfq(isd:ied,jsd:jed))
          allocate( corr_ps_day(isd:ied,jsd:jed))
          allocate( corr_ps_min(isd:ied,jsd:jed))
          allocate( corr_ps_6h(isd:ied,jsd:jed))
          allocate( corr_ps_3h(isd:ied,jsd:jed))
          allocate( corr_ps_1h(isd:ied,jsd:jed))
          allocate( corr_ps_5day(isd:ied,jsd:jed))

          allocate(corr_t_day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_t_5day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_t_6h(isd:ied,1:nlev,jsd:jed))
          allocate(corr_t_3h(isd:ied,1:nlev,jsd:jed))

          allocate(corr_to_day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_to_5day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_to_6h(isd:ied,1:nlev,jsd:jed))
          allocate(corr_to_3h(isd:ied,1:nlev,jsd:jed))


          allocate(corr_u_day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_u_5day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_u_6h(isd:ied,1:nlev,jsd:jed))
          allocate(corr_u_3h(isd:ied,1:nlev,jsd:jed))

          allocate(corr_v_day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_v_5day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_v_6h(isd:ied,1:nlev,jsd:jed))
          allocate(corr_v_3h(isd:ied,1:nlev,jsd:jed))

          allocate(corr_q_day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_q_5day(isd:ied,1:nlev,jsd:jed))
          allocate(corr_q_6h(isd:ied,1:nlev,jsd:jed))
          allocate(corr_q_3h(isd:ied,1:nlev,jsd:jed))

          corr_ps_hfq(:,:) = 0.0
          corr_ps_day(:,:) = 0.0
          corr_ps_6h(:,:) = 0.0
          corr_ps_3h(:,:) = 0.0
          corr_ps_1h(:,:) = 0.0
          corr_ps_5day(:,:) = 0.0
          corr_ps_min(:,:) = 0.0
!deallocate(corr_ps_hfq)
         corr_t_day(:,:,:) = 0.0
         corr_t_5day(:,:,:) = 0.0
         corr_t_6h(:,:,:) = 0.0
         corr_t_3h(:,:,:) = 0.0

         corr_to_day(:,:,:) = 0.0
         corr_to_5day(:,:,:) = 0.0
         corr_to_6h(:,:,:) = 0.0
         corr_to_3h(:,:,:) = 0.0

         corr_u_day(:,:,:) = 0.0
         corr_u_5day(:,:,:) = 0.0
         corr_u_6h(:,:,:) = 0.0
         corr_u_3h(:,:,:) = 0.0

         corr_v_day(:,:,:) = 0.0
         corr_v_5day(:,:,:) = 0.0
         corr_v_6h(:,:,:) = 0.0
         corr_v_3h(:,:,:) = 0.0

         corr_q_day(:,:,:) = 0.0
         corr_q_5day(:,:,:) = 0.0
         corr_q_6h(:,:,:) = 0.0
         corr_q_3h(:,:,:) = 0.0
! read obs lat lon
! read psfc t p loc
        file_obs1="/public/home/hou/wangj/obs/era5_atmos_obs_234.nc"
        ierr1=nf90_open(trim(file_obs1),nf90_nowrite,ncid2)
        ierr2=nf90_inq_varid(ncid2,'lat',varid)
        ierr3=nf90_get_var(ncid2,varid,lat)
        ierr4=nf90_inq_varid(ncid2,'lon',varid)
        ierr5=nf90_get_var(ncid2,varid,lon)
        ierr6 = nf90_close(ncid2)

        T_grid%x = lon(:)
        T_grid%y = lat(:)

! read u v loc
        file_obs2="/public/home/hou/wangj/obs/interp_u409_to_410_1.nc"

        file_obs3='/public/home/hou/wangj/obs/interp_v559_to_560_1.nc'

        ierr7=nf90_open(trim(file_obs2),nf90_nowrite,ncid2)
        ierr8=nf90_inq_varid(ncid2,'lat',varid)
        ierr9=nf90_get_var(ncid2,varid,latu)
        ierr10=nf90_inq_varid(ncid2,'lon',varid)
        ierr11=nf90_get_var(ncid2,varid,lonu)
        ierr12 = nf90_close(ncid2)

        T_grid%xu = lonu(:)
        T_grid%yu = latu(:)

        ierr13=nf90_open(trim(file_obs3),nf90_nowrite,ncid2)
        ierr14=nf90_inq_varid(ncid2,'lat',varid)
        ierr15=nf90_get_var(ncid2,varid,latv)
        ierr16=nf90_inq_varid(ncid2,'lon',varid)
        ierr17=nf90_get_var(ncid2,varid,lonv)
        ierr18 = nf90_close(ncid2)

        T_grid%xv = lonv(:)
        T_grid%yv = latv(:)

!        write(*,*) 'p',ierr1,ierr2,ierr3,ierr4,ierr5,ierr6
!        write(*,*) 'u',ierr7,ierr8,ierr9,ierr10,ierr11,ierr12 
!        write(*,*) 'v',ierr13,ierr14,ierr15,ierr16,ierr17,ierr18
!open(30,file="/home/lmkli/wangj/test/WRF371/v3_0509/WRFV3_case12_13_7_4_2/lonu.txt")
!write(30,*) T_grid%xu
!open(31,file="/home/lmkli/wangj/test/WRF371/v3_0509/WRFV3_case12_13_7_4_2/latu.txt")
!write(31,*) T_grid%yu
!open(32,file="/home/lmkli/wangj/test/WRF371/v3_0509/WRFV3_case12_13_7_4_2/lonv.txt")
!write(32,*) T_grid%xv
!open(33,file="/home/lmkli/wangj/test/WRF371/v3_0509/WRFV3_case12_13_7_4_2/latv.txt")
!write(33,*) T_grid%yv


      if(hea_day)     call read_data_day(grid,ens_day_ps0,isd, ied,jsd,jed )
      if(hea_6h)      call read_data_6h(grid,ens_6h_ps0,isd, ied,jsd,jed )
      if(hea_3h)      call read_data_3h(grid,ens_3h_ps0,isd, ied,jsd,jed )
      if(hea_1h)      call read_data_1h(grid,ens_1h_ps0,isd, ied,jsd,jed )
      if(hea_5day)      call read_data_5day(grid,ens_5day_ps0,isd, ied,jsd,jed )

      if(config_flags%hea_tday .or. config_flags%hea_t5day .or. config_flags%hea_t6h .or. config_flags%hea_t3h)  vari='t'
      if(config_flags%hea_t3h)    call read_data_t3h_tuv(grid,ens_3h_t0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_t6h)    call read_data_t6h_tuv(grid,ens_6h_t0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_tday)   call read_data_day_tuv(grid,ens_day_t0,isd,ied,jsd,jed,nk,vari) 
      if(config_flags%hea_t5day)  call read_data_5day_tuv(grid,ens_5day_t0,isd,ied,jsd,jed,nk,vari)

      if(config_flags%hea_uday .or. config_flags%hea_u5day .or.config_flags%hea_u6h .or. config_flags%hea_u3h)  vari='u'
      if(config_flags%hea_u3h)    call read_data_t3h_tuv(grid,ens_3h_u0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_u6h)    call read_data_t6h_tuv(grid,ens_6h_u0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_uday)   call read_data_day_tuv(grid,ens_day_u0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_u5day)  call read_data_5day_tuv(grid,ens_5day_u0,isd,ied,jsd,jed,nk,vari)

      if(config_flags%hea_vday .or. config_flags%hea_v5day .or.config_flags%hea_v6h .or. config_flags%hea_v3h)  vari='v'
      if(config_flags%hea_v3h)    call  read_data_t3h_tuv(grid,ens_3h_v0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_v6h)    call  read_data_t6h_tuv(grid,ens_6h_v0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_vday)   call read_data_day_tuv(grid,ens_day_v0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_v5day)  call read_data_5day_tuv(grid,ens_5day_v0,isd,ied,jsd,jed,nk,vari)

 if(config_flags%hea_q6h .or. config_flags%hea_q3h .or. &
        config_flags%hea_q5day .or. config_flags%hea_qday ) vari='q'
write(*,*) 'begread'

      if(config_flags%hea_q3h)  call read_data_t3h_tuv(grid,ens_3h_q0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_q6h)  call read_data_t6h_tuv(grid,ens_6h_q0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_qday)  call read_data_day_tuv(grid,ens_day_q0,isd,ied,jsd,jed,nk,vari)
      if(config_flags%hea_q5day)  call read_data_5day_tuv(grid,ens_5day_q0,isd,ied,jsd,jed,nk,vari)


      write(*,*) 'init_day1',atm_step,isd+10,jsd+10,ens_day_q0(1)%data(isd+10,jsd+10,1)
      write(*,*) 'init_day20',atm_step,isd+10,jsd+10,ens_day_q0(20)%data(isd+10,jsd+10,1)
      write(*,*) 'init_3h1',atm_step,isd+10,jsd+10,ens_3h_q0(1)%data(isd+10,jsd+10,1)
      write(*,*) 'init_3h20',atm_step,isd+10,jsd+10,ens_3h_q0(20)%data(isd+10,jsd+10,1)


        call ensemble_filter_atm_init(ensemble_size,heaada_enssize)

  end subroutine init_fv_ada



  subroutine fv_ada(grid,config_flags,Atm_Time,Atm_var_g,corr_ps_hfqo,corr_ps_mino,corr_ps_dayo,&
                        corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo, &
                        corr_t_dayo, corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo,&
                        corr_u_dayo, corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo,&
                        corr_v_dayo, corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo,&
                        corr_q_dayo, corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo,&
                        atm_step, isd, ied, jsd, jed, kms, kme, &
                        isc_flth, iec_flth, jsc_flth, jec_flth,z_top, &
                        hea_hfq_1,hea_hfq_min,hea_day,hea_5day,hea_6h,hea_3h,hea_1h, &
                        dt_const_ps_hfq,dt_const_ps_day,dt_const_ps_mon, &
                       ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,ens_hfq_pso04,&
                       ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,& 
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,& 
                        ens_hfq_pso17,ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,var_t,var_u)
! no ak bk
     TYPE ( domain ), INTENT(INOUT) :: grid
     TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
  TYPE(WRFU_Time), intent(in), target        :: Atm_time
  TYPE(assim_data_type), intent(inout) :: Atm_var_g
  real, dimension(:,:), allocatable,intent(inout) :: corr_ps_hfqo
  real, dimension(:,:), allocatable,intent(inout) :: corr_ps_mino
  real, dimension(:,:), allocatable,intent(inout) :: corr_ps_dayo,corr_psfc6ho,corr_psfc3ho,corr_psfc1ho,corr_psfc5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_t_dayo,corr_t6ho,corr_t3ho,corr_t1ho,corr_t5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_u_dayo,corr_u6ho,corr_u3ho,corr_u1ho,corr_u5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_v_dayo,corr_v6ho,corr_v3ho,corr_v1ho,corr_v5dayo
    real,dimension(:,:,:), allocatable,intent(inout) :: corr_q_dayo,corr_q6ho,corr_q3ho,corr_q1ho,corr_q5dayo
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso07,ens_hfq_pso08,ens_hfq_pso09
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso16,ens_hfq_pso17,ens_hfq_pso18
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso19,ens_hfq_pso20
   real,dimension(isd:ied, kms:kme,jsd:jed),intent(inout) :: var_t,var_u
       character(len=1) :: vari

  integer, intent(in)  :: atm_step
  integer, intent(in)  :: isc_flth, iec_flth, jsc_flth, jec_flth
  integer,intent(in)   :: isd, ied, jsd, jed, kms, kme,z_top !halo domain
  real ::  gridi,gridj,lat1,lon1
  real :: p_lev(6),levm(6),tc(6),td(6),wind(6),wind_dir(6),sp_hum(6)
  INTEGER                         :: yr, month, day, hr, minute, sec, rc
    integer :: m,n, k, i, j
    integer :: num_obs, num_obs_sice_up
    integer :: second0, day0
    integer :: is, ie, hs, he, vs, ve
    real    :: inflate
    integer :: i_ll, j_ll , k_ll !lulv
    integer :: i_t,j_t,k_t
    logical :: used
    logical,intent(in) :: hea_hfq_1,hea_hfq_min,hea_day,hea_5day,hea_6h,hea_3h,hea_1h
    real,intent(in) :: dt_const_ps_hfq,dt_const_ps_day,dt_const_ps_mon
!    integer, parameter :: ensemble_min_day=1440

!    ied  = min( ied,559 )
!    jed   = min( jed,409 )

     CALL WRFU_TimeGet( Atm_Time, YY=yr, dayOfYear=day, H=hr, M=minute,S=sec, rc=rc)

!write(*,*) 'ada_time',atm_step,yr,day,hr,minute

    CALL wrf_get_myproc( myproc )
    CALL wrf_get_nproc( nproc )

! ------- collect and update ps  ens
     if(hea_day)  call mk_ens_hf_ada_day2(config_flags,Atm_var_g,ens_day_ps0,atm_step,isd,ied,jsd,jed)
     if(hea_6h)   call mk_ens_ada_6h2(Atm_var_g,ens_6h_ps0,atm_step,isd,ied, jsd, jed)
     if(hea_3h)   call mk_ens_ada_3h2(Atm_var_g,ens_3h_ps0,atm_step,isd,ied, jsd, jed)
     if(hea_5day) call mk_ens_ada_5day2(Atm_var_g,ens_5day_ps0,atm_step,isd,ied,jsd,jed)

!  ------- collect and update t ens
  if(config_flags%hea_tday .or. config_flags%hea_t5day .or. config_flags%hea_t6h .or. config_flags%hea_t3h) vari='t'
  if(config_flags%hea_t3h)   call mk_ens_ada_3ht2(config_flags,Atm_var_g,ens_3h_t0,atm_step,isd,ied,jsd,jed,z_top,vari)
  if(config_flags%hea_t6h)   call mk_ens_ada_6ht2(config_flags,Atm_var_g,ens_6h_t0,atm_step,isd,ied,jsd,jed,z_top,vari)
  if(config_flags%hea_tday)  call mk_ens_ada_dayt2(config_flags,Atm_var_g,ens_day_t0,atm_step,isd,ied,jsd,jed,z_top,vari)
  if(config_flags%hea_t5day) call mk_ens_ada_5dayt2(config_flags,Atm_var_g,ens_5day_t0,atm_step,isd,ied,jsd,jed,z_top,vari)

!  ------- collect and update u ens
  if(config_flags%hea_u3h)   call mk_ens_ada_3hu2(config_flags,Atm_var_g,ens_3h_u0,atm_step,isd,ied,jsd,jed,z_top)
  if(config_flags%hea_u6h)   call mk_ens_ada_6hu2(config_flags,Atm_var_g,ens_6h_u0,atm_step,isd,ied,jsd,jed,z_top)
  if(config_flags%hea_uday)  call mk_ens_ada_dayu2(config_flags,Atm_var_g,ens_day_u0,atm_step,isd,ied,jsd,jed,z_top)
  if(config_flags%hea_u5day) call mk_ens_ada_5dayu2(config_flags,Atm_var_g,ens_5day_u0,atm_step,isd,ied,jsd,jed,z_top)

!  ------- collect and update v ens

  if(config_flags%hea_v3h)   call mk_ens_ada_3hv2(config_flags,Atm_var_g,ens_3h_v0,atm_step,isd,ied,jsd,jed,z_top)
  if(config_flags%hea_v6h)   call mk_ens_ada_6hv2(config_flags,Atm_var_g,ens_6h_v0,atm_step,isd,ied,jsd,jed,z_top)
  if(config_flags%hea_vday)  call mk_ens_ada_dayv2(config_flags,Atm_var_g,ens_day_v0,atm_step,isd,ied,jsd,jed,z_top)
  if(config_flags%hea_v5day) call mk_ens_ada_5dayv2(config_flags,Atm_var_g,ens_5day_v0,atm_step,isd,ied,jsd,jed,z_top)


!
   if (config_flags%debug_t) then !write ens1-20 for t.get_a_to_b_tuv(a,b1,b2...). pass the value of a to b1,b2...

        if(config_flags%hea_t3h) call get_a_to_b_tuv(ens_3h_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)

        if(config_flags%hea_t6h) call get_a_to_b_tuv(ens_6h_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)


        if(config_flags%hea_tday) call get_a_to_b_tuv(ens_day_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)
        if(config_flags%hea_t5day) call get_a_to_b_tuv(ens_5day_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)
   endif






        if(   atm_step==1 .or. mod(atm_step,360)==0   ) then
              call get_ncep_field(grid,config_flags,Atm_Time, Atmos_obs, assim_frequency,atm_step)
!               call get_obs_real(grid,config_flags,Atm_Time, Atmos_obs,assim_frequency,atm_step,gridi,gridj,p_lev,levm,tc,td,wind,wind_dir,sp_hum)


                if(hea_6h .or. config_flags%hea_t6h .or. config_flags%hea_u6h .or. config_flags%hea_v6h)  then


                     do j_t = jsd,jed
                     do i_t = isd,ied
                         corr_ps_6h(i_t,j_t)        =  Atm_var_g%ps(i_t,j_t)
                         ens_6h_ps%data(i_t,j_t,1) =  Atm_var_g%ps(i_t,j_t)
                     enddo
                     enddo

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_t_6h(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t)
                        ens_6h_t%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                        corr_u_6h(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        ens_6h_u%data(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        corr_v_6h(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                        ens_6h_v%data(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo


                   call ensemble_filter_atm_6h(grid,config_flags,ens_6h_t0(1:heaada_enssize),ens_6h_t, &
                   & ens_hfq_q(1:heaada_enssize),ens_6h_u0(1:ensemble_size),ens_6h_u,&
                   & ens_6h_v0(1:heaada_enssize),ens_6h_v, &
                   & ens_6h_ps0(1:heaada_enssize),ens_6h_ps,&
                   & Atm_var_g, Atmos_obs,isc_flth,iec_flth,jsc_flth,jec_flth,filter_halo, ni, nj, nk, &
                   &  T_grid, Atm_time,atm_step,isd, ied, jsd,jed) !no akbk

                     do j_ll = jsc_flth,jec_flth
                     do i_ll = isc_flth,iec_flth
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_ps_6h(i_ll,j_ll) = (ens_6h_ps%data(i_ll,j_ll,1)-corr_ps_6h(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        else
                         corr_ps_6h(i_ll,j_ll) =  (ens_6h_ps%data(i_ll,j_ll,1)-corr_ps_6h(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)*&
                                                & cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        endif

                        corr_psfc6ho(i_ll,j_ll) = corr_ps_6h(i_ll,j_ll)

                    enddo
                    enddo

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_t_6h(i_t,k_t,j_t) = (ens_6h_t%data(i_t,k_t,j_t)-corr_t_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_u_6h(i_t,k_t,j_t) = (ens_6h_u%data(i_t,k_t,j_t)-corr_u_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_v_6h(i_t,k_t,j_t) = (ens_6h_v%data(i_t,k_t,j_t)-corr_v_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         corr_t_6h(i_t,k_t,j_t)  = (ens_6h_t%data(i_t,k_t,j_t)-corr_t_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                         corr_u_6h(i_t,k_t,j_t)  = (ens_6h_u%data(i_t,k_t,j_t)-corr_u_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                         corr_v_6h(i_t,k_t,j_t)  =  (ens_6h_v%data(i_t,k_t,j_t)-corr_v_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        endif

                        corr_t6ho(i_t,k_t,j_t) = corr_t_6h(i_t,k_t,j_t)
                        corr_u6ho(i_t,k_t,j_t) = corr_u_6h(i_t,k_t,j_t)
                        corr_v6ho(i_t,k_t,j_t) = corr_v_6h(i_t,k_t,j_t)

                     enddo
                     enddo
                     enddo



          endif !hea_6h

                if(hea_3h .or. config_flags%hea_t3h .or. config_flags%hea_u3h .or. config_flags%hea_v3h)  then

!        write(*,*) 'running hea_3h ada_driver'
                     do j_t = jsd,jed
                     do i_t = isd,ied
                         corr_ps_3h(i_t,j_t)        =  Atm_var_g%ps(i_t,j_t)
                         ens_3h_ps%data(i_t,j_t,1) =  Atm_var_g%ps(i_t,j_t)
                     enddo
                     enddo


                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_t_3h(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t)
                        ens_3h_t%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                        corr_u_3h(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        ens_3h_u%data(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        corr_v_3h(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                        ens_3h_v%data(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo


                   call ensemble_filter_atm_3h(grid,config_flags,ens_3h_t0(1:heaada_enssize),ens_3h_t,&
                   & ens_hfq_q(1:heaada_enssize),ens_3h_u0(1:ensemble_size),ens_3h_u, &
                   & ens_3h_v0(1:heaada_enssize),ens_3h_v,&
                   & ens_3h_ps0(1:heaada_enssize),ens_3h_ps,&
                   & Atm_var_g,Atmos_obs,isc_flth,iec_flth,jsc_flth,jec_flth,filter_halo, ni, nj, nk, &
                   &  T_grid, Atm_time,atm_step,isd, ied, jsd,jed) !no akbk

                     do j_ll = jsc_flth,jec_flth
                     do i_ll = isc_flth,iec_flth
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_ps_3h(i_ll,j_ll) = (ens_3h_ps%data(i_ll,j_ll,1)-corr_ps_3h(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        else
                         corr_ps_3h(i_ll,j_ll) = (ens_3h_ps%data(i_ll,j_ll,1)-corr_ps_3h(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)*&
                                                & cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        endif
                        corr_psfc3ho(i_ll,j_ll) = corr_ps_3h(i_ll,j_ll)
                     enddo
                     enddo

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_t_3h(i_t,k_t,j_t) =  (ens_3h_t%data(i_t,k_t,j_t)-corr_t_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_u_3h(i_t,k_t,j_t) =  (ens_3h_u%data(i_t,k_t,j_t)-corr_u_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_v_3h(i_t,k_t,j_t) =  (ens_3h_v%data(i_t,k_t,j_t)-corr_v_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                          corr_t_3h(i_t,k_t,j_t)  = (ens_3h_t%data(i_t,k_t,j_t)-corr_t_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_u_3h(i_t,k_t,j_t)  = (ens_3h_u%data(i_t,k_t,j_t)-corr_u_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_v_3h(i_t,k_t,j_t)  = (ens_3h_v%data(i_t,k_t,j_t)-corr_v_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        endif

                        corr_t3ho(i_t,k_t,j_t) = corr_t_3h(i_t,k_t,j_t)
                        corr_u3ho(i_t,k_t,j_t) = corr_u_3h(i_t,k_t,j_t)
                        corr_v3ho(i_t,k_t,j_t) = corr_v_3h(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo


                endif !hea_3h

!                if(hea_1h)  then

!        write(*,*) 'running hea_1h ada_driver'
!                     do j_t = jsd,jed
!                     do i_t = isd,ied
!                         corr_ps_1h(i_t,j_t)        =  Atm_var_g%ps(i_t,j_t)
!                         ens_1h_ps%data(i_t,j_t,1) =  Atm_var_g%ps(i_t,j_t)
!                     enddo
!                     enddo
!             call  ensemble_filter_atm_1h(grid,config_flags,ens_hfq_t(1:heaada_enssize),ens_hfq_q(1:heaada_enssize),&
!             & ens_hfq_u(1:ensemble_size),ens_hfq_v(1:heaada_enssize),ens_1h_ps0(1:heaada_enssize),ens_1h_ps,&
!             &  Atm_var_g,Atmos_obs,isc_flth,iec_flth,jsc_flth,jec_flth,filter_halo, ni, nj, nk,&
!                   &  T_grid, Atm_time,atm_step,isd, ied, jsd,jed) !no akbk
!
!                     do j_ll = jsc_flth,jec_flth
!                     do i_ll = isc_flth,iec_flth
!                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
!                          corr_ps_1h(i_ll,j_ll) =  (ens_1h_ps%data(i_ll,j_ll,1)-corr_ps_1h(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
!                        else
!                         corr_ps_1h(i_ll,j_ll) = (ens_1h_ps%data(i_ll,j_ll,1)-corr_ps_1h(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)*&
!                                                &cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
!                        endif
!                        corr_psfc1ho(i_ll,j_ll) = corr_ps_1h(i_ll,j_ll)
!                     enddo
!                     enddo
!                endif !hea_1h




                if(hea_day .or. config_flags%hea_tday .or. config_flags%hea_uday .or. config_flags%hea_vday ) then
!write(*,*) 'running hea_day ada_driver'

                     do j_t = jsd,jed
                     do i_t = isd,ied
                         corr_ps_day(i_t,j_t) = Atm_var_g%ps(i_t,j_t)
                         ens_day_ps%data(i_t,j_t,1) = Atm_var_g%ps(i_t,j_t)
                         
                     enddo
                     enddo
!write(*,*) 'ada5'

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_t_day(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t) 
                        ens_day_t%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t) 
!                        corr_t_day(i_t,k_t,j_t) = var_t(i_t,k_t,j_t)
!                        ens_day_t%data(i_t,k_t,j_t)  = var_t(i_t,k_t,j_t)
!                        ens_day_u%data(i_t,k_t,j_t)  = var_u(i_t,k_t,j_t)
                        corr_u_day(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        ens_day_u%data(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        corr_v_day(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                        ens_day_v%data(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)

                     enddo
                     enddo
                     enddo
!write(*,*) 'ada6'



          call ensemble_filter_atm_day(grid,config_flags,ens_day_t0(1:heaada_enssize),ens_day_t,&
                        & ens_hfq_q(1:heaada_enssize), ens_day_u0(1:ensemble_size),ens_day_u,&
                        & ens_day_v0(1:heaada_enssize),ens_day_v, &
                 & ens_day_ps0(1:heaada_enssize),ens_day_ps,&
                 &  Atm_var_g, Atmos_obs, isc_flth,iec_flth,jsc_flth,jec_flth,filter_halo, ni, nj, nk, &
                   &  T_grid, Atm_time,atm_step,isd, ied, jsd,jed) !no akbk
!write(*,*) 'ada7'

!                     do j_ll = jsc_flth,jec_flth
!                     do i_ll = isc_flth,iec_flth
                do j_ll = jsd,jed
                do i_ll = isd,ied


                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_ps_day(i_ll,j_ll) = (ens_day_ps%data(i_ll,j_ll,1)-corr_ps_day(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        else
                         corr_ps_day(i_ll,j_ll) = (ens_day_ps%data(i_ll,j_ll,1)-corr_ps_day(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)*&
                                                & cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        endif

                        corr_ps_dayo(i_ll,j_ll) = corr_ps_day(i_ll,j_ll)

                    enddo
                    enddo
!write(*,*) 'ada8'
                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_t_day(i_t,k_t,j_t) =  (ens_day_t%data(i_t,k_t,j_t)-corr_t_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_u_day(i_t,k_t,j_t) =  (ens_day_u%data(i_t,k_t,j_t)-corr_u_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_v_day(i_t,k_t,j_t) =  (ens_day_v%data(i_t,k_t,j_t)-corr_v_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         corr_t_day(i_t,k_t,j_t) = (ens_day_t%data(i_t,k_t,j_t)-corr_t_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)

                         corr_u_day(i_t,k_t,j_t) = (ens_day_u%data(i_t,k_t,j_t)-corr_u_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)

                         corr_v_day(i_t,k_t,j_t) = (ens_day_v%data(i_t,k_t,j_t)-corr_v_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        endif

                        corr_t_dayo(i_t,k_t,j_t) = corr_t_day(i_t,k_t,j_t)
                        corr_u_dayo(i_t,k_t,j_t) = corr_u_day(i_t,k_t,j_t)
                        corr_v_dayo(i_t,k_t,j_t) = corr_v_day(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo
!write(*,*) 'ada9'


!        do j_ll = jsd,jed
!        do i_ll = isd,ied

!          write(*,'("ensend,corrend",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xf18.2,1xf18.2)')&
!           atm_step,myproc,i_ll,j_ll,T_grid%x(i_ll),T_grid%y(j_ll),ens_day_ps%data(i_ll,j_ll,1),corr_ps_day(i_ll,j_ll)

!        enddo
!        enddo


                endif !hea_day



                if(hea_hfq_1) then
write(*,*) 'running hea_hfq_1 ada_driver'

                     do j_t = jsd,jed
                     do i_t = isd,ied
                         corr_ps_hfq(i_t,j_t) = Atm_var_g%ps(i_t,j_t)
                         ens_hfq_ps%data(i_t,j_t,1) = Atm_var_g%ps(i_t,j_t)
                     enddo
                     enddo

                    call ensemble_filter_atm_hfq(grid,config_flags,ens_hfq_t(1:heaada_enssize),ens_hfq_q(1:heaada_enssize),&
                     & ens_hfq_u(1:ensemble_size),ens_hfq_v(1:heaada_enssize),ens_hfq_ps0(1:heaada_enssize),ens_hfq_ps,&
                     & Atm_var_g, Atmos_obs, isc_flth, iec_flth, jsc_flth,jec_flth,filter_halo, ni, nj, nk, &
                     &  T_grid, Atm_time,atm_step, & 
                     & isd, ied, jsd, jed) !no akbk


                     do j_ll = jsc_flth,jec_flth
                     do i_ll = isc_flth,iec_flth

!                        if(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG .le. 80) then !hudu to du

                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_ps_hfq(i_ll,j_ll) = (ens_hfq_ps%data(i_ll,j_ll,1) -corr_ps_hfq(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        else
                        corr_ps_hfq(i_ll,j_ll) = (ens_hfq_ps%data(i_ll,j_ll,1) -corr_ps_hfq(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        endif
                        
                        corr_ps_hfqo(i_ll,j_ll) = corr_ps_hfq(i_ll,j_ll)

                    enddo
                    enddo


             endif !hea_hfq
!         endif
!         endif 


        if(hea_5day .or. config_flags%hea_t5day .or. config_flags%hea_u5day  .or. config_flags%hea_v5day) then
!write(*,*) 'running 5day ada_driver'


                     do j_t = jsd,jed
                     do i_t = isd,ied
                         corr_ps_5day(i_t,j_t) = Atm_var_g%ps(i_t,j_t)
                         ens_5day_ps%data(i_t,j_t,1) = Atm_var_g%ps(i_t,j_t)
                     enddo
                     enddo


                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_t_5day(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t)
                        ens_5day_t%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                        corr_u_5day(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        ens_5day_u%data(i_t,k_t,j_t)  = Atm_var_g%u(i_t,k_t,j_t)
                        corr_v_5day(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                        ens_5day_v%data(i_t,k_t,j_t)  = Atm_var_g%v(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo

              call ensemble_filter_atm_5day(grid,config_flags,ens_5day_t0(1:heaada_enssize),ens_5day_t,&
                   & ens_hfq_q(1:heaada_enssize),ens_5day_u0(1:ensemble_size),ens_5day_u,&
                   & ens_5day_v0(1:heaada_enssize),ens_5day_v,ens_5day_ps0(1:heaada_enssize),ens_5day_ps,&
                   & Atm_var_g,Atmos_obs,isc_flth,iec_flth,jsc_flth,jec_flth,filter_halo, ni, nj, nk, &
                   &  T_grid, Atm_time,atm_step,isd, ied, jsd,jed) !no akbk

                     do j_ll = jsc_flth,jec_flth
                     do i_ll = isc_flth,iec_flth
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                     corr_ps_5day(i_ll,j_ll) =  (ens_5day_ps%data(i_ll,j_ll,1)-corr_ps_5day(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        else
                     corr_ps_5day(i_ll,j_ll) = (ens_5day_ps%data(i_ll,j_ll,1)-corr_ps_5day(i_ll,j_ll))*cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)*&
                                                & cos(Atm_var_g%y(i_ll,j_ll)*RAD_TO_DEG)
                        endif
                        corr_psfc5dayo(i_ll,j_ll) = corr_ps_5day(i_ll,j_ll)
                     enddo
                     enddo

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        if(Atm_var_g%y(i_ll,j_ll) .le. 80) then
                          corr_t_5day(i_t,k_t,j_t) = (ens_5day_t%data(i_t,k_t,j_t)-corr_t_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_u_5day(i_t,k_t,j_t)  = (ens_5day_u%data(i_t,k_t,j_t)-corr_u_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                          corr_v_5day(i_t,k_t,j_t)  = (ens_5day_v%data(i_t,k_t,j_t)-corr_v_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         corr_t_5day(i_t,k_t,j_t) = (ens_5day_t%data(i_t,k_t,j_t)-corr_t_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                         corr_u_5day(i_t,k_t,j_t) = (ens_5day_u%data(i_t,k_t,j_t)-corr_u_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                         corr_v_5day(i_t,k_t,j_t) = (ens_5day_v%data(i_t,k_t,j_t)-corr_v_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        endif
                        corr_t5dayo(i_t,k_t,j_t) = corr_t_5day(i_t,k_t,j_t)
                        corr_u5dayo(i_t,k_t,j_t) = corr_u_5day(i_t,k_t,j_t)
                        corr_v5dayo(i_t,k_t,j_t) = corr_v_5day(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo

                endif 



        endif !if(  mod(atm_step,360)==0   ) 

  end subroutine fv_ada  


     subroutine mda_o(grid,config_flags,Atm_Time,Atm_var_g,&
                                corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r,&
                                atm_step,isd, ied, jsd, jed, kms, kme,&
                                 its, ite, jts, jte, kts, kte, z_top,&
                                ids,  ide,  jds,  jde,  kds,  kde,&
                       ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,ens_hfq_pso04,&
                       ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,&
                        ens_hfq_pso17,ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20)

     TYPE ( domain ), INTENT(INOUT) :: grid
     TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
     TYPE(assim_data_type), intent(inout) :: Atm_var_g
  real,dimension(:,:,:), allocatable,intent(inout) :: corr_t_dayo_r,corr_t6ho_r,corr_t3ho_r,corr_t5dayo_r
  INTEGER,intent(in) :: atm_step,isd, ied, jsd, jed, kms, kme,&
                         its, ite, jts, jte, kts, kte, z_top, &
                        ids,  ide,  jds,  jde,  kds,  kde
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso07,ens_hfq_pso08,ens_hfq_pso09
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso16,ens_hfq_pso17,ens_hfq_pso18
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso19,ens_hfq_pso20
    TYPE(WRFU_Time), intent(in), target        ::  Atm_Time
   integer :: j_t,i_t,k_t
  real ::  gridi,gridj,lat1,lon1
  real :: p_lev(6),levm(6),tc(6),td(6),wind(6),wind_dir(6),sp_hum(6)
  real :: z_o2(6),errf(6)
  real,parameter :: dx=37.037
  real :: temp_q,temp_t,temp_u,temp_v

    if(config_flags%hea_ot5day)  call mk_ens_ada_5dayt2r(config_flags,Atm_var_g,ens_5day_t0,atm_step,isd,ied,jsd,jed,z_top) 
    if(config_flags%hea_otday)  call mk_ens_ada_dayt2r(config_flags,Atm_var_g,ens_day_t0,atm_step,isd,ied,jsd,jed,z_top)
    if(config_flags%hea_ot6h)  call mk_ens_ada_6ht2r(config_flags,Atm_var_g,ens_6h_t0,atm_step,isd,ied,jsd,jed,z_top)
    if(config_flags%hea_ot3h)  call mk_ens_ada_3ht2r(config_flags,Atm_var_g,ens_3h_t0,atm_step,isd,ied,jsd,jed,z_top)

    if (config_flags%debug_ot) then

        if(config_flags%hea_ot5day) call get_a_to_b_tuv(ens_5day_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                       ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)

        if(config_flags%hea_otday) call get_a_to_b_tuv(ens_day_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)

        if(config_flags%hea_ot6h) call get_a_to_b_tuv(ens_6h_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)

        if(config_flags%hea_ot3h) call get_a_to_b_tuv(ens_3h_t0,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,z_top)
!

   endif


   if(   atm_step==1 .or. mod(atm_step,360)==0   ) then
!GET REAL OBS
               call get_obs_real(grid,config_flags,Atm_Time, &
                Atmos_obs,assim_frequency,atm_step,gridi,gridj,p_lev,levm,tc,td,wind,wind_dir,sp_hum,lat1,lon1)


                call  get_z_o(grid%pb,grid%p_top,grid%p,&
                        gridi,gridj,p_lev,levm,isd, ied, jsd, jed, kms, kme,z_o2)

                call get_obs_errf(grid,config_flags,&
                        Atmos_obs,gridi,gridj,p_lev,levm,tc,td,wind,wind_dir,sp_hum,&
                        isd, ied, jsd, jed, kms, kme,z_o2,errf)


              if( config_flags%hea_ot6h   ) then
!
                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_to_6h(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t) !ok jifen
                        ens_6h_to%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo


             call ensemble_filter_atm_6h_r(grid,config_flags,ens_6h_q0(1:heaada_enssize),ens_6h_q,&
                   &ens_6h_t0(1:heaada_enssize),ens_6h_to,ens_6h_u0(1:ensemble_size),ens_6h_u,&
                   &ens_6h_v0(1:heaada_enssize),ens_6h_v,ens_6h_ps0(1:heaada_enssize),ens_6h_ps,&
                   &Atm_var_g,gridi, gridj, z_o2,errf,lat1,lon1,&
                   &  T_grid, Atm_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts,jte,kts,kte,&
                        ids,  ide,  jds,  jde,  kds,  kde, dx)

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied

                        if(Atm_var_g%y(i_t,j_t) .le. 80) then
                         temp_t = (ens_6h_to%data(i_t,k_t,j_t)-corr_to_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         temp_t = (ens_6h_to%data(i_t,k_t,j_t)-corr_to_6h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)

                        endif
                       corr_t6ho_r(i_t,k_t,j_t) =  temp_t

                     enddo
                     enddo
                     enddo


                endif !if( hea_t6h   ) then

              if( config_flags%hea_ot3h   ) then
!
                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_to_3h(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t) !ok jifen
                        ens_3h_to%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo

             call ensemble_filter_atm_3h_r(grid,config_flags,ens_3h_q0(1:heaada_enssize),ens_3h_q,&
                   &ens_3h_t0(1:heaada_enssize),ens_3h_to,ens_3h_u0(1:ensemble_size),ens_3h_u,&
                   &ens_3h_v0(1:heaada_enssize),ens_3h_v,ens_3h_ps0(1:heaada_enssize),ens_3h_ps,&
                   &Atm_var_g,gridi, gridj, z_o2,errf,lat1,lon1,&
                   &  T_grid, Atm_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts,jte,kts,kte,&
                        ids,  ide,  jds,  jde,  kds,  kde, dx)

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied

                        if(Atm_var_g%y(i_t,j_t) .le. 80) then
                         temp_t = (ens_3h_to%data(i_t,k_t,j_t)-corr_to_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         temp_t = (ens_3h_to%data(i_t,k_t,j_t)-corr_to_3h(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)

                        endif
                       corr_t3ho_r(i_t,k_t,j_t) =  temp_t

                     enddo
                     enddo
                     enddo

         endif ! hea_3h


              if( config_flags%hea_otday   ) then
!
                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_to_day(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t) 
                        ens_day_to%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo

             call ensemble_filter_atm_day_r(grid,config_flags,ens_day_q0(1:heaada_enssize),ens_day_q,&
                   &ens_day_t0(1:heaada_enssize),ens_day_to,ens_day_u0(1:ensemble_size),ens_day_u,&
                   &ens_day_v0(1:heaada_enssize),ens_day_v,ens_day_ps0(1:heaada_enssize),ens_day_ps,&
                   &Atm_var_g,gridi, gridj, z_o2,errf,lat1,lon1,&
                   &  T_grid, Atm_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts,jte,kts,kte,&
                        ids,  ide,  jds,  jde,  kds,  kde, dx)

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied

                        if(Atm_var_g%y(i_t,j_t) .le. 80) then
                         temp_t = (ens_day_to%data(i_t,k_t,j_t)-corr_to_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         temp_t = (ens_day_to%data(i_t,k_t,j_t)-corr_to_day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)

                        endif
                       corr_t_dayo_r(i_t,k_t,j_t) =  temp_t

                     enddo
                     enddo
                     enddo

         endif ! hea_day


              if( config_flags%hea_ot5day   ) then
!
                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied
                        corr_to_5day(i_t,k_t,j_t) = Atm_var_g%t(i_t,k_t,j_t)
                        ens_5day_to%data(i_t,k_t,j_t)  = Atm_var_g%t(i_t,k_t,j_t)
                     enddo
                     enddo
                     enddo

             call ensemble_filter_atm_5day_r(grid,config_flags,ens_5day_q0(1:heaada_enssize),ens_5day_q,&
                   &ens_5day_t0(1:heaada_enssize),ens_5day_to,ens_5day_u0(1:ensemble_size),ens_5day_u,&
                   &ens_5day_v0(1:heaada_enssize),ens_5day_v,ens_5day_ps0(1:heaada_enssize),ens_5day_ps,&
                   &Atm_var_g,gridi, gridj, z_o2,errf,lat1,lon1,&
                   &  T_grid, Atm_time,atm_step,&
                   & isd, ied, jsd, jed,kms, kme,its, ite, jts,jte,kts,kte,&
                        ids,  ide,  jds,  jde,  kds,  kde, dx)

                     do j_t = jsd,jed
                     do k_t = 1,z_top
                     do i_t = isd,ied

                        if(Atm_var_g%y(i_t,j_t) .le. 80) then
                         temp_t = (ens_5day_to%data(i_t,k_t,j_t)-corr_to_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)
                        else
                         temp_t = (ens_5day_to%data(i_t,k_t,j_t)-corr_to_5day(i_t,k_t,j_t))*cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)*&
                                                &cos(Atm_var_g%y(i_t,j_t)*RAD_TO_DEG)

                        endif
                       corr_t5dayo_r(i_t,k_t,j_t) =  temp_t

                     enddo
                     enddo
                     enddo

         endif ! hea_5day



   endif ! if(   atm_step==1 .or. mod(atm_step,360)==0   ) then




  end  subroutine mda_o



        subroutine mk_ens_hf_ada(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
! one mins intv
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)

         integer, intent(in)  :: it_step,isd,ied, jsd, jed
         integer  :: n, m,taup1, i, i_a, i_b,jj,ii  !lulv        


       if(it_step <= heaada_enssize) then
      
               do jj = jsd,jed
               do ii = isd,ied
                 ps_ens_tau(it_step)%data(ii,jj,1) = atmos_var%ps(ii,jj)
                enddo
                enddo
          else

           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1) 
                enddo
                enddo
          end do

                do jj = jsd,jed
                do ii = isd,ied
                   ps_ens_tau(heaada_enssize)%data(ii,jj,1) = atmos_var%ps(ii,jj)
                enddo
                enddo
       endif
        end subroutine mk_ens_hf_ada


        subroutine mk_ens_hf_ada_min(atmos_var, ps_ens_tau,it_step,isd,ied,jsd, jed)

         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize
        
         time_min_ens=heahfq_minintv*heaada_enssize
         real_heaada_enssize = REAL(heaada_enssize)

         i_a = mod(it_step,heahfq_minintv)
         if(i_a==0) i_a=heahfq_minintv

        if (it_step <= time_min_ens ) then !it_step<=400 mint
               do jj = jsd,jed
               do ii = isd,ied
                ps_ens_tmp(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

         if (i_a ==heahfq_minintv) then
           i_ensb = i_ensb + 1
           do m=1,heahfq_minintv

            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tau(i_ensb)%data(ii,jj,1) = ps_ens_tau(i_ensb)%data(ii,jj,1) &
                                                + ps_ens_tmp(m)%data(ii,jj,1)*1./real_heaada_enssize
            enddo
            enddo
           enddo
         endif !i_a <= heahfq_minintv
         endif 

 
        if (it_step > time_min_ens ) then

               do jj = jsd,jed
               do ii = isd,ied
                ps_ens_tmp(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo
 
        if( i_a == heahfq_minintv ) then      

               do m =1 , heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1) = ps_ens_tau(m+1)%data(ii,jj,1)
               enddo
               enddo
               enddo
  
         ps_ens_tau(heaada_enssize)%data(:,:,1) = 0.0
           do m=1,heahfq_minintv
           do jj = jsd,jed
           do ii = isd,ied
           ps_ens_tau(heaada_enssize)%data(ii,jj,1) = ps_ens_tau(heaada_enssize)%data(ii,jj,1) &
                                        + ps_ens_tmp(m)%data(ii,jj,1)*1./real_heaada_enssize
           enddo
           enddo
           enddo 

         endif !i_a == heahfq_minintv
        endif !it_step >= time_min_ens


        end subroutine mk_ens_hf_ada_min


!      subroutine mk_ens_hf_ada_day(atmos_var, ps_ens_tau,it_step,isd,ied,jsd,jed)
   subroutine mk_ens_hf_ada_day(atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed)

         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_day
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_day = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
        real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)

         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

!    write(*,*) 'ensemble_min_day,time_min_ens',ensemble_min_day,time_min_ens
!    write(*,*) 'it_step,i_a,ps_ens_tmp_min',it_step,i_a,ps_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

       if (it_step <= time_min_ens ) then !1440*20day
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) =ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

!     write(*,*) 'i_ensb,ps_ens_tmp_hour',i_ensb,ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1)

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv


        if (i_day==0) then
          i_ensb_day=i_ensb_day+1 !1~20 day
!          write(*,*) 'it_step,i_a,i_day,i_ensb_day',it_step,i_a,i_day,i_ensb_day
          ps_ens_tau(i_ensb_day)%data(:,:,1) =0.

            do m=1,24
            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(i_ensb_day)%data(ii,jj,1) = ps_ens_tau(i_ensb_day)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(m)%data(ii,jj,1)*1./24.
            enddo
            enddo
            enddo

!      write(*,*) 'it_step,i_ensb_day,ps_ens_tau',it_step,i_ensb_day,ps_ens_tau(i_ensb_day)%data(isd+10,jsd+10,1)

        endif
       endif

      if (it_step > time_min_ens ) then
       if(i_ensb==24) i_ensb=0

        if( i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!     write(*,*) 'i_a,i_day,i_ensb',it_step,i_a,i_day,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) = ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                                + ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo !m

             if (i_ensb ==24) i_ensb=0
         endif ! ia==60


       if( i_day == 0 ) then
!       write(*,*) 'i_a,i_day',it_step,i_a,i_day
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1)
                enddo
                enddo
          end do

          ps_ens_tau(heaada_enssize)%data(:,:,1) =0.

            do m=1,24
            do jj = jsd,jed
            do ii = isd,ied
        ps_ens_tau(heaada_enssize)%data(ii,jj,1) = ps_ens_tau(heaada_enssize)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(m)%data(ii,jj,1)*1./24.
            enddo
            enddo
            enddo

       endif !i_day==0
      endif !it_step >= time_min_ens

!      write(*,*) 'ps_ens_tau20',ps_ens_tau(heaada_enssize)%data(isd+10,jsd+10,1)

      end subroutine mk_ens_hf_ada_day

        subroutine  get_a_to_b(ps_ens_tau,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,&
!                        isc_flth,iec_flth, jsc_flth, jec_flth)
                        isd,ied,jsd, jed)
!        integer, intent(in)  :: isc_flth, iec_flth, jsc_flth, jec_flth

          integer, intent(in)  :: isd,ied,jsd, jed
         type(ada_field_type), intent(in) :: ps_ens_tau(:)
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso01,ens_hfq_pso02
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso03,ens_hfq_pso04
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso05,ens_hfq_pso06
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso07,ens_hfq_pso08
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso09,ens_hfq_pso10
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso11,ens_hfq_pso12
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso13,ens_hfq_pso14
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso15,ens_hfq_pso16
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso17,ens_hfq_pso18
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso19,ens_hfq_pso20
   integer :: i,j
!  write(*,*) 'get_a_to_b1',ps_ens_tau(1)%data(isd+10,jsd+10,1)

   do j= jsd, jed
   do i= isd,ied
        ens_hfq_pso01(i,1,j) = ps_ens_tau(1)%data(i,j,1)
        ens_hfq_pso02(i,1,j) = ps_ens_tau(2)%data(i,j,1)
        ens_hfq_pso03(i,1,j) = ps_ens_tau(3)%data(i,j,1)
        ens_hfq_pso04(i,1,j) = ps_ens_tau(4)%data(i,j,1)
        ens_hfq_pso05(i,1,j) = ps_ens_tau(5)%data(i,j,1)
        ens_hfq_pso06(i,1,j) = ps_ens_tau(6)%data(i,j,1)
        ens_hfq_pso07(i,1,j) = ps_ens_tau(7)%data(i,j,1)
        ens_hfq_pso08(i,1,j) = ps_ens_tau(8)%data(i,j,1)
        ens_hfq_pso09(i,1,j) = ps_ens_tau(9)%data(i,j,1)
        ens_hfq_pso10(i,1,j) = ps_ens_tau(10)%data(i,j,1)
        ens_hfq_pso11(i,1,j) = ps_ens_tau(11)%data(i,j,1)
        ens_hfq_pso12(i,1,j) = ps_ens_tau(12)%data(i,j,1)
        ens_hfq_pso13(i,1,j) = ps_ens_tau(13)%data(i,j,1)
        ens_hfq_pso14(i,1,j) = ps_ens_tau(14)%data(i,j,1)
        ens_hfq_pso15(i,1,j) = ps_ens_tau(15)%data(i,j,1)
        ens_hfq_pso16(i,1,j) = ps_ens_tau(16)%data(i,j,1)
        ens_hfq_pso17(i,1,j) = ps_ens_tau(17)%data(i,j,1)
        ens_hfq_pso18(i,1,j) = ps_ens_tau(18)%data(i,j,1)
        ens_hfq_pso19(i,1,j) = ps_ens_tau(19)%data(i,j,1)
        ens_hfq_pso20(i,1,j) = ps_ens_tau(20)%data(i,j,1)
   enddo
   enddo

!   write(*,*) 'get_a_to_b',ens_hfq_pso20(isd+10,1,jsd+10)
      end subroutine get_a_to_b

   subroutine get_a_to_b_tuv(ps_ens_tau,ens_hfq_pso01,ens_hfq_pso02,ens_hfq_pso03,&
                        ens_hfq_pso04,ens_hfq_pso05,ens_hfq_pso06,ens_hfq_pso07,ens_hfq_pso08,&
                        ens_hfq_pso09,ens_hfq_pso10,ens_hfq_pso11,ens_hfq_pso12,&
                        ens_hfq_pso13,ens_hfq_pso14,ens_hfq_pso15,ens_hfq_pso16,ens_hfq_pso17,&
                        ens_hfq_pso18,ens_hfq_pso19,ens_hfq_pso20,isd,ied,jsd,jed,nlev)

          integer, intent(in)  :: isd,ied,jsd, jed,nlev
         type(ada_field_type), intent(in) :: ps_ens_tau(:)
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso01,ens_hfq_pso02
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso03,ens_hfq_pso04
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso05,ens_hfq_pso06
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso07,ens_hfq_pso08
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso09,ens_hfq_pso10
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso11,ens_hfq_pso12
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso13,ens_hfq_pso14
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso15,ens_hfq_pso16
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso17,ens_hfq_pso18
   real, dimension(:,:,:), allocatable, intent(inout) :: ens_hfq_pso19,ens_hfq_pso20
   integer :: i,j,k

!  write(*,*) 'get_a_to_t',ps_ens_tau(1)%data(isd+10,1,jsd+10)
   do j= jsd, jed
   do k= 1,nlev
   do i= isd,ied
        ens_hfq_pso01(i,k,j) = ps_ens_tau(1)%data(i,k,j)
        ens_hfq_pso02(i,k,j) = ps_ens_tau(2)%data(i,k,j)
        ens_hfq_pso03(i,k,j) = ps_ens_tau(3)%data(i,k,j)
        ens_hfq_pso04(i,k,j) = ps_ens_tau(4)%data(i,k,j)
        ens_hfq_pso05(i,k,j) = ps_ens_tau(5)%data(i,k,j)
        ens_hfq_pso06(i,k,j) = ps_ens_tau(6)%data(i,k,j)
        ens_hfq_pso07(i,k,j) = ps_ens_tau(7)%data(i,k,j)
        ens_hfq_pso08(i,k,j) = ps_ens_tau(8)%data(i,k,j)
        ens_hfq_pso09(i,k,j) = ps_ens_tau(9)%data(i,k,j)
        ens_hfq_pso10(i,k,j) = ps_ens_tau(10)%data(i,k,j)
        ens_hfq_pso11(i,k,j) = ps_ens_tau(11)%data(i,k,j)
        ens_hfq_pso12(i,k,j) = ps_ens_tau(12)%data(i,k,j)
        ens_hfq_pso13(i,k,j) = ps_ens_tau(13)%data(i,k,j)
        ens_hfq_pso14(i,k,j) = ps_ens_tau(14)%data(i,k,j)
        ens_hfq_pso15(i,k,j) = ps_ens_tau(15)%data(i,k,j)
        ens_hfq_pso16(i,k,j) = ps_ens_tau(16)%data(i,k,j)
        ens_hfq_pso17(i,k,j) = ps_ens_tau(17)%data(i,k,j)
        ens_hfq_pso18(i,k,j) = ps_ens_tau(18)%data(i,k,j)
        ens_hfq_pso19(i,k,j) = ps_ens_tau(19)%data(i,k,j)
        ens_hfq_pso20(i,k,j) = ps_ens_tau(20)%data(i,k,j)
   enddo
   enddo
   enddo


   end subroutine get_a_to_b_tuv

      subroutine read_data_1h(grid,ens_1h_ps0,isd,ied, jsd, jed )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_1h_ps0(:)
              integer, intent(in)  :: isd,ied, jsd, jed
              integer :: i,j
   do j=jsd,jed
   do i=isd,ied

       ens_1h_ps0(1)%data(i,j,1)   =  grid%ens1h01(i,j)
       ens_1h_ps0(2)%data(i,j,1)   =  grid%ens1h02(i,j)
       ens_1h_ps0(3)%data(i,j,1)   =  grid%ens1h03(i,j)
       ens_1h_ps0(4)%data(i,j,1)   =  grid%ens1h04(i,j)
       ens_1h_ps0(5)%data(i,j,1)   =  grid%ens1h05(i,j)
       ens_1h_ps0(6)%data(i,j,1)   =  grid%ens1h06(i,j)
       ens_1h_ps0(7)%data(i,j,1)   =  grid%ens1h07(i,j)
       ens_1h_ps0(8)%data(i,j,1)   =  grid%ens1h08(i,j)
       ens_1h_ps0(9)%data(i,j,1)   =  grid%ens1h09(i,j)
       ens_1h_ps0(10)%data(i,j,1)  =  grid%ens1h10(i,j)
       ens_1h_ps0(11)%data(i,j,1)  =  grid%ens1h11(i,j)
       ens_1h_ps0(12)%data(i,j,1)  =  grid%ens1h12(i,j)
       ens_1h_ps0(13)%data(i,j,1)  =  grid%ens1h13(i,j)
       ens_1h_ps0(14)%data(i,j,1)  =  grid%ens1h14(i,j)
       ens_1h_ps0(15)%data(i,j,1)  =  grid%ens1h15(i,j)
       ens_1h_ps0(16)%data(i,j,1)  =  grid%ens1h16(i,j)
       ens_1h_ps0(17)%data(i,j,1)  =  grid%ens1h17(i,j)
       ens_1h_ps0(18)%data(i,j,1)  =  grid%ens1h18(i,j)
       ens_1h_ps0(19)%data(i,j,1)  =  grid%ens1h19(i,j)
       ens_1h_ps0(20)%data(i,j,1)  =  grid%ens1h20(i,j)

   enddo
   enddo

      end subroutine read_data_1h

      subroutine read_data_3h(grid,ens_3h_ps0,isd,ied, jsd, jed )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_3h_ps0(:)
              integer, intent(in)  :: isd,ied, jsd, jed
              integer :: i,j
   do j=jsd,jed
   do i=isd,ied
       ens_3h_ps0(1)%data(i,j,1)   =  grid%ens3h01(i,j)
       ens_3h_ps0(2)%data(i,j,1)   =  grid%ens3h02(i,j)
       ens_3h_ps0(3)%data(i,j,1)   =  grid%ens3h03(i,j)
       ens_3h_ps0(4)%data(i,j,1)   =  grid%ens3h04(i,j)
       ens_3h_ps0(5)%data(i,j,1)   =  grid%ens3h05(i,j)
       ens_3h_ps0(6)%data(i,j,1)   =  grid%ens3h06(i,j)
       ens_3h_ps0(7)%data(i,j,1)   =  grid%ens3h07(i,j)
       ens_3h_ps0(8)%data(i,j,1)   =  grid%ens3h08(i,j)
       ens_3h_ps0(9)%data(i,j,1)   =  grid%ens3h09(i,j)
       ens_3h_ps0(10)%data(i,j,1)  =  grid%ens3h10(i,j)
       ens_3h_ps0(11)%data(i,j,1)  =  grid%ens3h11(i,j)
       ens_3h_ps0(12)%data(i,j,1)  =  grid%ens3h12(i,j)
       ens_3h_ps0(13)%data(i,j,1)  =  grid%ens3h13(i,j)
       ens_3h_ps0(14)%data(i,j,1)  =  grid%ens3h14(i,j)
       ens_3h_ps0(15)%data(i,j,1)  =  grid%ens3h15(i,j)
       ens_3h_ps0(16)%data(i,j,1)  =  grid%ens3h16(i,j)
       ens_3h_ps0(17)%data(i,j,1)  =  grid%ens3h17(i,j)
       ens_3h_ps0(18)%data(i,j,1)  =  grid%ens3h18(i,j)
       ens_3h_ps0(19)%data(i,j,1)  =  grid%ens3h19(i,j)
       ens_3h_ps0(20)%data(i,j,1)  =  grid%ens3h20(i,j)

   enddo
   enddo

!  write(*,*) 'get_data_grid1',ens_3h_ps0(1)%data(isd+10,jsd+10,1)
!  write(*,*) 'get_data1_grid20',ens_3h_ps0(20)%data(isd+10,jsd+10,1)
end   subroutine read_data_3h


      subroutine read_data_5day(grid,ens_5day_ps0,isd,ied, jsd, jed )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_5day_ps0(:)
              integer, intent(in)  :: isd,ied, jsd, jed
              integer :: i,j
   do j=jsd,jed
   do i=isd,ied
       ens_5day_ps0(1)%data(i,j,1)   =  grid%ens5day1(i,j)
       ens_5day_ps0(2)%data(i,j,1)   =  grid%ens5day2(i,j)
       ens_5day_ps0(3)%data(i,j,1)   =  grid%ens5day3(i,j)
       ens_5day_ps0(4)%data(i,j,1)   =  grid%ens5day4(i,j)
       ens_5day_ps0(5)%data(i,j,1)   =  grid%ens5day5(i,j)
       ens_5day_ps0(6)%data(i,j,1)   =  grid%ens5day6(i,j)
       ens_5day_ps0(7)%data(i,j,1)   =  grid%ens5day7(i,j)
       ens_5day_ps0(8)%data(i,j,1)   =  grid%ens5day8(i,j)
       ens_5day_ps0(9)%data(i,j,1)   =  grid%ens5day9(i,j)
       ens_5day_ps0(10)%data(i,j,1)  =  grid%ens5day10(i,j)
       ens_5day_ps0(11)%data(i,j,1)  =  grid%ens5day11(i,j)
       ens_5day_ps0(12)%data(i,j,1)  =  grid%ens5day12(i,j)
       ens_5day_ps0(13)%data(i,j,1)  =  grid%ens5day13(i,j)
       ens_5day_ps0(14)%data(i,j,1)  =  grid%ens5day14(i,j)
       ens_5day_ps0(15)%data(i,j,1)  =  grid%ens5day15(i,j)
       ens_5day_ps0(16)%data(i,j,1)  =  grid%ens5day16(i,j)
       ens_5day_ps0(17)%data(i,j,1)  =  grid%ens5day17(i,j)
       ens_5day_ps0(18)%data(i,j,1)  =  grid%ens5day18(i,j)
       ens_5day_ps0(19)%data(i,j,1)  =  grid%ens5day19(i,j)
       ens_5day_ps0(20)%data(i,j,1)  =  grid%ens5day20(i,j)

   enddo
   enddo

!  write(*,*) 'get_data_grid5day1',ens_5day_ps0(1)%data(isd+10,jsd+10,1)
!  write(*,*) 'get_data1_grid5day20',ens_5day_ps0(20)%data(isd+10,jsd+10,1)
end   subroutine read_data_5day

      subroutine read_data_6h(grid,ens_6h_ps0,isd,ied, jsd, jed )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_6h_ps0(:)
              integer, intent(in)  :: isd,ied, jsd, jed
              integer :: i,j
   do j=jsd,jed
   do i=isd,ied
       ens_6h_ps0(1)%data(i,j,1)   =  grid%ens6h01(i,j)
       ens_6h_ps0(2)%data(i,j,1)   =  grid%ens6h02(i,j)
       ens_6h_ps0(3)%data(i,j,1)   =  grid%ens6h03(i,j)
       ens_6h_ps0(4)%data(i,j,1)   =  grid%ens6h04(i,j)
       ens_6h_ps0(5)%data(i,j,1)   =  grid%ens6h05(i,j)
       ens_6h_ps0(6)%data(i,j,1)   =  grid%ens6h06(i,j)
       ens_6h_ps0(7)%data(i,j,1)   =  grid%ens6h07(i,j)
       ens_6h_ps0(8)%data(i,j,1)   =  grid%ens6h08(i,j)
       ens_6h_ps0(9)%data(i,j,1)   =  grid%ens6h09(i,j)
       ens_6h_ps0(10)%data(i,j,1)  =  grid%ens6h10(i,j)
       ens_6h_ps0(11)%data(i,j,1)  =  grid%ens6h11(i,j)
       ens_6h_ps0(12)%data(i,j,1)  =  grid%ens6h12(i,j)
       ens_6h_ps0(13)%data(i,j,1)  =  grid%ens6h13(i,j)
       ens_6h_ps0(14)%data(i,j,1)  =  grid%ens6h14(i,j)
       ens_6h_ps0(15)%data(i,j,1)  =  grid%ens6h15(i,j)
       ens_6h_ps0(16)%data(i,j,1)  =  grid%ens6h16(i,j)
       ens_6h_ps0(17)%data(i,j,1)  =  grid%ens6h17(i,j)
       ens_6h_ps0(18)%data(i,j,1)  =  grid%ens6h18(i,j)
       ens_6h_ps0(19)%data(i,j,1)  =  grid%ens6h19(i,j)
       ens_6h_ps0(20)%data(i,j,1)  =  grid%ens6h20(i,j)

   enddo
   enddo
!  write(*,*) 'get_data_grid1',ens_6h_ps0(1)%data(isd+10,jsd+10,1)
!  write(*,*) 'get_data1_grid20',ens_6h_ps0(20)%data(isd+10,jsd+10,1)

      end subroutine read_data_6h

        subroutine read_data_day(grid,ens_day_ps0,isd,ied, jsd, jed )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_day_ps0(:)
              integer, intent(in)  :: isd,ied, jsd, jed
              integer :: i,j
!  write(*,*) 'read_data1_grid1',isc_flth+5,jsc_flth+5,grid%ENSDAY1(isd+10,jsd+10)
!  write(*,*) 'read_data1_grid20',isc_flth+5,jsc_flth+5,grid%ENSDAY20(isd+10,jsd+10)
   do j=jsd,jed
   do i=isd,ied
       ens_day_ps0(1)%data(i,j,1)   =  grid%ENSDAY1(i,j)
       ens_day_ps0(2)%data(i,j,1)   =  grid%ENSDAY2(i,j)
       ens_day_ps0(3)%data(i,j,1)   =  grid%ENSDAY3(i,j)
       ens_day_ps0(4)%data(i,j,1)   =  grid%ENSDAY4(i,j)
       ens_day_ps0(5)%data(i,j,1)   =  grid%ENSDAY5(i,j)
       ens_day_ps0(6)%data(i,j,1)   =  grid%ENSDAY6(i,j)
       ens_day_ps0(7)%data(i,j,1)   =  grid%ENSDAY7(i,j)
       ens_day_ps0(8)%data(i,j,1)   =  grid%ENSDAY8(i,j)
       ens_day_ps0(9)%data(i,j,1)   =  grid%ENSDAY9(i,j)
       ens_day_ps0(10)%data(i,j,1)  =  grid%ENSDAY10(i,j)
       ens_day_ps0(11)%data(i,j,1)  =  grid%ENSDAY11(i,j)
       ens_day_ps0(12)%data(i,j,1)  =  grid%ENSDAY12(i,j)
       ens_day_ps0(13)%data(i,j,1)  =  grid%ENSDAY13(i,j)
       ens_day_ps0(14)%data(i,j,1)  =  grid%ENSDAY14(i,j)
       ens_day_ps0(15)%data(i,j,1)  =  grid%ENSDAY15(i,j)
       ens_day_ps0(16)%data(i,j,1)  =  grid%ENSDAY16(i,j)
       ens_day_ps0(17)%data(i,j,1)  =  grid%ENSDAY17(i,j)
       ens_day_ps0(18)%data(i,j,1)  =  grid%ENSDAY18(i,j)
       ens_day_ps0(19)%data(i,j,1)  =  grid%ENSDAY19(i,j)
       ens_day_ps0(20)%data(i,j,1)  =  grid%ENSDAY20(i,j)

   enddo
   enddo
!  write(*,*) 'get_data_grid1',ens_day_ps0(1)%data(isd+10,jsd+10,1)
!  write(*,*) 'get_data1_grid20',ens_day_ps0(20)%data(isd+10,jsd+10,1)
        end subroutine read_data_day

  subroutine read_data_t3h_tuv(grid,ens_3h_t0,isd,ied, jsd, jed,nlev,vari )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_3h_t0(:)
              integer, intent(in)  :: isd,ied, jsd, jed,nlev
              integer :: i,j,k
              character(len=1) :: vari
        write(*,*) 'read_data_t3h_tuv'

  select case(vari)
     case('t')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied
        ens_3h_t0(1)%data(i,k,j)    =  grid%ENS3hT1(i,k,j)
        ens_3h_t0(2)%data(i,k,j)    =  grid%ENS3hT2(i,k,j)
        ens_3h_t0(3)%data(i,k,j)    =  grid%ENS3hT3(i,k,j)
        ens_3h_t0(4)%data(i,k,j)    =  grid%ENS3hT4(i,k,j)
        ens_3h_t0(5)%data(i,k,j)    =  grid%ENS3hT5(i,k,j)
        ens_3h_t0(6)%data(i,k,j)    =  grid%ENS3hT6(i,k,j)
        ens_3h_t0(7)%data(i,k,j)    =  grid%ENS3hT7(i,k,j)
        ens_3h_t0(8)%data(i,k,j)    =  grid%ENS3hT8(i,k,j)
        ens_3h_t0(9)%data(i,k,j)    =  grid%ENS3hT9(i,k,j)
        ens_3h_t0(10)%data(i,k,j)   =  grid%ENS3hT10(i,k,j)
        ens_3h_t0(11)%data(i,k,j)   =  grid%ENS3hT11(i,k,j)
        ens_3h_t0(12)%data(i,k,j)   =  grid%ENS3hT12(i,k,j)
        ens_3h_t0(13)%data(i,k,j)   =  grid%ENS3hT13(i,k,j)
        ens_3h_t0(14)%data(i,k,j)   =  grid%ENS3hT14(i,k,j)
        ens_3h_t0(15)%data(i,k,j)   =  grid%ENS3hT15(i,k,j)
        ens_3h_t0(16)%data(i,k,j)   =  grid%ENS3hT16(i,k,j)
        ens_3h_t0(17)%data(i,k,j)   =  grid%ENS3hT17(i,k,j)
        ens_3h_t0(18)%data(i,k,j)   =  grid%ENS3hT18(i,k,j)
        ens_3h_t0(19)%data(i,k,j)   =  grid%ENS3hT19(i,k,j)
        ens_3h_t0(20)%data(i,k,j)   =  grid%ENS3hT20(i,k,j)
       enddo
       enddo
       enddo

     case('u')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

        ens_3h_t0(1)%data(i,k,j)    =  grid%ENS3hU1(i,k,j)
        ens_3h_t0(2)%data(i,k,j)    =  grid%ENS3hU2(i,k,j)
        ens_3h_t0(3)%data(i,k,j)    =  grid%ENS3hU3(i,k,j)
        ens_3h_t0(4)%data(i,k,j)    =  grid%ENS3hU4(i,k,j)
        ens_3h_t0(5)%data(i,k,j)    =  grid%ENS3hU5(i,k,j)
        ens_3h_t0(6)%data(i,k,j)    =  grid%ENS3hU6(i,k,j)
        ens_3h_t0(7)%data(i,k,j)    =  grid%ENS3hU7(i,k,j)
        ens_3h_t0(8)%data(i,k,j)    =  grid%ENS3hU8(i,k,j)
        ens_3h_t0(9)%data(i,k,j)    =  grid%ENS3hU9(i,k,j)
        ens_3h_t0(10)%data(i,k,j)   =  grid%ENS3hU10(i,k,j)
        ens_3h_t0(11)%data(i,k,j)   =  grid%ENS3hU11(i,k,j)
        ens_3h_t0(12)%data(i,k,j)   =  grid%ENS3hU12(i,k,j)
        ens_3h_t0(13)%data(i,k,j)   =  grid%ENS3hU13(i,k,j)
        ens_3h_t0(14)%data(i,k,j)   =  grid%ENS3hU14(i,k,j)
        ens_3h_t0(15)%data(i,k,j)   =  grid%ENS3hU15(i,k,j)
        ens_3h_t0(16)%data(i,k,j)   =  grid%ENS3hU16(i,k,j)
        ens_3h_t0(17)%data(i,k,j)   =  grid%ENS3hU17(i,k,j)
        ens_3h_t0(18)%data(i,k,j)   =  grid%ENS3hU18(i,k,j)
        ens_3h_t0(19)%data(i,k,j)   =  grid%ENS3hU19(i,k,j)
        ens_3h_t0(20)%data(i,k,j)   =  grid%ENS3hU20(i,k,j)

           enddo
           enddo
           enddo


     case('v')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

        ens_3h_t0(1)%data(i,k,j)    =  grid%ENS3hV1(i,k,j)
        ens_3h_t0(2)%data(i,k,j)    =  grid%ENS3hV2(i,k,j)
        ens_3h_t0(3)%data(i,k,j)    =  grid%ENS3hV3(i,k,j)
        ens_3h_t0(4)%data(i,k,j)    =  grid%ENS3hV4(i,k,j)
        ens_3h_t0(5)%data(i,k,j)    =  grid%ENS3hV5(i,k,j)
        ens_3h_t0(6)%data(i,k,j)    =  grid%ENS3hV6(i,k,j)
        ens_3h_t0(7)%data(i,k,j)    =  grid%ENS3hV7(i,k,j)
        ens_3h_t0(8)%data(i,k,j)    =  grid%ENS3hV8(i,k,j)
        ens_3h_t0(9)%data(i,k,j)    =  grid%ENS3hV9(i,k,j)
        ens_3h_t0(10)%data(i,k,j)   =  grid%ENS3hV10(i,k,j)
        ens_3h_t0(11)%data(i,k,j)   =  grid%ENS3hV11(i,k,j)
        ens_3h_t0(12)%data(i,k,j)   =  grid%ENS3hV12(i,k,j)
        ens_3h_t0(13)%data(i,k,j)   =  grid%ENS3hV13(i,k,j)
        ens_3h_t0(14)%data(i,k,j)   =  grid%ENS3hV14(i,k,j)
        ens_3h_t0(15)%data(i,k,j)   =  grid%ENS3hV15(i,k,j)
        ens_3h_t0(16)%data(i,k,j)   =  grid%ENS3hV16(i,k,j)
        ens_3h_t0(17)%data(i,k,j)   =  grid%ENS3hV17(i,k,j)
        ens_3h_t0(18)%data(i,k,j)   =  grid%ENS3hV18(i,k,j)
        ens_3h_t0(19)%data(i,k,j)   =  grid%ENS3hV19(i,k,j)
        ens_3h_t0(20)%data(i,k,j)   =  grid%ENS3hV20(i,k,j)

           enddo
           enddo
           enddo


     case('q')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied
        ens_3h_t0(1)%data(i,k,j)    =  grid%ENS3HQ1(i,k,j)
        ens_3h_t0(2)%data(i,k,j)    =  grid%ENS3HQ2(i,k,j)
        ens_3h_t0(3)%data(i,k,j)    =  grid%ENS3HQ3(i,k,j)
        ens_3h_t0(4)%data(i,k,j)    =  grid%ENS3HQ4(i,k,j)
        ens_3h_t0(5)%data(i,k,j)    =  grid%ENS3HQ5(i,k,j)
        ens_3h_t0(6)%data(i,k,j)    =  grid%ENS3HQ6(i,k,j)
        ens_3h_t0(7)%data(i,k,j)    =  grid%ENS3HQ7(i,k,j)
        ens_3h_t0(8)%data(i,k,j)    =  grid%ENS3HQ8(i,k,j)
        ens_3h_t0(9)%data(i,k,j)    =  grid%ENS3HQ9(i,k,j)
        ens_3h_t0(10)%data(i,k,j)   =  grid%ENS3HQ10(i,k,j)
        ens_3h_t0(11)%data(i,k,j)   =  grid%ENS3HQ11(i,k,j)
        ens_3h_t0(12)%data(i,k,j)   =  grid%ENS3HQ12(i,k,j)
        ens_3h_t0(13)%data(i,k,j)   =  grid%ENS3HQ13(i,k,j)
        ens_3h_t0(14)%data(i,k,j)   =  grid%ENS3HQ14(i,k,j)
        ens_3h_t0(15)%data(i,k,j)   =  grid%ENS3HQ15(i,k,j)
        ens_3h_t0(16)%data(i,k,j)   =  grid%ENS3HQ16(i,k,j)
        ens_3h_t0(17)%data(i,k,j)   =  grid%ENS3HQ17(i,k,j)
        ens_3h_t0(18)%data(i,k,j)   =  grid%ENS3HQ18(i,k,j)
        ens_3h_t0(19)%data(i,k,j)   =  grid%ENS3HQ19(i,k,j)
        ens_3h_t0(20)%data(i,k,j)   =  grid%ENS3HQ20(i,k,j)
           enddo
           enddo
           enddo

!          write(*,*) 'get_t3h_grid1',ens_3h_t0(1)%data(isd+10,1,jsd+10),grid%ENS3HQ1(isd+10,1,jsd+10)
!          write(*,*) 'get_t3h_grid20',ens_3h_t0(20)%data(isd+10,1,jsd+10),grid%ENS3HQ20(isd+10,1,jsd+10)

   end select

  end subroutine read_data_t3h_tuv


  subroutine read_data_t6h_tuv(grid,ens_6h_t0,isd,ied, jsd, jed,nlev,vari )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_6h_t0(:)
              integer, intent(in)  :: isd,ied, jsd, jed,nlev
              integer :: i,j,k
              character(len=1) :: vari

write(*,*) 'read_6h'
  select case(vari)
     case('t')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

        ens_6h_t0(1)%data(i,k,j)    =  grid%ENS6hT1(i,k,j)
        ens_6h_t0(2)%data(i,k,j)    =  grid%ENS6hT2(i,k,j)
        ens_6h_t0(3)%data(i,k,j)    =  grid%ENS6hT3(i,k,j)
        ens_6h_t0(4)%data(i,k,j)    =  grid%ENS6hT4(i,k,j)
        ens_6h_t0(5)%data(i,k,j)    =  grid%ENS6hT5(i,k,j)
        ens_6h_t0(6)%data(i,k,j)    =  grid%ENS6hT6(i,k,j)
        ens_6h_t0(7)%data(i,k,j)    =  grid%ENS6hT7(i,k,j)
        ens_6h_t0(8)%data(i,k,j)    =  grid%ENS6hT8(i,k,j)
        ens_6h_t0(9)%data(i,k,j)    =  grid%ENS6hT9(i,k,j)
        ens_6h_t0(10)%data(i,k,j)   =  grid%ENS6hT10(i,k,j)
        ens_6h_t0(11)%data(i,k,j)   =  grid%ENS6hT11(i,k,j)
        ens_6h_t0(12)%data(i,k,j)   =  grid%ENS6hT12(i,k,j)
        ens_6h_t0(13)%data(i,k,j)   =  grid%ENS6hT13(i,k,j)
        ens_6h_t0(14)%data(i,k,j)   =  grid%ENS6hT14(i,k,j)
        ens_6h_t0(15)%data(i,k,j)   =  grid%ENS6hT15(i,k,j)
        ens_6h_t0(16)%data(i,k,j)   =  grid%ENS6hT16(i,k,j)
        ens_6h_t0(17)%data(i,k,j)   =  grid%ENS6hT17(i,k,j)
        ens_6h_t0(18)%data(i,k,j)   =  grid%ENS6hT18(i,k,j)
        ens_6h_t0(19)%data(i,k,j)   =  grid%ENS6hT19(i,k,j)
        ens_6h_t0(20)%data(i,k,j)   =  grid%ENS6hT20(i,k,j)

       enddo
       enddo
       enddo

!          write(*,*) 'get_t6h_grid1',ens_6h_t0(1)%data(isd+10,1,jsd+10)
!          write(*,*) 'get_t6h_grid20',ens_6h_t0(20)%data(isd+10,1,jsd+10)

     case('u')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

        ens_6h_t0(1)%data(i,k,j)    =  grid%ENS6hU1(i,k,j)
        ens_6h_t0(2)%data(i,k,j)    =  grid%ENS6hU2(i,k,j)
        ens_6h_t0(3)%data(i,k,j)    =  grid%ENS6hU3(i,k,j)
        ens_6h_t0(4)%data(i,k,j)    =  grid%ENS6hU4(i,k,j)
        ens_6h_t0(5)%data(i,k,j)    =  grid%ENS6hU5(i,k,j)
        ens_6h_t0(6)%data(i,k,j)    =  grid%ENS6hU6(i,k,j)
        ens_6h_t0(7)%data(i,k,j)    =  grid%ENS6hU7(i,k,j)
        ens_6h_t0(8)%data(i,k,j)    =  grid%ENS6hU8(i,k,j)
        ens_6h_t0(9)%data(i,k,j)    =  grid%ENS6hU9(i,k,j)
        ens_6h_t0(10)%data(i,k,j)   =  grid%ENS6hU10(i,k,j)
        ens_6h_t0(11)%data(i,k,j)   =  grid%ENS6hU11(i,k,j)
        ens_6h_t0(12)%data(i,k,j)   =  grid%ENS6hU12(i,k,j)
        ens_6h_t0(13)%data(i,k,j)   =  grid%ENS6hU13(i,k,j)
        ens_6h_t0(14)%data(i,k,j)   =  grid%ENS6hU14(i,k,j)
        ens_6h_t0(15)%data(i,k,j)   =  grid%ENS6hU15(i,k,j)
        ens_6h_t0(16)%data(i,k,j)   =  grid%ENS6hU16(i,k,j)
        ens_6h_t0(17)%data(i,k,j)   =  grid%ENS6hU17(i,k,j)
        ens_6h_t0(18)%data(i,k,j)   =  grid%ENS6hU18(i,k,j)
        ens_6h_t0(19)%data(i,k,j)   =  grid%ENS6hU19(i,k,j)
        ens_6h_t0(20)%data(i,k,j)   =  grid%ENS6hU20(i,k,j)

           enddo
           enddo
           enddo

!          write(*,*) 'get_u6h_grid1',ens_6h_t0(1)%data(isd+10,1,jsd+10)
!          write(*,*) 'get_u6h_grid20',ens_6h_t0(20)%data(isd+10,1,jsd+10)

     case('v')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

        ens_6h_t0(1)%data(i,k,j)    =  grid%ENS6hV1(i,k,j)
        ens_6h_t0(2)%data(i,k,j)    =  grid%ENS6hV2(i,k,j)
        ens_6h_t0(3)%data(i,k,j)    =  grid%ENS6hV3(i,k,j)
        ens_6h_t0(4)%data(i,k,j)    =  grid%ENS6hV4(i,k,j)
        ens_6h_t0(5)%data(i,k,j)    =  grid%ENS6hV5(i,k,j)
        ens_6h_t0(6)%data(i,k,j)    =  grid%ENS6hV6(i,k,j)
        ens_6h_t0(7)%data(i,k,j)    =  grid%ENS6hV7(i,k,j)
        ens_6h_t0(8)%data(i,k,j)    =  grid%ENS6hV8(i,k,j)
        ens_6h_t0(9)%data(i,k,j)    =  grid%ENS6hV9(i,k,j)
        ens_6h_t0(10)%data(i,k,j)   =  grid%ENS6hV10(i,k,j)
        ens_6h_t0(11)%data(i,k,j)   =  grid%ENS6hV11(i,k,j)
        ens_6h_t0(12)%data(i,k,j)   =  grid%ENS6hV12(i,k,j)
        ens_6h_t0(13)%data(i,k,j)   =  grid%ENS6hV13(i,k,j)
        ens_6h_t0(14)%data(i,k,j)   =  grid%ENS6hV14(i,k,j)
        ens_6h_t0(15)%data(i,k,j)   =  grid%ENS6hV15(i,k,j)
        ens_6h_t0(16)%data(i,k,j)   =  grid%ENS6hV16(i,k,j)
        ens_6h_t0(17)%data(i,k,j)   =  grid%ENS6hV17(i,k,j)
        ens_6h_t0(18)%data(i,k,j)   =  grid%ENS6hV18(i,k,j)
        ens_6h_t0(19)%data(i,k,j)   =  grid%ENS6hV19(i,k,j)
        ens_6h_t0(20)%data(i,k,j)   =  grid%ENS6hV20(i,k,j)

           enddo
           enddo
           enddo

!          write(*,*) 'get_v6h_grid1',ens_6h_t0(1)%data(isd+10,1,jsd+10)
!          write(*,*) 'get_v6h_grid20',ens_6h_t0(20)%data(isd+10,1,jsd+10)

     case('q')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

        ens_6h_t0(1)%data(i,k,j)    =  grid%ENS6HQ1(i,k,j)
        ens_6h_t0(2)%data(i,k,j)    =  grid%ENS6HQ2(i,k,j)
        ens_6h_t0(3)%data(i,k,j)    =  grid%ENS6HQ3(i,k,j)
        ens_6h_t0(4)%data(i,k,j)    =  grid%ENS6HQ4(i,k,j)
        ens_6h_t0(5)%data(i,k,j)    =  grid%ENS6HQ5(i,k,j)
        ens_6h_t0(6)%data(i,k,j)    =  grid%ENS6HQ6(i,k,j)
        ens_6h_t0(7)%data(i,k,j)    =  grid%ENS6HQ7(i,k,j)
        ens_6h_t0(8)%data(i,k,j)    =  grid%ENS6HQ8(i,k,j)
        ens_6h_t0(9)%data(i,k,j)    =  grid%ENS6HQ9(i,k,j)
        ens_6h_t0(10)%data(i,k,j)   =  grid%ENS6HQ10(i,k,j)
        ens_6h_t0(11)%data(i,k,j)   =  grid%ENS6HQ11(i,k,j)
        ens_6h_t0(12)%data(i,k,j)   =  grid%ENS6HQ12(i,k,j)
        ens_6h_t0(13)%data(i,k,j)   =  grid%ENS6HQ13(i,k,j)
        ens_6h_t0(14)%data(i,k,j)   =  grid%ENS6HQ14(i,k,j)
        ens_6h_t0(15)%data(i,k,j)   =  grid%ENS6HQ15(i,k,j)
        ens_6h_t0(16)%data(i,k,j)   =  grid%ENS6HQ16(i,k,j)
        ens_6h_t0(17)%data(i,k,j)   =  grid%ENS6HQ17(i,k,j)
        ens_6h_t0(18)%data(i,k,j)   =  grid%ENS6HQ18(i,k,j)
        ens_6h_t0(19)%data(i,k,j)   =  grid%ENS6HQ19(i,k,j)
        ens_6h_t0(20)%data(i,k,j)   =  grid%ENS6HQ20(i,k,j)

           enddo
           enddo
           enddo
!          write(*,*) 'get_t6h_grid1',ens_6h_t0(1)%data(isd+10,1,jsd+10),grid%ENS6HQ1(isd+10,1,jsd+10)
!          write(*,*) 'get_t6h_grid20',ens_6h_t0(20)%data(isd+10,1,jsd+10),grid%ENS6HQ20(isd+10,1,jsd+10)

   end select

  end subroutine read_data_t6h_tuv 

  subroutine read_data_day_tuv(grid,ens_day_t0,isd,ied, jsd, jed,nlev,vari )
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_day_t0(:)
              integer, intent(in)  :: isd,ied, jsd, jed,nlev
              integer :: i,j,k
              character(len=1) :: vari
  select case(vari)
     case('t')
!          write(*,*) 'read_tdata5_grid1',grid%ENS1DAYT1(isd+10,1,jsd+10),grid%ENS1DAYT20(isd+10,28,jsd+10)
!          write(*,*) 'read_tdata5_grid20',grid%ENS1DAYT20(isd+10,1,jsd+10),grid%ENS1DAYT20(isd+10,27,jsd+10)
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied
    
       ens_day_t0(1)%data(i,k,j)   =  grid%ENS1DAYT1(i,k,j)
       ens_day_t0(2)%data(i,k,j)   =  grid%ENS1DAYT2(i,k,j)
       ens_day_t0(3)%data(i,k,j)   =  grid%ENS1DAYT3(i,k,j)
       ens_day_t0(4)%data(i,k,j)   =  grid%ENS1DAYT4(i,k,j)
       ens_day_t0(5)%data(i,k,j)   =  grid%ENS1DAYT5(i,k,j)
       ens_day_t0(6)%data(i,k,j)   =  grid%ENS1DAYT6(i,k,j)
       ens_day_t0(7)%data(i,k,j)   =  grid%ENS1DAYT7(i,k,j)
       ens_day_t0(8)%data(i,k,j)   =  grid%ENS1DAYT8(i,k,j)
       ens_day_t0(9)%data(i,k,j)   =  grid%ENS1DAYT9(i,k,j)
       ens_day_t0(10)%data(i,k,j)  =  grid%ENS1DAYT10(i,k,j)
       ens_day_t0(11)%data(i,k,j)  =  grid%ENS1DAYT11(i,k,j)
       ens_day_t0(12)%data(i,k,j)  =  grid%ENS1DAYT12(i,k,j)
       ens_day_t0(13)%data(i,k,j)  =  grid%ENS1DAYT13(i,k,j)
       ens_day_t0(14)%data(i,k,j)  =  grid%ENS1DAYT14(i,k,j)
       ens_day_t0(15)%data(i,k,j)  =  grid%ENS1DAYT15(i,k,j)
       ens_day_t0(16)%data(i,k,j)  =  grid%ENS1DAYT16(i,k,j)
       ens_day_t0(17)%data(i,k,j)  =  grid%ENS1DAYT17(i,k,j)
       ens_day_t0(18)%data(i,k,j)  =  grid%ENS1DAYT18(i,k,j)
       ens_day_t0(19)%data(i,k,j)  =  grid%ENS1DAYT19(i,k,j)
       ens_day_t0(20)%data(i,k,j)  =  grid%ENS1DAYT20(i,k,j)
           enddo
           enddo
           enddo
!          write(*,*) 'get_tdata_grid1',ens_day_t0(1)%data(isd+10,1,jsd+10)
!          write(*,*) 'get_tdata1_grid20',ens_day_t0(20)%data(isd+10,1,jsd+10)

     case('u')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_day_t0(1)%data(i,k,j)   =  grid%ENS1DAYU1(i,k,j)
       ens_day_t0(2)%data(i,k,j)   =  grid%ENS1DAYU2(i,k,j)
       ens_day_t0(3)%data(i,k,j)   =  grid%ENS1DAYU3(i,k,j)
       ens_day_t0(4)%data(i,k,j)   =  grid%ENS1DAYU4(i,k,j)
       ens_day_t0(5)%data(i,k,j)   =  grid%ENS1DAYU5(i,k,j)
       ens_day_t0(6)%data(i,k,j)   =  grid%ENS1DAYU6(i,k,j)
       ens_day_t0(7)%data(i,k,j)   =  grid%ENS1DAYU7(i,k,j)
       ens_day_t0(8)%data(i,k,j)   =  grid%ENS1DAYU8(i,k,j)
       ens_day_t0(9)%data(i,k,j)   =  grid%ENS1DAYU9(i,k,j)
       ens_day_t0(10)%data(i,k,j)  =  grid%ENS1DAYU10(i,k,j)
       ens_day_t0(11)%data(i,k,j)  =  grid%ENS1DAYU11(i,k,j)
       ens_day_t0(12)%data(i,k,j)  =  grid%ENS1DAYU12(i,k,j)
       ens_day_t0(13)%data(i,k,j)  =  grid%ENS1DAYU13(i,k,j)
       ens_day_t0(14)%data(i,k,j)  =  grid%ENS1DAYU14(i,k,j)
       ens_day_t0(15)%data(i,k,j)  =  grid%ENS1DAYU15(i,k,j)
       ens_day_t0(16)%data(i,k,j)  =  grid%ENS1DAYU16(i,k,j)
       ens_day_t0(17)%data(i,k,j)  =  grid%ENS1DAYU17(i,k,j)
       ens_day_t0(18)%data(i,k,j)  =  grid%ENS1DAYU18(i,k,j)
       ens_day_t0(19)%data(i,k,j)  =  grid%ENS1DAYU19(i,k,j)
       ens_day_t0(20)%data(i,k,j)  =  grid%ENS1DAYU20(i,k,j)
           enddo
           enddo
           enddo
!          write(*,*) 'get_udata_grid1',ens_day_t0(1)%data(isd+10,1,jsd+10)
!          write(*,*) 'get_udata1_grid20',ens_day_t0(20)%data(isd+10,1,jsd+10)


     case('v')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_day_t0(1)%data(i,k,j)   =  grid%ENS1DAYV1(i,k,j)
       ens_day_t0(2)%data(i,k,j)   =  grid%ENS1DAYV2(i,k,j)
       ens_day_t0(3)%data(i,k,j)   =  grid%ENS1DAYV3(i,k,j)
       ens_day_t0(4)%data(i,k,j)   =  grid%ENS1DAYV4(i,k,j)
       ens_day_t0(5)%data(i,k,j)   =  grid%ENS1DAYV5(i,k,j)
       ens_day_t0(6)%data(i,k,j)   =  grid%ENS1DAYV6(i,k,j)
       ens_day_t0(7)%data(i,k,j)   =  grid%ENS1DAYV7(i,k,j)
       ens_day_t0(8)%data(i,k,j)   =  grid%ENS1DAYV8(i,k,j)
       ens_day_t0(9)%data(i,k,j)   =  grid%ENS1DAYV9(i,k,j)
       ens_day_t0(10)%data(i,k,j)  =  grid%ENS1DAYV10(i,k,j)
       ens_day_t0(11)%data(i,k,j)  =  grid%ENS1DAYV11(i,k,j)
       ens_day_t0(12)%data(i,k,j)  =  grid%ENS1DAYV12(i,k,j)
       ens_day_t0(13)%data(i,k,j)  =  grid%ENS1DAYV13(i,k,j)
       ens_day_t0(14)%data(i,k,j)  =  grid%ENS1DAYV14(i,k,j)
       ens_day_t0(15)%data(i,k,j)  =  grid%ENS1DAYV15(i,k,j)
       ens_day_t0(16)%data(i,k,j)  =  grid%ENS1DAYV16(i,k,j)
       ens_day_t0(17)%data(i,k,j)  =  grid%ENS1DAYV17(i,k,j)
       ens_day_t0(18)%data(i,k,j)  =  grid%ENS1DAYV18(i,k,j)
       ens_day_t0(19)%data(i,k,j)  =  grid%ENS1DAYV19(i,k,j)
       ens_day_t0(20)%data(i,k,j)  =  grid%ENS1DAYV20(i,k,j)
           enddo
           enddo
           enddo
     case('q')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_day_t0(1)%data(i,k,j)   =  grid%ENS1DAYQ1(i,k,j)
       ens_day_t0(2)%data(i,k,j)   =  grid%ENS1DAYQ2(i,k,j)
       ens_day_t0(3)%data(i,k,j)   =  grid%ENS1DAYQ3(i,k,j)
       ens_day_t0(4)%data(i,k,j)   =  grid%ENS1DAYQ4(i,k,j)
       ens_day_t0(5)%data(i,k,j)   =  grid%ENS1DAYQ5(i,k,j)
       ens_day_t0(6)%data(i,k,j)   =  grid%ENS1DAYQ6(i,k,j)
       ens_day_t0(7)%data(i,k,j)   =  grid%ENS1DAYQ7(i,k,j)
       ens_day_t0(8)%data(i,k,j)   =  grid%ENS1DAYQ8(i,k,j)
       ens_day_t0(9)%data(i,k,j)   =  grid%ENS1DAYQ9(i,k,j)
       ens_day_t0(10)%data(i,k,j)  =  grid%ENS1DAYQ10(i,k,j)
       ens_day_t0(11)%data(i,k,j)  =  grid%ENS1DAYQ11(i,k,j)
       ens_day_t0(12)%data(i,k,j)  =  grid%ENS1DAYQ12(i,k,j)
       ens_day_t0(13)%data(i,k,j)  =  grid%ENS1DAYQ13(i,k,j)
       ens_day_t0(14)%data(i,k,j)  =  grid%ENS1DAYQ14(i,k,j)
       ens_day_t0(15)%data(i,k,j)  =  grid%ENS1DAYQ15(i,k,j)
       ens_day_t0(16)%data(i,k,j)  =  grid%ENS1DAYQ16(i,k,j)
       ens_day_t0(17)%data(i,k,j)  =  grid%ENS1DAYQ17(i,k,j)
       ens_day_t0(18)%data(i,k,j)  =  grid%ENS1DAYQ18(i,k,j)
       ens_day_t0(19)%data(i,k,j)  =  grid%ENS1DAYQ19(i,k,j)
       ens_day_t0(20)%data(i,k,j)  =  grid%ENS1DAYQ20(i,k,j)
           enddo
           enddo
           enddo


   end select
        end subroutine read_data_day_tuv

        subroutine read_data_5day_tuv(grid,ens_nday_t0,isd,ied, jsd, jed,nlev,vari)
              TYPE ( domain ), INTENT(INOUT) :: grid
              type(ada_field_type), intent(inout) :: ens_nday_t0(:)
              integer, intent(in)  :: isd,ied, jsd, jed,nlev
              integer :: i,j,k
              character(len=1) :: vari
  select case(vari)
     case('t')

!          write(*,*) 'get_5tdata_grid1',grid%ENS5DAYT1(isd+10,1,jsd+10)
!          write(*,*) 'get_5tdata1_grid20',grid%ENS5DAYT20(isd+10,20,jsd+10)
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_nday_t0(1)%data(i,k,j)   =  grid%ENS5DAYT1(i,k,j)
       ens_nday_t0(2)%data(i,k,j)   =  grid%ENS5DAYT2(i,k,j)
       ens_nday_t0(3)%data(i,k,j)   =  grid%ENS5DAYT3(i,k,j)
       ens_nday_t0(4)%data(i,k,j)   =  grid%ENS5DAYT4(i,k,j)
       ens_nday_t0(5)%data(i,k,j)   =  grid%ENS5DAYT5(i,k,j)
       ens_nday_t0(6)%data(i,k,j)   =  grid%ENS5DAYT6(i,k,j)
       ens_nday_t0(7)%data(i,k,j)   =  grid%ENS5DAYT7(i,k,j)
       ens_nday_t0(8)%data(i,k,j)   =  grid%ENS5DAYT8(i,k,j)
       ens_nday_t0(9)%data(i,k,j)   =  grid%ENS5DAYT9(i,k,j)
       ens_nday_t0(10)%data(i,k,j)  =  grid%ENS5DAYT10(i,k,j)
       ens_nday_t0(11)%data(i,k,j)  =  grid%ENS5DAYT11(i,k,j)
       ens_nday_t0(12)%data(i,k,j)  =  grid%ENS5DAYT12(i,k,j)
       ens_nday_t0(13)%data(i,k,j)  =  grid%ENS5DAYT13(i,k,j)
       ens_nday_t0(14)%data(i,k,j)  =  grid%ENS5DAYT14(i,k,j)
       ens_nday_t0(15)%data(i,k,j)  =  grid%ENS5DAYT15(i,k,j)
       ens_nday_t0(16)%data(i,k,j)  =  grid%ENS5DAYT16(i,k,j)
       ens_nday_t0(17)%data(i,k,j)  =  grid%ENS5DAYT17(i,k,j)
       ens_nday_t0(18)%data(i,k,j)  =  grid%ENS5DAYT18(i,k,j)
       ens_nday_t0(19)%data(i,k,j)  =  grid%ENS5DAYT19(i,k,j)
       ens_nday_t0(20)%data(i,k,j)  =  grid%ENS5DAYT20(i,k,j)
           enddo
           enddo
           enddo
!          write(*,*) 'get_tdata_grid1',ens_nday_t0(1)%data(isd+10,1,jsd+10)
!          write(*,*) 'get_tdata1_grid20',ens_nday_t0(20)%data(isd+10,1,jsd+10)

     case('u')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_nday_t0(1)%data(i,k,j)   =  grid%ENS5DAYU1(i,k,j)
       ens_nday_t0(2)%data(i,k,j)   =  grid%ENS5DAYU2(i,k,j)
       ens_nday_t0(3)%data(i,k,j)   =  grid%ENS5DAYU3(i,k,j)
       ens_nday_t0(4)%data(i,k,j)   =  grid%ENS5DAYU4(i,k,j)
       ens_nday_t0(5)%data(i,k,j)   =  grid%ENS5DAYU5(i,k,j)
       ens_nday_t0(6)%data(i,k,j)   =  grid%ENS5DAYU6(i,k,j)
       ens_nday_t0(7)%data(i,k,j)   =  grid%ENS5DAYU7(i,k,j)
       ens_nday_t0(8)%data(i,k,j)   =  grid%ENS5DAYU8(i,k,j)
       ens_nday_t0(9)%data(i,k,j)   =  grid%ENS5DAYU9(i,k,j)
       ens_nday_t0(10)%data(i,k,j)  =  grid%ENS5DAYU10(i,k,j)
       ens_nday_t0(11)%data(i,k,j)  =  grid%ENS5DAYU11(i,k,j)
       ens_nday_t0(12)%data(i,k,j)  =  grid%ENS5DAYU12(i,k,j)
       ens_nday_t0(13)%data(i,k,j)  =  grid%ENS5DAYU13(i,k,j)
       ens_nday_t0(14)%data(i,k,j)  =  grid%ENS5DAYU14(i,k,j)
       ens_nday_t0(15)%data(i,k,j)  =  grid%ENS5DAYU15(i,k,j)
       ens_nday_t0(16)%data(i,k,j)  =  grid%ENS5DAYU16(i,k,j)
       ens_nday_t0(17)%data(i,k,j)  =  grid%ENS5DAYU17(i,k,j)
       ens_nday_t0(18)%data(i,k,j)  =  grid%ENS5DAYU18(i,k,j)
       ens_nday_t0(19)%data(i,k,j)  =  grid%ENS5DAYU19(i,k,j)
       ens_nday_t0(20)%data(i,k,j)  =  grid%ENS5DAYU20(i,k,j)
           enddo
           enddo
           enddo



     case('v')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_nday_t0(1)%data(i,k,j)   =  grid%ENS5DAYV1(i,k,j)
       ens_nday_t0(2)%data(i,k,j)   =  grid%ENS5DAYV2(i,k,j)
       ens_nday_t0(3)%data(i,k,j)   =  grid%ENS5DAYV3(i,k,j)
       ens_nday_t0(4)%data(i,k,j)   =  grid%ENS5DAYV4(i,k,j)
       ens_nday_t0(5)%data(i,k,j)   =  grid%ENS5DAYV5(i,k,j)
       ens_nday_t0(6)%data(i,k,j)   =  grid%ENS5DAYV6(i,k,j)
       ens_nday_t0(7)%data(i,k,j)   =  grid%ENS5DAYV7(i,k,j)
       ens_nday_t0(8)%data(i,k,j)   =  grid%ENS5DAYV8(i,k,j)
       ens_nday_t0(9)%data(i,k,j)   =  grid%ENS5DAYV9(i,k,j)
       ens_nday_t0(10)%data(i,k,j)  =  grid%ENS5DAYV10(i,k,j)
       ens_nday_t0(11)%data(i,k,j)  =  grid%ENS5DAYV11(i,k,j)
       ens_nday_t0(12)%data(i,k,j)  =  grid%ENS5DAYV12(i,k,j)
       ens_nday_t0(13)%data(i,k,j)  =  grid%ENS5DAYV13(i,k,j)
       ens_nday_t0(14)%data(i,k,j)  =  grid%ENS5DAYV14(i,k,j)
       ens_nday_t0(15)%data(i,k,j)  =  grid%ENS5DAYV15(i,k,j)
       ens_nday_t0(16)%data(i,k,j)  =  grid%ENS5DAYV16(i,k,j)
       ens_nday_t0(17)%data(i,k,j)  =  grid%ENS5DAYV17(i,k,j)
       ens_nday_t0(18)%data(i,k,j)  =  grid%ENS5DAYV18(i,k,j)
       ens_nday_t0(19)%data(i,k,j)  =  grid%ENS5DAYV19(i,k,j)
       ens_nday_t0(20)%data(i,k,j)  =  grid%ENS5DAYV20(i,k,j)
           enddo
           enddo
           enddo

    case('q')
       do k=1,nlev
       do j=jsd,jed
       do i=isd,ied

       ens_nday_t0(1)%data(i,k,j)   =  grid%ENS5DAYQ1(i,k,j)
       ens_nday_t0(2)%data(i,k,j)   =  grid%ENS5DAYQ2(i,k,j)
       ens_nday_t0(3)%data(i,k,j)   =  grid%ENS5DAYQ3(i,k,j)
       ens_nday_t0(4)%data(i,k,j)   =  grid%ENS5DAYQ4(i,k,j)
       ens_nday_t0(5)%data(i,k,j)   =  grid%ENS5DAYQ5(i,k,j)
       ens_nday_t0(6)%data(i,k,j)   =  grid%ENS5DAYQ6(i,k,j)
       ens_nday_t0(7)%data(i,k,j)   =  grid%ENS5DAYQ7(i,k,j)
       ens_nday_t0(8)%data(i,k,j)   =  grid%ENS5DAYQ8(i,k,j)
       ens_nday_t0(9)%data(i,k,j)   =  grid%ENS5DAYQ9(i,k,j)
       ens_nday_t0(10)%data(i,k,j)  =  grid%ENS5DAYQ10(i,k,j)
       ens_nday_t0(11)%data(i,k,j)  =  grid%ENS5DAYQ11(i,k,j)
       ens_nday_t0(12)%data(i,k,j)  =  grid%ENS5DAYQ12(i,k,j)
       ens_nday_t0(13)%data(i,k,j)  =  grid%ENS5DAYQ13(i,k,j)
       ens_nday_t0(14)%data(i,k,j)  =  grid%ENS5DAYQ14(i,k,j)
       ens_nday_t0(15)%data(i,k,j)  =  grid%ENS5DAYQ15(i,k,j)
       ens_nday_t0(16)%data(i,k,j)  =  grid%ENS5DAYQ16(i,k,j)
       ens_nday_t0(17)%data(i,k,j)  =  grid%ENS5DAYQ17(i,k,j)
       ens_nday_t0(18)%data(i,k,j)  =  grid%ENS5DAYQ18(i,k,j)
       ens_nday_t0(19)%data(i,k,j)  =  grid%ENS5DAYQ19(i,k,j)
       ens_nday_t0(20)%data(i,k,j)  =  grid%ENS5DAYQ20(i,k,j)
           enddo
           enddo
           enddo


   end select
        end subroutine read_data_5day_tuv


   subroutine mk_ens_hf_ada_day2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_day
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_day = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
        real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60


!     write(*,*) 'ensemble_min_day,time_min_ens',ensemble_min_day,time_min_ens,ps_ens_tmp_min(1)%data(isd+10,jsd+10,1)
!     write(*,*) 'it_step,i_a,ps_ens_tmp_min',it_step,i_a,ps_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

           

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) =ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

!     write(*,*) 'i_ensb,ps_ens_tmp_hour',i_ensb,ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1)

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day==0) then
!          write(*,*) 'it_step,i_a,i_day',it_step,i_a,i_day
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1)
                enddo
                enddo
          end do

          ps_ens_tau(heaada_enssize)%data(:,:,1) =0.

            do m=1,24
            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,jj,1) = ps_ens_tau(heaada_enssize)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(m)%data(ii,jj,1)*1./24.
            enddo
            enddo
            enddo

!      write(*,*) 'it_step,i_ensb_day,ps_ens_tau',it_step,i_ensb_day,ps_ens_tau(heaada_enssize)%data(isd+10,jsd+10,1)

        endif


!      write(*,*) 'ps_ens_tau20',ps_ens_tau(heaada_enssize)%data(isd+10,jsd+10,1)
      end subroutine mk_ens_hf_ada_day2


subroutine mk_ens_ada_6h(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize

        time_min_ens=360  !6h*60min
         i_a = mod(it_step,60)
         i_6h = mod(it_step,360)

         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

!    write(*,*) 'it_step,i_a,ps_ens_6h_min',it_step,i_a,ps_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

       if (it_step <= 7200 ) then !7200=20*360

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.        

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) =ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

!     write(*,*) 'i_ensb,ps_ens_tmp_hour',i_ensb,ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1)

            if (i_ensb ==6) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_6h==0) then
         i_ensb_6h=i_ensb_6h+1 !1~20 day
         ps_ens_tau(i_ensb_6h)%data(:,:,1) =0.
            do m=1,6
            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(i_ensb_6h)%data(ii,jj,1) = ps_ens_tau(i_ensb_6h)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(m)%data(ii,jj,1)*1./6.
            enddo
            enddo
            enddo

!      write(*,*) 'it_step,i_ensb_6h,ps_ens_tau',it_step,i_ensb_6h,ps_ens_tau(i_ensb_6h)%data(isd+10,jsd+10,1)
        endif
       endif

end subroutine mk_ens_ada_6h


subroutine mk_ens_ada_1h(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_1h
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize

         time_min_ens=60  !3h*60min
         i_a = mod(it_step,60)
         i_1h = mod(it_step,60)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

!    write(*,*) 'it_step,i_a,ps_ens_1h_min',it_step,i_a,ps_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

       if (it_step <= 1200 ) then !1200=20*1*60

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tau(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tau(i_ensb)%data(ii,jj,1) = ps_ens_tau(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

! write(*,*) 'it_step,i_ensb_1h,ps_ens_tau',it_step,i_ensb,ps_ens_tau(i_ensb)%data(isd+10,jsd+10,1)
    endif
    endif

end subroutine mk_ens_ada_1h


 subroutine mk_ens_ada_5day2(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_day,i_5day
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_day = 0
         integer,  save ::  i_ensb_5day = 0
         real ::  intv,real_heaada_enssize
         integer, parameter :: ensemble_min_day=1440
 
        time_min_ens=ensemble_min_day*heaada_enssize !1440*20
        real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         i_5day = mod(it_step,7200) !1440*5
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) '5dayit_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) = ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo


            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day == 0 ) then

         i_ensb_day  = i_ensb_day+1
         ps_ens_tmp_day(i_ensb_day)%data(:,:,1) =0.      

           do m=1,24
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_day(i_ensb_day)%data(ii,jj,1) = ps_ens_tmp_day(i_ensb_day)%data(ii,jj,1) +&
                                        ps_ens_tmp_hour(m)%data(ii,jj,1)*1./24.
            enddo
            enddo
           enddo

         if(i_ensb_day == 5) i_ensb_day=0

        endif


      if(i_5day==0) then


           do m = 1, heaada_enssize-1
                write(*,*) 'm',ps_ens_tau(m)%data(isd+10,jsd+10,1)
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1)
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,1) = 0. 

            do m=1,5
            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,jj,1) = ps_ens_tau(heaada_enssize)%data(ii,jj,1) + &
                                              ps_ens_tmp_day(m)%data(ii,jj,1)*1./5.
            enddo
            enddo
            enddo
        endif

 end subroutine  mk_ens_ada_5day2

subroutine mk_ens_ada_1h2(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_3h = 0
         real ::  intv,real_heaada_enssize

         i_a  = mod(it_step,60)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

        if (i_a == 60 ) then

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(20)%data(ii,jj,1) = ps_ens_tmp_hour(20)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1)
                enddo
                enddo
           end do


          ps_ens_tau(heaada_enssize)%data(:,:,1) =0.

            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,jj,1) = ps_ens_tau(heaada_enssize)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(20)%data(ii,jj,1)
            enddo
            enddo
       write(*,*) 'it_step,ps_ens_tau',it_step,ps_ens_tau(20)%data(isd+10,jsd+10,1)
        endif

end subroutine mk_ens_ada_1h2



subroutine mk_ens_ada_3h2(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_3h = 0
         real ::  intv,real_heaada_enssize

         i_a  = mod(it_step,60)
         i_3h = mod(it_step,180)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) = ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

            write(*,*)  'i_ensb,ps_ens_tmp_hour',i_ensb,ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1)

            if (i_ensb ==3) i_ensb=0
         endif !i_a <= heahfq_minintv


        if (i_3h==0) then

          write(*,*) 'it_step,i_a,i_3h',it_step,i_a,i_3h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1)
                enddo
                enddo
           end do

          ps_ens_tau(heaada_enssize)%data(:,:,1) =0.

            do m=1,3
            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,jj,1) = ps_ens_tau(heaada_enssize)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(m)%data(ii,jj,1)*1./3.
            enddo
            enddo
            enddo
           write(*,*) 'it_step,i_ensb_3h,ps_ens_tau',it_step,i_ensb_3h,ps_ens_tau(20)%data(isd+10,jsd+10,1)
        endif

end subroutine mk_ens_ada_3h2


subroutine mk_ens_ada_6h2(atmos_var, ps_ens_tau,it_step,isd,ied, jsd, jed)
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed
!         integer,intent(in) :: ensemble_min_day
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize

         time_min_ens=360  !6h*60min
         i_a = mod(it_step,60)
         i_6h = mod(it_step,360)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do ii = isd,ied
                    ps_ens_tmp_min(i_a)%data(ii,jj,1)  = atmos_var%ps(ii,jj)
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           ps_ens_tmp_hour(i_ensb)%data(:,:,1) =0.

           do m=1,60
            do jj = jsd,jed
            do ii = isd,ied
                 ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) = ps_ens_tmp_hour(i_ensb)%data(ii,jj,1) +&
                                        ps_ens_tmp_min(m)%data(ii,jj,1)*1./60.
            enddo
            enddo
           enddo

!            write(*,*)  'i_ensb,ps_ens_tmp_hour',i_ensb,ps_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1)

            if (i_ensb ==6) i_ensb=0
         endif !i_a <= heahfq_minintv


        if (i_6h==0) then

!          write(*,*) 'it_step,i_a,i_6h',it_step,i_a,i_6h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,jj,1)=ps_ens_tau(m+1)%data(ii,jj,1)
                enddo
                enddo
           end do

          ps_ens_tau(heaada_enssize)%data(:,:,1) =0.

            do m=1,6
            do jj = jsd,jed
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,jj,1) =  ps_ens_tau(heaada_enssize)%data(ii,jj,1) + &
                                              ps_ens_tmp_hour(m)%data(ii,jj,1)*1./6.
            enddo
            enddo
            enddo

!           write(*,*) 'it_step,i_ensb_6h,ps_ens_tau',it_step,i_ensb_6h,ps_ens_tau(20)%data(isd+10,jsd+10,1)
        endif

end subroutine mk_ens_ada_6h2

   subroutine mk_ens_ada_3ht2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top,vari)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440
        character(len=1),intent(in) :: vari

         i_a  = mod(it_step,60)
         i_3h = mod(it_step,180)
         if(i_a==0) i_a=60

  select case(vari)
     case('t')
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo
   end select

!     write(*,*) 'tensemble_3h_day,time_min_ens',ensemble_min_day,time_min_ens,t_ens_tmp_min(1)%data(isd+10,jsd+10,1)
!     write(*,*) 'tit_step,i_a,t_ens_tmp_3h',it_step,i_a,t_ens_tmp_min(i_a)%data(isd+10,jsd+10,20)


         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo
!     write(*,*)  'ti_ensb,t_ens_tmp_3h',i_ensb,t_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
!                t_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)

            if (i_ensb == 3) i_ensb=0
         endif ! gather hour

        if (i_3h==0) then

!          write(*,*) 'it_step,i_a,i_3h',it_step,i_a,i_3h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,3
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                               t_ens_tmp_hour(m)%data(ii,jj,kk)*1./3.

            enddo
            enddo
            enddo
            enddo
!            write(*,*) 'it_3hstep,ps_ens_tau',it_step,ps_ens_tau(20)%data(isd+10,1,jsd+10)
         endif  !end i_3h ==0



   end subroutine mk_ens_ada_3ht2



   subroutine mk_ens_ada_6ht2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top,vari)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440
        character(len=1),intent(in) :: vari
!         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
!         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_6h = mod(it_step,360)
         if(i_a==0) i_a=60

  select case(vari)
     case('t')
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo
   end select
!     write(*,*) 'tensemble_6h_day,time_min_ens',ensemble_min_day,time_min_ens,t_ens_tmp_min(1)%data(isd+10,jsd+10,1)
!     write(*,*) 'tit_step,i_a,t_ens_tmp_6h',it_step,i_a,t_ens_tmp_min(i_a)%data(isd+10,jsd+10,20)

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo
!     write(*,*) 'ti_ensb,t_ens_tmp_6h',i_ensb,t_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
!                t_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)

            if (i_ensb == 6) i_ensb=0
         endif ! gather hour

        if (i_6h==0) then
!          write(*,*) 'it_6hstep,i_a,i_6h',it_step,i_a,i_6h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,6
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                                t_ens_tmp_hour(m)%data(ii,jj,kk)*1./6.

            enddo
            enddo
            enddo
            enddo
!           write(*,*) 'it_6hstep,i_ensb_6h,ps_ens_tau',it_step,i_ensb_6h,ps_ens_tau(20)%data(isd+10,1,jsd+10)
         endif  !end i_6h


   end subroutine mk_ens_ada_6ht2


   subroutine mk_ens_ada_3hu2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         i_a  = mod(it_step,60)
         i_3h = mod(it_step,180)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    u_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%u(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           u_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        u_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb == 3) i_ensb=0
         endif ! gather hour


        if (i_3h==0) then

!          write(*,*) 'it_step,i_a,i_3h',it_step,i_a,i_3h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

              ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,3
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                               u_ens_tmp_hour(m)%data(ii,jj,kk)*1./3.

            enddo
            enddo
            enddo
            enddo
!            write(*,*) 'it_3hustep,ps_ens_tau',it_step,ps_ens_tau(20)%data(isd+10,1,jsd+10)
         endif  !end i_3h ==0


   end subroutine mk_ens_ada_3hu2

   subroutine mk_ens_ada_3hv2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         i_a  = mod(it_step,60)
         i_3h = mod(it_step,180)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    v_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%v(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           v_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        v_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb == 3) i_ensb=0
         endif ! gather hour

        if (i_3h==0) then

!          write(*,*) 'it_step,i_a,i_3h',it_step,i_a,i_3h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

              ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,3
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                               v_ens_tmp_hour(m)%data(ii,jj,kk)*1./3.

            enddo
            enddo
            enddo
            enddo
!            write(*,*) 'it_3hvstep,ps_ens_tau',it_step,ps_ens_tau(20)%data(isd+10,1,jsd+10)
         endif  !end i_3h ==0


   end subroutine   mk_ens_ada_3hv2


    subroutine mk_ens_ada_6hu2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         i_a = mod(it_step,60)
         i_6h = mod(it_step,360)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    u_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%u(ii,kk,jj)
               enddo
               enddo
               enddo

!     write(*,*) 'uense_6h,time_min_ens',ensemble_min_day,time_min_ens,u_ens_tmp_min(1)%data(isd+10,jsd+10,1)
!     write(*,*) 'uit_step,i_a,u_ens_6h',it_step,i_a,u_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

! 2 gather t/u/v for hour 1-24
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           u_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        u_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo
            if (i_ensb ==6) i_ensb=0
         endif !gather hour

        if (i_6h==0) then
!          write(*,*) 'uit_6hstep,i_a,i_6h',it_step,i_a,i_6h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.
            do m=1,6
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                                u_ens_tmp_hour(m)%data(ii,jj,kk)*1./6.

            enddo
            enddo
            enddo
            enddo

!           write(*,*) 'it_6hstep,i_ensb_6h,u_tau',it_step,i_ensb_6h,ps_ens_tau(20)%data(isd+10,1,jsd+10)
         endif  !end i_6h


   end subroutine mk_ens_ada_6hu2

   subroutine mk_ens_ada_6hv2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440
         i_a = mod(it_step,60)
         i_6h = mod(it_step,360)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    v_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%v(ii,kk,jj)
               enddo
               enddo
               enddo


! 2 gather t/u/v for hour 1-24
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
           v_ens_tmp_hour(i_ensb)%data(:,:,:) =0.
           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        v_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb ==6) i_ensb=0
         endif !gather hour

        if (i_6h==0) then
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.
            do m=1,6
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                                v_ens_tmp_hour(m)%data(ii,jj,kk)*1./6.

            enddo
            enddo
            enddo
            enddo
         endif  !end i_6h


   end subroutine  mk_ens_ada_6hv2


   subroutine mk_ens_ada_dayt2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top,vari)
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


!   write(*,*) 'z_top',z_top

         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60

! 1 gather  t/u/v for min 1-60
  select case(vari)
     case('t')
     write(*,*) 't_collect'
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo


   end select


! 2 gather t/u/v for hour 1-24
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

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
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              t_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo

        endif


      write(*,*)'tps_ens_tau20',ps_ens_tau(heaada_enssize)%data(isd+10,1,jsd+10),ps_ens_tau(heaada_enssize)%data(isd+10,27,jsd+10)
    end subroutine mk_ens_ada_dayt2


   subroutine mk_ens_ada_dayu2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
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


         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    u_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%u(ii,kk,jj)
               enddo
               enddo
               enddo

     write(*,*) 'uensemble_min_day,time_min_ens',ensemble_min_day,time_min_ens,u_ens_tmp_min(1)%data(isd+10,jsd+10,1)
     write(*,*) 'uit_step,i_a,u_ens_tmp_min',it_step,i_a,u_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           u_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        u_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

     write(*,*) 'ui_ensb,u_ens_tmp_hour',i_ensb,u_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
                u_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day==0) then
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
          ps_ens_tau(heaada_enssize)%data(:,:,:) =0.

            do m=1,24
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              u_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo

        write(*,*) 'uit_step,u_ens_tau',it_step,ps_ens_tau(heaada_enssize)%data(isd+10,1,jsd+10),&
                ps_ens_tau(heaada_enssize)%data(isd+10,27,jsd+10)

        endif
 write(*,*)'ups_ens_tau20',ps_ens_tau(heaada_enssize)%data(isd+10,1,jsd+10),ps_ens_tau(heaada_enssize)%data(isd+10,27,jsd+10)

    end subroutine  mk_ens_ada_dayu2


   subroutine mk_ens_ada_dayv2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
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

         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    v_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%u(ii,kk,jj)
               enddo
               enddo
               enddo
     write(*,*) 'vensemble_min_day,time_min_ens',ensemble_min_day,time_min_ens,v_ens_tmp_min(1)%data(isd+10,jsd+10,1)
     write(*,*) 'vit_step,i_a,u_ens_tmp_min',it_step,i_a,v_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
            v_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        v_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo
     write(*,*) 'vi_ensb,v_ens_tmp_hour',i_ensb,v_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
                v_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day==0) then
          write(*,*) 'vit_step,i_a,i_day',it_step,i_a,i_day
           do m = 1, heaada_enssize-1
               do kk = 1,z_top
               do jj = jsd,jed
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)= ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
          end do

          ps_ens_tau(heaada_enssize)%data(:,:,:) =0.

            do m=1,24
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              v_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo

        write(*,*) 'vit_step,v_ens_tau',it_step,ps_ens_tau(heaada_enssize)%data(isd+10,1,jsd+10),&
                ps_ens_tau(heaada_enssize)%data(isd+10,27,jsd+10)

        endif

   end subroutine  mk_ens_ada_dayv2

   subroutine mk_ens_ada_5dayt2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top,vari)
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
        integer  :: i_5day

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         i_5day = mod(it_step,7200) !1440*5
         if(i_a==0) i_a=60

   select case(vari)
     case('t')
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo
   end select


     write(*,*) '5tensemble_min_day,time_min_ens',ensemble_min_day,time_min_ens,t_ens_tmp_min(1)%data(isd+10,jsd+10,1)
     write(*,*) '5tit_step,i_a,t_ens_tmp_min',it_step,i_a,t_ens_tmp_min(i_a)%data(isd+10,jsd+10,1)

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) '5tdayit_step,i_a,i_ensb',it_step,i_a,i_ensb
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
            enddo
           enddo
           
     write(*,*) '5ti_ensb,t_ens_tmp_hour',i_ensb,t_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
                t_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)
  
            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day == 0 ) then

         i_ensb_day  = i_ensb_day+1
         t_ens_tmp_day(i_ensb_day)%data(:,:,1) = 0.
           do m=1,24
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) = t_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) +&
                                        t_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
           enddo

         if(i_ensb_day == 5) i_ensb_day=0

        endif

      if(i_5day==0) then
           do m = 1, heaada_enssize-1
!                write(*,*) 'm',m,ps_ens_tau(m)%data(isd+10,1,jsd+10)
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,5
!                write(*,*) 'mm',ps_ens_tmp_day(m)%data(isd+10,jsd+10,1)
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              t_ens_tmp_day(m)%data(ii,jj,kk)*1./5.
            enddo
            enddo
            enddo
            enddo
           write(*,*)  '5tit_step,i_ensb_day,ps_ens_tau',it_step,i_ensb_day,ps_ens_tau(20)%data(isd+10,1,jsd+10)
        endif

  end subroutine mk_ens_ada_5dayt2

  subroutine mk_ens_ada_5dayu2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
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
        integer  :: i_5day

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         i_5day = mod(it_step,7200) !1440*5
         if(i_a==0) i_a=60
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    u_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%u(ii,kk,jj)

               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) '5udayit_step,i_a,i_ensb',it_step,i_a,i_ensb
           u_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = u_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        u_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
            enddo
           enddo
     write(*,*) '5ui_ensb,u_ens_tmp_hour',i_ensb,u_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
                v_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)
            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day == 0 ) then

         i_ensb_day  = i_ensb_day+1
         u_ens_tmp_day(i_ensb_day)%data(:,:,1) = 0.
           do m=1,24
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 u_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) = u_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) +&
                                        u_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo
         if(i_ensb_day == 5) i_ensb_day=0

        endif


      if(i_5day==0) then
           do m = 1, heaada_enssize-1
!                write(*,*) 'm',m,ps_ens_tau(m)%data(isd+10,1,jsd+10)
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,5
!                write(*,*) 'mm',ps_ens_tmp_day(m)%data(isd+10,jsd+10,1)
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              u_ens_tmp_day(m)%data(ii,jj,kk)*1./5.
            enddo
            enddo
            enddo
            enddo
           write(*,*) '5uit_step,i_ensb_5uday,ps_ens_tau',it_step,i_ensb_day,ps_ens_tau(20)%data(isd+10,1,jsd+10)
        endif



  end subroutine mk_ens_ada_5dayu2

  subroutine mk_ens_ada_5dayv2(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)

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
        integer  :: i_5day

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         i_5day = mod(it_step,7200) !1440*5
         if(i_a==0) i_a=60
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    v_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%v(ii,kk,jj)

               enddo
               enddo
               enddo


         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) '5vdayit_step,i_a,i_ensb',it_step,i_a,i_ensb
            v_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = v_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        v_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
            enddo
           enddo

     write(*,*) '5vi_ensb,v_ens_tmp_hour',i_ensb,v_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,1),&
                v_ens_tmp_hour(i_ensb)%data(isd+10,jsd+10,27)
            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv


        if (i_day == 0 ) then

         i_ensb_day  = i_ensb_day+1
         v_ens_tmp_day(i_ensb_day)%data(:,:,1) = 0.
           do m=1,24
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 v_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) = v_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) +&
                                        v_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo
         if(i_ensb_day == 5) i_ensb_day=0

        endif

      if(i_5day==0) then
           do m = 1, heaada_enssize-1
!                write(*,*) 'm',m,ps_ens_tau(m)%data(isd+10,1,jsd+10)
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do
               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,5
!                write(*,*) 'mm',ps_ens_tmp_day(m)%data(isd+10,jsd+10,1)
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              v_ens_tmp_day(m)%data(ii,jj,kk)*1./5.
            enddo
            enddo
            enddo
            enddo
           write(*,*) '5vit_step,i_ensb_5vday,ps_ens_tau',it_step,i_ensb_day,ps_ens_tau(20)%data(isd+10,1,jsd+10)
        endif

  end subroutine  mk_ens_ada_5dayv2

  subroutine mk_ens_ada_6hqr(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

        i_a = mod(it_step,60)
         i_6h = mod(it_step,360)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    q_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%q(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!   write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb

          q_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        q_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb == 6) i_ensb=0
         endif !

        if (i_6h==0) then
!          write(*,*) 'it_6hqstep,i_a,i_6h',it_step,i_a,i_6h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do


               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,6
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                                q_ens_tmp_hour(m)%data(ii,jj,kk)*1./6.

            enddo
            enddo
            enddo
            enddo
!   write(*,*)
!   'it_6hstep,i_ensb_6h,ps_ens_tau',it_step,i_ensb_6h,ps_ens_tau(20)%data(isd+10,1,jsd+10)
         endif  !end i_6h


   end subroutine mk_ens_ada_6hqr


  subroutine mk_ens_ada_5dayqr(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)

         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_day
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_day = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440
        integer  :: i_5day

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         i_5day = mod(it_step,7200) !1440*5
         if(i_a==0) i_a=60
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    q_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%q(ii,kk,jj)
               enddo
               enddo
               enddo
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) '5tdayit_step,i_a,i_ensb',it_step,i_a,i_ensb
           q_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        q_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
            enddo
           enddo

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day == 0 ) then

         i_ensb_day  = i_ensb_day+1
         q_ens_tmp_day(i_ensb_day)%data(:,:,1) = 0.
           do m=1,24
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 q_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) = q_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) +&
                                        q_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
           enddo
         if(i_ensb_day == 5) i_ensb_day=0

        endif

      if(i_5day==0) then
           do m = 1, heaada_enssize-1
!                write(*,*) 'm',m,ps_ens_tau(m)%data(isd+10,1,jsd+10)
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do
               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,5
!                write(*,*) 'mm',ps_ens_tmp_day(m)%data(isd+10,jsd+10,1)
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              q_ens_tmp_day(m)%data(ii,jj,kk)*1./5.
            enddo
            enddo
            enddo
            enddo
        endif

 end subroutine mk_ens_ada_5dayqr

  subroutine mk_ens_ada_dayqr(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_day
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_day = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440
         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    q_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%q(ii,kk,jj)
               enddo
               enddo
               enddo
         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        q_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv
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

          ps_ens_tau(heaada_enssize)%data(:,:,:) =0.
            do m=1,24
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              q_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo

        endif



  end subroutine mk_ens_ada_dayqr


  subroutine mk_ens_ada_3hqr(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         i_a  = mod(it_step,60)
         i_3h = mod(it_step,180)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    q_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%q(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           q_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = q_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        q_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb == 3) i_ensb=0
         endif ! gather hour

        if (i_3h==0) then

!          write(*,*) 'it_step,i_a,i_3h',it_step,i_a,i_3h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do
               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,3
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                               q_ens_tmp_hour(m)%data(ii,jj,kk)*1./3.

            enddo
            enddo
            enddo
            enddo
     endif  !end i_3h ==0

  end subroutine mk_ens_ada_3hqr


   subroutine mk_ens_ada_5dayt2r(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
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
        integer  :: i_5day

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         i_5day = mod(it_step,7200) !1440*5
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
            write(*,*) '5tdayit_step,i_a,i_ensb',it_step,i_a,i_ensb
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
            enddo
           enddo

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv
        if (i_day == 0 ) then

         i_ensb_day  = i_ensb_day+1
         t_ens_tmp_day(i_ensb_day)%data(:,:,1) = 0.
           do m=1,24
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) = t_ens_tmp_day(i_ensb_day)%data(ii,jj,kk) +&
                                        t_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
           enddo

         if(i_ensb_day == 5) i_ensb_day=0

        endif

      if(i_5day==0) then
           do m = 1, heaada_enssize-1
!                write(*,*) 'm',m,ps_ens_tau(m)%data(isd+10,1,jsd+10)
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,5
!                write(*,*) 'mm',ps_ens_tmp_day(m)%data(isd+10,jsd+10,1)
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              t_ens_tmp_day(m)%data(ii,jj,kk)*1./5.
            enddo
            enddo
            enddo
            enddo
        endif
    end subroutine mk_ens_ada_5dayt2r

    subroutine mk_ens_ada_dayt2r(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
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
         time_min_ens=ensemble_min_day*heaada_enssize !1440*20
         real_heaada_enssize = REAL(heaada_enssize) !20

         i_a = mod(it_step,60)
         i_day = mod(it_step,1440)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb ==24) i_ensb=0
         endif !i_a <= heahfq_minintv

        if (i_day==0) then
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
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                              t_ens_tmp_hour(m)%data(ii,jj,kk)*1./24.
            enddo
            enddo
            enddo
            enddo

        endif

    end subroutine mk_ens_ada_dayt2r

   subroutine mk_ens_ada_6ht2r(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_6h
         integer,  save ::  i_ensb = 0
         integer,  save ::  i_ensb_6h = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         i_a = mod(it_step,60)
         i_6h = mod(it_step,360)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo

            if (i_ensb == 6) i_ensb=0
         endif ! gather hour

        if (i_6h==0) then
!          write(*,*) 'it_6hstep,i_a,i_6h',it_step,i_a,i_6h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do

               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,6
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied

                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                                t_ens_tmp_hour(m)%data(ii,jj,kk)*1./6.

            enddo
            enddo
            enddo
            enddo
        endif
   end subroutine mk_ens_ada_6ht2r

 subroutine mk_ens_ada_3ht2r(config_flags,atmos_var,ps_ens_tau,it_step,isd,ied,jsd,jed,z_top)
         TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
         type(assim_data_type), intent(in) :: atmos_var
         type(ada_field_type), intent(inout) :: ps_ens_tau(:)
         integer, intent(in)  :: it_step,isd,ied, jsd, jed,z_top
         integer  :: n, m,taup1, i, i_a, i_b,ii,jj,kk,time_min_ens,i_3h
         integer,  save ::  i_ensb = 0
         real ::  intv,real_heaada_enssize
        integer, parameter :: ensemble_min_day=1440

         i_a  = mod(it_step,60)
         i_3h = mod(it_step,180)
         if(i_a==0) i_a=60

               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                    t_ens_tmp_min(i_a)%data(ii,jj,kk)  = atmos_var%t(ii,kk,jj)
               enddo
               enddo
               enddo

         if (i_a == 60 ) then
           i_ensb = i_ensb + 1 !1~24 hour
!            write(*,*) 'it_step,i_a,i_ensb',it_step,i_a,i_ensb
           t_ens_tmp_hour(i_ensb)%data(:,:,:) =0.

           do m=1,60
            do kk = 1,z_top
            do jj = jsd,jed
            do ii = isd,ied
                 t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) = t_ens_tmp_hour(i_ensb)%data(ii,jj,kk) +&
                                        t_ens_tmp_min(m)%data(ii,jj,kk)*1./60.
            enddo
            enddo
           enddo
           enddo
            if (i_ensb == 3) i_ensb=0
         endif ! gather hour

        if (i_3h==0) then

!          write(*,*) 'it_step,i_a,i_3h',it_step,i_a,i_3h
           do m = 1, heaada_enssize-1
               do jj = jsd,jed
               do kk = 1,z_top
               do ii = isd,ied
                  ps_ens_tau(m)%data(ii,kk,jj)=ps_ens_tau(m+1)%data(ii,kk,jj)
                enddo
                enddo
                enddo
           end do
               ps_ens_tau(heaada_enssize)%data(:,:,:) = 0.

            do m=1,3
            do jj = jsd,jed
            do kk = 1,z_top
            do ii = isd,ied
                ps_ens_tau(heaada_enssize)%data(ii,kk,jj) = ps_ens_tau(heaada_enssize)%data(ii,kk,jj) + &
                                               t_ens_tmp_hour(m)%data(ii,jj,kk)*1./3.

            enddo
            enddo
            enddo
            enddo
        endif

 end subroutine mk_ens_ada_3ht2r


 subroutine get_z_o(pb,p_top,p,gridi,gridj,p_lev,levm,& 
        isd, ied, jsd, jed, kms, kme,z_o2)
! calculate eta level for obs_station
    REAL,   INTENT(IN) ::  pb( isd:ied, kms:kme, jsd:jed ) 
    REAL,   INTENT(IN) ::  p_top
    REAL,   INTENT(IN) :: p( isd:ied, kms:kme, jsd:jed )  
    real,intent(in)  :: gridi,gridj,p_lev(6),levm(6)
      INTEGER , INTENT(IN) :: isd, ied, jsd, jed, kms, kme
        real :: dx_o,dy_o,i_o,j_o
        real :: pbase_o_md(1:kme),pp_o_md(1:kme)
        real :: p_obs
        real,intent(out) :: z_o2(6)
        integer :: k,n,z_o
        real,parameter :: staggerx_t = 0.5
        real,parameter :: staggery_t = 0.5
        integer :: nsta

  nsta = 6
        z_o2(:) = 999.0

    do  n = 1,nsta
        p_obs = p_lev(n)



        if (p_obs .eq. 999.0 .or. p_obs .lt. 5000) then
                z_o = 999.0                
                z_o2(n) = 999.0
                goto 208


         elseif ((isd <= gridi .AND. gridi <= ied) .AND. (jsd <= gridj .AND. gridj <= jed)) then

 
                dx_o = gridi - staggerx_t -FLOAT(int( gridi - staggerx_t))
                dy_o = gridj - staggery_t -FLOAT(int( gridj - staggery_t))
                i_o  = MIN0(MAX0(1,int(gridi-staggerx_t)),ied-1)
                j_o  = MIN0(MAX0(1,int(gridj-staggery_t)),jed-1)

                do k = 1, kme-1
                   pbase_o_md(k) = 0.01 *( (1- dy_o) * ( &
                        (1.- dx_o)*pb(i_o,k,j_o) +  dx_o *pb(i_o+1,k,j_o) ) +  &
                         dy_o* ( (1.-dx_o)*pb(i_o,k,j_o+1) +   &
                                  dx_o *pb(i_o+1,k,j_o+1) ) )

                   pp_o_md(k) = 0.01 *( (1- dy_o) * ( &
                           (1.- dx_o)*p(i_o,k,j_o) +  dx_o *p(i_o+1,k,j_o) ) +  &
                           dy_o* ( (1.-dx_o)*p(i_o,k,j_o+1) +   &
                                  dx_o *p(i_o+1,k,j_o+1) ) )                
                enddo

                do k = 1, kme-2
                     z_o = k



                        if ( (0.01*p_obs <= pbase_o_md(k)+pp_o_md(k)) .and. &
                                        ( 0.01*p_obs)>= pbase_o_md(k+1)+pp_o_md(k+1) ) then 

                   z_o2(n) = real(z_o+1)-  &
                        ( ( 0.01*p_obs - pbase_o_md(z_o+1)- pp_o_md(z_o+1)) / &
                        (pbase_o_md(z_o)+pp_o_md(z_o)-pbase_o_md(z_o+1)-pp_o_md(z_o+1)) )
                   z_o2(n)=max(z_o2(n),1.0)

                   endif
                enddo

        endif
208     continue


     enddo !n=6


 end subroutine  get_z_o

 subroutine get_obs_errf(grid,config_flags, &
                        Atmos_obs,gridi,gridj,p_lev,levm,tc,td,wind,wind_dir,sp_hum,&
         isd, ied, jsd, jed, kms, kme,z_o2,errf)

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    type(assim_data_type), intent(inout) :: Atmos_obs
    real,intent(out)  :: gridi,gridj
    real, intent(inout) :: p_lev(6),levm(6),tc(6),td(6),wind(6),wind_dir(6),sp_hum(6)
    integer, intent(in) :: isd, ied, jsd, jed, kms, kme
    real :: z_o2(6)
    real , intent(out)  :: errf(6)
    real ::  t_obs_loc,p_obs
    integer :: nsta,n
   real,parameter :: r_d = 287.04
   real,parameter :: cp  = 1004.67
   real :: obs_pt
    integer :: i,j,k
    integer :: i1, j1, k1, i2, j2, k2
    real :: dx, dy,dz
    real :: w1, w2, w3, w4, w5, w6, w7, w8 

   nsta = 6


   obs_pt = 999. !potential temperature

   do n=1,nsta

    t_obs_loc = tc(n)
    p_obs = p_lev(n)
    dz    = z_o2(n)



           if (p_obs .eq. 999.0 .or. p_obs .lt. 5000 .or. t_obs_loc .eq. 999. ) then
                z_o2(n) = 999.0
                goto 208


            elseif ((isd <= gridi .AND. gridi <= ied) .AND. (jsd <= gridj .AND. gridj <= jed)) then


              obs_pt =  t_obs_loc * (100000.0 / p_obs)**(r_d / cp) - 300



            i1 = floor(gridi)
            j1 = floor(gridj)
            k1 = floor(z_o2(n))
            i2 = i1 + 1
            j2 = j1 + 1
            k2 = k1 + 1

            dx = gridi - real(i1)
            dy = gridj - real(j1)
            dz = z_o2(n) - real(k1)
            w1 = (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
            w2 = dx * (1.0 - dy) * (1.0 - dz)
            w3 = (1.0 - dx) * dy * (1.0 - dz)
            w4 = dx * dy * (1.0 - dz)
            w5 = (1.0 - dx) * (1.0 - dy) * dz
            w6 = dx * (1.0 - dy) * dz
            w7 = (1.0 - dx) * dy * dz
            w8 = dx * dy * dz

            errf(n) =   w1 * Atmos_obs%t(i1, k1, j1) + w2 * Atmos_obs%t(i2, k1,j1) + &
                        w3 * Atmos_obs%t(i1, k1, k1) + w4 * Atmos_obs%t(i2,k1,j2) + &
                        w5 * Atmos_obs%t(i1, k2, k2) + w6 * Atmos_obs%t(i2,k2,j1) + &
                        w7 * Atmos_obs%t(i1, k2, k2) + w8 * Atmos_obs%t(i2,k2,j2) - obs_pt


         endif

208     continue
    enddo !nsta

 end subroutine get_obs_errf








end module module_ada_driver_fv_mod
