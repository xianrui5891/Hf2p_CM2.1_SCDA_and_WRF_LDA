module module_ada_core_mod

  ! ADA modules
  USE module_utility
  use module_ada_types_mod, only : ada_grid_type, atmos_obs_type, assim_data_type
!  use bgrid_vert_mod, only: vert_grid_type, compute_pres_full
!  use ada_driver_fv_mod, only : assim_frequency !lulv
     USE module_domain, ONLY : domain
  USE netcdf
  USE module_configure, ONLY : grid_config_rec_type
  implicit none

!#include <netcdf.inc>
!include '/home/PublicSoftware/intel2017-software/netcdf-4.4.1.1-parallel/include/netcdf.inc'

  private
  public ::   get_ncep_field,get_obs_real

  integer, parameter :: MAX_ATM_OBS = 1 !x1y 15000
  integer, parameter :: MAX_FILES = 10 
  integer, parameter :: MAX_OBSLVS = 24

  ! ADA_CORE_NML namelist variables
  real, private :: max_misfit = 5.0 !< this is used to inflate observation errors where the difference from the first guess is large
  real, private :: ass_start_lat = -90.0
  real, private :: ass_end_lat = 90.0 ! snz add to set obs domain
!  integer :: atm_data_freq = 6 !lulv
  real, private, allocatable, dimension(:) :: x_grid, y_grid
    
  integer, private :: isc, iec, jsc, jec, isd, ied, jsd, jed !< indices for local domain on model grid
  integer, private :: isg, ieg, jsg, jeg
  integer, private :: nk, iy0_prev, idx_prev

  logical :: ada_core_initialized = .false.
  logical :: atm_obs_initialized = .false.
  logical :: use_temp_as_obs = .false.
  logical :: use_ps_as_obs = .false.
  logical :: use_uv_as_obs = .false.
  logical :: use_q_as_obs = .false.

  type(atmos_obs_type), target, private  :: atmos_observations(MAX_ATM_OBS)  
!  type(time_type) , dimension(0:100), public :: time_window !< time window for DROP, MOORING and SATELLITE data respectively
 
contains


  subroutine open_profile_dataset()
    ! Placeholder routine
  end subroutine open_profile_dataset



!lulv comment  subroutine get_ncep_field(model_time, atm_var)
!  subroutine get_ncep_field(model_time, atm_var, assim_frequency,T_grid)
   subroutine get_ncep_field(grid,config_flags,model_time, atm_var, assim_frequency,atm_step)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    TYPE(WRFU_Time), intent(in), target        ::  model_time
    type(assim_data_type), intent(inout) :: atm_var
    integer,intent(in) :: atm_step
!    type(ada_grid_type), intent(inout) :: T_grid
    integer, intent(in) :: assim_frequency !lulv add
    
!    integer :: unit, stdout_unit, istat !< mpp_open unit number
    integer :: ndim, nvar, natt, ntime,ncid2,ncid !< Variables for netcdf variable sized
    integer ::  nlev !< variable holds grid information
    integer :: i, i_m,ii,jj,zz !< counting/loop varaibles
    integer :: time_idx,obs_time_idx !time_idx is model time_idx; obs_time_idx is used for read obs.nc
    integer,save :: time_idx0 
    ! Date variables
    integer :: iy0, in0, id0, ih0, im0, is0 !< date variables
    INTEGER                         :: yr, day, hr, minute, sec, rc, hr3
    real :: hr2
    integer, dimension(12) :: n_days !< num days in month

    character(len=128) :: file_name
    character(len=32) :: axisname

        integer ierr1,ierr2,ierr3,ierr4,varid,ierr5,ierr6,ierr7,ierr8
        integer ierr9,ierr10,ierr11,ierr12,ierr13,ierr14,ierr15,ierr16
        integer,parameter :: nlon=560, nlat=410 !,ntime1=1000
    real, pointer, dimension(:,:) :: ps => NULL()
    real, pointer, dimension(:,:) :: phis => NULL()
    real, pointer, dimension(:,:,:) ::  temp => NULL(),uobs => NULL(),vobs =>NULL()
    character(len=200) :: file_obs1,file_obs2,file_obs3,file_obs
    character(len=200) :: path_file,file_t_name,file_t_obs
    character(len=200) :: file_ps_name,file_ps_obs,file_u_name,file_u_obs,file_v_name,file_v_obs

! sth for num2str and delete left space
    character(5) num_obs_time_idx1
    character(len=5) num_obs_time_idx2

 
        allocate(ps(1:nlon,1:nlat))
        allocate(phis(1:nlon,1:nlat))
        allocate(temp(1:nlon,1:nlat,1:27))
        allocate(uobs(1:nlon,1:nlat,1:27))
        allocate(vobs(1:nlon,1:nlat,1:27))

!file_obs1="/home/lmkli/wangj/test/WRF371/v3_0509/fortran_readnc/era5_atmos_obs_1.nc"
!file_obs2="/home/lmkli/wangj/test/WRF371/v3_0509/fortran_readnc/era5_atmos_obs_2.nc"
!file_obs3="/home/lmkli/wangj/test/WRF371/v3_0509/fortran_readnc/era5_atmos_obs_3.nc"
!file_obs1="/home/lmkli/wangj/test/WRF371/v3_0509/fortran_readnc/era5_atmos_obs_234.nc"
path_file="/public/home/hou/wangj/obs/obs_all/"
file_t_name='interp_t559_to_t560_'
file_ps_name='interp_psfc559_to_560_'
file_u_name='interp_u409_to_410_'
file_v_name='interp_v559_to_560_'


        CALL WRFU_TimeGet( model_time, YY=yr, dayOfYear=day, H=hr, M=minute,S=sec,rc=rc)
      iy0=yr
!    stdout_unit = stdout()
!w

    if(assim_frequency==6) then!6hourly
        hr2 = hr/6
        hr2 = floor(hr2) 
        hr3 = int(hr2) 
        if(hr == 0) time_idx = (day-1)*4 + 1
        if(hr /=0)  time_idx = (day-1)*4 + hr3 + 1
        if (day == 366 .and. hr==23) then
                iy0_prev = yr
                idx_prev = time_idx
        endif 
        if( time_idx ==1 .and. iy0_prev/=0 .and. yr>iy0_prev) then
                time_idx = idx_prev +1
        endif

        elseif(assim_frequency==1) then  ! 1 hourly
                if(hr == 0) time_idx = (day-1)*24 + 1
                if(hr /=0)  time_idx = (day-1)*24 + hr +1
                if (day == 366 .and. hr==23) then
                        iy0_prev = yr
                        idx_prev = time_idx
                endif
        if( time_idx ==1 .and. iy0_prev/=0 .and. yr>iy0_prev) then
                time_idx = idx_prev +1
        endif
    endif


!   write(*,*) 'obs_1',atm_step,time_idx,yr,day,hr,minute


    if(atm_step==360) then
         time_idx0=time_idx
         obs_time_idx = 2
         time_idx = 2
    endif

      if(atm_step/=360) then
          obs_time_idx = (time_idx-time_idx0)+2
          time_idx = (time_idx-time_idx0)+2
      endif



!   write(*,*) 'obs_2',atm_step,time_idx0,yr,day,hr,minute,obs_time_idx,time_idx

        if (assim_frequency==6 .or. assim_frequency==1) then 
!! 1000/6=166.*  2000/6=333  if obs interval is 1 hourly
!! now obs interval is 6 hourly
        !        if ( time_idx .gt. 1 .and. time_idx .lt. 166 )     file_obs = file_obs1
        !        if ( time_idx .gt. 167 .and. time_idx .lt. 333 )   file_obs = file_obs2
        !        if ( time_idx .gt. 334 )                           file_obs = file_obs3
              if ( time_idx >= 1 .and. time_idx <= 500 )      file_obs =file_obs1
                if ( time_idx >= 501 .and. time_idx <= 1000 )  then
                         file_obs =file_obs2
                         obs_time_idx = time_idx-500  !defined by file_obs2 time variable
                endif
                if ( time_idx >= 1001 .and. time_idx <= 1500 )  then
                        file_obs =file_obs3
                        obs_time_idx = time_idx-1000  
                ENDIF
!                if ( time_idx .gt. 1501 .and. time_idx .lt. 2000 )  then
!                        file_obs =file_obs4
!                        obs_time_idx = time_idx-1500
!                ENDIF
!                if ( time_idx .gt. 2001 .and. time_idx .lt. 2500 )  then
!                        file_obs =file_obs5
!                        obs_time_idx = time_idx-2000
!                ENDIF

!                if ( time_idx .gt. 2501 )        then
!                         file_obs =file_obs6                
!                          obs_time_idx = time_idx-2500
!                endif

        endif


!        if (assim_frequency==1) then

!              if ( time_idx >= 1 .and. time_idx <= 500 )      file_obs=file_obs1
!                if ( time_idx >= 501 .and. time_idx <= 1000 )  then
!                         file_obs =file_obs2
!                         obs_time_idx = time_idx-500  !defined by file_obs2 timevariable
!                endif
!                if ( time_idx >= 1001 .and. time_idx <= 1500 )  then
!                        file_obs =file_obs3
!                        obs_time_idx = time_idx-1000
!                ENDIF


!        endif

 
!        obs_time_idx = (time_idx-time_idx0)+1           !used for restart 

        if(atm_step==1) obs_time_idx=1


!         do ii=1,nlon
!         do jj=1,nlat 
!          atm_var%ps(ii,jj)=0.0
!          atm_var%phis(ii,jj)=0.0
!          ps(ii,jj)=0.0
!          phis(ii,jj) = 0.0
!        enddo
!        enddo


          atm_var%u_2(:,:,:)=0.0
          atm_var%t_2(:,:,:)=0.0

          atm_var%ps(:,:)=0.0
          atm_var%phis(:,:)=0.0
          ps(:,:)=0.0
          phis(:,:) = 0.0

        atm_var%t(:,:,:) = 0.        
        atm_var%u(:,:,:) = 0.
        atm_var%v(:,:,:) = 0.
        temp(:,:,:) = 0.
        uobs(:,:,:) = 0.
        vobs(:,:,:) = 0.

!
! change num into character 
!
     write(num_obs_time_idx1,'(i4)') obs_time_idx
     num_obs_time_idx2=adjustl(num_obs_time_idx1)


     file_ps_obs = trim(path_file)//trim(file_ps_name)//trim(num_obs_time_idx2)//'.nc'
     file_t_obs=trim(path_file)//trim(file_t_name)//trim(num_obs_time_idx2)//'.nc' 
     file_u_obs=trim(path_file)//trim(file_u_name)//trim(num_obs_time_idx2)//'.nc'
     file_v_obs=trim(path_file)//trim(file_v_name)//trim(num_obs_time_idx2)//'.nc'

! read psfc

        ierr1=nf90_open(trim(file_ps_obs),nf90_nowrite,ncid)
        ierr2=nf90_inq_varid(ncid,'PSFC',varid)
        ierr3=nf90_get_var(ncid,varid,ps)
        ierr4 = nf90_close(ncid)

!         do ii=1,nlon
         do jj=1,nlat
         do ii=1,nlon
                 atm_var%ps(ii,jj) = ps(ii,jj)
                 atm_var%u_2(ii,1,jj) = ps(ii,jj)
         enddo
         enddo

!    write(*,*) 'obs_value1',atm_var%ps(560,410),atm_var%ps(560,405),atm_var%ps(550,410)
!    write(*,*) 'obs_value2',atm_var%ps(560,410),atm_var%ps(560,405),atm_var%ps(550,410)

! read t
        ierr5=nf90_open(trim(file_t_obs),nf90_nowrite,ncid)
        ierr6=nf90_inq_varid(ncid,'T',varid)
        ierr7=nf90_get_var(ncid,varid,temp)
        ierr8 = nf90_close(ncid)

        do zz=1,27
        do jj=1,nlat
        do ii=1,nlon        
        atm_var%t(ii,zz,jj) = temp(ii,jj,zz)
        enddo
        enddo
        enddo

!    write(*,*) 'obs_4',atm_step,obs_time_idx,time_idx,ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7,ierr8
!    write(*,*) 'file_t_obs',file_t_obs,file_ps_obs
!         write(*,*) 'file_u_obs',file_u_obs
!read u
        ierr9=nf90_open(trim(file_u_obs),nf90_nowrite,ncid)
        ierr10=nf90_inq_varid(ncid,'U',varid)
        ierr11=nf90_get_var(ncid,varid,uobs)
        ierr12 = nf90_close(ncid)

!    write(*,*) 'obs_5',atm_step,obs_time_idx,ierr9,ierr10,ierr11,ierr12
    


        do zz=1,27
        do jj=1,nlat
        do ii=1,nlon
        atm_var%u(ii,zz,jj) = uobs(ii,jj,zz)
        enddo
        enddo
        enddo
!    write(*,*) 'u_value1',atm_var%u(560,410,1),atm_var%u(560,405,1),atm_var%u(550,410,1)
!    write(*,*) 'u_value2',atm_var%u(560,410,1),atm_var%u(560,405,1),atm_var%u(550,410,1)


!read v
        ierr13=nf90_open(trim(file_v_obs),nf90_nowrite,ncid)
        ierr14=nf90_inq_varid(ncid,'V',varid)
        ierr15=nf90_get_var(ncid,varid,vobs)
        ierr16 = nf90_close(ncid)

        do zz=1,27
        do jj=1,nlat
        do ii=1,nlon
        atm_var%v(ii,zz,jj) = vobs(ii,jj,zz)
        atm_var%t_2(ii,zz,jj) = vobs(ii,jj,zz)
        enddo
        enddo
        enddo

!        write(*,*) 'file_v_obs',file_v_obs

!    write(*,*) 'v_value1',atm_var%v(560,410,1),atm_var%v(560,405,1),atm_var%v(550,410,1)
!    write(*,*) 'v_value2',atm_var%v(560,410,1),atm_var%v(560,405,1),atm_var%v(550,410,1)
!    write(*,*) 'obs_6',atm_step,obs_time_idx,ierr13,ierr14,ierr15,ierr16

      deallocate(ps)
      deallocate(phis)
      deallocate(temp)      
      deallocate(uobs)
      deallocate(vobs)
  end subroutine get_ncep_field


  subroutine get_obs_real(grid,config_flags,model_time, atm_var,assim_frequency,atm_step, gridi,gridj,p_lev,levm,&
                          tc,td,wind,wind_dir,sp_hum,lat1,lon1)
    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    TYPE(WRFU_Time), intent(in), target        ::  model_time
    type(assim_data_type), intent(inout) :: atm_var
    integer,intent(in) :: atm_step
    integer, intent(in) :: assim_frequency !lulv add
    integer :: ndim, nvar, natt, ntime,ncid2,ncid !< Variables for netcdfvariable sized
    integer ::  nlev !< variable holds grid information
    integer :: i, i_m,ii,jj,zz !< counting/loop varaibles
    integer :: time_idx,obs_time_idx !time_idx is model time_idx; obs_time_idxis used for read obs.nc
    integer,save :: time_idx0
    integer :: iy0, in0, id0, ih0, im0, is0 !< date variables
    INTEGER                         :: yr, day, hr, minute, sec, rc, hr3
    real :: hr2
    real,intent(out)  :: gridi,gridj,lat1,lon1
    real, intent(inout) :: p_lev(6),levm(6),tc(6),td(6),wind(6),wind_dir(6),sp_hum(6)
    integer, dimension(12) :: n_days !< num days in month

    character(len=128) :: file_name
    character(len=32) :: axisname
        integer ierr1,ierr2,ierr3,ierr4,varid,ierr5,ierr6,ierr7,ierr8
        integer ierr9,ierr10,ierr11,ierr12,ierr13,ierr14,ierr15,ierr16,open_num
        integer,parameter :: nlon=560, nlat=410 !,ntime1=1000
    real, pointer, dimension(:,:) :: ps => NULL()
    real, pointer, dimension(:,:) :: phis => NULL()
    real, pointer, dimension(:,:,:) ::  temp => NULL(),uobs => NULL(),vobs=>NULL()
    character(len=200) :: file_obs1,file_obs2,file_obs3,file_obs
    character(len=200) :: path_file,file_t_name,file_t_obs
    character(len=200) :: file_obs_r_name
    character(len=10) :: char_date
    character(5) num_obs_time_idx1
    character(len=5) num_obs_time_idx2
        real :: es(6), e(6), t(6), td_t(6), t_t(6) 
       CALL WRFU_TimeGet( model_time, YY=yr, dayOfYear=day, H=hr,M=minute,S=sec,rc=rc)

    if(assim_frequency==6) then!6hourly
        hr2 = hr/6
        hr2 = floor(hr2)
        hr3 = int(hr2)
        if(hr == 0) time_idx = (day-1)*4 + 1
        if(hr /=0)  time_idx = (day-1)*4 + hr3 + 1
        if (day == 366 .and. hr==23) then
                iy0_prev = yr
                idx_prev = time_idx
        endif
        if( time_idx ==1 .and. iy0_prev/=0 .and. yr>iy0_prev) then
                time_idx = idx_prev +1
        endif

        elseif(assim_frequency==1) then  ! 1 hourly
                if(hr == 0) time_idx = (day-1)*24 + 1
                if(hr /=0)  time_idx = (day-1)*24 + hr +1
                if (day == 366 .and. hr==23) then
                        iy0_prev = yr
                        idx_prev = time_idx
                endif
        if( time_idx ==1 .and. iy0_prev/=0 .and. yr>iy0_prev) then
                time_idx = idx_prev +1
        endif
    endif

    if(atm_step==360) then
         time_idx0=time_idx
         obs_time_idx = 2
         time_idx = 2
    endif

      if(atm_step/=360) then
          obs_time_idx = (time_idx-time_idx0)+2
          time_idx = (time_idx-time_idx0)+2
      endif

       if (assim_frequency==6 .or. assim_frequency==1) then
              if ( time_idx >= 1 .and. time_idx <= 500 )      file_obs =file_obs1
                if ( time_idx >= 501 .and. time_idx <= 1000 )  then
                         file_obs =file_obs2
                         obs_time_idx = time_idx-500  !defined by file_obs2 time variable
                endif
                if ( time_idx >= 1001 .and. time_idx <= 1500 )  then
                        file_obs =file_obs3
                        obs_time_idx = time_idx-1000
                ENDIF
        endif

        if(atm_step==1) obs_time_idx=1
     write(num_obs_time_idx1,'(i4)') obs_time_idx
     num_obs_time_idx2=adjustl(num_obs_time_idx1)

 
    path_file = "/public/home/hou/wangj/WRFV3_mda/run/"
    file_obs_r_name = trim(path_file)//'OBS_1_'//trim(num_obs_time_idx2)//'.txt' 

!    open_num = 177+obs_time_idx
        open_num = 177

      open(unit=open_num,file=file_obs_r_name,form='formatted',status='old')
      read(open_num,201) char_date
 201  FORMAT(1x,a10)
      read(open_num,202) lat1,lon1,gridj,gridi
 202  FORMAT(1x,4(f9.4,1x))


      write(*,*) 'nameobs',file_obs_r_name
        write(*,*) char_date
        write(*,*)  lat1,lon1,gridj,gridi
 

    DO i=1,6
      read(open_num,203) p_lev(i),levm(i),tc(i),td(i),wind(i),wind_dir(i)
 203 FORMAT(1x,6(f11.3,1x))


!  write(*,*) p_lev(i),levm(i),tc(i),td(i),wind(i),wind_dir(i)

! ! ------- cal specific_humidity ------------------------------------
        t(i) = tc(i) + 273.15! convert Celsius to Kelvin
        td_t(i) = td(i)
       ! Saturation vapor pressure calculation for temperature t
        es(6) = 6.112 * exp((17.67 * t(i)) / (t(i) + 243.5))
       ! Saturation vapor pressure calculation for dew point temperature td_t
        e(i) = 6.112 * exp((17.67 * td_t(i)) / (td_t(i) + 243.5))
        sp_hum(i) = 0.622 * e(i) / (p_lev(i) - 0.378 * e(i) )
! ! ------- end cal specific_humidity  wangjj------------------------------------

     END DO


      close(unit=open_num)

  end subroutine get_obs_real



   
end module module_ada_core_mod
