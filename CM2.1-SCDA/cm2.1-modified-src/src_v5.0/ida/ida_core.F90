module ida_core_mod

  use fms_mod, only : file_exist,read_data
  use mpp_mod, only : mpp_error, FATAL, NOTE, mpp_sum, stdout, mpp_sync_self
  use mpp_mod, only : mpp_pe, mpp_root_pe
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, &
       domain2d, mpp_get_global_domain, mpp_update_domains
  use time_manager_mod, only : time_type, operator( <= ), operator( - ), &
       operator( > ), operator ( < ),  set_time, set_date, get_date, get_time
!  use get_cal_time_mod, only : get_cal_time
  use axis_utils_mod, only : frac_index  

  use ida_types_mod, only : sice_obs_type, ida_grid_type

  use ada_types_mod, only : ada_grid_type

  use horiz_interp_type_mod, only: horiz_interp_type
  use horiz_interp_bilinear_mod, only : horiz_interp_bilinear_new

  implicit none

  real, private :: max_misfit = 5.0 ! this is used to inflate observation errors where
                                    ! the difference from the first guess is large
  real :: sice_obs_lats = 30.0 ! snz add to set obs domain
    
  integer, parameter, private :: max_sice_obs=100000
  integer, parameter, private :: max_files=10 
  integer :: max_obslvs = 5 ! for vd test

  real, parameter :: pai = 3.1415926
  real, parameter :: rad = pai/180.0

  type(sice_obs_type), target, private  :: sice_observations(max_sice_obs)  

  integer, private :: isc, iec, jsc, jec, isd, ied, jsd, jed  ! indices for local domain on model grid
  integer, private :: isg, ieg, jsg, jeg
  integer, private :: nk
  integer, private :: num_sice_obs

  ! time window for DROP, MOORING and SATELLITE data respectively

  type(time_type) , dimension(0:100), public :: time_window

  type(horiz_interp_type) :: Interp

  real, allocatable, dimension(:, :) :: x_grid, y_grid
  real :: lon_out(1, 1), lat_out(1, 1)

  logical :: use_cn_as_obs = .false.
  logical :: use_sst_as_obs = .false.
  logical :: use_iuv_as_obs = .false.
  logical :: use_sice_obs = .false.
  
  public :: copy_sice_obs, ida_core_init, get_sice_obs_up, get_sice_obs_dn, &
            put_sice_obs_dn
  
  contains

  subroutine init_sice_obs(time_s, time_e, localize)  

    use fms_mod, only : open_namelist_file,close_file,check_nml_error
    use mpp_io_mod, only : mpp_open, MPP_ASCII, MPP_RDONLY, MPP_MULTI, MPP_SINGLE
    use mpp_domains_mod, only : mpp_global_field
    
    type(time_type), intent(in) :: time_s, time_e
    logical, intent(in), optional :: localize

    integer :: cn_window = 1, iuv_window = 0, unknown_window = 30
    logical :: cn_obs, sst_obs, iuv_obs
    
    integer :: i,j, obs_variable
    integer :: CN_ID = 1, UI_ID = 2, VI_ID = 3, HI_ID = 4
    
    type obs_entry_type
       character(len=128) :: filename
       character(len=16)  :: file_type
    end type obs_entry_type

    namelist /sice_obs_nml/ cn_window, iuv_window, sice_obs_lats, &
                            cn_obs, sst_obs, iuv_obs

    character(len=128) :: input_files(max_files) = ''
    integer :: nfiles, filetype(max_files), ioun, io_status, ierr,&
                unit, nrecs, n
    character(len=256) :: record
    type(obs_entry_type) :: tbl_entry

    ioun = open_namelist_file()
    read(ioun,nml=sice_obs_nml,iostat = io_status)
    ierr = check_nml_error(io_status,'sice_obs_nml')
    call close_file(ioun)    

    if(sst_obs .or. iuv_obs) use_sice_obs = .true.
    use_cn_as_obs = cn_obs
    use_sst_as_obs = sst_obs
    use_iuv_as_obs = iuv_obs

    ! time window for DROP, MOORING and SATELLITE data respectively
    ! will be available from namelist

    time_window(:) = set_time(0,unknown_window)
    time_window(cn_window+1:cn_window+1) = set_time(0,cn_window)
    time_window(iuv_window:iuv_window+10) = set_time(0,iuv_window)

    ! get local indices for Model grid

    return

  end subroutine init_sice_obs

  subroutine get_sice_obs_dn(model_time, ens_size, sice_obs, num)

    use mpp_io_mod, only : mpp_open, mpp_close, mpp_read, &
        mpp_get_fields, mpp_get_axes, mpp_get_atts, mpp_get_info, &
        MPP_RDONLY, MPP_NETCDF, MPP_OVERWR, MPP_APPEND, &
        axistype, fieldtype, mpp_get_axis_data, MPP_MULTI, &
        MPP_SINGLE, mpp_get_times, mpp_copy_meta
    use axis_utils_mod, only : get_axis_cart, get_axis_bounds

    implicit none

#include <netcdf.inc>

    type(time_type), intent(in) :: model_time
    type(sice_obs_type), dimension(:), intent(inout) :: sice_obs
    integer, intent(in) :: ens_size
    integer, intent(inout) :: num

    integer :: k,yr,mon,day,hr,min,sec, kk, k_lev
    type(time_type) :: obs_time, tdiff

    real, dimension(:), allocatable   :: lon_src, lat_src
    type(axistype), dimension(:), allocatable, target :: axes
    type(axistype), pointer :: lon_axis, lat_axis, z_axis, t_axis
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_time, field_3d, field_2d
    integer   :: ndim,nvar,natt,ntime
    real      :: missing_value = 1.0e+20
    real, allocatable :: lon(:), lat(:)
    real, allocatable :: obs_cn(:,:),obs_u(:,:),obs_v(:,:), obs_t(:,:)
    integer :: nlon, nlat, nt
    integer :: iy0,in0,id0,ih0,im0,is0, n_days(12)
    integer :: ncid, rcode, start(3), nread(3), dims(3), id_cn, id_u, id_v, id_t

    integer::  unit, unit_out, i,j,n, nv, ii, jj, ii0, jj0, time_idx, i_m
    character(len=128) :: file_name, time_units, sst_filename
    character(len=32) :: fldname, axisname, anal_fldname

    data n_days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/ 

    file_name = 'INPUT/sice_obs.nc'
    call mpp_open(unit, trim(file_name), MPP_RDONLY, &
                  MPP_NETCDF, threading=MPP_MULTI, fileset = MPP_SINGLE)
    call mpp_get_info(unit,ndim,nvar,natt,ntime)
    allocate(axes(ndim))
    call mpp_get_axes(unit,axes)
    do i=1,ndim
      call mpp_get_atts(axes(i),name=axisname)
      select case (trim(axisname))
      case ('xt')
        lon_axis => axes(i)
      case ('yt')
        lat_axis => axes(i)
      case ('time')
        t_axis => axes(i)
      end select
    end do
    call mpp_get_atts(lon_axis,len=nlon)
    call mpp_get_atts(lat_axis,len=nlat)
    allocate(lon(nlon), lat(nlat))
    allocate(obs_cn(nlon,nlat), obs_u(nlon,nlat), obs_v(nlon,nlat), obs_t(nlon,nlat))

    call get_date(model_time, iy0,in0,id0,ih0,im0,is0)

!! monthly data use one line below
    time_idx = (iy0-1976)*12+in0
!! otherwise use codes below (daily or weekly?)

!    if(in0 == 1)then
!      time_idx = (iy0-1976)*365 + id0
!    else
!      time_idx = 0
!      do i_m = 1, in0-1
!        time_idx = time_idx + n_days(i_m)
!      end do
!      time_idx = (iy0-1976)*365 + time_idx + id0
!    end if

    if(mpp_pe() == mpp_root_pe())print*,'time_idx for sice=',time_idx

    rcode = nf_open(trim(file_name), NF_NOWRITE, ncid)
    if(rcode /= 0)print*,'rcode.nf_open=',rcode
    rcode = nf_inq_varid(ncid, 'EXT', id_cn)
    if(rcode /= 0)print*,'rcode.nf_inq_varid.cn=',rcode
    rcode = nf_inq_varid(ncid, 'UI', id_u)
    if(rcode /= 0)print*,'rcode.nf_inq_varid.iu=',rcode
    rcode = nf_inq_varid(ncid, 'VI', id_v)
    if(rcode /= 0)print*,'rcode.nf_inq_varid.iv=',rcode
    rcode = nf_inq_vardimid(ncid, id_cn, dims)
    if(rcode /= 0)print*,'rcode.nf_inq_vardimid=',rcode
    rcode = nf_inq_dimlen(ncid, dims(1), nlon)
    if(rcode /= 0)print*,'rcode.nf_inq_dimlen.lon=',rcode
    rcode = nf_inq_dimlen(ncid, dims(2), nlat)
    if(rcode /= 0)print*,'rcode.nf_inq_dimlen.lat=',rcode
    rcode = nf_inq_dimlen(ncid, dims(3), nt)
    if(rcode /= 0)print*,'rcode.nf_inq_dimlen.nt=',rcode

    if(nlon /= 360 .or. nlat /= 200) print*,'nlon,nlat in cn=',nlon,nlat

    start = 1
    nread(1) = nlon; nread(2) = nlat

    start(3) = time_idx
    nread(3) = 1
    rcode = nf_get_vara_double(ncid, id_cn, start, nread, obs_cn)
    if(rcode /= 0)print*,'rcode.sst=',rcode

    rcode = nf_get_vara_double(ncid, id_u, start, nread, obs_u)
    if(rcode /= 0)print*,'rcode.sst=',rcode

    rcode = nf_get_vara_double(ncid, id_v, start, nread, obs_v)
    if(rcode /= 0)print*,'rcode.sst=',rcode

! read sst
    if(in0 == 1)then
      time_idx = (iy0-1976)*365 + id0
    else
      time_idx = 0
      do i_m = 1, in0-1
        time_idx = time_idx + n_days(i_m)
      end do
      time_idx = (iy0-1976)*365 + time_idx + id0
    end if

    if(mpp_pe() == mpp_root_pe())print*,'time_idx for sst=',time_idx

    sst_filename = "INPUT/sst_daily.nc"
    rcode = nf_open(trim(sst_filename), NF_NOWRITE, ncid)
    if(rcode /= 0)print*,'rcode.nf_open=',rcode
    rcode = nf_inq_varid(ncid, 'SST', id_t)
    if(rcode /= 0)print*,'rcode.nf_inq_varid=',rcode
    rcode = nf_inq_vardimid(ncid, id_t, dims)
    if(rcode /= 0)print*,'rcode.nf_inq_vardimid=',rcode
    rcode = nf_inq_dimlen(ncid, dims(1), nlon)
    if(rcode /= 0)print*,'rcode.nf_inq_dimlen.lon=',rcode
    rcode = nf_inq_dimlen(ncid, dims(2), nlat)
    if(rcode /= 0)print*,'rcode.nf_inq_dimlen.lat=',rcode
    rcode = nf_inq_dimlen(ncid, dims(3), nt)
    if(rcode /= 0)print*,'rcode.nf_inq_dimlen.nt=',rcode

    if(nlon /= 360 .or. nlat /= 200) print*,'nlon,nlat in sst=',nlon,nlat

    start = 1
    nread(1) = nlon; nread(2) = nlat

    start(3) = time_idx
    nread(3) = 1
    rcode = nf_get_vara_double(ncid, id_t, start, nread, obs_t)
    if(rcode /= 0)print*,'rcode.sst=',rcode
! read sst done

    obs_time = set_date(iy0,in0,id0,ih0,im0,is0)

    num = 0
    do j = 1, nlat
    do i = 1, nlon

      if ( (obs_cn(i,j) >= 0.0 .and. obs_cn(i,j) <= 1.0) .and. &
           (abs(y_grid(i,j)) > sice_obs_lats) .and. abs(obs_t(i,j)) < 100.0 ) then
        num = num + 1
        sice_observations(num)%lon = x_grid(i,j)
	sice_observations(num)%lat = y_grid(i,j) 
	sice_observations(num)%flag = .true.
	if (sice_observations(num)%lat < 64.0) then
	  sice_observations(num)%i_index_dn = frac_index(sice_observations(num)%lon, x_grid(:,1))
	  sice_observations(num)%j_index_dn = frac_index(sice_observations(num)%lat, y_grid(90,:))
        else
          lon_out(1,1) = sice_observations(num)%lon*rad
          lat_out(1,1) = sice_observations(num)%lat*rad   
          call horiz_interp_bilinear_new (Interp, x_grid*rad, y_grid*rad, &
                 lon_out, lat_out)
          if (Interp%wti(1,1,2) < 1.0) then
	    sice_observations(num)%i_index_dn =Interp%i_lon(1,1,1) +&
	                                    Interp%wti(1,1,2)
	  else
	    sice_observations(num)%i_index_dn =Interp%i_lon(1,1,2)
	  endif
	  if (Interp%wtj(1,1,2) < 1.0) then
	    sice_observations(num)%j_index_dn =Interp%j_lat(1,1,1) +&
	                                    Interp%wtj(1,1,2)
	  else
	    sice_observations(num)%j_index_dn =Interp%j_lat(1,1,2)
	  endif
        end if
	if (use_cn_as_obs) then
	  sice_observations(num)%cn = obs_cn(i,j)
	  sice_observations(num)%flag_cn = .true.
        end if
	if (use_iuv_as_obs) then
	  sice_observations(num)%ui = obs_u(i,j)
	  sice_observations(num)%vi = obs_v(i,j)
	  sice_observations(num)%flag_uv = .true.
        end if
	if (use_sst_as_obs) then
          sice_observations(num)%ti = obs_t(i,j)
          sice_observations(num)%flag_ti = .true.
        end if
        if (sice_observations(num)%cn == 1.0 .and. &
            sice_observations(num)%ti > -1.5) then
!          print*,'lon,lat,sst=',sice_observations(num)%lon,&
!               sice_observations(num)%lat,sice_observations(num)%ti
          sice_observations(num)%ti = -1.76 
        end if
        if (sice_observations(num)%cn == 0.0 .and. &
            sice_observations(num)%ti < -1.5) then
!          print*,'lon,lat,sst=',sice_observations(num)%lon,&
!               sice_observations(num)%lat,sice_observations(num)%ti
          sice_observations(num)%ti = -1.0
        end if
        sice_observations(num)%time = obs_time
        sice_observations(num)%ens_size = ens_size

        if (ASSOCIATED(sice_observations(num)%enso_ti)) then
          DEALLOCATE(sice_observations(num)%enso_ti)
          NULLIFY(sice_observations(num)%enso_ti)
        end if
        ALLOCATE(sice_observations(num)%enso_ti(ens_size))
        sice_observations(num)%enso_ti(:) = 0.0
        if (ASSOCIATED(sice_observations(num)%enso_ui)) then
          DEALLOCATE(sice_observations(num)%enso_ui)
          NULLIFY(sice_observations(num)%enso_ui)
        end if
        ALLOCATE(sice_observations(num)%enso_ui(ens_size))
        sice_observations(num)%enso_ui(:) = 0.0
        if (ASSOCIATED(sice_observations(num)%enso_vi)) then
          DEALLOCATE(sice_observations(num)%enso_vi)
          NULLIFY(sice_observations(num)%enso_vi)
        end if
        ALLOCATE(sice_observations(num)%enso_vi(ens_size))
        sice_observations(num)%enso_vi(:) = 0.0

        call copy_sice_obs(sice_observations(num:num),sice_obs(num:num))

      end if

    end do
    end do
        
    CALL mpp_close(unit)
    return

  end subroutine get_sice_obs_dn

  subroutine put_sice_obs_dn(sice_obs_dn, num_dn)

    type(sice_obs_type), dimension(:), intent(in) :: sice_obs_dn
    integer, intent(in) :: num_dn

    integer :: i

    do i = 1, num_dn
      call copy_sice_obs(sice_obs_dn(i:i), sice_observations(i:i))
    end do
    num_sice_obs = num_dn

  end subroutine put_sice_obs_dn

  subroutine get_sice_obs_up(model_time, FV_grid, sice_obs_up, num_up)

    type(time_type), intent(in) :: model_time
    type(ada_grid_type), intent(in) :: FV_grid
    type(sice_obs_type), dimension(:), intent(inout) :: sice_obs_up
    integer, intent(inout) :: num_up

    real :: lon_up

    integer :: i

    num_up = 0
    do i = 1, num_sice_obs
      if (sice_observations(i)%flag) then
        num_up = num_up + 1
	lon_up = sice_observations(i)%lon
	if (lon_up < 0.0) lon_up = lon_up + 360.0
	if (lon_up > 360.0) lon_up = lon_up - 360.0
        sice_observations(i)%lon = lon_up
        sice_observations(i)%i_index_up = frac_index(lon_up, FV_grid%x(:))
	sice_observations(i)%j_index_up = frac_index(sice_observations(i)%lat, FV_grid%y(:))

	call copy_sice_obs(sice_observations(i:i), sice_obs_up(num_up:num_up))
      end if
    end do

    return

  end subroutine get_sice_obs_up

  subroutine ida_core_init(sice_grid, time_s, time_e, localize)

    use fms_mod, only : open_namelist_file, check_nml_error, close_file

    type(ida_grid_type), intent(in) :: sice_grid
    logical, intent(in), optional :: localize
    type(time_type), intent(in) :: time_s, time_e

    integer :: ioun, ierr, io_status, i, j
    namelist /ida_core_nml/ max_misfit

    ioun = open_namelist_file()
    read(ioun,nml=ida_core_nml,iostat = io_status)
    ierr = check_nml_error(io_status,'ida_core_nml')
    call close_file(ioun)

    allocate(x_grid(sice_grid%ni,sice_grid%nj), y_grid(sice_grid%ni,sice_grid%nj))
    x_grid(:,:) = sice_grid%x_ocn(:,:)
    y_grid(:,:) = sice_grid%y_ocn(:,:)

    do j = 1, sice_grid%nj
    do i = 1, sice_grid%ni
      if (x_grid(i,j) < 80.5) x_grid(i,j) = x_grid(i,j) + 360.0
    end do
    end do

    call init_sice_obs(time_s, time_e, localize)

  end subroutine ida_core_init

  subroutine copy_sice_obs(obs_in, obs_out)

    type(sice_obs_type), dimension(:), intent(in) :: obs_in
    type(sice_obs_type), dimension(:), intent(inout) :: obs_out

    integer :: n

    if (size(obs_in) .ne. size(obs_out)) call mpp_error(FATAL)

    do n=1,size(obs_in)
       Obs_out(n)%lon = Obs_in(n)%lon
       Obs_out(n)%lat = Obs_in(n)%lat
       Obs_out(n)%time = Obs_in(n)%time
       Obs_out(n)%ens_size = Obs_in(n)%ens_size
       Obs_out(n)%i_index_up = Obs_in(n)%i_index_up
       Obs_out(n)%j_index_up = Obs_in(n)%j_index_up
       Obs_out(n)%i_index_dn = Obs_in(n)%i_index_dn
       Obs_out(n)%j_index_dn = Obs_in(n)%j_index_dn
       Obs_out(n)%i_index_up_uv = Obs_in(n)%i_index_up_uv
       Obs_out(n)%j_index_up_uv = Obs_in(n)%j_index_up_uv
       Obs_out(n)%i_index_dn_uv = Obs_in(n)%i_index_dn_uv
       Obs_out(n)%j_index_dn_uv = Obs_in(n)%j_index_dn_uv
       Obs_out(n)%cn = Obs_in(n)%cn
       Obs_out(n)%ui = Obs_in(n)%ui
       Obs_out(n)%vi = Obs_in(n)%vi
       Obs_out(n)%ti = Obs_in(n)%ti
       Obs_out(n)%flag = Obs_in(n)%flag
       Obs_out(n)%flag_cn = Obs_in(n)%flag_cn
       Obs_out(n)%flag_uv = Obs_in(n)%flag_uv
       Obs_out(n)%flag_ti = Obs_in(n)%flag_ti

       if (ASSOCIATED(Obs_out(n)%enso_ti)) then
         DEALLOCATE(Obs_out(n)%enso_ti)
	 NULLIFY(Obs_out(n)%enso_ti)
       end if
       ALLOCATE(Obs_out(n)%enso_ti(Obs_in(n)%ens_size))
       Obs_out(n)%enso_ti(:) = Obs_in(n)%enso_ti(:)
       if (ASSOCIATED(Obs_out(n)%enso_ui)) then
         DEALLOCATE(Obs_out(n)%enso_ui)
	 NULLIFY(Obs_out(n)%enso_ui)
       end if
       ALLOCATE(Obs_out(n)%enso_ui(Obs_in(n)%ens_size))
       Obs_out(n)%enso_ui(:) = Obs_in(n)%enso_ui(:)
       if (ASSOCIATED(Obs_out(n)%enso_vi)) then
         DEALLOCATE(Obs_out(n)%enso_vi)
	 NULLIFY(Obs_out(n)%enso_vi)
       end if
       ALLOCATE(Obs_out(n)%enso_vi(Obs_in(n)%ens_size))
       Obs_out(n)%enso_vi(:) = Obs_in(n)%enso_vi(:)
    end do

  end subroutine copy_sice_obs

end module ida_core_mod
