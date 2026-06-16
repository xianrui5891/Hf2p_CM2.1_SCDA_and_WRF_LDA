module module_assim_tools_mod
  ! A variety of operations required by assimilation.

!  use fms_mod, only : file_exist, open_namelist_file, check_nml_error, write_version_number, close_file
!  use fms_mod, only : error_mesg, WARNING, FATAL
!  use mpp_mod, only : mpp_pe,  stdlog, stdout
USE module_loc_and_dist_mod, only : loc_type
USE module_configure, ONLY : grid_config_rec_type
  logical :: module_initialized = .false.

  private
  public assim_tools_init, obs_increment, update_from_obs_inc, obs_increment_prf_eta_hyb,&
       & update_from_obs_inc_eta_hyb, update_from_obs_inc_prf_hyb, update_from_inc_bt_hyb, obs_increment_ps,&
       & obs_increment_prf_eta_hyb_stn, obs_increment_hfq,update_from_obs_inc_k !lulv add *stn

  !---- namelist ASSIM_TOOLS_NML with default values

  real :: cor_cutoff = 0.05
  real :: obs_diff_cutoff = 100.0

!  namelist /assim_tools_nml/ cor_cutoff, obs_diff_cutoff
! wangj  if  use above line,remember to add it

  !---- module name and version number
  character(len=*), parameter :: MODULE_NAME = 'assim_tools_mod'
  character(len=*), parameter :: VERS_NUM = '$Id: assim_tools.F90,v 1.1.2.1 2013/09/04 14:05:20 Seth.Underwood Exp $'
  INTEGER :: debug_level
   INCLUDE 'mpif.h'

contains

!  subroutine assim_tools_init()

!    integer :: unit, istat, stdlog_unit, stdout_unit

!    if ( module_initialized ) return

    ! Read namelist for run time control
!    if ( file_exist('input.nml') ) then
!       unit = open_namelist_file()
!       read(UNIT=unit, NML=assim_tools_nml, IOSTAT=istat)
!       call close_file(unit)
!    else
       ! Set istat to an arbitrary positive number if input.nml does not exist
!       istat = 100
!    end if

!    if ( check_nml_error(istat, 'assim_tools_nml') < 0 ) then
!       call error_mesg('assim_tools_mod::assim_tools_init', 'ASSIM_TOOLS_NML not found in input.nml,  Using defaults.', WARNING)
!    end if


!    stdout_unit = stdout()
!    stdlog_unit = stdlog()
!    call write_version_number(VERS_NUM, MODULE_NAME)
!    write (UNIT=stdlog_unit, NML=assim_tools_nml)
!    write (UNIT=stdout_unit, NML=assim_tools_nml)

!    module_initialized = .true.
!  end subroutine assim_tools_init

  subroutine obs_increment_prf_eta_hyb(ens, ens_size, obs, obs_var, obs_inc_eakf, obs_var_oi, obs_inc_oi, std_oi)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens
    real, intent(in) :: obs, obs_var, obs_var_oi
    real, intent(inout), dimension(ens_size) :: obs_inc_eakf, obs_inc_oi
    real, intent(in) :: std_oi

    real :: a, cov, obsmmean
    real :: mean, new_cov, new_mean
    real :: var_oi_bg

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    ! Default return values for obs_inc_eakf and obs_inc_oi
    obs_inc_eakf = 0.0
    obs_inc_oi = 0.0

    var_oi_bg = std_oi * std_oi

    ! Need to compute the prior covariance and mean; copy to 2d for now but
    ! get rid of this step at some point; just fits interface
    mean = sum(ens) / ens_size

    cov = 0.0
    do ie=1, ens_size
       cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
    end do
    cov = cov / ens_size

    if ( cov == 0.0 ) return

    if ( cov < var_oi_bg ) then
       new_cov = (cov *obs_var)/(cov + obs_var)
      
       obsmmean = obs - mean
       if (abs(obsmmean) <= obs_diff_cutoff) then
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
       elseif (obsmmean > obs_diff_cutoff) then
          new_mean = new_cov * (obs_var*mean + cov*(mean+obs_diff_cutoff))/(cov*obs_var)
       else
          new_mean = new_cov * (obs_var*mean + cov*(mean-obs_diff_cutoff))/(cov*obs_var)
       end if

       a = sqrt(new_cov/cov)
       obs_inc_eakf = a * (ens - mean) + new_mean - ens

       if ( var_oi_bg > 1.e-8 ) then
          new_cov = (var_oi_bg * obs_var_oi)/(var_oi_bg + obs_var_oi)
      
          if (abs(obsmmean) <= obs_diff_cutoff) then
             new_mean = new_cov * (obs_var_oi*mean + var_oi_bg*obs)/(cov*obs_var)
          elseif (obsmmean > obs_diff_cutoff) then
             new_mean = new_cov * (obs_var_oi*mean + var_oi_bg*(mean+obs_diff_cutoff))/(var_oi_bg*obs_var_oi)
          else
             new_mean = new_cov * (obs_var_oi*mean + var_oi_bg*(mean-obs_diff_cutoff))/(var_oi_bg*obs_var_oi)
          end if

          a = sqrt(new_cov/var_oi_bg)
          obs_inc_oi = a * (ens - mean) + new_mean - ens
       else
          obs_inc_oi = 0.0
       end if
    else
       if ( abs(cov) > 1.e-4 ) then
          new_cov = (cov *obs_var)/(cov + obs_var)
      
          obsmmean = obs - mean
          if (abs(obsmmean) <= obs_diff_cutoff) then
             new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          elseif (obsmmean > obs_diff_cutoff) then
             new_mean = new_cov * (obs_var*mean + cov*(mean+obs_diff_cutoff))/(cov*obs_var)
          else
             new_mean = new_cov * (obs_var*mean + cov*(mean-obs_diff_cutoff))/(cov*obs_var)
          end if

          a = sqrt(new_cov/cov)

          obs_inc_eakf = a * (ens - mean) + new_mean - ens
          obs_inc_oi = 0.0
       else
          obs_inc_eakf = 0.0

          obs_inc_oi = 0.0
       end if
    end if


  end subroutine obs_increment_prf_eta_hyb

!lulv add subroutine obs_increment_prf_eta_hyb_stn below
  subroutine obs_increment_prf_eta_hyb_stn(ensb,ens, ens_size, obs, obs_var, obs_inc_eakf, obs_var_oi, obs_inc_oi, std_oi)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens, ensb
    real, intent(in) :: obs, obs_var, obs_var_oi
    real, intent(inout), dimension(ens_size) :: obs_inc_eakf, obs_inc_oi
    real, intent(in) :: std_oi

    real :: a, cov, obsmmean
    real :: mean, new_cov, new_mean, meanb
    real :: var_oi_bg

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    ! Default return values for obs_inc_eakf and obs_inc_oi
    obs_inc_eakf = 0.0
    obs_inc_oi = 0.0

    var_oi_bg = std_oi * std_oi

    ! Need to compute the prior covariance and mean; copy to 2d for now but
    ! get rid of this step at some point; just fits interface
    mean = sum(ens) / ens_size
    meanb = sum(ensb) / ens_size

    cov = 0.0
    do ie=1, ens_size
       cov = cov + (ensb(ie)-meanb)*(ensb(ie)-meanb)
    end do
    cov = cov / ens_size

    if ( cov == 0.0 ) return

    if ( cov < var_oi_bg ) then
       new_cov = (cov *obs_var)/(cov + obs_var)
      
       obsmmean = obs - mean
       if (abs(obsmmean) <= obs_diff_cutoff) then
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
       elseif (obsmmean > obs_diff_cutoff) then
          new_mean = new_cov * (obs_var*mean + cov*(mean+obs_diff_cutoff))/(cov*obs_var)
       else
          new_mean = new_cov * (obs_var*mean + cov*(mean-obs_diff_cutoff))/(cov*obs_var)
       end if

       a = sqrt(new_cov/cov)
       obs_inc_eakf = a * (ens - mean) + new_mean - ens

       if ( var_oi_bg > 1.e-8 ) then
          new_cov = (var_oi_bg * obs_var_oi)/(var_oi_bg + obs_var_oi)
      
          if (abs(obsmmean) <= obs_diff_cutoff) then
             new_mean = new_cov * (obs_var_oi*mean + var_oi_bg*obs)/(cov*obs_var)
          elseif (obsmmean > obs_diff_cutoff) then
             new_mean = new_cov * (obs_var_oi*mean + var_oi_bg*(mean+obs_diff_cutoff))/(var_oi_bg*obs_var_oi)
          else
             new_mean = new_cov * (obs_var_oi*mean + var_oi_bg*(mean-obs_diff_cutoff))/(var_oi_bg*obs_var_oi)
          end if

          a = sqrt(new_cov/var_oi_bg)
          obs_inc_oi = a * (ens - mean) + new_mean - ens
       else
          obs_inc_oi = 0.0
       end if
    else
       if ( abs(cov) > 1.e-4 ) then
          new_cov = (cov *obs_var)/(cov + obs_var)
      
          obsmmean = obs - mean
          if (abs(obsmmean) <= obs_diff_cutoff) then
             new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          elseif (obsmmean > obs_diff_cutoff) then
             new_mean = new_cov * (obs_var*mean + cov*(mean+obs_diff_cutoff))/(cov*obs_var)
          else
             new_mean = new_cov * (obs_var*mean + cov*(mean-obs_diff_cutoff))/(cov*obs_var)
          end if

          a = sqrt(new_cov/cov)
          obs_inc_eakf = a * (ens - mean) + new_mean - ens

          obs_inc_oi = 0.0
       else
          obs_inc_eakf = 0.0

          obs_inc_oi = 0.0
       end if
    end if
    
  end subroutine obs_increment_prf_eta_hyb_stn
!lulv add subroutine obs_increment_prf_eta_hyb_stn above

  subroutine obs_increment(ens, ens_size, obs, obs_var, obs_inc_eakf, obs_inc_oi, std_oi)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens
    real, intent(in) :: obs, obs_var
    real, intent(inout), dimension(ens_size) :: obs_inc_eakf, obs_inc_oi
    real, intent(in) ::  std_oi

    real :: a, cov
    real :: mean, new_cov, new_mean
    real :: var_oi

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi = std_oi * std_oi

    ! Need to compute the prior covariance and mean; copy to 2d for now but
    ! get rid of this step at some point; just fits interface
    mean = sum(ens) / ens_size

    cov = 0.0

    do ie=1, ens_size
       cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
    end do

    cov = cov / ens_size

    if ( cov >= 1.e-8 ) then
       if ( cov < var_oi ) then
          new_cov = (cov *obs_var)/(cov + obs_var)
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          a = sqrt(new_cov/cov)
          obs_inc_eakf = a * (ens - mean) + new_mean - ens

          new_cov = (var_oi * obs_var)/(var_oi + obs_var)
          new_mean = new_cov * (obs_var*mean + var_oi*obs)/(var_oi*obs_var)
          a = sqrt(new_cov/var_oi)
          obs_inc_oi = a * (ens - mean) + new_mean - ens
       else
          new_cov = (cov *obs_var)/(cov + obs_var)
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          a = sqrt(new_cov/cov)
          obs_inc_eakf = a * (ens - mean) + new_mean - ens
          !if(mpp_pe()==122) print*,'in obs_increment,mean,cov,obs,obs_var=',mean,cov,obs,obs_var !lulv 20201010
          !if(mpp_pe()==122) print*,'in obs_increment,ens,ens-mean,obs_inc_eakf=',ens,ens-mean,obs_inc_eakf !lulv 20201010

!!$          do ie=1, ens_size
!!$             if ( obs_inc_eakf(ie) > 100.0 .or. obs_inc_eakf(ie) < -100.0 ) then
!!$                print*,'obs_inc_eakf=',obs_inc_eakf(ie),'new_mean=',new_mean,'new_cov=',new_cov,'cov=',cov
!!$                print*,'mean=',mean,'ens=',ens(ie)
!!$             end if
!!$          end do

          obs_inc_oi = 0.0
       end if
    else
       obs_inc_eakf = 0.0
       obs_inc_oi = 0.0
    end if
  end subroutine obs_increment

!lulv  subroutine obs_increment_hfq(ens, ens_size, obs, obs_var, obs_inc_eakf, obs_inc_oi, std_oi)
  subroutine obs_increment_hfq(config_flags,ensb, ens, ens_size, obs, obs_var, obs_inc_eakf, obs_inc_oi, std_oi,i_o,j_o,atm_step)
        TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens, ensb !lulv 20200825 add ensb
    real, intent(in) :: obs, obs_var
    real, intent(inout), dimension(ens_size) :: obs_inc_eakf, obs_inc_oi
    real, intent(in) ::  std_oi
    integer, intent(in) :: i_o,j_o,atm_step !wangj
    real,dimension(ens_size)  :: diffa_m
    INTEGER myproc

    real :: a, cov
    real :: mean, meanb,new_cov, new_mean,ens_sum,ensb_sum
    real :: var_oi

    integer :: ie

   CALL wrf_get_myproc( myproc )
!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi = std_oi * std_oi

    ! Need to compute the prior covariance and mean; copy to 2d for now but
    ! get rid of this step at some point; just fits interface

   ens_sum = 0.0
   ensb_sum = 0.0
   diffa_m = 0.0
   mean = 0.0
   meanb = 0.0
   new_cov = 0.0 
   new_mean = 0.0

    do ie=1, ens_size
     ens_sum = ens_sum+ens(ie)
     ensb_sum = ensb_sum+ensb(ie)
    end do  
     mean = ens_sum/ens_size
     meanb =ensb_sum/ens_size

!    mean = sum(ens) / ens_size
!    meanb = sum(ensb) / ens_size !lulv 20200825


                do ie =1,ens_size
                    diffa_m(ie) = ens(ie) - mean
                enddo

!IF(diffa_m(1)/=0) then
!      write(*,'("diffa_m1no0",1xi5,1xi3,1xi3,1xi3,20(1xf16.2),1xf16.2)')&
!                     atm_step,myproc,i_o,j_o,ens(:),mean
!endIF

    cov = 0.0

    do ie=1, ens_size
    !lulv   cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
       cov = cov + (ensb(ie)-meanb)*(ensb(ie)-meanb) !lulv 20200825
    end do

    cov = cov / ens_size


    if ( cov >= 1.e-8 ) then
       if ( cov < var_oi ) then
          new_cov = (cov *obs_var)/(cov + obs_var)
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          a = sqrt(new_cov/cov)
          obs_inc_eakf = a * (ens - mean) + new_mean - ens

          new_cov = (var_oi * obs_var)/(var_oi + obs_var)
          new_mean = new_cov * (obs_var*mean + var_oi*obs)/(var_oi*obs_var)
          a = sqrt(new_cov/var_oi)
          obs_inc_oi = a * (ens - mean) + new_mean - ens
       else
          new_cov = (cov *obs_var)/(cov + obs_var)
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          a = sqrt(new_cov/cov)

          do ie=1, ens_size    
!          obs_inc_eakf(ie) = a * (ens(ie) - mean) + new_mean - ens(ie) !wangj
           obs_inc_eakf(ie) =  new_mean - ens(ie)
          enddo

          obs_inc_oi = 0.0
       end if
    else
       obs_inc_eakf = 0.0
       obs_inc_oi = 0.0
    end if



  end subroutine obs_increment_hfq

  subroutine obs_increment_ps(ens, ens_size, obs, obs_var, obs_inc_eakf, obs_inc_oi, std_oi, cov_inflate)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens
    real, intent(in) :: obs, obs_var
    real, intent(inout), dimension(ens_size) :: obs_inc_eakf, obs_inc_oi
    real, intent(in) ::  std_oi, cov_inflate

    real :: a, cov
    real :: mean, new_cov, new_mean
    real :: var_oi

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi = std_oi * std_oi

    ! Need to compute the prior covariance and mean; copy to 2d for now but
    ! get rid of this step at some point; just fits interface
    mean = sum(ens) / ens_size

    cov = 0.0

    do ie=1, ens_size
       cov = cov + (ens(ie)-mean)*(ens(ie)-mean)
    end do

    cov = cov / ens_size

    cov = cov * cov_inflate

    if ( cov >= 1.e-8 ) then
       if ( cov < var_oi ) then
          new_cov = (cov *obs_var)/(cov + obs_var)
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          a = sqrt(new_cov/cov)
          obs_inc_eakf = a * (ens - mean) + new_mean - ens

          new_cov = (var_oi * obs_var)/(var_oi + obs_var)
          new_mean = new_cov * (obs_var*mean + var_oi*obs)/(var_oi*obs_var)
          a = sqrt(new_cov/var_oi)
          obs_inc_oi = a * (ens - mean) + new_mean - ens
       else
          new_cov = (cov *obs_var)/(cov + obs_var)
          new_mean = new_cov * (obs_var*mean + cov*obs)/(cov*obs_var)
          a = sqrt(new_cov/cov)
          obs_inc_eakf = a * (ens - mean) + new_mean - ens

!!$          do ie=1, ens_size
!!$             if ( obs_inc_eakf(ie) > 100.0 .or. obs_inc_eakf(ie) < -100.0 ) then
!!$                print*,'obs_inc_eakf=',obs_inc_eakf(ie),'new_mean=',new_mean,'new_cov=',new_cov,'cov=',cov
!!$                print*,'mean=',mean,'ens=',ens(ie)
!!$             end if
!!$          end do

          obs_inc_oi = 0.0
       end if
    else
       obs_inc_eakf = 0.0
       obs_inc_oi = 0.0
    end if
  end subroutine obs_increment_ps

  subroutine obs_increment_eta_oi(ens, ens_size, obs, obs_var, obs_inc_oi, std_bg)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: ens
    real, intent(in) :: obs, obs_var
    real, intent(inout), dimension(ens_size) :: obs_inc_oi
    real, intent(in) :: std_bg

    real :: a, mean, new_cov, new_mean
    real :: var_oi

!    if ( .not.module_initialized ) call assim_tools_init()

    ! Need to compute the prior covariance and mean; copy to 2d for now but
    ! get rid of this step at some point; just fits interface
    mean = sum(ens) / ens_size

    var_oi = std_bg*std_bg

    if ( std_bg > 1.0e-5 ) then
       new_cov = (var_oi * obs_var)/(var_oi + obs_var)
       new_mean = new_cov * (obs_var*mean + var_oi*obs)/(var_oi*obs_var)
       a = sqrt(new_cov/var_oi)
       obs_inc_oi = a * (ens - mean) + new_mean - ens
    else
       obs_inc_oi = 0.0
    end if
  end subroutine obs_increment_eta_oi

  subroutine update_from_obs_inc_prf_hyb(obs, obs_inc_eakf, obs_inc_oi, state, ens_size, new_state,&
       & cov_factor, cor_oi, std_oi_o, std_oi_g, ass_method, ass_variable, flag_hyb)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs, obs_inc_eakf, obs_inc_oi, state
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor, cor_oi, std_oi_o, std_oi_g
    integer, intent(in) :: ass_method, ass_variable, flag_hyb

    real, dimension(ens_size) :: new_state_oi, new_state_eakf
    real :: var_oi_o, var_oi_g, mean_s, mean_o, cv_s, cv_o, cv, cr, cv_oi
    real :: std_o, std_g, r_ens, r_oi, total_r, total_r_inv
    real :: cor_cut

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi_o = std_oi_o*std_oi_o
    var_oi_g = std_oi_g*std_oi_g

    ! Compute statistics up to the second moment
    !::sdu:: temporary setting of cor_cut until an actual value is obtained.
    cor_cut = cor_cutoff
    if ( ass_variable == 2 ) cor_cut = cor_cutoff*0.5

    mean_s = sum(state) / ens_size
    mean_o = sum(obs) / ens_size

    cv_s = 0.0
    cv_o = 0.0
    cv   = 0.0

    do ie=1, ens_size
       cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
       cv_o = cv_o + (obs(ie) - mean_o)*(obs(ie) - mean_o)
       cv   = cv + (state(ie) - mean_s)*(obs(ie) - mean_o)
    end do

    cv_s = cv_s / ens_size
    cv_o = cv_o / ens_size
    cv   = cv / ens_size
!print*,'cv_s,cv_o,cv=',cv_s,cv_o,cv !lulv 2020-04-08

    if ( cv_s /= 0.0 .and. cv_o /= 0.0 ) then
       std_o = sqrt(cv_o)
       std_g = sqrt(cv_s)
       cr   = cv /(std_g*std_o)
    else
       std_o = 0.0
       std_g = 0.0
       cr = 0.0
    end if

!if (flag_hyb == 3361*6*33) then
!  print*,'cv_s,cv_o,cv,cr=',cv_s,cv_o,cv,cr,'cov_factor=',cov_factor
!end if

    new_state_eakf = 0.0
    new_state_oi = 0.0

    if ( std_o > std_oi_o ) then
       if ( abs(cr) >= cor_cut .and. cv_o > 1.e-4 ) then
          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
       else
          new_state_eakf = 0.0
       end if
       new_state = new_state_eakf
    else
       if ( abs(cr) >= cor_cut .and. cv_o > 1.e-4 ) then
          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
       else
          new_state_eakf = 0.0
       end if

!!$       cv_oi   = std_oi_g * cor_oi * std_oi_o ! snz test!!!!!!!!
       cv_oi   = std_oi_g * cr * std_oi_o

!!$       if ( cor_oi >= cor_cut .and. var_oi_o > 1.e-8 ) then ! snz test!!!!!!!!
       if ( abs(cr) >= cor_cut .and. var_oi_o > 1.e-4 .and. flag_hyb == 0 ) then
          new_state_oi = (cov_factor * cv_oi / var_oi_o) * obs_inc_oi
       else
          new_state_oi = 0.0
       end if

       if ( std_o > 1.e-4 ) then
          r_ens = std_g/std_o
       else
          r_ens = 0.0
       end if

       if ( std_oi_o > 1.e-4 ) then
          r_oi = std_oi_g/std_oi_o
       else
          r_oi = 0.0
       end if

       total_r = r_ens + r_oi

       if ( total_r > 0.0 ) then
          total_r_inv = 1.0/total_r
       else
          total_r_inv = 0.0
       end if

       if ( sum(new_state_oi) /= 0.0 ) then
          new_state =  r_ens*total_r_inv * new_state_eakf + r_oi*total_r_inv * new_state_oi
       else
          new_state = new_state_eakf

!!$          if ( r_oi > r_ens ) then
!!$             new_state =  r_ens*total_r_inv * new_state_eakf + r_oi*total_r_inv * new_state_oi
!!$          end if
       end if
    end if

!   print*,'obs_inc_eakf,cov_factor,cv,cv_o,new_state_eakf=',&
!      & obs_inc_eakf(1),cov_factor,cv,cv_o,new_state_eakf(1) !lulv 2020-04-10
!!$    new_state = new_state_oi ! for test purpose only use snz_oi (sigma_0) (snz)
  end subroutine update_from_obs_inc_prf_hyb

  subroutine update_from_inc_bt_hyb(obs, obs_inc_eakf, state, ens_size, new_state, cov_factor, std_oi_o, std_oi_g)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs, obs_inc_eakf, state
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor, std_oi_o, std_oi_g

    real, dimension(ens_size) :: new_state_oi, new_state_eakf
    real :: mean_s, mean_o, cv_s, cv_o, cv, cr, cv_oi
    real :: std_o, std_g, r_ens, r_oi, total_r, total_r_inv
    real :: cor_cut, var_oi_o, var_oi_g

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi_o = std_oi_o*std_oi_o
    var_oi_g = std_oi_g*std_oi_g

    ! Compute statistics up to the second moment
    cor_cut = cor_cutoff*0.5

    mean_s = sum(state) / ens_size
    mean_o = sum(obs) / ens_size

    cv_s = 0.0
    cv_o = 0.0
    cv   = 0.0

    do ie=1, ens_size
       cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
       cv_o = cv_o + (obs(ie) - mean_o)*(obs(ie) - mean_o)
       cv   = cv + (state(ie) - mean_s)*(obs(ie) - mean_o)
    end do

    cv_s = cv_s / ens_size
    cv_o = cv_o / ens_size
    cv   = cv / ens_size
    if ( cv_s /= 0.0 .and. cv_o /= 0.0 ) then
       std_o = sqrt(cv_o)
       std_g = sqrt(cv_s)
       cr   = cv /(std_g*std_o)
    else
       std_o = 0.0
       std_g = 0.0
       cr = 0.0
    end if

    new_state_eakf = 0.0
    new_state_oi = 0.0

    if ( std_o > std_oi_o ) then
       if ( abs(cr) >= cor_cut .and. cv_o > 1.e-8) then
          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
       else
          new_state_eakf = 0.0
       end if

       new_state = new_state_eakf
    else
       if ( abs(cr) >= cor_cut .and. cv_o > 1.e-8 ) then
          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
       else
          new_state_eakf = 0.0
       end if

       cv_oi   = std_oi_g * cr * std_oi_o ! snz test!!!!!!!!

       if ( abs(cr) >= cor_cut .and. var_oi_o > 1.e-8 ) then ! snz test!!!!!!!!!!
          new_state_oi = (cov_factor * cv_oi / var_oi_o) * obs_inc_eakf
       else
          new_state_oi = 0.0
       end if

       if ( std_o > 1.e-4 ) then
          r_ens = std_g/std_o
       else
          r_ens = 0.0
       end if

       if ( std_oi_o > 1.e-4 ) then
          r_oi = std_oi_g/std_oi_o
       else
          r_oi = 0.0
       end if

       total_r = r_ens + r_oi

       if ( total_r > 0.0 ) then
          total_r_inv = 1.0/total_r
       else
          total_r_inv = 0.0
       end if

       if ( sum(new_state_oi) /= 0.0 ) then
          new_state =  r_ens*total_r_inv * new_state_eakf + r_oi*total_r_inv * new_state_oi
       else
          new_state = new_state_eakf
       end if
    end if
  end subroutine update_from_inc_bt_hyb

  subroutine update_from_obs_inc_eta_hyb(obs, obs_inc_eakf, obs_inc_oi, state, ens_size, new_state,&
       & cov_factor, cor_oi, std_oi_o, std_oi_g, ass_method, ass_variable)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs, obs_inc_eakf, obs_inc_oi, state
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor, cor_oi, std_oi_o, std_oi_g
    integer, intent(in) :: ass_method, ass_variable

    real, dimension(ens_size) :: new_state_oi, new_state_eakf
    real :: var_oi_o, var_oi_g, mean_s, mean_o, cv_s, cv_o, cv, cr, cv_oi
    real :: std_o, std_g, r_ens, r_oi, total_r, total_r_inv
    real :: cor_cut

    integer :: ie

!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi_o = std_oi_o*std_oi_o
    var_oi_g = std_oi_g*std_oi_g

    ! Compute statistics up to the second moment
    !::sdu:: temporary setting of cor_cut until an actual value is obtained.
    cor_cut = cor_cutoff
    if ( ass_variable == 2 ) cor_cut = cor_cutoff*0.01

    mean_s = sum(state) / ens_size
    mean_o = sum(obs) / ens_size

    cv_s = 0.0
    cv_o = 0.0
    cv   = 0.0

    do ie=1, ens_size
       cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
       cv_o = cv_o + (obs(ie) - mean_o)*(obs(ie) - mean_o)
       cv   = cv + (state(ie) - mean_s)*(obs(ie) - mean_o)
    end do

    cv_s = cv_s / ens_size
    cv_o = cv_o / ens_size
    cv   = cv / ens_size
    if ( cv_s /= 0.0 .and. cv_o /= 0.0 ) then
       std_g = sqrt(cv_s)
       std_o = sqrt(cv_o)
       cr   = cv /(std_g*std_o)
    else
       std_g = 0.0
       std_o = 0.0
       cr = 0.0
    end if

    if ( std_g > std_oi_g ) then
       if ( abs(cr) >= cor_cut .and. cv_o > 1.e-8 ) then
          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
       else
          new_state_eakf = 0.0
       end if
       new_state = new_state_eakf
    else
       if ( abs(cr) >= cor_cut .and. cv_o > 1.e-8 ) then
          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
       else
          new_state_eakf = 0.0
       end if

       cv_oi   = std_oi_g * cr * std_oi_o !!!!snz test!

       if ( abs(cr) >= cor_cut .and. var_oi_o > 1.e-8 ) then
          new_state_oi = (cov_factor * cv_oi / var_oi_o) * obs_inc_oi
       else
          new_state_oi = 0.0
       end if

       if ( std_o > 1.e-4 ) then
          r_ens = std_g/std_o
       else
          r_ens = 0.0
       end if

       if ( std_oi_o > 1.e-4 ) then
          r_oi = std_oi_g/std_oi_o
       else
          r_oi = 0.0
       end if

       total_r = r_ens + r_oi

       if ( total_r > 0.0 ) then
          total_r_inv = 1.0/total_r
       else
          total_r_inv = 0.0
       end if

       new_state =  r_ens*total_r_inv * new_state_eakf + r_oi*total_r_inv * new_state_oi
    end if
  end subroutine update_from_obs_inc_eta_hyb

  subroutine update_from_obs_inc(config_flags,obs, obs_inc_eakf, obs_inc_oi, state, ens_size, new_state, cov_factor,&
       & cor_oi, std_oi_o, std_oi_g, ass_method,ass_variable,i_o,j_o,i_m,j_m,atm_step,&
        & isc,iec,jsc,jec,obs_loc,vari)
      TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs, obs_inc_eakf, obs_inc_oi, state
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor, cor_oi, std_oi_o, std_oi_g
    integer, intent(in) :: ass_method,ass_variable,i_o,j_o,i_m,j_m,isc,iec,jsc,jec

    real, dimension(ens_size) :: new_state_oi, new_state_eakf
    real :: var_oi_o, var_oi_g, mean_s, mean_o, cv_s, cv_o, cv, cr
    real :: cor_cut
INTEGER myproc
integer, intent(in) :: atm_step
    integer :: ie
    type(loc_type),intent(in)   ::  obs_loc
    character(len=1),intent(in) :: vari

!    if ( .not.module_initialized ) call assim_tools_init()
 CALL wrf_get_myproc( myproc )

    var_oi_o = std_oi_o*std_oi_o
    var_oi_g = std_oi_g*std_oi_g

    ! Compute statistics up to the second moment

    select case ( ass_variable )
    case (1)
       cor_cut = cor_cutoff
    case (2)
       cor_cut = cor_cutoff*0.5
    case (3)
       cor_cut = cor_cutoff*4.0
    case (4)
       cor_cut = cor_cutoff*2.0
    case default
       !::sdu:: Possibly uninitialized.  Temporary value.
       cor_cut = cor_cutoff
    end select

    mean_s = sum(state) / ens_size
    mean_o = sum(obs) / ens_size

        if (config_flags%debug_mda) then 
          write(*,*) 'atm_step,obs,i_o,j_o,obs_loc%lon,obs_loc%lat,i_m,j_m,mean_s,mean_o',&
                vari,atm_step,obs,i_o,j_o,obs_loc%lon,obs_loc%lat,i_m,j_m,mean_s,mean_o
     endif

    cv_s = 0.0
    cv_o = 0.0
    cv   = 0.0

    do ie=1, ens_size
       cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
       cv_o = cv_o + (obs(ie) - mean_o)*(obs(ie) - mean_o)
       cv   = cv + (state(ie) - mean_s)*(obs(ie) - mean_o)
    end do

    cv_s = cv_s / ens_size
    cv_o = cv_o / ens_size
    cv   = cv / ens_size


    if ( cv_s /= 0.0 .and. cv_o /= 0.0 ) then
       cr   = cv /(sqrt(cv_s)*sqrt(cv_o))
    else
       cr = 0.0
    end if

    if ( abs(cr) >= cor_cut .and. cv_o > 1.e-8 ) then
       if ( cv_o < var_oi_o ) then
          if ( std_oi_o == 0.0 ) then
             new_state_oi = 0.0
          else
!!$             cv = cor_oi * std_oi_g*std_oi_o
             cv = cr * std_oi_g*std_oi_o
             new_state_oi = (cov_factor * cv / var_oi_o) * obs_inc_oi
          end if

          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
          new_state =  cv_o/(cv_o+var_oi_o) * new_state_eakf + var_oi_o/(cv_o+var_oi_o) * new_state_oi
       else
          new_state = (cov_factor * cv / cv_o) * obs_inc_eakf

        if (config_flags%debug_mda) then 
           write(*,'("cov_factor",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2)')&
              atm_step,myproc,i_o,j_o,obs_loc%lon,obs_loc%lat,i_m,j_m,cov_factor

           write(*,'("cv,cv_o",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2,1xf16.2)')&
              atm_step,myproc,i_o,j_o,obs_loc%lon,obs_loc%lat,i_m,j_m,cv,cv_o

           write(*,'("obs,state20,obs_inc_eakfu,new_state",1xi5,1xi3,1xi3,1xi3,1xf6.2,1xf6.2,1xi3,1xi3,1xf16.2,1xf16.2,1xf16.2,1xf16.2)')&
              atm_step,myproc,i_o,j_o,obs_loc%lon,obs_loc%lat,i_m,j_m,obs(20),state(20),obs_inc_eakf(20),new_state(20)
  endif
       end if
    else
       new_state = 0.0
    end if
  end subroutine update_from_obs_inc

  subroutine update_from_obs_inc_eta_oi(obs_inc_oi, ens_size, new_state, cov_factor, cor_oi, std_oi_o, std_oi_g)
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs_inc_oi
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor, cor_oi, std_oi_g, std_oi_o

    real :: var_oi_g, var_oi_o, cv
    real :: cor_cut

!    if ( .not.module_initialized ) call assim_tools_init()

    var_oi_o = std_oi_o*std_oi_o
    var_oi_g = std_oi_g*std_oi_g

    !::sdu:: cor_cut is uninitialized!!!!!
    if ( cor_oi > cor_cut*1.0e-2 ) then
       cv = cor_oi * std_oi_g*std_oi_o
       new_state = (cov_factor * cv / var_oi_o) * obs_inc_oi
    else
       new_state = 0.0
    end if
  end subroutine update_from_obs_inc_eta_oi

  subroutine update_from_obs_inc_k(config_flags,obs, obs_inc_eakf, obs_inc_oi,state, ens_size, new_state, cov_factor,&
       & cor_oi, std_oi_o, std_oi_g,ass_method,ass_variable,i,j,n,k,atm_step,&
        & isc,iec,jsc,jec,obs_loc)
      TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    integer, intent(in) :: ens_size
    real, intent(in), dimension(ens_size) :: obs, obs_inc_eakf, obs_inc_oi,state
    real, intent(inout), dimension(ens_size) :: new_state
    real, intent(in) :: cov_factor, cor_oi, std_oi_o, std_oi_g
    integer, intent(in) :: ass_method,ass_variable,i,j,n,k,isc,iec,jsc,jec
    real, dimension(ens_size) :: new_state_oi, new_state_eakf
    real :: var_oi_o, var_oi_g, mean_s, mean_o, cv_s, cv_o, cv, cr
    real :: cor_cut
INTEGER myproc
integer, intent(in) :: atm_step
    integer :: ie
    type(loc_type),intent(in)   ::  obs_loc


    var_oi_o = std_oi_o*std_oi_o
    var_oi_g = std_oi_g*std_oi_g

    ! Compute statistics up to the second moment

    select case ( ass_variable )
    case (1)
       cor_cut = cor_cutoff
    case (2)
       cor_cut = cor_cutoff*0.5
    case (3)
       cor_cut = cor_cutoff*4.0
    case (4)
       cor_cut = cor_cutoff*2.0
    case default
       !::sdu:: Possibly uninitialized.  Temporary value.
       cor_cut = cor_cutoff
    end select

    mean_s = sum(state) / ens_size
    mean_o = sum(obs) / ens_size

    cv_s = 0.0
    cv_o = 0.0
    cv   = 0.0


    do ie=1, ens_size
       cv_s = cv_s + (state(ie) - mean_s)*(state(ie) - mean_s)
       cv_o = cv_o + (obs(ie) - mean_o)*(obs(ie) - mean_o)
       cv   = cv + (state(ie) - mean_s)*(obs(ie) - mean_o)
    end do

    cv_s = cv_s / ens_size
    cv_o = cv_o / ens_size
    cv   = cv / ens_size


    if ( cv_s /= 0.0 .and. cv_o /= 0.0 ) then
       cr   = cv /(sqrt(cv_s)*sqrt(cv_o))
    else
       cr = 0.0
    end if

    if ( abs(cr) >= cor_cut .and. cv_o > 1.e-8 ) then
       if ( cv_o < var_oi_o ) then
          if ( std_oi_o == 0.0 ) then
             new_state_oi = 0.0
          else
!!$             cv = cor_oi * std_oi_g*std_oi_o
             cv = cr * std_oi_g*std_oi_o
             new_state_oi = (cov_factor * cv / var_oi_o) * obs_inc_oi
          end if

          new_state_eakf = (cov_factor * cv / cv_o) * obs_inc_eakf
          new_state =  cv_o/(cv_o+var_oi_o) * new_state_eakf + var_oi_o/(cv_o+var_oi_o) * new_state_oi
       else
          new_state = (cov_factor * cv / cv_o) * obs_inc_eakf
       end if
    else


       new_state = 0.0
    end if




 end subroutine update_from_obs_inc_k

end module module_assim_tools_mod
