module module_loc_and_dist_mod

  ! Two d spherical location and distance module. Inputs are latitude and
  ! longitude in degrees???

!  use fms_mod, only : error_mesg, FATAL
!  use constants_mod, only : DEG_TO_RAD

  private
  public loc_type, get_dist, set_loc, get_loc

  type loc_type
     real :: lon, lat
  end type loc_type

contains
  ! computes pseudo-distance between two locations; just relative distance
  ! matters so can compute Euclidean distance in 3d rather than on surface
  ! of sphere
  real function get_dist(a, b,atm_step)
    type (loc_type), intent(in) :: a, b
        integer, intent(in)  :: atm_step
    real :: xa, xb, ya, yb, za, zb
     real,parameter :: PI      = 3.14159265358979323846
     real,parameter :: DEG_TO_RAD=PI/180.

    xa = cos(a%lon * DEG_TO_RAD) * cos(a%lat * DEG_TO_RAD)
    xb = cos(b%lon * DEG_TO_RAD) * cos(b%lat * DEG_TO_RAD)
    ya = sin(a%lon * DEG_TO_RAD) * cos(a%lat * DEG_TO_RAD)
    yb = sin(b%lon * DEG_TO_RAD) * cos(b%lat * DEG_TO_RAD)
    za = sin(a%lat * DEG_TO_RAD)
    zb = sin(b%lat * DEG_TO_RAD)

!write(*,*) 'testdist',atm_step,DEG_TO_RAD,a%lon * DEG_TO_RAD,cos(a%lon * DEG_TO_RAD)
!write(*,*) 'xa,xb,ya,yb,za,zb',atm_step,a%lon,a%lat,b%lon,b%lat,xa,xb,ya,yb,za,zb

    get_dist = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) + (za-zb)*(za-zb)

!write(*,*) 'get_dist',get_dist

  end function get_dist

end module module_loc_and_dist_mod
