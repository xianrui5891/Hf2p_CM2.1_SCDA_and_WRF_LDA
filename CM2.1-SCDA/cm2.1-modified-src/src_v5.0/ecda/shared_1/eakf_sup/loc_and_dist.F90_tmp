module loc_and_dist_mod

  ! Two d spherical location and distance module. Inputs are latitude and
  ! longitude in degrees???

  use fms_mod, only : error_mesg, FATAL
  use constants_mod, only : DEG_TO_RAD

  private
  public loc_type, get_dist, set_loc, get_loc

  type loc_type
     real :: lon, lat
  end type loc_type

contains

! Given a location type and a double precision value between 0 and 1
! puts this value into the location.
  subroutine set_loc(loc, lon, lat)
    type (loc_type), intent(out) :: loc
    real, intent(in) :: lon, lat

    if ( lon < 0.0 .or. lon > 360.0 .or. lat < -90.0 .or. lat > 90.0 ) then
       call error_mesg('loc_and_dist_mod::set_loc', 'Value of lon or lat is out range', FATAL)
    end if
    loc%lon = lon
    loc%lat = lat
  end subroutine set_loc

  ! Given a location type, return lon and lat
  subroutine get_loc(loc, lon, lat)
    type(loc_type), intent(in) :: loc
    real, intent(out) :: lon, lat

    lon = loc%lon
    lat = loc%lat
  end subroutine get_loc

  ! computes pseudo-distance between two locations; just relative distance
  ! matters so can compute Euclidean distance in 3d rather than on surface
  ! of sphere
  real function get_dist(a, b)
    type (loc_type), intent(in) :: a, b

    real :: xa, xb, ya, yb, za, zb

    xa = cos(a%lon * DEG_TO_RAD) * cos(a%lat * DEG_TO_RAD)
    xb = cos(b%lon * DEG_TO_RAD) * cos(b%lat * DEG_TO_RAD)
    ya = sin(a%lon * DEG_TO_RAD) * cos(a%lat * DEG_TO_RAD)
    yb = sin(b%lon * DEG_TO_RAD) * cos(b%lat * DEG_TO_RAD)
    za = sin(a%lat * DEG_TO_RAD)
    zb = sin(b%lat * DEG_TO_RAD)

    get_dist = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) + (za-zb)*(za-zb)
  end function get_dist
end module loc_and_dist_mod
