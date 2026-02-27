module cov_cutoff_mod
  ! Computes a covariance cutoff function from Gaspari and Cohn (their eqn. 4.10)
  ! QJRMS, 125, 723-757.
  !
  ! z_in is the distance while c is the cutoff distance. For distances greater
  ! than 2c, the cov_factor returned goes to 0.
  implicit none

contains

  real function comp_cov_factor(z_in, c)
    real, intent(in) :: z_in, c

    real :: z, r
    real :: aa, bb

    aa = 0.083333333333333
    bb = 1.666666666666667

    z = z_in
    r = z / c

    if ( z >= 2*c ) then
       comp_cov_factor = 0.0
    else if ( z >= c .and. z < 2*c ) then
       comp_cov_factor = aa*r*r*r*r*r - 0.5*r*r*r*r + 0.625*r*r*r + bb*r*r - 5.*r + 4. - (2.*c)/(3.*z)
    else
       comp_cov_factor = -0.25*r*r*r*r*r + 0.5*r*r*r*r + 0.625*r*r*r - bb*r*r + 1.
    end if
  end function comp_cov_factor
end module cov_cutoff_mod
