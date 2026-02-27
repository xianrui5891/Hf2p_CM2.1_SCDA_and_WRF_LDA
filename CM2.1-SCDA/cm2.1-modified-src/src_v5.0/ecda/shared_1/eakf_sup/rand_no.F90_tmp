module rand_no_mod
  ! $Id$
  ! Platform independent random number generator from
  ! Numerical Recipies
  ! Mark Webb July 1999

  implicit none

  public gau0

contains

  real function ran0(idum)
    integer, intent (INOUT) :: idum

    integer :: IA,IM,IQ,IR,k
    real :: AM

    IA=16807
    IM=2147483647
    AM=1.0/IM
    IQ=127773
    IR=2836

    if ( idum.eq.0 ) then
       print*,'ran0', 'idum=0, ZERO seed not allowed.'
    end if

    k = idum/IQ
    idum = ia*(idum-k*iq)-ir*k
    if ( idum.lt.0 ) idum = idum+im
    ran0 = am*idum
  end function ran0

  real function gau0()

    real :: a, b, w0
    integer :: idum1=100
    integer :: idum2=200

    w0 = 1.0
    do while ( w0 >= 1.0 )
       a = 2.0*ran0(idum1) - 1.0
       b = 2.0*ran0(idum1) - 1.0
       w0 = a*a + b*b
    end do
    w0 = sqrt((-2.0*alog(w0))/w0)
    gau0 = a*w0
  end function gau0
end module rand_no_mod
