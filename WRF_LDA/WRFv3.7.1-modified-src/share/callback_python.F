module python_interface

!<description>
!This file is intended to make the interface/connection
!with the python
!
!Xianrui Zhu 2025/6/2
!

public python_da


contains

subroutine python_da(current_timestr, x_size, y_size, z_size, u, v, t, q, xlat ,xlon)
    implicit none
    character (len=*), intent(in)           :: current_timestr
    integer, intent(in)                     :: x_size, y_size, z_size
    real, dimension(x_size, y_size, z_size), intent(inout)   :: u
    real, dimension(x_size, y_size, z_size), intent(inout)   :: v
    real, dimension(x_size, y_size, z_size), intent(inout)   :: t
    real, dimension(x_size, z_size), intent(inout)   :: q
    real, dimension(x_size, z_size), intent(inout)   :: xlat
    real, dimension(x_size, z_size), intent(inout)   :: xlon
    !f2py intent(callback, hide) python_foo
    external python_foo

    call python_foo(current_timestr, x_size, y_size, z_size, u, v, t, q, xlat, xlon)

end subroutine 

end module python_interface