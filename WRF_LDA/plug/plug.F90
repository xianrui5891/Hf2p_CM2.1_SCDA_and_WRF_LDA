module plug
use main

public call_wrf_main

contains
    subroutine call_wrf_main
        !f2py intent(callback, hide) python_foo
        external python_foo
        call wrf
    end subroutine call_wrf_main

end module plug