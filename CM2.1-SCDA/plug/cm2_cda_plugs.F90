module cm2_cda_plugs

    use cm2_cda_mainsubs, only: tool_cm2_cda_maininit,&
        tool_atmos_step_pre, tool_atmos_step,&
        tool_ocean_step_pre, tool_ocean_step,&
        tool_restart_step, tool_main_end


    public :: cm2_cda_maininit,&
            atmos_step_pre,&
            atmos_step,&
            ocean_step_pre,&
            ocean_step,&
            restart_step,&
            main_end

contains

    subroutine cm2_cda_maininit(num_cpld_calls_o, num_atmos_calls_o)
        integer,intent(out):: num_cpld_calls_o, num_atmos_calls_o
        integer :: num_cpld_calls_t, num_atmos_calls_t
!        integer :: nlna = 144, nlta = 90 ! for global atmosphere bottom dimensions
!        real(kind=8), intent(inout) :: U_bot(144, 90), V_bot(144, 90), T_bot(144, 90)
!        integer :: nlno = 360, nlto = 200 ! for global ocean top dimensions
!        real(kind=8), intent(inout) :: SSU(360, 200), SSV(360, 200), SST(360, 200)

        call tool_cm2_cda_maininit(num_cpld_calls_t, num_atmos_calls_t)

        num_cpld_calls_o = num_cpld_calls_t
        num_atmos_calls_o = num_atmos_calls_t

    end subroutine cm2_cda_maininit

    subroutine atmos_step_pre(nc)
!        integer :: nlno = 360, nlto = 200 ! for global ocean top dimensions
!        real(kind=8), intent(inout) :: SSU(360, 200), SSV(360, 200), SST(360, 200)
        integer, intent(in) :: nc

        call tool_atmos_step_pre(nc)

    end subroutine atmos_step_pre

    subroutine atmos_step(U_bot, V_bot, T_bot, P_bot, nc, na, nu_scda)
!        integer :: nlna = 144, nlta = 90 ! for global atmosphere bottom dimensions
        real(kind=8), intent(inout) :: U_bot(144, 90), V_bot(144, 90), T_bot(144, 90), P_bot(144, 90)
        integer, intent(in) :: nc, na, nu_scda

        print*,'in cm2_plug nc, na=', nc, na
        call tool_atmos_step(U_bot, V_bot, T_bot, P_bot, nc, na, nu_scda)

    end subroutine atmos_step

    subroutine ocean_step_pre()

        call tool_ocean_step_pre()

    end subroutine ocean_step_pre

    subroutine ocean_step(SSU, SSV, SST, SSZ, nc, na, nu_scda)
!        integer :: nlno = 360, nlto = 200 ! for global ocean top dimensions
        real(kind=8), intent(inout) :: SSU(360, 200), SSV(360, 200), SST(360, 200), SSZ(360, 200)
        integer, intent(in) :: nc, na, nu_scda

        call tool_ocean_step(SSU, SSV, SST, SSZ, nc, na, nu_scda)

    end subroutine ocean_step

    subroutine restart_step(nc)
        integer, intent(in) :: nc

        call tool_restart_step(nc)

    end subroutine restart_step

    subroutine main_end()

        call tool_main_end()

    end subroutine main_end

end module cm2_cda_plugs
