module fftw_module
!        include "fftw3.f90"
        include "fftw3.f"
        integer ( kind = 8 ) plan_backward
        integer ( kind = 8 ) plan_forward
        !private::plan_backward,plan_forward
contains
subroutine FFT_FORWARD(inarr,outarr)
        double complex,dimension(:)::inarr,outarr
        integer::N
        N=size(inarr)
        call dfftw_plan_dft_1d (plan_forward,N,inarr,outarr,FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan_forward)
        call dfftw_destroy_plan (plan_forward)
end subroutine

subroutine FFT_BACKWARD(inarr,outarr)
        double complex,dimension(:)::inarr,outarr
        integer::N
        N=size(inarr)
        call dfftw_plan_dft_1d (plan_backward,N,inarr,outarr,FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute (plan_backward)
        call dfftw_destroy_plan (plan_backward)
end subroutine
end module

