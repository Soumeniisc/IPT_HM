!----------------------------------------------------------------------------
!  date:   May 23, 2017
!
!----------------------------------------------------------------------------
module IPT_MOD
  
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use FFT_MOD
  use imaginary_time_frequency_FFTs ! this is the module form arijit
  implicit none
  
  private
  public :: eq_ipt                     
            
contains
  
  subroutine eq_ipt(parm_,dmft_)
    type(parm), intent(in)       :: parm_
    type(dmft), intent(inout)    :: dmft_
    integer                      :: i, j
    double precision             :: w, tau
    
    do i=1, parm_%N_tau
       dmft_%self_energy%matsubara_t(i)= parm_%U_i**2*dmft_%weiss_green%matsubara_t(i)*dmft_%weiss_green%matsubara_t(parm_%N_tau-i+1)*dmft_%weiss_green%matsubara_t(i)
    end do
    
    ! Fourier transformation
    ! self_energy(t) -> self_energy(w)
    call fft_t2w(parm_,dmft_%self_energy%matsubara_t,dmft_%self_energy%matsubara_w)
    !call tau2wn_mid_point_taus(dmft_%self_energy%matsubara_t,dmft_%self_energy%matsubara_w,parm_%beta) ! by arijit
    
    ! solve the impurity Dyson equation
    ! G(w)=[G0(w)^{-1}-self_energy(w)]^{-1}
    do j=1, parm_%N_tau
       dmft_%local_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)/(1.d0-dmft_%weiss_green%matsubara_w(j)*dmft_%self_energy%matsubara_w(j))
    end do
    
  end subroutine eq_ipt
  
  

end module IPT_MOD
