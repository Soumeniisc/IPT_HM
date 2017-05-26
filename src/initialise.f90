program init_test
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use FFT_MOD
  use imaginary_time_frequency_FFTs

  implicit none
  type(parm) :: parm_
  type(dmft) :: dmft_


    double precision                :: ek, e_max, de, w, tau
    integer                         :: i, j, k,N_e
    call set_parm(parm_)
    call construct_dmft(parm_,dmft_)
    ! initialize G0(w)
    open (unit=13, file='Gf_wn', status='unknown', action='write')
    
    if (parm_%dos=='semicircular') then
       e_max=2.d0
       N_e = 1000
       de=e_max/dble(N_e)
       dmft_%weiss_green%matsubara_w(j)=0.d0
       do j=1, parm_%N_tau
          w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
          do k=2, 2*N_e
             ek=dble(k-N_e-1)*de
             dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+de*sqrt(4.d0-ek**2)/(2.d0*pi)/(xj*w-ek)
          end do
          write(13,*) w, real(dmft_%weiss_green%matsubara_w(j)), Imag(dmft_%weiss_green%matsubara_w(j))
       end do
    end if
    close(13)
    
    ! Fourier transformation: G0(w) -> G0(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)-parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    !call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t)
    call wn2tau_mid_point_taus(dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t, parm_%beta) ! by arijit
    open (unit=13, file='Gf_tau', status='unknown', action='write')
    do i=1, parm_%N_tau
       tau=dble(i-0.5)*parm_%dtau
       dmft_%weiss_green%matsubara_t(i)=dmft_%weiss_green%matsubara_t(i) - 0.5d0 + 0.25d0*parm_%e2*tau*(parm_%beta-tau)
       write(13,*) tau, real(dmft_%weiss_green%matsubara_t(i)), Imag(dmft_%weiss_green%matsubara_t(i))
       !write(13,*) tau, dmft_%weiss_green%matsubara_t(j)
    end do
   close(13)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)+parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    print*, "it's working fine"
end program init_test
