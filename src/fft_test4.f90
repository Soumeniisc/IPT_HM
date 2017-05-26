! here we are testing inverse fourier transform 
program fft_test
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use FFT_MOD
  use imaginary_time_frequency_FFTs

  implicit none
  type(parm) :: parm_
  type(dmft) :: dmft_
  integer    :: j
  real       :: w,e,nf, tau
  double complex,allocatable,dimension(:)::ftaus

  call set_parm(parm_)
  call construct_dmft(parm_,dmft_)

  e = -0.1
  nf = 1.d0/(1.d0+exp(dble(e*parm_%beta)))

  open (unit=7, file='G_wn', action='write') 
  do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j)=1.d0/(xj*w - e)
       write (7,'(f15.10,f15.10,f15.10)') w , Real(dmft_%local_green%matsubara_w(j)), Imag(dmft_%local_green%matsubara_w(j))
       dmft_%local_green%matsubara_w(j)=1.d0/(xj*w - e) - 1.d0/(xj*w)
  end do
  close (7)

  !call fft_w2t(parm_,dmft_%local_green%matsubara_w,dmft_%local_green%matsubara_t)

  open (unit=7, file='G_tau_input', action='write')
  do j=1, parm_%N_tau
  	 
	tau = (j - 0.5)*parm_%dtau
        dmft_%local_green%matsubara_t(j) = -exp(e*(parm_%beta-tau)) * nf
	write (7,'(f15.10,f15.10)') tau , dmft_%local_green%matsubara_t(j)
  end do
  close(7)
  
  allocate(ftaus(parm_%N_tau))
  !ftaus = dcmplx(dmft_%local_green%matsubara_t+0.5,0.0) 
  dmft_%local_green%matsubara_t = dmft_%local_green%matsubara_t + 0.5
  !print*, "doing inverse fourier transform"
  call fft_t2w(parm_,dmft_%local_green%matsubara_t,dmft_%local_green%matsubara_w)
  !call tau2wn_mid_point_taus(dmft_%local_green%matsubara_t,dmft_%local_green%matsubara_w,parm_%beta)
  !deallocate(ftaus)
  print*, "inverse fourier transform is done"

  open (unit=7, file='G_wn_out', action='write')
  do j=1, parm_%N_tau
       !print*, j, dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j) = dmft_%local_green%matsubara_w(j) + 1.d0/dcmplx(0.d0,w)
       write (7,*) w , Real(dmft_%local_green%matsubara_w(j)), "    ", Imag(dmft_%local_green%matsubara_w(j)) 
  end do
  close (7)

end program fft_test