program fft_test
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use FFT_MOD

  implicit none
  type(parm) :: parm_
  type(dmft) :: dmft_
  integer    :: j
  real       :: w

  call set_parm(parm_)
  call construct_dmft(parm_,dmft_)

  open (unit=7, file='G_wn', action='write') 
  do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j)=1.d0/((xj*w - 0.1)*(xj*w - 0.1))
       write (7,'(f15.10,f15.10,f15.10)') w , Real(dmft_%local_green%matsubara_w(j)), Imag(dmft_%local_green%matsubara_w(j))
       dmft_%local_green%matsubara_w(j) = dmft_%local_green%matsubara_w(j) - 1.d0/(xj*w)
  end do
  close (7)

  call fft_w2t(parm_,dmft_%local_green%matsubara_w,dmft_%local_green%matsubara_t)

  open (unit=7, file='G_tau_free', action='write')
  do j=1, parm_%N_tau+1
  	write (7,'(f15.10,f15.10)') (j-1)*parm_%dtau , dmft_%local_green%matsubara_t(j) - 0.5
        dmft_%local_green%matsubara_t(j) = dmft_%local_green%matsubara_t(j) - 0.5
  end do
  close (7)
  
  call fft_t2w(parm_,dmft_%local_green%matsubara_t,dmft_%local_green%matsubara_w)

  open (unit=7, file='G_wn_out', action='write')
  do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       write (7,'(f15.10,f15.10,f15.10)') w , Real(dmft_%local_green%matsubara_w(j)), "    ", Imag(dmft_%local_green%matsubara_w(j))
  end do
  close (7)
end program fft_test