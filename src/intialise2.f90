subroutine initialize_eq_green(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: ek, e_max, de, w, tau
    integer                         :: i, j, k
    
    ! initialize G0(w)
    if (parm_%dos=='semicircular') then
       e_max=2.d0
       de=e_max/dble(parm_%N_e)
       dmft_%weiss_green%matsubara_w(j)=0.d0
       do j=1, parm_%N_tau
          w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
          do k=2, 2*parm_%N_e
             ek=dble(k-parm_%N_e-1)*de
             dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+de*sqrt(4.d0-ek**2)/(2.d0*pi)/(xj*w-ek)
          end do
       end do
    end if
    
    ! Fourier transformation: G0(w) -> G0(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)-parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t)
    do i=1, parm_%N_tau+1
       tau=dble(i-1)*parm_%dtau
       dmft_%weiss_green%matsubara_t(i)=dmft_%weiss_green%matsubara_t(i)-0.5d0+0.25d0*parm_%e2*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)+parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    
  end subroutine initialize_eq_green

