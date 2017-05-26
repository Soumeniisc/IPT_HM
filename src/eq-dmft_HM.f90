!----------------------------------------------------------------------------
!  data:  May 23, 2017
!  update: Its modified by soumen on 05/01/17
!
!  update: read green function is added, print green function is added also added when to read the green function .init veriable
!  update: high frequency ecponents are calculated here
! 
!
!----------------------------------------------------------------------------
module EQ_DMFT_MOD
  
  use CONST_MOD
  use PARM_MOD
  use DMFT_MOD
  use GREEN_MOD
  use IPT_MOD
  use FFT_MOD
  use INTEGRAL_MOD
  !use imaginary_time_frequency_FFTs
  implicit none
  
  private
  public :: start_eq_dmft
  
contains
  
  subroutine start_eq_dmft(parm_,dmft_)
    type(parm), intent(in)        :: parm_
    type(dmft), intent(inout)     :: dmft_
    integer                       :: i, j
    double precision              :: w, tau
    
    dmft_%converged=.false.
    dmft_%G0_diff=1.d0
    dmft_%time_step=1
    dmft_%iteration=1

    if(parm_%init.eq.0) then
	print*, "non interacting grren finction is being set as  host green function"
    	!call initialize_eq_green(parm_,dmft_)
    else 
	print*, "reading eq green function"
	call read_eq_green(parm_,dmft_)
    end if
	
    ! equilibrium DMFT self-consistency loop
    do while (.not. dmft_%converged .and. dmft_%iteration<=parm_%N_iter)
       call eq_ipt(parm_,dmft_)
       call eq_dmft_self_consistency(parm_,dmft_)
       dmft_%iteration=dmft_%iteration+1
       if (dmft_%G0_diff<=parm_%tolerance) then
          dmft_%converged=.true.
          write (*,'(A)') 'Equilibrium DMFT is converged'
       end if
       if (dmft_%iteration==parm_%N_iter+1 .and. dmft_%G0_diff>parm_%tolerance) then
          write (*,'(A)') 'Equilibrium DMFT is NOT converged !!!!!'
       end if
    end do
    
    ! Fourier transformation: G(w) -> G(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j)=dmft_%local_green%matsubara_w(j)-1.d0/(xj*w)!-(parm_%e2+0.25d0*parm_%U_i**2)/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%local_green%matsubara_w,dmft_%local_green%matsubara_t)
    !call wn2tau_mid_point_taus(dmft_%local_green%matsubara_w,dmft_%local_green%matsubara_t, parm_%beta) ! by arijit
    do i=1, parm_%N_tau
       tau=dble(i-0.5)*parm_%dtau
       dmft_%local_green%matsubara_t(i)=dmft_%local_green%matsubara_t(i)-0.5d0!+0.25d0*(parm_%e2+0.25d0*parm_%U_i**2)*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%local_green%matsubara_w(j)=dmft_%local_green%matsubara_w(j)+1.d0/(xj*w)!+(parm_%e2+0.25d0*parm_%U_i**2)/(xj*w)/(xj*w)/(xj*w)
    end do
    
    call measure_density(parm_,dmft_)
    call measure_double_occupancy(parm_,dmft_)
    call measure_kinetic_energy(parm_,dmft_)
    call print_green_function(parm_,dmft_)
  end subroutine start_eq_dmft
  
  
  subroutine eq_dmft_self_consistency(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i, j, k
    double precision                :: w, tau, ah(2)

    ah(:) = 0.d0
    if (parm_%dos=='semicircular') then
       ! solve the lattice Dyson equation: G0(w)=1/[i*w-G(w)]
       do j=1, parm_%N_tau
          w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
          dmft_%weiss_green%matsubara_w(j)=1.d0/(xj*w-dmft_%local_green%matsubara_w(j))
       end do
    end if
    
    ! Fourier transformation: G0(w) -> G0(t)
    !call calculate_coeffi(parm_,100,dmft_%weiss_green,ah)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)!-(parm_%e2+0.25d0*parm_%U_i**2)/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green_new%matsubara_t)
    !call wn2tau_mid_point_taus(dmft_%weiss_green%matsubara_w,dmft_%weiss_green_new%matsubara_t, parm_%beta) ! by arijit
    do i=1, parm_%N_tau
       tau=dble(i-0.5)*parm_%dtau
       dmft_%weiss_green_new%matsubara_t(i)=dmft_%weiss_green_new%matsubara_t(i)-0.5d0!+ 0.25d0*(parm_%e2+0.25d0*parm_%U_i**2)*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)!+(parm_%e2+0.25d0*parm_%U_i**2)/(xj*w)/(xj*w)/(xj*w)
    end do
    
    ! evaluate |G0_{new}-G0_{old}|
    dmft_%G0_diff=0.d0
    do i=1, parm_%N_tau
       dmft_%G0_diff=dmft_%G0_diff+abs(dmft_%weiss_green_new%matsubara_t(i)-dmft_%weiss_green%matsubara_t(i))
    end do
    write (*,'(A,I4)') '  Iteration #', dmft_%iteration !, ah(1), ah(2)
    write (*,'(A,f16.14)') ' |G0_new-G0_old| =', dmft_%G0_diff
    ! G0_{old} <= G0_{new}
    dmft_%weiss_green%matsubara_t(:)= (1-parm_%mix)*dmft_%weiss_green_new%matsubara_t(:)+ parm_%mix*dmft_%weiss_green%matsubara_t(:)

    open (unit=7, file='density_', action='write', position='append')
    write(7,*) dmft_%iteration, -dmft_%weiss_green%matsubara_t(parm_%N_tau)
    close(7)
  end subroutine eq_dmft_self_consistency
  
  
  subroutine initialize_eq_green(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: ek, e_max, de, w, tau
    integer                         :: i, j, k,N_e
    
    ! initialize G0(w)
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
       end do
    end if
    
    ! Fourier transformation: G0(w) -> G0(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)-parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t)
    !call wn2tau_mid_point_taus(dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t, parm_%beta) ! by arijit
    do i=1, parm_%N_tau
       tau=dble(i-1)*parm_%dtau
       dmft_%weiss_green%matsubara_t(i)=dmft_%weiss_green%matsubara_t(i)-0.5d0+0.25d0*parm_%e2*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)+parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    
  end subroutine initialize_eq_green

  subroutine read_eq_green(parm_,dmft_)

    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                	    :: i,j
    double precision                :: w, re, im, tau
    !open(unit=24,file='sig1upA.dat',status='unknown')
    open (unit=13, file='Gf_old', status='unknown', action='read')
    do j=1, parm_%N_tau
       read (13,*) w, re, im
       dmft_%weiss_green%matsubara_w(j)= 1/(xj*w  - re - xj*im)
       !dmft_%weiss_green%matsubara_w(j)= re + im
   end do
     
    ! Fourier transformation: G0(w) -> G0(t)
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)-1.d0/(xj*w)-parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    call fft_w2t(parm_,dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t)
    !call wn2tau_mid_point_taus(dmft_%weiss_green%matsubara_w,dmft_%weiss_green%matsubara_t, parm_%beta) ! by arijit
    do i=1, parm_%N_tau
       tau=dble(i-0.5)*parm_%dtau
       dmft_%weiss_green%matsubara_t(i)=dmft_%weiss_green%matsubara_t(i)-0.5d0+0.25d0*parm_%e2*tau*(parm_%beta-tau)
    end do
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       dmft_%weiss_green%matsubara_w(j)=dmft_%weiss_green%matsubara_w(j)+1.d0/(xj*w)+parm_%e2/(xj*w)/(xj*w)/(xj*w)
    end do
    
    close (13)
  end subroutine read_eq_green
  
  
  subroutine measure_density(parm_,dmft_)
    !  n(t)=<c^+(t)c(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    double precision                :: slop
    integer                         :: npoint
    
    npoint = 1 ! no of point for liner fit
    slop = (dmft_%local_green%matsubara_t(parm_%N_tau) - dmft_%local_green%matsubara_t(parm_%N_tau-npoint))/dble(npoint)
    dmft_%density = -(dmft_%local_green%matsubara_t(parm_%N_tau) + slop/2.0)
    open (unit=7, file='density', status='old', action='write', position='append')
    write (7,'(f15.10,f15.10,f15.10)') parm_%U_i, dmft_%density 
    close (7)
    
  end subroutine measure_density
  
  
  subroutine measure_double_occupancy(parm_,dmft_)
    !  d(t)=<n_up(t)*n_do(t)>
    type(parm), intent(in)        :: parm_
    type(dmft), intent(inout)     :: dmft_
    integer                       :: i
    double precision, allocatable :: SxG(:)
    
    ! d(t)=n_up(t)*n_do(t)-1/U*\int dtau self_energy^{M}(beta-tau)*G^{M}(tau)
    if (parm_%U_i==0.d0) then
       dmft_%double_occu = (-dmft_%local_green%matsubara_t(parm_%N_tau))**2 ! density^2
    else
       allocate(SxG(parm_%N_tau))
       do i=1, parm_%N_tau+1
          SxG(i)=dmft_%self_energy%matsubara_t(parm_%N_tau-i+2)*dmft_%local_green%matsubara_t(i)
       end do
       dmft_%double_occu=dmft_%density**2 -1.d0/parm_%U_i*parm_%dtau*trapezoid_d(SxG,1,parm_%N_tau)
       deallocate(SxG)
    end if
    
    !open (unit=7, file='double-occupancy', status='old', action='write', position='append')
    !write (7,'(f15.10,f15.10)') t, double_occu
    !close (7)   
  end subroutine measure_double_occupancy
  
  
  subroutine measure_kinetic_energy(parm_,dmft_)
    !  E_{kin}(t)=2\sum_{k} E(k)<c_{k}^+(t)c_{k}(t)>
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                         :: i
    double precision, allocatable   :: GxG(:)
    
    if (parm_%dos=='semicircular') then
       ! E_{kin}(0)=-2*\int_0^{beta} dtau G^M(tau)*G^M(beta-tau)
       allocate(GxG(parm_%N_tau))
       do i=1, parm_%N_tau
          GxG(i)=dmft_%local_green%matsubara_t(i)*dmft_%local_green%matsubara_t(parm_%N_tau-i+1)
       end do
       dmft_%kinetic_energy =-2.d0*parm_%dtau*trapezoid_d(GxG,1,parm_%N_tau)
       deallocate(GxG)
    end if
    
    !open (unit=7, file='kinetic-energy', status='old', action='write', position='append')
    !write (7,'(f15.10,f15.10)') t, kinetic_energy
    !close (7)    
  end subroutine measure_kinetic_energy

   subroutine print_green_function(parm_,dmft_)
    type(parm), intent(in)          :: parm_
    type(dmft), intent(inout)       :: dmft_
    integer                	    :: j
    double precision                :: w, tau
    
    
    open (unit=12, file='Gf', status='old', action='write')
    open (unit=13, file='G0', status='old', action='write')
    open (unit=14, file='Sig', status='old', action='write')
    do j=1, parm_%N_tau
       w=dble(2*j-parm_%N_tau-1)*pi/parm_%beta
       !write (12,'(f15.10,f15.10,f15.10)') w, REALPART(dmft_%local_green%matsubara_w(j)), IMAGPART(dmft_%local_green%matsubara_w(j))
       write (12,*) w, REALPART(dmft_%local_green%matsubara_w(j)), IMAGPART(dmft_%local_green%matsubara_w(j))
       write (13,*) w, REALPART(dmft_%weiss_green%matsubara_w(j)), IMAGPART(dmft_%weiss_green%matsubara_w(j))
       write (14,*) w, REALPART(dmft_%self_energy%matsubara_w(j)), IMAGPART(dmft_%self_energy%matsubara_w(j))
    end do    
    close (12)
    close (13)
    close (14)
    
    open (unit=15, file='Gf_tau', status='old', action='write')
    open (unit=16, file='G0_tau', status='old', action='write')
    open (unit=17, file='Sig_tau', status='old', action='write')
    do j=1, parm_%N_tau
        tau=dble(j-0.5)*parm_%dtau
        write (15,*) tau, real(dmft_%local_green%matsubara_t(j))!, IMAGPART(dmft_%local_green%matsubara_t(j))
        write (16,*) tau, real(dmft_%weiss_green%matsubara_t(j))!, IMAGPART(dmft_%weiss_green%matsubara_t(j))
        write (17,*) tau, real(dmft_%self_energy%matsubara_t(j))!, IMAGPART(dmft_%self_energy%matsubara_t(j))
    end do    
    close (15)
    close (16)
    close (17)

    open (unit=17,file="observables", status="replace", action="write")
    write(17,*) "# U_i   density   double_occu  kinetic_energy interaction_enegy total_energy"
    write(17,*) parm_%U_i, dmft_%density, dmft_%double_occu, dmft_%kinetic_energy, dmft_%double_occu*parm_%U_i, dmft_%kinetic_energy + dmft_%double_occu*parm_%U_i
    close(17)
  end subroutine print_green_function
  
  
end module EQ_DMFT_MOD
