module imaginary_time_frequency_FFTs
	use const_mod
	use fftw_module
	implicit none
contains
subroutine wn2tau_left_point_taus(fwns,ftaus,beta)
	implicit none
	double precision,intent(in)::beta
	double complex,dimension(:)::fwns,ftaus
	double complex::phase_fac
	double precision::phase
	integer::N,k
	N=size(fwns)
	call FFT_FORWARD(fwns,ftaus)
	do k=1,N
		phase=PI*(k-1)*(1-dble(N))/dble(N)
		phase_fac=exp(dcmplx(0.0d0,-phase))
		ftaus(k)=(ftaus(k)*phase_fac)/beta
	enddo
end subroutine
subroutine tau2wn_left_point_taus(ftaus,fwns,beta)
	implicit none
	double precision,intent(in)::beta
	double complex,dimension(:)::fwns,ftaus
	integer::N,k,j
	double precision::phase
	double complex::phase_facs
	double complex,dimension(:),allocatable::ftil_taus
	N=size(ftaus)
	allocate(ftil_taus(N))
	do k=1,N
		phase=PI*(k-1)*(1-dble(N))/dble(N)
		phase_facs=exp(dcmplx(0.0d0,phase))
		ftil_taus(k)=ftaus(k)*phase_facs
	enddo
	call FFT_BACKWARD(ftil_taus,fwns)

	do j=1,N
		fwns(j)=fwns(j)*(beta/dble(N))
	enddo
	deallocate(ftil_taus)
end subroutine



subroutine tau2wn_mid_point_taus(ftaus,fwns,beta) !***
	implicit none
	double precision,intent(in)::beta
	double complex,dimension(:)::fwns,ftaus
	double complex::phase_fac
	double precision::phase
	integer::N,k,j
	double complex,dimension(:),allocatable::ftil_taus
	N=size(fwns)
	allocate(ftil_taus(N))
	do k=1,N
		phase=(k-1)*PI*((1-N)/dble(N))
		phase_fac=exp(dcmplx(0.0d0,phase))
		ftil_taus(k)=ftaus(k)*phase_fac
	enddo
	call FFT_BACKWARD(ftil_taus,fwns)
	
	do j=1,N
		phase=(j-1)*PI/dble(N)+0.5d0*PI*(1-N)/dble(N)
		phase_fac=exp(dcmplx(0.0d0,phase))
		fwns(j)=(beta/dble(N))*fwns(j)*phase_fac
	enddo
	deallocate(ftil_taus)
end subroutine

subroutine wn2tau_mid_point_taus(fwns,ftaus,beta) !***
	implicit none
	double precision,intent(in)::beta
	double complex,dimension(:)::fwns,ftaus
	double complex::phase_fac
	double precision::phase
	integer::N,k,j
	double complex,dimension(:),allocatable::ftil_wns
	N=size(fwns)
	allocate(ftil_wns(N))
	do j=1,N
		phase=(j-1)*PI/dble(N)
		phase_fac=exp(dcmplx(0.0d0,-phase))
		ftil_wns(j)=fwns(j)*phase_fac
	enddo
	call FFT_FORWARD(ftil_wns,ftaus)
	
	do k=1,N
		phase=(k-1)*PI*((1-N)/dble(N))+0.5d0*PI*(1/dble(N)-1)
		phase_fac=exp(dcmplx(0.0d0,-phase))
		ftaus(k)=(1.0d0/beta)*ftaus(k)*phase_fac
	enddo
	deallocate(ftil_wns)
end subroutine


end module
