!----------------------------------------------------------------------------
!
!  date:   May 23, 2017
!
!----------------------------------------------------------------------------
module PARM_MOD
  
  use CONST_MOD
  implicit none
  
  private
  public :: set_parm
  
  type, public :: parm
     character(100)     :: dos              ! density of states
     double precision   :: e2               ! <e^2>=\int de D(e)*e^2
     double precision   :: temperature      ! temperature of the initial thermal state
     double precision   :: beta             ! inverse temperature (beta=1/temperature)
     double precision   :: U_i              ! initial value of the interaction parameter U
     integer            :: N_tau            ! # of imaginary-time steps
     integer            :: N_iter           ! # of maximum iterations of DMFT self-consistency
     double precision   :: dtau             ! dtau=beta/N_tau
     double precision   :: tolerance        ! tolerance for DMFT convergence
     character(100)     :: solver           ! impurity solver
     integer            :: init             ! # 0 when G0 is non-interacting GF. 1 when G0 is interacting GF of previous U value
     double precision   :: delta            ! # ionic potetial
     double precision   :: mix              ! # G = (1-mix)*G_new + G except first iteration

  end type parm
  
contains
  
  subroutine set_parm(parm_)
    type(parm), intent(inout) :: parm_
    character(100)            :: arg
    integer                   :: n_arg
    
    call getarg(1,parm_%dos)
    if (parm_%dos=='semicircular') then
       parm_%e2=0.d0
    end if

    call getarg(2,arg)
    read(arg,*) parm_%beta
    parm_%temperature=1.d0/parm_%beta
    call getarg(3,arg)
    read(arg,*) parm_%U_i
    call getarg(4,arg)
    read(arg,*) parm_%N_tau
    parm_%dtau=parm_%beta/dble(parm_%N_tau)
    call getarg(5,arg)
    read(arg,*) parm_%N_iter
    call getarg(6,arg)
    read(arg,*) parm_%tolerance
    call getarg(7,parm_%solver)
    print*, parm_%solver,"------------------------"
    call getarg(8,arg)
    read(arg,*) parm_%init
    print*, parm_%init,"------init------------------"
    call getarg(9,arg)
    read(arg,*) parm_%mix
    !read(arg,*) parm_%init
    !parm_%init = 1
    call getarg(10,arg)
    read(arg,*) parm_%delta
    print*,"mixing--------------",parm_%mix

    open(unit=13, file='PARAMS', status='unknown', action='write')
    write(13,*) "dos:  ", parm_%dos
    write(13,*) "beta:  ",parm_%beta
    write(13,*) "U:  ", parm_%U_i
    write(13,*) "N_tau:  ",parm_%N_tau
    write(13,*) "N_iter:  ", parm_%N_iter
    write(13,*) "tolarence:  ",parm_%tolerance
    write(13,*) "solver:  ", parm_%solver
    write(13,*) "init:  ",parm_%init
    write(13,*) "mix:  ",parm_%mix
    write(13,*) "delta:  ",parm_%delta
    close(13)

  end subroutine set_parm
  
  
 
  
end module PARM_MOD