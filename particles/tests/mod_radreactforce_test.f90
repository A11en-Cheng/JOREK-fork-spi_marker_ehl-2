module mod_radreactforce_test 
use fruit
use mod_particle_types, only: particle_kinetic_relativistic
use mod_particle_types, only: particle_gc_relativistic
implicit none
private
public :: run_fruit_radreactforce
!> Variables --------------------------------------
integer(kind=1),parameter :: q=-1
real*8,parameter          :: tol_1dn6=1.d-6
real*8,parameter          :: tol_1dn5=1.d-5
real*8,parameter          :: predefval=457450.d0
real*8,parameter          :: scale_factor=4.d0
real*8,parameter          :: energy=1e6
real*8,parameter          :: pitch_parallel=1.d0
real*8,parameter          :: pitch_half=5d-1
real*8,parameter          :: dt=1.d-3
real*8,parameter          :: radreactforce_par=0.d0
real*8,parameter          :: radreactforce_pitchhalf_scaled=1.6d1
real*8,parameter          :: Bnorm=1.d0
real*8,dimension(3),parameter :: B=(/0.d0,0.d0,1.d0/)
real*8,dimension(3),parameter :: x=(/1.d0,0.d0,0.d0/)
type(particle_kinetic_relativistic) :: prt
type(particle_gc_relativistic)      :: gc
real*8 :: mass
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_radreactforce
  implicit none
  write(*,'(/A)') "  ... setting-up: radreactforce test"
  call setup
  write(*,'(/A)') "  ... running: radreactforce test"
  call run_test_case(test_radreactforce_parallel,'test_radreactforce_parallel')
  call run_test_case(test_radreactforce_Bscaling,'test_radreactforce_Bscaling')
  call run_test_case(test_radreactforce_integration,'test_radreactforce_integration')
  write(*,'(/A)') "  ... tearing-down: radreactforce test"
  call teardown
end subroutine run_fruit_radreactforce

!> Set-ups and tear-downs -------------------------
subroutine setup()
  use constants, only: MASS_ELECTRON,ATOMIC_MASS_UNIT
  implicit none
  mass = MASS_ELECTRON/ATOMIC_MASS_UNIT
  prt%q = q; gc%q = q; prt%x = x; gc%x = x;
  prt%p = 0.d0; gc%p = 0.d0;
end subroutine setup

subroutine teardown()
  implicit none
  mass = 0.d0
  prt%q = 0; gc%q = 0; prt%x = 0.d0; gc%x = 0.d0;
  prt%p = 0.d0; gc%p = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------
!> Check radiation reaction force is 0 for purely parallel motion
subroutine test_radreactforce_parallel()
  use mod_radreactforce, only: radreactforce_kinetic
  use mod_radreactforce, only: radreactforce_gc
  implicit none
  real*8              :: p0_gc,p_gc
  real*8,dimension(3) :: p0_prt
  !> initial condition
  call momentum_prt_and_gc(energy,pitch_parallel,mass,Bnorm,prt,gc)
  p0_prt = prt%p; p0_gc = gc2momentum(gc,mass,Bnorm)
  !> execute one step
  call radreactforce_kinetic(B,dt,mass,prt)
  call radreactforce_gc(B,dt,mass,gc)
  p_gc = gc2momentum(gc,mass,Bnorm)
  !> check error
  call assert_equals(radreactforce_par,norm2(prt%p-p0_prt)/norm2(p0_prt),&
  tol_1dn6,'Error parallel radreact force kinetic: force not 0!')
  call assert_equals(radreactforce_par,abs((p_gc-p0_gc)/p0_gc),tol_1dn6,&
  'Error parallel radreact force gc: force not 0!')
  !> restore zero momentum
  prt%p = 0d0; gc%p = 0d0;
end subroutine test_radreactforce_parallel

!> Check B^-2 scaling of radiation reaction force
subroutine test_radreactforce_Bscaling()
  use mod_radreactforce, only: radreactforce_kinetic
  use mod_radreactforce, only: radreactforce_gc
  implicit none
  real*8              :: p0_gc,ptemp_gc
  real*8,dimension(3) :: p0_prt,ptemp_prt
  call momentum_prt_and_gc(energy,pitch_half,mass,Bnorm,prt,gc)
  p0_prt = prt%p; p0_gc = gc%p(1);
  !> perform one step
  call radreactforce_kinetic(B,dt,mass,prt)
  call radreactforce_gc(B,dt,mass,gc)
  ptemp_prt = prt%p; ptemp_gc = gc%p(1);
  !> compute new solution with scaled magnetic field
  call momentum_prt_and_gc(energy,pitch_half,mass,scale_factor*Bnorm,prt,gc)
  call radreactforce_kinetic(scale_factor*B,dt,mass,prt)
  call radreactforce_gc(scale_factor*B,dt,mass,gc) 
  !> check solution
  call assert_equals(radreactforce_pitchhalf_scaled,&
  norm2(prt%p-p0_prt)/norm2(ptemp_prt-p0_prt),tol_1dn6,&
  'Error radreact force kinetic w scaled B and pitch_half: incorrect scaling!')
  call assert_equals(radreactforce_pitchhalf_scaled,&
  (gc%p(1)-p0_gc)/(ptemp_gc-p0_gc),tol_1dn6,&
  'Error radreact force kinetic w scaled B and pitch_half: incorrect scaling!')
  !> restore zero
  prt%p = 0.d0; gc%p = 0.d0;
end subroutine test_radreactforce_Bscaling

!> Running a proper physics test (bump-on-tail) is too expensive so we only
!> check that the operator has not changed since the physics tests were done.
!> We do this by comparing the result to a precomputed value (evaluated
!> when the operator was confirmed to pass the physics test). Also
!> check that prt == gc (when the field is uniform)
subroutine test_radreactforce_integration()
  use mod_radreactforce, only: radreactforce_kinetic
  use mod_radreactforce, only: radreactforce_gc
  use mod_radreactforce, only: radreactforce_gc_rhs
  implicit none
  real*8              :: p0_gc,p_gc
  real*8,dimension(2) :: gc_p_rhs
  real*8,dimension(3) :: p0_prt
  !> compute momentum before integration step
  call momentum_prt_and_gc(energy,pitch_half,mass,Bnorm,prt,gc)
  call radreactforce_kinetic(B,dt,mass,prt)
  call radreactforce_gc(B,dt,mass,gc)
  p0_prt = prt%p; p0_gc = gc2momentum(gc,mass,Bnorm);
  !> compute an integration step
  call momentum_prt_and_gc(energy,pitch_half,mass,Bnorm,prt,gc)
  call radreactforce_gc_rhs(B,mass,int(gc%q),gc%p,gc_p_rhs)
  gc%p = gc%p + dt * gc_p_rhs; p_gc = gc2momentum(gc,mass,Bnorm);
  !> checks
  call assert_equals(0.d0,(p0_gc-norm2(p0_prt))/p0_gc,tol_1dn6,&
  'Error radreactforce integration: different kinetic and gc momenta!')
  call assert_equals(0.d0,(p0_gc-p_gc)/p0_gc,tol_1dn6,&
  'Error radreactforce integration: different GC integration and integration RHS!')
  call assert_equals(0.d0,(p0_gc-predefval)/p0_gc,tol_1dn5,&
  'Error radreactforce integration: different GC momentum than expected value!')
  !> restore zeros
  prt%p = 0.d0; gc%p = 0.d0;
end subroutine test_radreactforce_integration

!> Tools ------------------------------------------
!> single function to initialize particle and gc momentum coordinates
!> from energy and pitch angle
subroutine momentum_prt_and_gc(energy,pitch,mass,Bnorm,prt,gc)
  use constants, only: EL_CHG,SPEED_OF_LIGHT,ATOMIC_MASS_UNIT
  implicit none
  real*8,intent(in) :: energy,pitch,mass,Bnorm
  type(particle_kinetic_relativistic),intent(inout) :: prt
  type(particle_gc_relativistic),intent(inout)      :: gc
  real*8 :: p
  p = sqrt( (energy * EL_CHG / ( mass * SPEED_OF_LIGHT**2 * ATOMIC_MASS_UNIT) + &
      1.d0)**2 - 1.d0)*mass*SPEED_OF_LIGHT
  prt%p = (/sqrt(1.d0 - pitch**2) * p, pitch * p, 0.d0/)
  gc%p  = (/pitch * p, ( 1.d0 - pitch**2 ) * p**2 / (2.d0 * mass * Bnorm )/)
end subroutine momentum_prt_and_gc

!> Evaluate momentum (norm) from gc momentum coordinates
function gc2momentum(gc,mass,Bnorm) result(p)
  implicit none
  type(particle_gc_relativistic),intent(in) :: gc
  real*8,intent(in) :: mass, Bnorm !< mass [amu], B-field magnitude [T]
  real*8 :: p
  p = sqrt(gc%p(1)**2 + 2.d0 * mass * Bnorm * gc%p(2) )
end function gc2momentum
!> ------------------------------------------------
end module mod_radreactforce_test
