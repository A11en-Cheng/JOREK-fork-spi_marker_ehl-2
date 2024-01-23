module mod_pusher_boris_spec_test
use fruit
use constants,          only: TWOPI,EL_CHG,ATOMIC_MASS_UNIT
use mod_particle_types, only: particle_kinetic_leapfrog
implicit none
private
public :: run_fruit_pusher_boris_spec
!> Variables --------------------------------------
!> TODO: get accurate v^(-1/2), Delzanno, see JCP(2013) pusher_test
integer(kind=1),parameter       :: q=1
real*8,parameter                :: tol_1den1=1.d-1
real*8,parameter                :: tol_1den3=1.d-3
real*8,parameter                :: tol_1den5=1.d-5
real*8,parameter                :: tol_1den7=1.d-7
real*8,parameter                :: time_half_gyro=5.d-1
real*8,parameter                :: timestep_base=1.d-2
real*8,parameter                :: mass=EL_CHG/(TWOPI*ATOMIC_MASS_UNIT)
real*8,dimension(3),parameter   :: x_init=(/1.d0,0.d0,0.d0/)
real*8,dimension(3),parameter   :: v_init=(/0.d0,TWOPI,0.d0/)
real*8,dimension(3),parameter   :: B=(/0.d0,0.d0,1.d0/)
real*8,dimension(3),parameter   :: E=(/0.d0,0.d0,0.d0/)
type(particle_kinetic_leapfrog) :: prt
real*8 :: gyro_pulse
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_pusher_boris_spec
  implicit none
  write(*,'(/A)') "  ... setting-up: pusher boris spec"
  call setup
  write(*,'(/A)') "  ... running: pusher boris spec"
  call run_test_case(test_single_gyro_orbit,'test_single_gyro_orbit')
  call run_test_case(test_half_gyro_orbit_convergence,'test_half_gyro_orbit_convergence')
  write(*,'(/A)') "  ... tearing-down: pusher boris spec"
  call teardown
end subroutine run_fruit_pusher_boris_spec

!> Set-ups and teardown
subroutine setup()
  implicit none
  prt%x = x_init; prt%v = v_init; prt%q = q;
  gyro_pulse =(real(q,kind=8)*EL_CHG*B(3))/(ATOMIC_MASS_UNIT*mass)
  write(*,*) gyro_pulse/TWOPI
end subroutine setup

subroutine teardown()
  implicit none
  prt%x = 0.d0; prt%v = 0.d0; prt%q = 0; gyro_pulse = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------
!> Test the gyromotion of a single particle, orbiting with radius 1 in 1s
subroutine test_single_gyro_orbit()
  implicit none
  integer,parameter :: n_tests_loc=4
  real*8,parameter  :: time_loc=2.5d-1
  real*8,dimension(n_tests_loc),parameter :: tol_x_loc=(/1.d-1,1.d-3,1.d-1,1.d-7/)
  real*8,dimension(n_tests_loc),parameter :: tol_y_loc=(/1.d-1,1.d-1,1.d-1,1.d-7/)
  real*8,dimension(n_tests_loc),parameter :: x_sol_loc=(/0.d0,-1.d0,0.d0,1.d0/)
  real*8,dimension(n_tests_loc),parameter :: y_sol_loc=(/1.d0,0.d0,-1.d0,0.d0/)
  integer :: ii
  real*8,dimension(2) :: xy_sol,vxvy_sol
  character(len=70) :: message_x,message_y
  do ii=1,n_tests_loc 
    call restore_particle_kinetic_leapfrog(x_init,v_init,prt)
    call integrate_orbit(timestep_base,time_loc*real(ii,kind=8),mass,B,E,prt)
    call analytical_gyro_orbit_xy(time_loc*real(ii,kind=8),gyro_pulse,&
    x_init(1:2),v_init(1:2),xy_sol,vxvy_sol)
    write(message_x,'(A,I0,A)') 'Error single gryo-orbit test N#: ',ii,' wrong x-coord!'
    write(message_y,'(A,I0,A)') 'Error single gyro-orbit test N#: ',ii,' wrong y-coord!'
    call assert_equals(xy_sol(1),prt%x(1),tol_x_loc(ii),message_x)
    call assert_equals(xy_sol(2),prt%x(2),tol_y_loc(ii),message_y)
  enddo
  call restore_particle_kinetic_leapfrog(x_init,v_init,prt)
end subroutine test_single_gyro_orbit

!> TODO: we still have linear convergence for the y-component
!> (which is also the quickly-varying one)
subroutine test_half_gyro_orbit_convergence()
  implicit none
  integer,parameter :: n_tests_loc=5
  real*8,dimension(n_tests_loc),parameter :: tol_x_loc=(/1.d-3,1.d-5,1.d-7,1.d-7,1.d-7/)
  real*8,dimension(n_tests_loc),parameter :: tol_y_loc=(/1.d-1,1.d-2,1.d-3,1.d-4,1.d-5/)
  integer :: ii
  real*8,dimension(2) :: xy_sol,vxvy_sol
  character(len=70) :: message_x,message_y
  call analytical_gyro_orbit_xy(time_half_gyro,gyro_pulse,x_init(1:2),v_init(1:2),xy_sol,vxvy_sol)
  do ii=1,n_tests_loc
    call restore_particle_kinetic_leapfrog(x_init,v_init,prt)
    call integrate_orbit(10.d0**(-(ii+1)),time_half_gyro,mass,B,E,prt)
    write(message_x,'(A,I0,A)') 'Error half orbit convergence test N#: ',ii,' x-coord not converged!'
    write(message_y,'(A,I0,A)') 'Error half orbit convergence test N#: ',ii,' y-coord not converged!'
    call assert_equals(xy_sol(1),prt%x(1),tol_x_loc(ii),trim(message_x))
    call assert_equals(xy_sol(2),prt%x(2),tol_y_loc(ii),trim(message_y))
  enddo
  call restore_particle_kinetic_leapfrog(x_init,v_init,prt)
end subroutine test_half_gyro_orbit_convergence

!> Tools ------------------------------------------
subroutine restore_particle_kinetic_leapfrog(x,v,particle)
  implicit none
  real*8,dimension(3),intent(in) :: x,v
  type(particle_kinetic_leapfrog),intent(inout) :: particle
  particle%x = x; particle%v = v;
end subroutine restore_particle_kinetic_leapfrog

subroutine integrate_orbit(timestep,time,mass_loc,B_loc,E_loc,particle)
  use mod_boris, only: boris_push_cartesian
  implicit none
  real*8,intent(in) :: timestep,time,mass_loc
  real*8,dimension(3),intent(in) :: E_loc,B_loc
  type(particle_kinetic_leapfrog),intent(inout) :: particle
  integer :: ii
  do ii=1,nint(time/timestep)
    call boris_push_cartesian(particle,mass_loc,E_loc,B_loc,timestep)
  enddo
end subroutine integrate_orbit

!> analytical solution to the gyromotion
subroutine analytical_gyro_orbit_xy(time,pulse,xyinit,vxvyinit,xy,vxvy)
  implicit none
  real*8,intent(in) :: time,pulse
  real*8,dimension(2),intent(in)  :: xyinit,vxvyinit
  real*8,dimension(2),intent(out) :: xy,vxvy
  real*8 :: cpt,spt
  cpt = cos(pulse*time); spt=sin(pulse*time)
  xy(1) = xyinit(1) + (vxvyinit(1)*spt - vxvyinit(2)*(cpt-1.d0))/pulse
  xy(2) = xyinit(2) + (vxvyinit(1)*(cpt-1.d0) + vxvyinit(2)*spt)/pulse
  vxvy(1) = vxvyinit(1)*cpt + vxvyinit(2)*spt
  vxvy(2) = -vxvyinit(1)*spt + vxvyinit(2)*cpt
end subroutine analytical_gyro_orbit_xy

!> ------------------------------------------------
end module mod_pusher_boris_spec_test
