module mod_pusher_boris_penning_spec_test
use fruit
implicit none
private
public :: run_fruit_pusher_boris_penning_spec
!> Variables --------------------------------------
integer,parameter                   :: n_tests=2
integer,parameter                   :: message_len=100
real*8,parameter                    :: time_sol=0.d0
real*8,parameter                    :: expect_sol=0.d0
real*8,dimension(n_tests),parameter :: dts=(/1.d-4,1.d-5/)
real*8,dimension(n_tests),parameter :: tols=(/1.1d-4,1.1d-6/)
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_pusher_boris_penning_spec
  implicit none
  write(*,'(/A)') "  ... setting-up: pusher boris penning spec"
  write(*,'(/A)') "  ... running: pusher boris penning spec"
  call run_test_case(test_penning_case_cartesian,'test_penning_case_cartesian')
  write(*,'(/A)') "  ... tearing-down: pusher_boris_penning_spec"
end subroutine run_fruit_pusher_boris_penning_spec

!> Tests ------------------------------------------
subroutine test_penning_case_cartesian()
  use mod_particle_types, only: particle_base
  use mod_particle_types, only: particle_kinetic_leapfrog
  use mod_boris,          only: boris_initial_half_step_backwards_XYZ
  use mod_boris,          only: boris_push_cartesian
  use mod_penning_case,   only: case_penning_cartesian
  implicit none
  class(particle_base),allocatable :: particle
  type(case_penning_cartesian)     :: case
  integer                    :: ii,jj,n_steps
  real*8                     :: err
  character(len=message_len) :: message
  do jj=1,n_tests
    !> initialize
    n_steps = nint(case%time_end/dts(jj))
    allocate(particle_kinetic_leapfrog::particle)
    call case%initialize(particle)
    !> run test case
    select type (p=>particle)
    type is (particle_kinetic_leapfrog)
      call assert_equals(1,int(p%q,kind=4),&
      'Error penning case cartesian test: unexpected particle charge!')
      call boris_initial_half_step_backwards_XYZ(p,case%mass,&
      case%E(p%x,time_sol),case%B(p%x,time_sol),dts(jj))
      !> compute orbit
      do ii=1,n_steps
        call boris_push_cartesian(p,case%mass,case%E(p%x,time_sol),& 
        case%B(p%x,time_sol),dts(jj))
      enddo
      err = 9.99d9; err = case%calc_error(p) !< compute error
    end select
    !> checks
    write(message,'(A,F0.16,A)') 'Error penning case cartesian test: error > ',tols(jj),'!'
    call assert_equals(expect_sol,err,tols(jj),trim(message))
    !> cleanup
    if(allocated(particle)) deallocate(particle)
  enddo
end subroutine test_penning_case_cartesian

subroutine test_penning_case_cylindrical()
  use mod_particle_types, only: particle_kinetic_leapfrog
  use mod_boris,          only: boris_initial_half_step_backwards_RZPhi
  use mod_boris,          only: boris_push_cylindrical
  use mod_penning_case,   only: case_penning_cylindrical
  implicit none
  type(particle_kinetic_leapfrog) :: particle
  type(case_penning_cylindrical)  :: case
  integer                    :: ii,jj,n_steps
  real*8                     :: err
  character(len=message_len) :: message
  do jj=1,n_tests
    !> initialize
    n_steps = nint(case%time_end/dts(jj))
    call case%initialize(particle)
    call assert_equals(1,int(particle%q,kind=4),&
    'Error penning case cylindrical test: unexpected particle charge!')
    call boris_initial_half_step_backwards_RZPhi(particle,case%mass,&
    case%E(particle%x,time_sol),case%B(particle%x,time_sol),dts(jj))
    !> compute orbit
    do ii=1,n_steps
      call boris_push_cylindrical(particle,case%mass,case%E(particle%x,time_sol),&
      case%B(particle%x,time_sol),dts(jj))
    enddo
    err = 9.99d9; err = case%calc_error(particle) !< compute error
    !> checks
    write(message,'(A,F0.16,A)') 'Error penning case cylindrical test: error > ',tols(jj),'!'
    call assert_equals(expect_sol,err,tols(jj),trim(message))
  enddo
end subroutine test_penning_case_cylindrical

!> ------------------------------------------------
end module mod_pusher_boris_penning_spec_test
