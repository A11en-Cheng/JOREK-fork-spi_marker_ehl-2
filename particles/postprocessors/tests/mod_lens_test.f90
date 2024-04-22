!> mod_lens_test contains all variable and procedures required
!> for testing the lens class. Due to the fact that lens is an 
!> abstract class, the pinhole_lens class instead
module mod_lens_test
  use fruit
  implicit none

  private
  public :: run_fruit_lens

!> Variables --------------------------------------------------
integer,parameter :: n_x=3
real*8,parameter :: tol_r8=5.d-16
real*8,dimension(n_x),parameter :: center_uppbnd=(/4.5d3,4.67d2,-2.34d-1/)
real*8,dimension(n_x),parameter :: center_lowbnd=(/-1.d2,2.45d1,-2.3d2/)
real*8,dimension(n_x) :: center_sol

!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> basket running the set-up, test and teard-down procedures
subroutine run_fruit_lens()
  implicit none
  write(*,'(/A)') "  ... setting-up: lens tests"
  call setup
  write(*,'(/A)') "  ... running: lens tests"
  call run_test_case(test_lens_allocation_deallocation,&
  'test_lens_allocation_deallocation')
  write(*,'(/A)') "  ... tearing-down: lens tests"
  call teardown
end subroutine run_fruit_lens

!> Set-up and teard-down---------------------------------------
!> set-up test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> initialize center array
  call gnu_rng_interval(n_x,center_lowbnd,center_uppbnd,center_sol)
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  center_sol = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------------------
!> test the lens allocation and deallocation
subroutine test_lens_allocation_deallocation()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_pinhole_lens,        only: pinhole_lens
  implicit none
  !> variables
  type(pinhole_lens) :: lens_loc
  !> allocate lens initialised to zero
  call lens_loc%allocate_lens(n_x)
  call assert_equals(n_x,lens_loc%n_x,&
  "Error lens allocation to zero: n_x mismatch!")
  call assert_equals_allocatable_arrays(n_x,0.d0,lens_loc%center,&
  "Error lens allocation to zero: center array not allocated!")
  !> deallocate lens
  call lens_loc%deallocate_lens()
  call assert_equals(0,lens_loc%n_x,&
  "Error lens deallocation: n_x not to zero!")
  call assert_false(allocated(lens_loc%center),&
  "Error lens deallocation: center not deallocated!")
  !> allocate lens with initialisation
  call lens_loc%allocate_lens(n_x,center_sol)
  call assert_equals(n_x,lens_loc%n_x,&
  "Error lens allocation to value: n_x mismatch!")
  call assert_equals_allocatable_arrays(n_x,center_sol,lens_loc%center,&
  tol_r8,"Error lens allocation to value: center array")
  !> deallocate lens
  call lens_loc%deallocate_lens()
  call assert_equals(0,lens_loc%n_x,&
  "Error lens deallocation: n_x not to zero!")
  call assert_false(allocated(lens_loc%center),&
  "Error lens deallocation: center not deallocated!")
end subroutine test_lens_allocation_deallocation

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_lens_test
