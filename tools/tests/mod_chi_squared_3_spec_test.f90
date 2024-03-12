module mod_chi_squared_3_spec_test
use fruit
use mod_pcg32_rng
use mod_sobseq_rng
use mod_sampling
implicit none

private
public :: run_fruit_chi_squared_3_spec

contains

!> Fruit basket -----------------------------------------
subroutine run_fruit_chi_squared_3_spec
  implicit none
  write(*,'(/A)') "  ... setting-up: chi squared 3 spec tests"
  write(*,'(/A)') "  ... running: chi squared 3 spec tests"
  call run_test_case(test_chi_squared_mean,'test_chi_squared_mean')
  call run_test_case(test_chi_squared_mean_sobseq,'test_chi_squared_mean_sobseq')
  write(*,'(/A)') "  ... tearing-down: chi squared 3 spec tests"
end subroutine run_fruit_chi_squared_3_spec 

!> Tests ------------------------------------------------
subroutine test_chi_squared_mean
  type(pcg32_rng) :: rng
  integer :: ifail, i
  integer, parameter :: n_tries = 1000
  real*8 :: x(n_tries)
  call rng%initialize(n_dims=1, seed=1231789264, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  do i=1,n_tries
    call rng%next(x(i:i))
    x(i) = sample_chi_squared_3(x(i))
  end do

  call assert_equals(sum(x)/n_tries, 3.d0, 0.2d0, "Mean should be 3")
end subroutine test_chi_squared_mean

!> Test whether the minimum and maximum values are between 0 and 1
subroutine test_chi_squared_mean_sobseq
  type(sobseq_rng) :: rng
  integer :: ifail, i
  integer, parameter :: n_tries = 1000
  real*8 :: x(n_tries)
  call rng%initialize(n_dims=1, seed=0, n_streams=1, i_stream=1, ierr=ifail)
  call assert_equals(0, ifail)
  do i=1,n_tries
    call rng%next(x(i:i))
    x(i) = sample_chi_squared_3(x(i))
  end do

  call assert_equals(sum(x)/n_tries, 3.d0, 0.1d0, "Mean should be 3")
end subroutine test_chi_squared_mean_sobseq

!> ------------------------------------------------------
end module mod_chi_squared_3_spec_test
