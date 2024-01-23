!> Module testing our quicksort implementation
module mod_quicksort_test
use fruit
implicit none
private
public :: run_fruit_quicksort
!> Variables --------------------------------------
integer,parameter :: message_len=50
integer,parameter :: n_samples=100
real*8,parameter  :: test_value=1d0
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_quicksort
  implicit none
  write(*,'(/A)') "  ... setting-up: quicksort "
  write(*,'(/A)') "  ... running: quicksort "
  call run_test_case(test_quicksort_1,'test_quicksort_1')
  call run_test_case(test_quicksort_array,'test_quicksort_array')
  write(*,'(/A)') "  ... tearing-down: quicksort"
end subroutine run_fruit_quicksort

!> Tests ------------------------------------------
subroutine test_quicksort_1()
  use mod_quicksort, only: quicksort
  implicit none
  real*8,dimension(1) :: sample
  sample(1) = test_value; call quicksort(sample);
  call assert_equals(test_value,sample(1),&
  'Error quicksort 1 test: incorrect value for sorted array!')
end subroutine test_quicksort_1

subroutine test_quicksort_array()
  use mod_quicksort, only: quicksort
  implicit none
  integer :: ii
  real*8,dimension(n_samples) :: samples
  character(len=message_len)  :: message
  samples = 0d0
  do ii=1,n_samples
    samples(ii) = real((n_samples-ii)**2,kind=8)
  enddo
  call quicksort(samples)
  do ii=1,n_samples
    write(message,'(A,I0,A)') 'Error quicksort array test: ',ii,&
    'sample mismatch!'
    call assert_equals(real((ii-1)**2,kind=8),samples(ii),trim(message))
  enddo
end subroutine test_quicksort_array

!> ------------------------------------------------
end module mod_quicksort_test
