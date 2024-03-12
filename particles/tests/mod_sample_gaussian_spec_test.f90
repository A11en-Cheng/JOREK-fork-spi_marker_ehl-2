module mod_sample_gaussian_spec_test
use fruit
use mod_pcg32_rng
implicit none
private
public :: run_fruit_sample_gaussian_spec
!> Variables --------------------------------------
type(pcg32_rng)           :: rng
integer,parameter         :: n_tries=1000000
integer,parameter         :: ndims=1
integer,parameter         :: seed=1231789264
integer,parameter         :: nstreams=1
integer,parameter         :: istream=1
real*8,parameter          :: tol_mean=1.d-3
real*8,parameter          :: tol_var=5.d-4
integer                   :: ifail
real*8,dimension(n_tries) :: x
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_sample_gaussian_spec()
  implicit none
  write(*,'(/A)') "  ... setting-up: sample gaussian spec"
  call setup
  write(*,'(/A)') "  ... running: sample gaussian spec"
  call run_test_case(test_gaussian_mean_and_variance,'test_gaussian_mean_and_variance')
  write(*,'(/A)') "  ... tearing-down: sample gaussian spec"
  call teardown
end subroutine run_fruit_sample_gaussian_spec

!> Set-ups and tear-downs -------------------------
subroutine setup()
  implicit none
  call rng%initialize(n_dims=ndims,seed=seed,n_streams=nstreams,&
  i_stream=istream,ierr=ifail)
  x = 0.d0
  call assert_equals(0,ifail,'PCG32 RNG not initialised!')
end subroutine setup

subroutine teardown()
  implicit none
  x = 0.d0
end subroutine teardown

!> Tests ------------------------------------------
subroutine test_gaussian_mean_and_variance()
  use mod_sampling, only: sample_gaussian
  implicit none
  integer :: ii
  do ii=1,n_tries
    call rng%next(x(ii:ii))
    x(ii) = sample_gaussian(x(ii))
  enddo
  call assert_equals(sum(x)/real(n_tries,kind=8),0.d0,tol_mean,&
  'Error gaussian sampling: mean should be 0!')
  call assert_equals(dot_product(x,x)/real(n_tries,kind=8),1.d0,&
  tol_var,'Error gaussian sampling: variance should be 1!')
end subroutine test_gaussian_mean_and_variance
!> ------------------------------------------------
end module mod_sample_gaussian_spec_test
