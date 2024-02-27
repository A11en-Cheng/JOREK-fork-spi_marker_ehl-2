module mod_rootfinding_spec_test
use fruit
implicit none
private
public :: run_fruit_rootfinding_spec
!> Variables --------------------------------------
integer,parameter         :: n_tests=8
real*8,parameter          :: tol=1d-8
real*8,parameter          :: divide=1.d1
real*8,parameter          :: delta=8.d-2
real*8,dimension(n_tests) :: y,ref
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_rootfinding_spec
  implicit none
  write(*,'(/A)') "  ... setting-up: root finding spec"
  call setup
  write(*,'(/A)') "  ... running: root finding spec"
  call run_test_case(test_newtons_method_f,'test_newtons_method_f')
  call run_test_case(test_halleys_method_f,'test_halleys_method_f')
  write(*,'(/A)') "  ... tearing-down: root finding spec"
  call teardown
end subroutine run_fruit_rootfinding_spec

!> Set-ups and tear-down --------------------------
subroutine setup()
  implicit none
  integer :: ii
  y   = (/(real(ii,kind=8)/divide,ii=1,n_tests)/)
  ref = (/(invf(y(ii)), ii=1,n_tests)/)
end subroutine setup

subroutine teardown
  implicit none
  y = 0.d0; ref = 0.d0;
end subroutine teardown

!> Tests ------------------------------------------
subroutine test_newtons_method_f()
  use mod_rootfinding, only: newtons_method
  implicit none
  integer :: ii,ierr
  real*8  :: x
  do ii=1,n_tests
    !> test approaching from above and below
    ierr=0; call newtons_method(f,df,y(ii),ref(ii)-delta,x,ierr)
    call assert_equals(ierr,0,'Error root finding newtons from below: solution not found!')
    call assert_equals(ref(ii),x,tol,'Error root finding newtons from below: incorrect root!')
    ierr=0; call newtons_method(f,df,y(ii),ref(ii)+delta,x,ierr)
    call assert_equals(ierr,0,'Error root finding newtons from above: solution not found!')
    call assert_equals(ref(ii),x,tol,'Error root finding newtons from above: incorrect root!')
  enddo
end subroutine test_newtons_method_f

subroutine test_halleys_method_f()
  use mod_rootfinding, only: halleys_method
  implicit none
  integer :: ii,ierr
  real*8  :: x
  do ii=1,n_tests
    !> test approaching from above and below
    ierr=0; call halleys_method(f,df,ddf,y(ii),ref(ii)-delta,x,ierr)
    call assert_equals(ierr,0,'Error root finding halleys from below: solution not found!')
    call assert_equals(ref(ii),x,tol,'Error root finding halleys from below: incorrect root!')
    ierr=0; call halleys_method(f,df,ddf,y(ii),ref(ii)+delta,x,ierr)
    call assert_equals(ierr,0,'Error root finding halleys from above: solution not found!')
    call assert_equals(ref(ii),x,tol,'Error root finding halleys from above: incorrect root!')
  enddo
end subroutine test_halleys_method_f

!> Tools ------------------------------------------
!> function, set to be tanh
pure function f(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: f
  f = tanh(x)
end function f

!> inverse of the function, set to be atanh
pure function invf(y)
  implicit none
  real*8,intent(in) :: y
  real*8            :: invf
  invf = atanh(y)
end function invf

!> first derivative of the function, 
!> set to be dtanh/dx
pure function df(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: df
  df = 1.d0 - tanh(x)*tanh(x)
end function df

!> second derivative of the function, 
!> set to be d^2tanh/dx^2
pure function ddf(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: ddf
  ddf = -2.d0*tanh(x)/cosh(x)/cosh(x)
end function ddf
!> ------------------------------------------------
end module mod_rootfinding_spec_test
