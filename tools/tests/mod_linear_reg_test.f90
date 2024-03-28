!> mod_linear_reg_test implements unit tests for 
!> linear regression methods
module mod_linear_reg_test
use fruit
use mod_dynamic_array_tools, only: allocate_check
use mod_dynamic_array_tools, only: deallocate_check
implicit none

private
public :: run_fruit_linear_reg 

!> Variables ---------------------------------------
integer,parameter :: N_points=1000  !< number of points
integer,parameter :: N_coeff=2      !< number of coefficients
real*8,parameter  :: tol=5.d-14     !< tolerance for checking results
real*8,parameter  :: tol_rse=1.d-13 !< expected residual std error
real*8,dimension(2),parameter       :: x_interval=(/-2.d1,1.d2/) !< x interval
!> exact line coefficients: slope,intersect
real*8,dimension(N_coeff),parameter :: line_coeff=(/-2.d-1,2.3d1/)
!> interval for the random noise
real*8,dimension(2),parameter :: noise_interval=(-1.d1,3.d0)
real*8,dimension(:),allocatable :: x               !< abscissas
real*8,dimension(:),allocatable :: y_exact,y_noise !< exact and noisy values

contains

!> Fruit basket ------------------------------------
!> fruit baskets: executes the set-up, runs the tests and
!> tears-down the linear regression features
subroutine run_fruit_linear_reg()
  implicit none

  !> execute set-up -> test -> tear-down
  write(*,'(/A)') "  ... setting-up: linear regression tests"
  call setup
  write(*,'(/A)') "  ... running: linear regression tests"
  call run_test_case(linear_ref_exact_y_test,'linear_ref_exact_y_test')
  call run_test_case(linear_reg_noise_y_test,'linear_reg_noise_y_test')
  write(*,'(/A)') "  ... tearing-down: linear regression tests"
  call teardown

end subroutine run_fruit_linear_reg

!> Set-up and tear-down-----------------------------
!> set-up the linear regression test features
subroutine setup()
  implicit none
  
  !> allocate arrays
  call allocate_check(n_points,x)
  call allocate_check(n_points,y_exact)
  call allocate_check(n_points,y_noise)

  !> compute abscissa
  call random_number(x)
  x = x_interval(1)+(x_interval(2)-x_interval(1))*x
  !> compute exact and noisy values
  y_exact = line_coeff(2) + line_coeff(1)*x
  call random_number(y_noise)
  y_noise = y_exact + noise_interval(1) + &
  (noise_interval(2)-noise_interval(1))*y_noise

end subroutine setup

!> tear-down the linear regression test features
subroutine teardown()
  implicit none
  !> deallocate arrays
  call deallocate_check(x)
  call deallocate_check(y_exact)
  call deallocate_check(y_noise)
end subroutine teardown

!> Tests -------------------------------------------

!> linear_reg_exact_y_test test the linear regression using 
!> the exact values
subroutine linear_ref_exact_y_test()
  use mod_linear_reg, only: linear_regression
  implicit none
  !> variables
  real*8,dimension(N_coeff) :: lin_reg_coeff

  !> compute linear regression
  call linear_regression(N_points,x,y_exact,lin_reg_coeff)
  !> tests
  call assert_equals(line_coeff/line_coeff,lin_reg_coeff/line_coeff,N_coeff,tol,&
  "Error linear regression exact: regression and line coefficients do not match!")

end subroutine linear_ref_exact_y_test

!> linear_reg_noise_y_test test residual std error w.r.t. a reference value
subroutine linear_reg_noise_y_test()
  use mod_linear_reg, only: linear_regression
  implicit none
  !> variables
  real*8 :: residual_std_err_reg
  real*8,dimension(N_coeff) :: lin_reg_coeff
  real*8,dimension(N_points) :: y_reg

  !> compute the regression
  call linear_regression(N_points,x,y_noise,lin_reg_coeff)
  !> compute the residual standard error
  y_reg = lin_reg_coeff(2) + lin_reg_coeff(1)*x
  residual_std_err_reg = sum((y_noise-y_reg)*(y_noise-y_reg))
  residual_std_err_reg = sqrt(residual_std_err_reg/(N_points-2))
  !> check
  call assert_true(residual_std_err_reg.lt.tol_rse,&
  "Error linear regression noise: the regression residual std error overflows!")

end subroutine linear_reg_noise_y_test

end module mod_linear_reg_test
