!> mod_linear_reg implement a basic examples of linear regression method
module mod_linear_reg
  implicit none
  private
  public :: linear_regression

!> Interfaces ---------------------------------------------------
  interface linear_regression
    module procedure linear_regression_1d
  end interface
contains

!> Procedures ---------------------------------------------------
!> compute the linear regression for for 1D data
!> inputs:
!>   N: (integer) number of points
!>   x: (double)(N) abscissas of the measurements
!>   y: (double)(N) measurements
!> outputs:
!>   coeff: (real8)(2) coefficients of the linear regression
!>          1: m slope of the regression line
!>          2: b: y-intercept of the regression line
subroutine linear_regression_1d(N,x,y,coeff)
  implicit none
  !> inputs
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x,y
  !> outputs
  real*8,dimension(2),intent(out) :: coeff
  !> variables
  integer :: ii
  real*8 :: sum_x,sum_y,sum_x2,sum_xy

  !> compute averages
  sum_x = sum(x); sum_x2 = sum(x*x);
  sum_xy = sum(x*y); sum_y = sum(y); 

  !> compute coefficients
  coeff = (/N*sum_xy-sum_x*sum_y,&
  sum_y*sum_x2-sum_x*sum_xy/)
  coeff = coeff/(N*sum_x2 - sum_x*sum_x)
end subroutine linear_regression_1d

end module mod_linear_reg
