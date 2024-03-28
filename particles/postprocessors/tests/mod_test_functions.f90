!> mod_test_functions contains some functions and their
!> analytical integrals used for testing
module mod_test_functions
implicit none

private
public :: sin2x,int_sin2x
public :: expxsin2x,int_expxsin2x

!> Interfaces ----------------------------------------------------------------
interface sin2x
  module procedure sin2x_serial,sin2x_vector
end interface

interface int_sin2x
  module procedure int_sin2x_serial,int_sin2x_vector
end interface

interface expxsin2x
  module procedure expxsin2x_serial,expxsin2x_vector
end interface

interface int_expxsin2x
  module procedure int_expxsin2x_serial,int_expxsin2x_vector
end interface

contains

!> sin(x)^2 ------------------------------------------------------------------

!> sin^2(x) function
function sin2x_serial(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: sin2x_serial
  sin2x_serial = sin(x)*sin(x)
end function sin2x_serial

!> integral of the sin^2(x) function
function int_sin2x_serial(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: int_sin2x_serial
  int_sin2x_serial = 5.d-1*(x-sin(x)*cos(x))
end function int_sin2x_serial

!> sin^2(x) function
function sin2x_vector(N,x)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x
  real*8,dimension(N) :: sin2x_vector
  sin2x_vector = sin(x)*sin(x)
end function sin2x_vector

!> integral of the sin^2(x) function
function int_sin2x_vector(N,x)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x
  real*8,dimension(N) :: int_sin2x_vector
  int_sin2x_vector = 5.d-1*(x-sin(x)*cos(x))
end function int_sin2x_vector

!> exp(x)*sin(x)^2
function expxsin2x_serial(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: expxsin2x_serial
  expxsin2x_serial = exp(x)*sin(x)*sin(x)
end function expxsin2x_serial

!> integral of the exp(x)*sin(x)^2 function
function int_expxsin2x_serial(x)
  implicit none
  real*8,intent(in) :: x
  real*8 :: int_expxsin2x_serial
  int_expxsin2x_serial = -1.d-1*exp(x)*(2.d0*sin(2.d0*x)+cos(2.d0*x)-5.d0)
end function int_expxsin2x_serial

!> exp(x)*sin(x)^2 function
function expxsin2x_vector(N,x)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x
  real*8,dimension(N) :: expxsin2x_vector
  expxsin2x_vector = exp(x)*sin(x)*sin(x)
end function expxsin2x_vector

!> integral of the exp(x)*sin(x)^2 function
function int_expxsin2x_vector(N,x)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(N),intent(in) :: x
  real*8,dimension(N) :: int_expxsin2x_vector
  int_expxsin2x_vector = -1.d-1*exp(x)*(2.d0*sin(2.d0*x)+cos(2.d0*x)-5.d0)
end function int_expxsin2x_vector

!> ---------------------------------------------------------------------------

end module mod_test_functions
