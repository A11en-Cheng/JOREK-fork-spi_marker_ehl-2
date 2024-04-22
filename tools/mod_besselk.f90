! The module boost_besselk contains all interfaces
! and procedures for importing and bining the 
! boost function cyl_bessel_k computing
! the modified bessel function of the second kind
! with fractional order.
! WARNING: cyl_bessel_k(nu,x) accept only scalal
!          inputs
module mod_besselk
implicit none

private
public ::  f_besselk,besselk
#ifdef UNIT_TESTS
public :: handle_float_exceptions
#endif

! binding interface to the cyl_bessel_k
interface 
  function besselk_cpp(nu,x) bind(C,name="besselk_cpp")
    use iso_c_binding, only: c_double
    implicit none
    real(c_double) :: besselk_cpp
    real(c_double) :: nu,x
  end function besselk_cpp
end interface

! interface for the different procedure in the module
interface f_besselk    
  module procedure f_besselk_single_cpp
end interface f_besselk

interface besselk
    module procedure besselk_single_cpp 
    module procedure besselk_x_array_cpp
    module procedure besselk_nu_array_cpp
    module procedure besselk_x_nu_array_cpp
end interface besselk

contains

! besselk is a specialization of the cyl_bessel_k
! function to double datatypes
! inputs:
!   nu: (real8) bessel function fractional order
!   x:  (real8) value at which the bessel function
!       is computed
! outputs:
!   bknu: (real8) value of the bessel functions
function f_besselk_single_cpp(nu,x) result(bknu)
  implicit none 
  real*8, intent(in) :: nu,x
  ! outputs
  real*8 :: bknu
  bknu = handle_float_exceptions(besselk_cpp(nu,x))
end function f_besselk_single_cpp

! besselk is a specialization of the cyl_bessel_k
! function to double datatypes
! inputs:
!   nu: (real8) bessel function fractional order
!   x:  (real8) value at which the bessel function
!       is computed
! outputs:
!   bknu: (real8) value of the bessel functions
subroutine besselk_single_cpp(nu,x,bknu)
  implicit none 
  real*8, intent(in) :: nu,x
  ! outputs
  real*8,intent(out) :: bknu
  bknu = handle_float_exceptions(besselk_cpp(nu,x))
end subroutine besselk_single_cpp

! besselk_x_array computes the modified bessel
! function of the second kind for an array of x
! inputs:
!   Nx: (integer) number of x values
!   nu: (real8) bessel function fractional order
!   x:  (real8)(Nx) array of x values
! outputs:
!   bknu: (real8)(Nx) array of bessel functions
subroutine besselk_x_array_cpp(Nx,nu,x,bknu)
  implicit none
  ! inputs
  integer,intent(in)              :: Nx
  real*8,intent(in)               :: nu
  real*8,dimension(Nx),intent(in) :: x
  ! outputs
  real*8,dimension(Nx),intent(out) :: bknu
  ! variables
  integer :: ii

  do ii=1,Nx
    bknu(ii) = handle_float_exceptions(besselk_cpp(nu,x(ii)))
  enddo
end subroutine besselk_x_array_cpp

! besselk_nu_array computes the modified bessel
! function of the second kind for an array of nu
! inputs:
!   Nnu: (integer) number of fractional orders
!   nu:  (real)(Nu) array of fractional orders
!   x:   (real8) value at which the bessel function
!        is computed
! outputs:
!   bknu: (real8)(Nu) array of bessel functions
subroutine besselk_nu_array_cpp(Nnu,nu,x,bknu)
  implicit none
  ! inputs
  integer,intent(in)               :: Nnu
  real*8,dimension(Nnu),intent(in) :: nu
  real*8,intent(in)                :: x
  !outputs
  real*8,dimension(Nnu),intent(out) :: bknu
  ! variables
  integer :: ii

  do ii=1,Nnu
    bknu(ii) = handle_float_exceptions(besselk_cpp(nu(ii),x))
  enddo
end subroutine besselk_nu_array_cpp

! besselk_x_nu_array computes the modified bessel
! function of the second kind for an array of x
! and nu being the x the first index
! inputs:
!   Nx:  (integer) number of x values
!   Nnu: (integer) number of fractional orders
!   nu:  (real)(Nu) array of fractional orders
!   x:   (real8)(Nx) array of x values
! outputs:
!   bknu: (real8)(Nx,Nu) array of bessel functions
subroutine besselk_x_nu_array_cpp(Nx,Nnu,nu,x,bknu)
  implicit none
  ! inputs
  integer,intent(in)               :: Nx,Nnu
  real*8,dimension(Nnu),intent(in) :: nu
  real*8,dimension(Nx),intent(in)  :: x
  ! outputs
  real*8,dimension(Nx,Nnu),intent(out) :: bknu
  ! variables
  integer :: ii,jj
  do jj=1,Nnu
    do ii=1,Nx
      bknu(ii,jj) = handle_float_exceptions(besselk_cpp(nu(jj),x(ii)))
    enddo
  enddo
end subroutine besselk_x_nu_array_cpp

! method used for handling bessel function exceptions
! inputs:
!   val:  (real8) value to be checked
! outpus:
!   val: (real8) val if tests pass 0 otherwise
function handle_float_exceptions(val_in) result(val)
  implicit none
  !> inputs
  real*8,intent(in) :: val_in
  !> outputs
  real*8 :: val
  !> checks
  val = val_in
  if(isnan(val).or.(abs(val).gt.huge(0d0)).or.(abs(val).lt.tiny(0d0))) then
    val = 0d0; return;
  endif
end function handle_float_exceptions

end module mod_besselk
