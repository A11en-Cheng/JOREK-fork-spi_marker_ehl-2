module mod_hermite_birkhoff_test
use fruit
implicit none
private
public :: run_fruit_hermite_birkhoff
!> Fruit basket -----------------------------------
integer,parameter :: n_tests=10
real*8,parameter  :: t0=0d0
real*8,parameter  :: t1=1d0
real*8,parameter  :: ftol=3d-3
real*8,parameter  :: dftol=1d-2
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_hermite_birkhoff
  implicit none
  write(*,'(/A)') "  ... setting-up:hermite birkhoff "
  write(*,'(/A)') "  ... running: hermite birkhoff"
  call run_test_case(test_interp_hermite_birkhoff_cos,&
  'test_interp_hermite_birkhoff_cos')
  write(*,'(/A)') "  ... tearing-down: hermite birkhoff"
end subroutine run_fruit_hermite_birkhoff

!> Tests ------------------------------------------
subroutine test_interp_hermite_birkhoff_cos
  use mod_hermite_birkhoff, only: HB_interp,HB_interp_dt
  implicit none
  real*8, dimension(1) :: y0, y1, y, dy0, dy1
  integer :: ii
  real*8 :: t
  y0 = cos(t0); y1 = cos(t1); dy0 = -sin(t0); dy1 = -sin(t1);

  do ii=1,n_tests
    t = real(ii-1,kind=8)/real(n_tests-1,kind=8) ! include endpoints
    call HB_interp(t0, t1, 1, y0, y1, dy0, dy1, t, y)
    call assert_equals(cos(t), y(1), ftol, &
    "Error hermite-birkhoff cosinus interp: value mismatch!")
    call HB_interp_dt(t0, t1, 1, y0, y1, dy0, dy1, t, y)
    call assert_equals(-sin(t), y(1), dftol, &
    "Error hermite-birkhoff cosinus interp: derivative value mismatch!")
  end do
end subroutine test_interp_hermite_birkhoff_cos
!> ------------------------------------------------
end module  mod_hermite_birkhoff_test
