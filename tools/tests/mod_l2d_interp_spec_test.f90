!> This module contains some testcases for the linear 2d interpolation
module mod_l2d_interp_spec_test
use fruit
use mod_interp_splinear
implicit none
private
public :: run_fruit_l2d_interp_spec
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_l2d_interp_spec
  implicit none
  write(*,'(/A)') "  ... setting-up: l2d interp spec"
  write(*,'(/A)') "  ... running: l2d interp spec"
  call run_test_case(test_l2d_interp_linear,'test_l2d_interp_linear')
  write(*,'(/A)') "  ... tearing-down: l2d interp spec"
end subroutine run_fruit_l2d_interp_spec

!> Tests ------------------------------------------
subroutine test_l2d_interp_linear
  integer, parameter :: nx = 2, ny = 3
  integer :: ix, iy
  real*8, parameter :: ax = 2.5, ay = 7.5 !< Slopes in x and y
  real*8 :: tx(nx)
  real*8 :: ty(ny)
  real*8 :: f(nx,ny)
  real*8 :: x, y
  real*8, parameter :: tolerance = 1.d-30

  ! Prepare function and nodes
  do ix=1,nx
    tx(ix) = real(ix,8)
    do iy=1,ny
      ty(iy) = real(iy,8)
      f(ix,iy) = tx(ix)*ax + ty(iy)*ay
    enddo
  enddo

  ! Test a few cases
  x = 1.5d0; y = 1.5d0
  call assert_equals(L2Dinterp(tx,ty,f,x,y), x*ax + y*ay, tolerance, "In domain")

  ! Outside the domain
  x = 0.0d0; y = 0.0d0
  ! extrapolation is no longer implemented
  call assert_equals(L2Dinterp(tx,ty,f,x,y), min(max(x,minval(tx)),maxval(tx))*ax +&
  + min(max(x,minval(ty)),maxval(ty))*ay, tolerance, "Outside domain")

  ! On one of the nodes
  x = 1.5d0; y = 2.0d0
  call assert_equals(L2Dinterp(tx,ty,f,x,y), x*ax + y*ay, tolerance, "On one node")

  ! On both of the nodes
  x = 2.d0; y = 2.d0
  call assert_equals(L2Dinterp(tx,ty,f,x,y), x*ax + y*ay, tolerance, "On both nodes")
end subroutine test_l2d_interp_linear
!> ------------------------------------------------
end module mod_l2d_interp_spec_test 
