!> mod_pinhole_lens_test contains all variable and procedures required
!> for testing the pinhole lens class.
module mod_pinhole_lens_test
  use fruit
  implicit none

  private
  public :: run_fruit_pinhole_lens

!> Variables --------------------------------------------------
character(len=19),parameter :: input_file='pinhole_lens_inputs'
integer,parameter :: read_unit=43
integer,parameter :: n_x=3
integer,parameter :: n_samples=12453
real*8,parameter  :: tol_r8=5.d-16
real*8,dimension(n_samples),parameter :: ones_r8=1.d0
real*8,dimension(n_x),parameter :: center_uppbnd=(/4.5d3,4.67d2,-2.34d-1/)
real*8,dimension(n_x),parameter :: center_lowbnd=(/-1.d2,2.45d1,-2.3d2/)
real*8,dimension(n_x)           :: center_sol
real*8,dimension(n_x,n_samples) :: x_sol

!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> basket running the set-up, test and teard-down procedures
subroutine run_fruit_pinhole_lens()
  implicit none
  write(*,'(/A)') "  ... setting-up: lens tests"
  call setup
  call write_pinhole_lens_inputs
  write(*,'(/A)') "  ... running: lens tests"
  call run_test_case(test_pinhole_lens_inputs,'test_pinhole_lens_inputs')
  call run_test_case(test_pinhole_lens_sampling,'test_pinhole_lens_sampling')
  call run_test_case(test_pinhole_initialisation,'test_pinhole_initialisation')
  call run_test_case(test_pinhole_lens_pdf,'test_pinhole_lens_pdf')
  write(*,'(/A)') "  ... tearing-down: lens tests"
  call teardown
end subroutine run_fruit_pinhole_lens

!> Set-up and teard-down---------------------------------------
!> set-up test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  integer :: ii
  !> initialize center array
  call gnu_rng_interval(n_x,center_lowbnd,center_uppbnd,center_sol)
  do ii=1,n_samples
    x_sol(:,ii) = center_sol
  enddo
end subroutine setup

!> write input file for testing
subroutine write_pinhole_lens_inputs()
  implicit none
  integer :: ifail
  open(read_unit,file=input_file,status='unknown',action='write',iostat=ifail)
  write(read_unit,'(/A)') '&pinhole_in'
  write(read_unit,'(/A,I10)') 'n_x = ',n_x
  write(read_unit,'(/A)') '/'
  write(read_unit,'(/A)') '&center_in'
  write(read_unit,'(/A,F28.16)') 'pinhole_center(1) = ', center_sol(1)
  write(read_unit,'(/A,F28.16)') 'pinhole_center(2) = ', center_sol(2)
  write(read_unit,'(/A,F28.16)') 'pinhole_center(3) = ', center_sol(3)
  write(read_unit,'(/A)') '/'
  close(read_unit)
end subroutine write_pinhole_lens_inputs

!> tear-down test features
subroutine teardown()
  implicit none
  center_sol = 0.d0; x_sol = 0.d0;
  call system("rm "//input_file)
end subroutine teardown

!> Tests ------------------------------------------------------
!> test the input file readre
subroutine test_pinhole_lens_inputs()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> variables
  type(pinhole_lens) :: pinhole
  integer :: rank,ifail
  integer,dimension(2) :: n_inputs
  integer,dimension(:),allocatable :: int_param
  real*8,dimension(:),allocatable  :: real_param
  !> initialisations
  rank = 0
  !> read data
  open(read_unit,file=input_file,status='old',action='read',iostat=ifail)
  call pinhole%read_lens_inputs(rank,read_unit,n_inputs,int_param,real_param)
  close(read_unit)
  !> checks
  call assert_equals((/1,3/),n_inputs,2,&
  "Error pinhole lens read inputs: N# inputs mismatch!")
  call assert_equals((/n_x/),int_param,1,&
  "Error pinhole lens read inputs: integer parameters mismatch!")
  call assert_equals(center_sol,real_param,n_x,tol_r8,&
  "Error pinhole lens read inputs: real parameters mismatch!")
  !> cleanup
  if(allocated(int_param)) deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine test_pinhole_lens_inputs

!> Test pinhole initialisation
subroutine test_pinhole_initialisation()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  type(pinhole_lens) :: pinhole
  call pinhole%init_pinhole(n_x,center_sol)
  call assert_equals(n_x,pinhole%n_x,&
  "Error pinhole lens initialisation: n_x mismatch!")
  call assert_true(allocated(pinhole%center),&
  "Error pinhole lens initialisation: center not allocated!")
  if(allocated(pinhole%center)) &
  call assert_equals(center_sol,pinhole%center,n_x,tol_r8,&
  "Error pinhole lens initialisation: center mismatch!")
  call pinhole%deallocate_lens
end subroutine test_pinhole_initialisation
 
!> Test pinhole sampling routine
subroutine test_pinhole_lens_sampling()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  type(pinhole_lens)  :: pinhole
  real*8,dimension(n_x,n_samples) :: x_loc
  call pinhole%init_pinhole(n_x,center_sol)
  call pinhole%sampling(n_samples,x_loc)
  call assert_equals(x_sol,x_loc,n_x,n_samples,tol_r8,&
  "Error pinhole lens sampling: positions mismatch!")
  call pinhole%deallocate_lens
end subroutine test_pinhole_lens_sampling

!> Test pinhole pdf routine
subroutine test_pinhole_lens_pdf()
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  type(pinhole_lens)  :: pinhole
  real*8,dimension(n_samples) :: pdf_loc
  call pinhole%init_pinhole(n_x,center_sol)
  call pinhole%pdf(n_samples,x_sol,pdf_loc)
  call assert_equals(ones_r8,pdf_loc,n_samples,tol_r8,&
  "Error pinhole lens sampling: positions mismatch!")
  call pinhole%deallocate_lens
end subroutine test_pinhole_lens_pdf

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_pinhole_lens_test
