! The mod_tiles_test module contains variables and procedures used
! for testing the good functionality of the mod_tiles modules
module mod_tiles_test
use fruit
implicit none

private
public :: run_fruit_tiles

! Module variables ------------------------------------------------
integer,parameter              :: N_rows=100 !< number of tile rows
integer,parameter              :: N_cols=100 !< number of tile columns
real*8,parameter               :: tol_real8=1.d-16 !< tolerance for error check
integer,dimension(2),parameter :: interval_int_1d=(/-1000,1000/)
real*8,dimension(2),parameter  :: interval_real8_1d=(/-1.d2,1.d2/)
! define interger/double 1d and 2d test arrays
integer :: N_rows_new,N_cols_new
integer :: N_rows_data,N_cols_data,N_rows_data_fail,N_cols_data_fail
integer :: offset_rows,offset_cols,offset_rows_fail,offset_cols_fail
integer,dimension(:),allocatable   :: data_int_1d
integer,dimension(:,:),allocatable :: data_int_2d
real*8,dimension(:),allocatable    :: data_real8_1d
real*8,dimension(:,:),allocatable  :: data_real8_2d

contains

! Tests basketes --------------------------------------------------

! run_fruit_tiles executes the tiles test set-up, tear-down
! and run the tests
subroutine run_fruit_tiles()
  implicit none
  
  ! execute setup -> tests -> teardown
  write(*,'(/A)') "  ... setting-up: tiles tests"
  call setup 
  write(*,'(/A)') "  ... running: tiles tests"
  call run_test_case(test_alloc_dealloc_noinit,'test_alloc_dealloc_noinit')
  call run_test_case(test_alloc_dealloc_init,'test_alloc_dealloc_init')
  call run_test_case(test_tile_resize,'test_tile_resize')
  call run_test_case(test_tile_resize_offsets,'test_tile_resize_offsets')
  call run_test_case(test_tile_resize_fail,'test_tile_resize_fail')
  call run_test_case(test_tile_resize_fail_offsets,'test_tile_resize_fail_offsets')
  write(*,'(/A)') "  ... tearing-down: tiles tests"
  call teardown

end subroutine run_fruit_tiles

! Set-up and tear-down --------------------------------------------

! setup initiliases the module variables
subroutine setup()
  use mod_dynamic_array_tools, only: allocate_check
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none

  ! initialise variables
  N_rows_new = floor(0.75*N_rows) 
  N_cols_new = floor(0.8*N_cols)
  ! set N_data variables
  N_rows_data      = floor(0.8*N_rows_new)
  N_rows_data_fail = N_rows_new + floor(0.2*N_rows_new)
  N_cols_data      = floor(0.65*N_cols_new)
  N_cols_data_fail = N_cols_new + floor(0.5*N_cols_new)
  ! set offset
  offset_rows      = floor(0.1*N_rows_new)
  offset_rows_fail = floor(0.3*N_rows_new)
  offset_cols      = floor(0.3*N_cols_new)
  offset_cols_fail = floor(0.45*N_cols_new)
  ! check for good defintions
  call assert_true((N_rows_new.gt.0).and.(N_rows_new.lt.N_rows),"Error in setting new number of rows")
  call assert_true((N_cols_new.gt.0).and.(N_cols_new.lt.N_cols),"Error in setting new number of columns")
  call assert_true(N_rows_data.gt.0,"Error N# data rows must be >0")  
  call assert_true(N_rows_data_fail.gt.0,"Error N# data rows fail must be >0")
  call assert_true(N_cols_data.gt.0,"Error N# data cols must be >0")
  call assert_true(N_cols_data_fail.gt.0,"Error N# data cols fail must be >0")
  call assert_true(offset_rows.gt.0,"Error offset rows must be >0")
  call assert_true(offset_rows_fail.gt.0,"Error offset rows fail must be >0")
  call assert_true(offset_cols.gt.0,"Error offset columns must be >0")
  call assert_true(offset_cols_fail.gt.0,"Error offset columns fail must be >0")
  ! just check the _fail and the offsets
  call assert_true(N_rows_data.le.N_rows_new,"Error N# data rows larger than N# rows new")
  call assert_true(N_rows_data_fail.gt.N_rows_new,"Error N# data rows fail smaller than N# rows new")
  call assert_true(N_cols_data.le.N_cols_new,"Error N# data columns larger than N# columns new")
  call assert_true(N_rows_data_fail.gt.N_rows_new,"Error N# data columns fails smaller than N# columns new")
  call assert_true((offset_rows+N_rows_data).le.N_rows_new,&
  "Error N# data rows + offset rows larger than N# rows new")
  call assert_true((offset_rows_fail+N_rows_data).gt.N_rows_new,&
  "Error N# data rows + offset rows fail smaller than N# rows new")
  call assert_true((offset_cols+N_cols_data).le.N_cols_new,&
  "Error N# data columns + offset columns larger than N# columns new")
  call assert_true((offset_cols_fail+N_cols_data).gt.N_cols_new,&
  "Error N# data columns + offset columns fail smaller than N# rows new")

  ! allocate integer test arrays
  call allocate_check(N_rows,data_int_1d)
  call allocate_check(N_rows,N_cols,data_int_2d)
  ! allocate double test arryas
  call allocate_check(N_rows,data_real8_1d)
  call allocate_check(N_rows,N_cols,data_real8_2d)

  ! generate integer random number within range
  call gnu_rng_interval(N_rows,interval_int_1d,data_int_1d)
  call gnu_rng_interval(N_rows,N_cols,interval_int_1d,data_int_2d)
  ! generate double random number within range
  call gnu_rng_interval(N_rows,interval_real8_1d,data_real8_1d)
  call gnu_rng_interval(N_rows,N_cols,interval_real8_1d,data_real8_2d)

end subroutine setup

! teardown cleans-up the module variables
subroutine teardown()
  use mod_dynamic_array_tools, only: deallocate_check
  implicit none

  ! deallocate integer test arrays
  call deallocate_check(data_int_1d)
  call deallocate_check(data_int_2d)
  ! deallocate double test arrays
  call deallocate_check(data_real8_1d)
  call deallocate_check(data_real8_2d)

end subroutine teardown

! Tests -----------------------------------------------------------

! test allocation and deallocation of every type of tiles
! initialization is set to zero
subroutine test_alloc_dealloc_noinit()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  integer :: ierr
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d
  integer,dimension(N_rows)        :: int_zero_array_1d
  integer,dimension(N_rows,N_cols) :: int_zero_array_2d
  real*8,dimension(N_rows)         :: real8_zero_array_1d 
  real*8,dimension(N_rows,N_cols)  :: real8_zero_array_2d

  ! init all arrays to zero
  ierr = 0
  int_zero_array_1d = 0
  int_zero_array_2d = 0
  real8_zero_array_1d = 0.d0
  real8_zero_array_2d = 0.d0

  ! check allocation and allocation for each tile type
  !> int_tile_2d
  call int_tile_1d%allocate_tile(N_rows,ierr)
  call assert_equals(N_rows,int_tile_1d%N_rows,&
  "Error: allocation of tile interger-1D N_rows does not match!")
  call assert_equals(int_zero_array_1d,int_tile_1d%data_array,N_rows,&
  "Error: allocation of tile interger-1D failed!")
  call int_tile_1d%deallocate_tile(ierr)
  call assert_equals(0,int_tile_1d%N_rows,&
  "Error: deallocation of tile interger-1D N_rows is not 0!")
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-1D failed!")
  !> int_tile_2d
  call int_tile_2d%allocate_tile(N_rows,N_cols,ierr)
  call assert_equals(N_rows,int_tile_2d%N_rows,&
  "Error: allocation of tile interger-2D N_rows does not match!")
  call assert_equals(N_cols,int_tile_2d%N_cols,&
  "Error: allocation of tile interger-1D N_cols does not match!")
  call assert_equals(int_zero_array_2d,int_tile_2d%data_array,N_rows,&
  N_cols,"Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_2d%deallocate_tile(ierr)
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-2D failed!")
  call assert_equals(0,int_tile_2d%N_rows,&
  "Error: deallocation of tile interger-2D N_rows is not 0!")
  call assert_equals(0,int_tile_2d%N_cols,&
  "Error: deallocation of tile interger-1D N_cols is not 0!")
  !> real8_tile_1d
  call real8_tile_1d%allocate_tile(N_rows,ierr)
  call assert_equals(N_rows,real8_tile_1d%N_rows,&
  "Error: allocation of tile double-1D N_rows does not match!")
  call assert_equals(real8_zero_array_1d,real8_tile_1d%data_array,&
  N_rows,tol_real8,"Error: allocation of tile double-1D failed!")
  call real8_tile_1d%deallocate_tile(ierr)
  call assert_equals(0,real8_tile_1d%N_rows,&
  "Error: deallocation of tile interger-1D N_rows is not 0!")
  call assert_false(allocated(real8_tile_1d%data_array),&
  "Error: deallocation of tile double-1D failed!")
  !> real8_tile_2d
  call real8_tile_2d%allocate_tile(N_rows,N_cols,ierr)
  call assert_equals(N_rows,real8_tile_2d%N_rows,&
  "Error: allocation of tile double-2D N_rows does not match!")
  call assert_equals(N_cols,real8_tile_2d%N_cols,&
  "Error: allocation of tile double-2D N_cols does not match!")
  call assert_equals(real8_zero_array_2d,real8_tile_2d%data_array,N_rows,N_cols,&
  tol_real8,"Error: allocation and double to data of tile double-1D failed!")
  call real8_tile_2d%deallocate_tile(ierr)
  call assert_false(allocated(real8_tile_2d%data_array),&
  "Error: deallocation of tile double-2D failed!")
  call assert_equals(0,real8_tile_2d%N_rows,&
  "Error: deallocation of tile double-2D N_rows is not 0!")
  call assert_equals(0,real8_tile_2d%N_cols,&
  "Error: deallocation of tile double-2D N_cols is not 0!")

  call assert_equals(0,ierr,"Error: tile allocation, an allocation is skipped")

end subroutine test_alloc_dealloc_noinit

!> test allocation and deallocation of tiles with value initialisation
subroutine test_alloc_dealloc_init()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  integer :: ierr
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d

  ! initialisation
  ierr = 0

  ! check allocation with data initialisation and deallocation
  !> int_tile_1d
  call int_tile_1d%allocate_tile(N_rows,ierr,data_int_1d)
  call assert_equals(N_rows,int_tile_1d%N_rows,&
  "Error: allocation and init to data of tile interger-1D N_rows does not match!")
  call assert_equals(data_int_1d,int_tile_1d%data_array,N_rows,&
  "Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_1d%deallocate_tile(ierr)
  call assert_equals(0,int_tile_1d%N_rows,&
  "Error: deallocation of tile interger-1D N_rows is not 0!")
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-1D failed!")
  !> int_tile_2d
  call int_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call assert_equals(N_rows,int_tile_2d%N_rows,&
  "Error: allocation and init to data of tile interger-2D N_rows does not match!")
  call assert_equals(N_cols,int_tile_2d%N_cols,&
  "Error: allocation and init to data of tile interger-1D N_cols does not match!")
  call assert_equals(data_int_2d,int_tile_2d%data_array,N_rows,&
  N_cols,"Error: allocation and init to data of tile interger-1D failed!")
  call int_tile_2d%deallocate_tile(ierr)
  call assert_false(allocated(int_tile_2d%data_array),&
  "Error: deallocation of tile integer-2D failed!")
  call assert_equals(0,int_tile_2d%N_rows,&
  "Error: deallocation of tile interger-2D N_rows is not 0!")
  call assert_equals(0,int_tile_2d%N_cols,&
  "Error: deallocation of tile interger-1D N_cols is not 0!")
  !> real8_tile_1d
  call real8_tile_1d%allocate_tile(N_rows,ierr,data_real8_1d)
  call assert_equals(N_rows,real8_tile_1d%N_rows,&
  "Error: allocation and init to data of tile double-1D N_rows does not match!")
  call assert_equals(data_real8_1d,real8_tile_1d%data_array,N_rows,&
  tol_real8,"Error: allocation and init to data of tile double-1D failed!")
  call real8_tile_1d%deallocate_tile(ierr)
  call assert_equals(0,real8_tile_1d%N_rows,&
  "Error: deallocation of tile interger-1D N_rows is not 0!")
  call assert_false(allocated(real8_tile_1d%data_array),&
  "Error: deallocation of tile double-1D failed!")
  !> real8_tile_2d
  call real8_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)
  call assert_equals(N_rows,real8_tile_2d%N_rows,&
  "Error: allocation and init to data of tile double-2D N_rows does not match!")
  call assert_equals(N_cols,real8_tile_2d%N_cols,&
  "Error: allocation and init to data of tile double-2D N_cols does not match!")
  call assert_equals(data_real8_2d,real8_tile_2d%data_array,N_rows,N_cols,&
  tol_real8,"Error: allocation and double to data of tile double-1D failed!")
  call real8_tile_2d%deallocate_tile(ierr)
  call assert_false(allocated(real8_tile_2d%data_array),&
  "Error: deallocation of tile double-2D failed!")
  call assert_equals(0,real8_tile_2d%N_rows,&
  "Error: deallocation of tile double-2D N_rows is not 0!")
  call assert_equals(0,real8_tile_2d%N_cols,&
  "Error: deallocation of tile double-2D N_cols is not 0!")

  call assert_equals(0,ierr,"Error: tile allocation, an allocation is skipped")

end subroutine test_alloc_dealloc_init

!> test the ability to resize the tile array
subroutine test_tile_resize()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d
  integer :: ierr
  integer,dimension(N_rows_new)            :: resize_data_int_1d
  integer,dimension(N_rows_new,N_cols_new) :: resize_data_int_2d
  real*8,dimension(N_rows_new)             :: resize_data_real8_1d
  real*8,dimension(N_rows_new,N_cols_new)  :: resize_data_real8_2d

  ! initialise resized data array
  ierr = 0
  !> resized integer 1D
  resize_data_int_1d = 0
  resize_data_int_1d(1:N_rows_data) = data_int_1d(1:N_rows_data)
  !> resized interger 2D
  resize_data_int_2d = 0
  resize_data_int_2d(1:N_rows_data,1:N_cols_data) = &
  data_int_2d(1:N_rows_data,1:N_cols_data)
  !> resized double 1D
  resize_data_real8_1d = 0.d0
  resize_data_real8_1d(1:N_rows_data) = data_real8_1d(1:N_rows_data)
  !> resized double 2D
  resize_data_real8_2d = 0.d0
  resize_data_real8_2d(1:N_rows_data,1:N_cols_data) = &
  data_real8_2d(1:N_rows_data,1:N_cols_data) 

  ! initialise tile array
  call int_tile_1d%allocate_tile(N_rows,ierr,data_int_1d)
  call int_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call real8_tile_1d%allocate_tile(N_rows,ierr,data_real8_1d)
  call real8_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)

  ! resize tiles
  call int_tile_1d%resize_tile(N_rows_new,N_rows_data,ierr)
  call int_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,ierr)
  call real8_tile_1d%resize_tile(N_rows_new,N_rows_data,ierr)
  call real8_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,ierr)

  ! check procedures
  call assert_equals(0,ierr,"Error: a resize tile routine internally failed")
  call assert_equals(N_rows_new,int_tile_1d%N_rows,&
  "Error: resize tile arre integer-1D tile N_rows does not match N_rows_new")
  call assert_equals(resize_data_int_1d,int_tile_1d%data_array,N_rows_new,&
  "Error: resize tile array interger-1D failed!")
  call assert_equals(N_rows_new,int_tile_2d%N_rows,&
  "Error: resize tile arre integer-2D tile N_rows does not match N_rows_new")
  call assert_equals(N_cols_new,int_tile_2d%N_cols,&
  "Error: resize tile arre integer-2D tile N_cols does not match N_cols_new")
  call assert_equals(resize_data_int_2d,int_tile_2d%data_array,N_rows_new,&
  N_cols_new,"Error: resize tile array interger-2D failed!")
  call assert_equals(N_rows_new,real8_tile_1d%N_rows,&
  "Error: resize tile arre double-1D tile N_rows does not match N_rows_new")
  call assert_equals(resize_data_real8_1d,real8_tile_1d%data_array,N_rows_new,&
  tol_real8,"Error: resize tile array double-1D failed!")
  call assert_equals(N_rows_new,real8_tile_2d%N_rows,&
  "Error: resize tile arre double-2D tile N_rows does not match N_rows_new")
  call assert_equals(N_cols_new,real8_tile_2d%N_cols,&
  "Error: resize tile arre double-2D tile N_cols does not match N_cols_new")
  call assert_equals(resize_data_real8_2d,real8_tile_2d%data_array,N_rows_new,&
  N_cols_new,tol_real8,"Error: resize tile array double-2D failed!")

  ! deallocate tiles
  call int_tile_1d%deallocate_tile(ierr)
  call int_tile_2d%deallocate_tile(ierr)
  call real8_tile_1d%deallocate_tile(ierr)
  call real8_tile_2d%deallocate_tile(ierr)
  
end subroutine test_tile_resize

!> test the ability to resize the tile array using the offsets
subroutine test_tile_resize_offsets()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d
  integer :: ierr
  integer,dimension(N_rows_new)            :: resize_data_int_1d
  integer,dimension(N_rows_new,N_cols_new) :: resize_data_int_2d
  real*8,dimension(N_rows_new)             :: resize_data_real8_1d
  real*8,dimension(N_rows_new,N_cols_new)  :: resize_data_real8_2d

  ! initialise resized data array
  ierr = 0
  !> resized integer 1D
  resize_data_int_1d = 0
  resize_data_int_1d(offset_rows+1:offset_rows+N_rows_data) = &
  data_int_1d(offset_rows+1:offset_rows+N_rows_data)
  !> resized interger 2D
  resize_data_int_2d = 0
  resize_data_int_2d(offset_rows+1:offset_rows+N_rows_data,offset_cols+1:offset_cols+N_cols_data) = &
  data_int_2d(offset_rows+1:offset_rows+N_rows_data,offset_cols+1:offset_cols+N_cols_data)
  !> resized double 1D
  resize_data_real8_1d = 0.d0
  resize_data_real8_1d(offset_rows+1:offset_rows+N_rows_data) = &
  data_real8_1d(offset_rows+1:offset_rows+N_rows_data)
  !> resized double 2D
  resize_data_real8_2d = 0.d0
  resize_data_real8_2d(offset_rows+1:offset_rows+N_rows_data,offset_cols+1:offset_cols+N_cols_data) = &
  data_real8_2d(offset_rows+1:offset_rows+N_rows_data,offset_cols+1:offset_cols+N_cols_data) 

  ! initialise tile array
  call int_tile_1d%allocate_tile(N_rows,ierr,data_int_1d)
  call int_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call real8_tile_1d%allocate_tile(N_rows,ierr,data_real8_1d)
  call real8_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)

  ! resize tiles
  call int_tile_1d%resize_tile(N_rows_new,N_rows_data,ierr,offset_rows)
  call int_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,ierr,&
  offset_rows,offset_cols)
  call real8_tile_1d%resize_tile(N_rows_new,N_rows_data,ierr,offset_rows)
  call real8_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,&
  ierr,offset_rows,offset_cols)

  ! check procedures
  call assert_equals(0,ierr,"Error: a resize tile routine internally failed")
  call assert_equals(N_rows_new,int_tile_1d%N_rows,&
  "Error: resize tile arre integer-1D tile N_rows does not match N_rows_new")
  call assert_equals(resize_data_int_1d,int_tile_1d%data_array,N_rows_new,&
  "Error: resize tile array interger-1D failed!")
  call assert_equals(N_rows_new,int_tile_2d%N_rows,&
  "Error: resize tile arre integer-2D tile N_rows does not match N_rows_new")
  call assert_equals(N_cols_new,int_tile_2d%N_cols,&
  "Error: resize tile arre integer-2D tile N_cols does not match N_cols_new")
  call assert_equals(resize_data_int_2d,int_tile_2d%data_array,N_rows_new,&
  N_cols_new,"Error: resize tile array interger-2D failed!")
  call assert_equals(N_rows_new,real8_tile_1d%N_rows,&
  "Error: resize tile arre double-1D tile N_rows does not match N_rows_new")
  call assert_equals(resize_data_real8_1d,real8_tile_1d%data_array,N_rows_new,&
  tol_real8,"Error: resize tile array double-1D failed!")
  call assert_equals(N_rows_new,real8_tile_2d%N_rows,&
  "Error: resize tile arre double-2D tile N_rows does not match N_rows_new")
  call assert_equals(N_cols_new,real8_tile_2d%N_cols,&
  "Error: resize tile arre double-1D tile N_cols does not match N_cols_new")
  call assert_equals(resize_data_real8_2d,real8_tile_2d%data_array,N_rows_new,&
  N_cols_new,tol_real8,"Error: resize tile array double-2D failed!")


  ! deallocate tiles
  call int_tile_1d%deallocate_tile(ierr)
  call int_tile_2d%deallocate_tile(ierr)
  call real8_tile_1d%deallocate_tile(ierr)
  call real8_tile_2d%deallocate_tile(ierr)
  
end subroutine test_tile_resize_offsets

!> test the fail of tile_resize if data size larger new tile size
subroutine test_tile_resize_fail()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d,int_tile_2d_2
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d,real8_tile_2d_2
  integer :: ierr

  ! initialise tile array
  call int_tile_1d%allocate_tile(N_rows,ierr,data_int_1d)
  call int_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call int_tile_2d_2%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call real8_tile_1d%allocate_tile(N_rows,ierr,data_real8_1d)
  call real8_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)
  call real8_tile_2d_2%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)

  ! resize tiles and check for failures
  ierr = 0
  call int_tile_1d%resize_tile(N_rows_new,N_rows_data_fail,ierr)
  call assert_equals(1,ierr,"Error: resize tile array integer-1D no failure for data overflow")
  ierr = 0
  call int_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data_fail,N_cols_data,ierr)
  call assert_equals(1,ierr,"Error: resize tile array integer-2D no failure for rows data overflow")
  ierr = 0
  call int_tile_2d_2%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data_fail,ierr)
  call assert_equals(1,ierr,"Error: resize tile array integer-2D no failure for columns data overflow")
  ierr = 0
  call real8_tile_1d%resize_tile(N_rows_new,N_rows_data_fail,ierr)
  call assert_equals(1,ierr,"Error: resize tile array double-1D no failure for data overflow")
  ierr = 0
  call real8_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data_fail,N_cols_data,ierr)
  call assert_equals(1,ierr,"Error: resize tile array double-2D no failure for row data overflow")
  ierr = 0
  call real8_tile_2d_2%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data_fail,ierr)
  call assert_equals(1,ierr,"Error: resize tile array double-2D no failure for column data overflow")

  ! deallocate tiles
  call int_tile_1d%deallocate_tile(ierr)
  call int_tile_2d%deallocate_tile(ierr)
  call int_tile_2d_2%deallocate_tile(ierr)
  call real8_tile_1d%deallocate_tile(ierr)
  call real8_tile_2d%deallocate_tile(ierr)
  call real8_tile_2d_2%deallocate_tile(ierr)
  
end subroutine test_tile_resize_fail

!> test the resize tile safeguards in case of offset overflow
subroutine test_tile_resize_fail_offsets()
  use mod_tiles,only: tile_int_1d,tile_int_2d
  use mod_tiles,only: tile_real8_1d,tile_real8_2d
  implicit none

  ! variables
  integer :: ierr
  type(tile_int_1d)   :: int_tile_1d
  type(tile_int_2d)   :: int_tile_2d,int_tile_2d_2
  type(tile_real8_1d) :: real8_tile_1d
  type(tile_real8_2d) :: real8_tile_2d,real8_tile_2d_2

  ! initialise tile array
  call int_tile_1d%allocate_tile(N_rows,ierr,data_int_1d)
  call int_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call int_tile_2d_2%allocate_tile(N_rows,N_cols,ierr,data_int_2d)
  call real8_tile_1d%allocate_tile(N_rows,ierr,data_real8_1d)
  call real8_tile_2d%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)
  call real8_tile_2d_2%allocate_tile(N_rows,N_cols,ierr,data_real8_2d)


  ! resize tiles and check for failures
  ierr = 0
  call int_tile_1d%resize_tile(N_rows_new,N_rows_data,ierr,offset_rows_fail)
  call assert_equals(1,ierr,"Error: resize tile array integer-1D no failure for offset overflow")
  ierr = 0
  call int_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,ierr,&
  offset_rows_fail,offset_cols)
  call assert_equals(1,ierr,"Error: resize tile array integer-2D no failure for row offset overflow")
  ierr = 0
  call int_tile_2d_2%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,ierr,&
  offset_rows,offset_cols_fail)
  call assert_equals(1,ierr,"Error: resize tile array integer-2D no failure for column offset overflow")
  ierr = 0
  call real8_tile_1d%resize_tile(N_rows_new,N_rows_data,ierr,offset_rows_fail)
  call assert_equals(1,ierr,"Error: resize tile array double-1D no failure for offset overflow")
  ierr = 0
  call real8_tile_2d%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,&
  ierr,offset_rows_fail,offset_cols)
  call assert_equals(1,ierr,"Error: resize tile array double-2D no failure for row offset overflow")
  ierr = 0
  call real8_tile_2d_2%resize_tile(N_rows_new,N_cols_new,N_rows_data,N_cols_data,&
  ierr,offset_rows,offset_cols_fail)
  call assert_equals(1,ierr,"Error: resize tile array double-2D no failure for column offset overflow")

  ! deallocate tiles
  call int_tile_1d%deallocate_tile(ierr)
  call int_tile_2d%deallocate_tile(ierr)
  call int_tile_2d_2%deallocate_tile(ierr)
  call real8_tile_1d%deallocate_tile(ierr)
  call real8_tile_2d%deallocate_tile(ierr)
  call real8_tile_2d_2%deallocate_tile(ierr)
 
end subroutine test_tile_resize_fail_offsets

end module mod_tiles_test
