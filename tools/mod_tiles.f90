! The module mod_tiles implement the tile types required by the
! mod_tile_arrays data structure. A tiled array is basically 
! an array of arrays, called tiles, such that each tile can have 
! an independet size and indexing. The idea is to split an input array 
! in multiple tiles and resize only a small amount of tiles when needed. 
! In addition, it is expected that tiled_array would allow some degree 
! of manual memory management if needed.
! WARNING: the mod_tiled_array module tries to be as less polymorphic 
!          as possible for avoiding incompatibilities with OpenMP and
!          offloading
module mod_tiles
implicit none
private 
public :: tile_int_1d, tile_int_2d
public :: tile_real8_1d, tile_real8_2d

! Types definitions -----------------------------------------------

! definition of tile types integer_1D, integer_2D, 
! double_1D and double_2D tile types
! WARNING: data in the tile must be stored with increasing index
! Attributes are:
!   N_rows:     (integer) number of rows of the data_array
!   N_cols:     (integer) number fo columns of the data array
!   data_array: (integer/real8)(N_rows/N_rows*N_columns) tile array
!               in which data are stored 
type,abstract :: tile_base
  integer :: N_rows
end type tile_base
!> define 1D intger tile array
type,extends(tile_base) :: tile_int_1d
  integer,dimension(:),allocatable :: data_array
contains
  !procedure(alloc_int_tile_1d),deferred :: allocate_tile => allocate_int_tile_1d
  procedure,pass(tile) :: allocate_tile   => allocate_int_tile_1d
  procedure,pass(tile) :: deallocate_tile => deallocate_int_tile_1d
  procedure,pass(tile) :: resize_tile     => resize_int_tile_1d 
end type tile_int_1d

!> define 2D integer tile array
type,extends(tile_base) :: tile_int_2d
  integer :: N_cols
  integer,dimension(:,:),allocatable :: data_array
contains
  procedure,pass(tile) :: allocate_tile   => allocate_int_tile_2d
  procedure,pass(tile) :: deallocate_tile => deallocate_int_tile_2d
  procedure,pass(tile) :: resize_tile     => resize_int_tile_2d
end type tile_int_2d

!> define 1D double array
type,extends(tile_base) :: tile_real8_1d 
  real*8,dimension(:),allocatable :: data_array
contains
  procedure,pass(tile) :: allocate_tile   => allocate_real8_tile_1d
  procedure,pass(tile) :: deallocate_tile => deallocate_real8_tile_1d
  procedure,pass(tile) :: resize_tile     => resize_real8_tile_1d
end type tile_real8_1d

!> define 2D double array
type,extends(tile_base) :: tile_real8_2d
  integer :: N_cols
  real*8,dimension(:,:),allocatable :: data_array
contains
  procedure,pass(tile) :: allocate_tile   => allocate_real8_tile_2d
  procedure,pass(tile) :: deallocate_tile => deallocate_real8_tile_2d
  procedure,pass(tile) :: resize_tile     => resize_real8_tile_2d
end type tile_real8_2d

contains

! Allocate tile ---------------------------------------------------

! The allocate_int_tile_1d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_int_1d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   ierr:   (integer) if 1 and error occurred
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile: (tile_int_1d) the allocated tile
!   ierr: (integer) if 1 and error occurred
subroutine allocate_int_tile_1d(tile,N_rows,ierr,dat)
  implicit none
  ! inputs-outputs
  integer,intent(inout)            :: ierr
  class(tile_int_1d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows
  integer,dimension(N_rows),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) then
     tile%N_rows = N_rows; allocate(tile%data_array(N_rows));
  else
     ierr = 1
     return
  endif
  if(present(dat)) then
    tile%data_array = dat
  else
    tile%data_array = 0
  endif
end subroutine allocate_int_tile_1d

! The allocate_int_tile_2d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_int_2d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   N_cols: (integer) number of tile columns
!   ierr:   (integer) if 1 and error occurred
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile: (tile_int_2d) the allocated tile
!   ierr: (integer) if 1 and error occurred
subroutine allocate_int_tile_2d(tile,N_rows,N_cols,ierr,dat)
!  implicit none
  ! inputs-outputs
  integer,intent(inout)            :: ierr 
  class(tile_int_2d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows,N_cols
  integer,dimension(N_rows,N_cols),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) then
    tile%N_rows = N_rows; tile%N_cols = N_cols
    allocate(tile%data_array(N_rows,N_cols))
  else
    ierr = 1
    return
  endif
  if(present(dat)) then
    tile%data_array = dat
  else
    tile%data_array = 0
  endif

end subroutine allocate_int_tile_2d

! The allocate_real8_tile_1d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_real8_1d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   ierr:   (integer) if 1 and error occurred
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile: (tile_real8_1d) the allocated tile
!   ierr: (integer) if 1 and error occurred
subroutine allocate_real8_tile_1d(tile,N_rows,ierr,dat)
  implicit none
  ! inputs-outputs
  integer,intent(inout)              :: ierr
  class(tile_real8_1d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows
  real*8,dimension(N_rows),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) then
    tile%N_rows = N_rows; allocate(tile%data_array(N_rows))
  else
    ierr = 1
    return
  endif
  if(present(dat)) then
    tile%data_array = dat
  else
    tile%data_array = 0.d0
  endif
end subroutine allocate_real8_tile_1d

! The allocate_real8_tile_2d procedure allocate an integer 1d tile
! If present, tile array is initilized to dat, to 0 otherwise
! inputs:
!   tile:   (tile_real8_2d) the tile to be allocated
!   N_rows: (integer) number of tile rows
!   N_cols: (integer) number of tile columns
!   ierr:   (integer) if 1 and error occurred
!   dat:    (double)(N_rows,N_cols)(optional) data array for initialisation
! outputs:
!   tile: (tile_real8_2d) the allocated tile
!   ierr: (integer) if 1 and error occurred
subroutine allocate_real8_tile_2d(tile,N_rows,N_cols,ierr,dat)
  implicit none
  ! inputs-outputs
  integer,intent(inout)              :: ierr
  class(tile_real8_2d),intent(inout) :: tile
  ! inputs
  integer,intent(in) :: N_rows,N_cols
  real*8,dimension(N_rows,N_cols),intent(in),optional :: dat

  if(.not.allocated(tile%data_array)) then
    tile%N_rows = N_rows; tile%N_cols = N_cols;
    allocate(tile%data_array(N_rows,N_cols))
  else
    ierr = 1
    return
  endif
  if(present(dat)) then
     tile%data_array = dat
  else
     tile%data_array = 0.d0
  endif

end subroutine allocate_real8_tile_2d

! Deallocate tile -------------------------------------------------

! Deallocate tile_int_1d data array
! inputs:
!   tile: (tile_int_1d) the allocated tile
!   ierr: (integer) if 1 and error occurred
! outputs:
!   tile: (tile_int_1d) the deallocated tile
!   ierr: (integer) if 1 and error occurred
subroutine deallocate_int_tile_1d(tile,ierr)
  implicit none
  ! inputs-outputs
  class(tile_int_1d),intent(inout) :: tile
  integer,intent(inout)            :: ierr
  
  if(allocated(tile%data_array)) then
    tile%N_rows = 0; deallocate(tile%data_array);
  else
    ierr = 1
  endif

end subroutine deallocate_int_tile_1d

! Deallocate tile_int_2d data array
! inputs:
!   tile: (tile_int_2d) the allocated tile
!   ierr: (integer) if 1 and error occurred
! outputs:
!   tile: (tile_int_2d) the deallocated tile
!   ierr: (integer) if 1 and error occurred
subroutine deallocate_int_tile_2d(tile,ierr)
  implicit none
  ! inputs-outputs
  class(tile_int_2d),intent(inout) :: tile
  integer,intent(inout)            :: ierr
  
  if(allocated(tile%data_array)) then
    tile%N_rows = 0; tile%N_cols = 0;
    deallocate(tile%data_array)
  else
    ierr = 1
  endif

end subroutine deallocate_int_tile_2d

! Deallocate tile_int_1d data array
! inputs:
!   tile: (tile_real8_1d) the allocated tile
!   ierr: (integer) if 1 and error occurred
! outputs:
!   tile: (tile_real8_1d) the deallocated tile
!   ierr: (integer) if 1 and error occurred
subroutine deallocate_real8_tile_1d(tile,ierr)
  implicit none
  ! inputs-outputs
  class(tile_real8_1d),intent(inout) :: tile
  integer,intent(inout)              :: ierr
  
  if(allocated(tile%data_array)) then
    tile%N_rows = 0;  deallocate(tile%data_array);
  else
    ierr = 1
  endif

end subroutine deallocate_real8_tile_1d

! Deallocate tile_int_1d data array
! inputs:
!   tile: (tile_real8_2d) the allocated tile
!   ierr: (integer) if 1 and error occurred
! outputs:
!   tile: (tile_real8_2d) the deallocated tile
!   ierr: (integer) if 1 and error occurred
subroutine deallocate_real8_tile_2d(tile,ierr)
  implicit none
  ! inputs-outputs
  class(tile_real8_2d),intent(inout) :: tile
  integer,intent(inout)              :: ierr
  
  if(allocated(tile%data_array)) then
    tile%N_rows = 0; tile%N_cols = 0;
    deallocate(tile%data_array)
  else
    ierr = 1
  endif

end subroutine deallocate_real8_tile_2d

! Resize tile -----------------------------------------------------

! resize_int_tile_1d resizes a integer 1d data_array of a tile
! preserving the data within the offset and offset+N_data
! data are stored in offset:offset+N_data
! inputs:
!   tile:       (tile_int_1d) the tile to resize
!   N_rows_new: (integer) desired size tile
!   N_data:     (integer) number of stored data
!   ierr:       (integer) if 1 and error occurred
!   offset_in:  (integer)(optional) initial index for storing data
! outputs:
!   tile: (tile_int_1d) the resized tile
!   ierr: (integer) if 1 and error occurred
subroutine resize_int_tile_1d(tile,N_rows_new,N_data,ierr,offset_in)
  implicit none
  ! inputs-outputs:
  class(tile_int_1d),intent(inout) :: tile
  integer,intent(inout)            :: ierr
  ! inputs:
  integer,intent(in)          :: N_rows_new,N_data
  integer,intent(in),optional :: offset_in
  ! variables
  integer                   :: offset
  integer,dimension(N_data) :: tmp_data !< temporary data array

  ! initialization
  if(present(offset_in)) then
    offset = offset_in
  else
    offset = 0
  endif
  if((offset+N_data).gt.tile%N_rows) offset = max(0,tile%N_rows-N_data)
  if((offset+N_data).gt.N_rows_new) then
    ierr = 1
    return
  endif

  ! resizing
  tmp_data = tile%data_array(offset+1:offset+N_data)
  call tile%deallocate_tile(ierr)
  call tile%allocate_tile(N_rows_new,ierr)
  tile%data_array(offset+1:offset+N_data) = tmp_data
  
end subroutine resize_int_tile_1d

! resize_int_tile_2d resizes a integer 2d data_array of a tile
! preserving the data vwithin the offset and offset+N_data
! data are stored in offset:offset+N_data
! inputs:
!   tile:           (tile_int_2d) the tile to resize
!   N_rows_new:     (integer) desired number of rows
!   N_cols_new:     (integer) desired number of columns
!   N_data_rows:    (integer) number of stored data rows
!   N_data_cols:    (integer) number of stored data columns
!   ierr:           (integer) if 1 and error occurred
!   offset_rows_in: (integer)(optional) initial row index for storing data
!   offset_cols_in: (integer)(optional) initial columns index for storing data
! outputs:
!   tile: (tile_int_2d) the resized tile
!   ierr: (integer) if 1 and error occurred
subroutine resize_int_tile_2d(tile,N_rows_new,N_cols_new,N_data_rows,N_data_cols,&
ierr,offset_rows_in,offset_cols_in)
  implicit none
  ! inputs-outputs:
  class(tile_int_2d),intent(inout) :: tile
  integer,intent(inout)            :: ierr
  ! inputs:
  integer,intent(in) :: N_rows_new,N_cols_new
  integer,intent(in) :: N_data_rows,N_data_cols
  integer,intent(in),optional :: offset_rows_in,offset_cols_in
  ! variables
  integer :: offset_rows,offset_cols
  integer,dimension(N_data_rows,N_data_cols) :: tmp_data !< temporary data array

  ! initialization
  if(present(offset_rows_in)) then
    offset_rows = offset_rows_in
  else
    offset_rows = 0
  endif
  if(present(offset_cols_in)) then
    offset_cols = offset_cols_in
  else
    offset_cols = 0
  endif
  if((offset_rows+N_data_rows).gt.tile%N_rows) offset_rows = max(0,tile%N_rows-N_data_rows)
  if((offset_cols+N_data_cols).gt.tile%N_cols) offset_cols = max(0,tile%N_cols-N_data_cols)
  if((offset_rows+N_data_rows).gt.N_rows_new) then
    ierr = 1
    return
  endif
  if((offset_cols+N_data_cols).gt.N_cols_new) then
    ierr = 1
    return
  endif

  ! resizing
  tmp_data = &
  tile%data_array(offset_rows+1:offset_rows+N_data_rows,offset_cols+1:offset_cols+N_data_cols)
  call tile%deallocate_tile(ierr)
  call tile%allocate_tile(N_rows_new,N_cols_new,ierr)
  tile%data_array(offset_rows+1:offset_rows+N_data_rows,offset_cols+1:offset_cols+N_data_cols) = &
  tmp_data
  
end subroutine resize_int_tile_2d

! resize_real8_1d resizes a double 1d data_array of a tile
! preserving the data within the offset and offset+N_data
! data are stored in offset:offset+N_data
! inputs:
!   tile:       (tile_real8_1d) the tile to resize
!   N_rows_new: (integer) desired size tile
!   N_data:     (integer) number of stored data
!   ierr:       (integer) if 1 and error occurred
!   offset_in:  (integer)(optional) initial index for storing data
! outputs:
!   tile: (tile_real8_1d) the resized tile
!   ierr: (integer) if 1 and error occurred
subroutine resize_real8_tile_1d(tile,N_rows_new,N_data,ierr,offset_in)
  implicit none
  ! inputs-outputs:
  class(tile_real8_1d),intent(inout) :: tile
  integer,intent(inout)              :: ierr
  ! inputs:
  integer,intent(in)          :: N_rows_new,N_data
  integer,intent(in),optional :: offset_in
  ! variables
  integer :: offset
  real*8,dimension(N_data) :: tmp_data !< temporary data array

  ! initialization
  if(present(offset_in)) then
    offset = offset_in
  else
    offset = 0
  endif
  if((offset+N_data).gt.tile%N_rows) offset = max(0,tile%N_rows-N_data)
  if((offset+N_data).gt.N_rows_new) then
    ierr = 1
    return
  endif

  ! resizing
  tmp_data = tile%data_array(offset+1:offset+N_data)
  call tile%deallocate_tile(ierr)
  call tile%allocate_tile(N_rows_new,ierr)
  tile%data_array(offset+1:offset+N_data) = tmp_data
  
end subroutine resize_real8_tile_1d

! resize_real8_2d resizes a double 2d data_array of a tile
! preserving the data vwithin the offset and offset+N_data
! data are stored in offset:offset+N_data
! inputs:
!   tile:           (tile_real8_2d) the tile to resize
!   N_rows_new:     (integer) desired number of rows
!   N_cols_new:     (integer) desired number of columns
!   N_data_rows:    (integer) number of stored data rows
!   N_data_cols:    (integer) number of stored data columns
!   ierr:           (integer) if 1 and error occurred
!   offset_rows_in: (integer)(optional) initial row index for storing data
!   offset_cols_in: (integer)(optional) initial columns index for storing dat
! outputs:
!   tile: (tile_real8_2d) the resized tile
!   ierr: (integer) if 1 and error occurred
subroutine resize_real8_tile_2d(tile,N_rows_new,N_cols_new,N_data_rows,N_data_cols,&
ierr,offset_rows_in,offset_cols_in)
  implicit none
  ! inputs-outputs:
  class(tile_real8_2d),intent(inout) :: tile
  integer,intent(inout)              :: ierr 
  ! inputs:
  integer,intent(in)          :: N_rows_new,N_cols_new
  integer,intent(in)          :: N_data_rows,N_data_cols
  integer,intent(in),optional :: offset_rows_in,offset_cols_in
  ! variables
  integer :: offset_rows,offset_cols
  real*8,dimension(N_data_rows,N_data_cols) :: tmp_data !< temporary data array

  ! initialization
  if(present(offset_rows_in)) then
    offset_rows = offset_rows_in
  else
    offset_rows = 0
  endif
  if(present(offset_cols_in)) then
    offset_cols = offset_cols_in
  else
    offset_cols = 0
  endif
  if((offset_rows+N_data_rows).gt.tile%N_rows) offset_rows = max(0,tile%N_rows-N_data_rows)
  if((offset_cols+N_data_cols).gt.tile%N_cols) offset_cols = max(0,tile%N_cols-N_data_cols)
  if((offset_rows+N_data_rows).gt.N_rows_new) then
    ierr = 1
    return
  endif
  if((offset_cols+N_data_cols).gt.N_cols_new) then
    ierr = 1
    return
  endif

  ! resizing
  tmp_data = &
  tile%data_array(offset_rows+1:offset_rows+N_data_rows,offset_cols+1:offset_cols+N_data_cols)
  call tile%deallocate_tile(ierr)
  call tile%allocate_tile(N_rows_new,N_cols_new,ierr)
  tile%data_array(offset_rows+1:offset_rows+N_data_rows,offset_cols+1:offset_cols+N_data_cols) = &
  tmp_data
  
end subroutine resize_real8_tile_2d

!------------------------------------------------------------------

end module mod_tiles
