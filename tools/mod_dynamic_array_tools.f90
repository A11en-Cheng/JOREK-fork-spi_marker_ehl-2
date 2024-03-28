! mod_dynamic_array_tools contains tool procedures for allocating
! and deallocating dynamic arrays with allocation checks
module mod_dynamic_array_tools
implicit none

private
public :: allocate_check,deallocate_check

! procedures for checking and allocating arrays
interface allocate_check
  module procedure allocate_check_integer
  module procedure allocate_check_integer_2d
  module procedure allocate_check_double
  module procedure allocate_check_double_2d
end interface

! procedures for checking and deallocating arrays
interface deallocate_check
  module procedure deallocate_check_integer
  module procedure deallocate_check_integer_2d
  module procedure deallocate_check_double
  module procedure deallocate_check_double_2d
end interface

contains

! Memory procedures ------------------------------------------------

! allocate integer array if not allocated
subroutine allocate_check_integer(N,array)
  implicit none
  integer,intent(in) :: N
  integer,dimension(:),allocatable,intent(inout) :: array
  if(.not.allocated(array)) allocate(array(N))
end subroutine allocate_check_integer

! allocate integer 2D-integer array if not allocated
subroutine allocate_check_integer_2d(N1,N2,array)
  implicit none
  integer,intent(in) :: N1,N2
  integer,dimension(:,:),allocatable,intent(inout) :: array
  if(.not.allocated(array)) allocate(array(N1,N2))
end subroutine allocate_check_integer_2d

! allocate double array if not allocated
subroutine allocate_check_double(N,array)
  implicit none
  integer,intent(in) :: N
  real*8,dimension(:),allocatable,intent(inout) :: array
  if(.not.allocated(array)) allocate(array(N))
end subroutine allocate_check_double

! allocate double 2D-array if not allocated
subroutine allocate_check_double_2d(N1,N2,array)
  implicit none
  integer,intent(in) :: N1,N2
  real*8,dimension(:,:),allocatable,intent(inout) :: array
  if(.not.allocated(array)) allocate(array(N1,N2))
end subroutine allocate_check_double_2d

! deallocate integer array if allocated
subroutine deallocate_check_integer(array)
  implicit none
  integer,dimension(:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_integer

! deallocate 2D-integer array if allocated
subroutine deallocate_check_integer_2d(array)
  implicit none
  integer,dimension(:,:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_integer_2d

! deallocate double array if allocated
subroutine deallocate_check_double(array)
  implicit none
  real*8,dimension(:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_double

! deallocate double 2D array if allocated
subroutine deallocate_check_double_2d(array)
  implicit none
  real*8,dimension(:,:),allocatable,intent(inout) :: array
  if(allocated(array)) deallocate(array)
end subroutine deallocate_check_double_2d

end module mod_dynamic_array_tools
