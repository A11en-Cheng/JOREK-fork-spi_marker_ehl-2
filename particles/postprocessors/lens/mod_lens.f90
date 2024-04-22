!> mod_lens is the base dataype containing variables and
!> procedure used for defining camera lenses
module mod_lens
implicit none

private
public :: lens

!> Variables ------------------------------------------------------
type,abstract :: lens
  integer :: n_x
  real*8,dimension(:),allocatable :: center !< lens center in cartesian coord.
  contains
  procedure,pass(lens_inout) :: allocate_lens
  procedure,pass(lens_inout) :: deallocate_lens
  procedure(ret_n_lens_inps_int),pass(lens_in),deferred   :: return_n_lens_inputs
  procedure(read_lens_inps_int),pass(lens_inout),deferred :: read_lens_inputs
  procedure(lens_sampling_int),pass(lens_inout),deferred  :: sampling
  procedure(lens_pdf_int),pass(lens_inout),deferred       :: pdf
end type lens

!> Interfaces -----------------------------------------------------
interface
  !> interface of the procedure for reading the lens inputs
  !> inputs:
  !>   lens_inout: (lens) lens class
  !>   my_id:      (integer) MPI rank
  !>   r_unit:     (integer) read unit id
  !>   n_inputs:   (integer)(2) number of integer and real inputs
  !> outputs:
  !>   lens_inout: (lens) lens class
  !>   n_inputs:   (integer)(2) number of integer and real inputs
  !>   int_param:  (integer)(:) integer input parameters
  !>   real_param: (real8)(:) real8 input parameters
  subroutine read_lens_inps_int(lens_inout,my_id,r_unit,&
  n_inputs,int_param,real_param)
    IMPORT :: lens
    implicit none
    !> inputs-outputs:
    class(lens),intent(inout)          :: lens_inout
    integer,dimension(2),intent(inout) :: n_inputs
    !> inputs:
    integer,intent(in) :: r_unit,my_id
    !> outputs:
    integer,dimension(:),allocatable,intent(out) :: int_param
    real*8,dimension(:),allocatable,intent(out)  :: real_param
  end subroutine read_lens_inps_int

  !> interface of the function returning the number of
  !> integer and real8 inputs
  !> inputs:
  !>   lens_in:  (lens) lens class
  !> outputs:
  !>   n_inputs: (integer)(2) number of integer and real inputs
  function ret_n_lens_inps_int(lens_in) result(n_inputs)
    IMPORT :: lens
    implicit none
    !> inputs:
    class(lens),intent(in) :: lens_in
    !> outputs:
    integer,dimension(2) :: n_inputs
  end function ret_n_lens_inps_int

  !> interface for the lens sampling procedure
  !> inputs:
  !>   lens_inout: (lens) lens class
  !>   n_samples:  (integer) number of samples to generate
  !> outputs:
  !>   x_pos:      (n_x,n_samples) positions on the lens
  subroutine lens_sampling_int(lens_inout,n_samples,x_pos)
    IMPORT :: lens
    implicit none
    !> inputs-outputs
    class(lens),intent(inout) :: lens_inout
    !> inputs:
    integer,intent(in)        :: n_samples
    !> outputs:
    real*8,dimension(lens_inout%n_x,n_samples),intent(out) :: x_pos
  end subroutine lens_sampling_int

  !> interface for the computation of the lens pdf
  !> inputs:
  !>   lens_inout: (lens) lens class
  !>   n_samples:  (integer) number of samples to generate
  !>   x_pos:     (n_x,n_samples) positions on the lens
  !> outputs:
  !>   pdf:       (n_x,n_samples) positions on the lens
  subroutine lens_pdf_int(lens_inout,n_samples,x_pos,pdf)
    IMPORT :: lens
    implicit none
    !> inputs-outputs
    class(lens),intent(inout) :: lens_inout
    !> inputs
    integer,intent(in) :: n_samples
    real*8,dimension(lens_inout%n_x,n_samples),intent(in) :: x_pos
    !> outputs
    real*8,dimension(n_samples),intent(out) :: pdf
  end subroutine lens_pdf_int  
end interface

contains
!> Procedures -----------------------------------------------------

!> allocate lens basic type
!> inputs:
!>   lens_inout: (lens) lens class to be allocated
!>   n_x:        (integer) n_x number of coefficients for position
!>   center:     (real8)(n_x) center for initialisation, default: 0
!> outputs:
!>   lens_inout: (lens) allocated lens
subroutine allocate_lens(lens_inout,n_x,center)
  implicit none
  !> inputs-outputs
  class(lens),intent(inout) :: lens_inout
  !> inputs
  integer :: n_x
  real*8,dimension(n_x),intent(in),optional :: center
  if(allocated(lens_inout%center)) then
    if(lens_inout%n_x.ne.n_x) deallocate(lens_inout%center)
  endif
  if(.not.allocated(lens_inout%center)) then
    allocate(lens_inout%center(n_x)); lens_inout%n_x = n_x;
  endif
  if(present(center)) then
    lens_inout%center = center
  else
    lens_inout%center = 0.d0
  endif
end subroutine allocate_lens

!> deallocate lens basic type and reset arrays
!>   lens_inout: (lens) lens class to be deallocated
!> outputs:
!>   lens_inout: (lens) deallocated lens
subroutine deallocate_lens(lens_inout)
  implicit none
  !> inputs-outputs
  class(lens),intent(inout) :: lens_inout
  if(allocated(lens_inout%center)) deallocate(lens_inout%center)
  lens_inout%n_x = 0
end subroutine deallocate_lens
!>-----------------------------------------------------------------
end module mod_lens
