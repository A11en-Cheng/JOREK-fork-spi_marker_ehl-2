!> the mod_pinhole_lens contains variables and procedures defining
!> the pinhole lens model
module mod_pinhole_lens
use mod_lens, only: lens
implicit none

private
public :: pinhole_lens

!> Variables ----------------------------------------------------------
type,extends(lens) :: pinhole_lens
  contains
  procedure,pass(lens_inout) :: init_pinhole
  procedure,pass(lens_in)    :: return_n_lens_inputs => return_n_pinhole_inputs
  procedure,pass(lens_inout) :: read_lens_inputs => read_pinhole_inputs
  procedure,pass(lens_inout) :: sampling => pinhole_sampling
  procedure,pass(lens_inout) :: pdf      => pinhole_pdf
end type pinhole_lens

!> Interfaces ---------------------------------------------------------
contains

!> Procedures ---------------------------------------------------------

!> read the number pinhole lens inputs
!> procedure for reading the pinhole lens input parameters
!> inputs:
!>   lens_inout: (pinhole_lens) pinhole lens
!>   my_id:      (integer) mpi rank
!>   r_unit:     (integer) read unit id
!>   n_inputs:   (integer) number of integer and real input arrays
!> outputs:
!>   lens_inout: (pinhole_lens) pinhole lens
!>   n_inputs:   (integer) number of integer and real input arrays
!>   int_param:  (integer)(:) inputs parameters
!>   real_param: (real8)(:) integer parameters
subroutine read_pinhole_inputs(lens_inout,my_id,r_unit,n_inputs,&
int_param,real_param)
  implicit none
  !> inputs-outputs:
  class(pinhole_lens),intent(inout)  :: lens_inout
  integer,dimension(2),intent(inout) :: n_inputs
  !> inputs:
  integer,intent(in) :: my_id,r_unit
  !> outputs:
  integer,dimension(:),allocatable,intent(out) :: int_param
  real*8,dimension(:),allocatable,intent(out)  :: real_param
  !> variables:
  integer :: n_x
  real*8,dimension(:),allocatable :: pinhole_center
  !> initialisations and definitions
  namelist /pinhole_in/ n_x
  namelist /center_in/ pinhole_center
  !> read data
  if(my_id.eq.0) then
    n_inputs = lens_inout%return_n_lens_inputs()
    if(allocated(int_param)) deallocate(int_param)
    allocate(int_param(n_inputs(1)))
    read(r_unit,pinhole_in); 
    int_param(1) = n_x; n_inputs(2) = n_x;
    allocate(pinhole_center(n_x))
    read(r_unit,center_in)
    if(allocated(real_param)) deallocate(real_param)
    allocate(real_param(n_x)); real_param = pinhole_center;
  endif
  !> cleanup
  if(allocated(pinhole_center)) deallocate(pinhole_center)
end subroutine read_pinhole_inputs

!> return the number of integer and real parameters
!> inputs:
!>   lens_in:  (pinhole_lens) pinhole lens
!> outputs:
!>   n_inputs: (integer) number of integer and real input 
!>             arrays, -1 means that the number is defined
!>             by the inputs
function return_n_pinhole_inputs(lens_in) result(n_inputs)
  implicit none
  !> inputs:
  class(pinhole_lens),intent(in) :: lens_in
  !> outputs:
  integer,dimension(2) :: n_inputs
  n_inputs = (/1,-1/)
end function return_n_pinhole_inputs

!> initialise the pinhole camera
!> inputs:
!>   lens_inout:     (pinhole_lens) pinhole lens to be initialised
!>   n_x:            (integer) size of the pinhole position array
!>   pinhole_center: (n_x) position of the pinhole cartesian coord.
!> outputs:
!>   lens_inout: (pinhole_lens) initialised pinhole lens
subroutine init_pinhole(lens_inout,n_x,pinhole_center)
  implicit none
  !> inputs-outputs
  class(pinhole_lens),intent(inout) :: lens_inout
  !> inputs
  integer,intent(in) :: n_x
  real*8,dimension(n_x),intent(in) :: pinhole_center
  call lens_inout%allocate_lens(n_x,pinhole_center)
end subroutine init_pinhole

!> return a set of samples on the pinhole lens (the pinhole position)
!> inputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   n_samples:  (integer) number of samples to be generated
!> outputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   x_pos:      (n_x,n_samples) arrays of the sampled lens points
subroutine pinhole_sampling(lens_inout,n_samples,x_pos)
  implicit none
  !> inputs-outputs
  class(pinhole_lens),intent(inout) :: lens_inout
  !> inputs
  integer,intent(in) :: n_samples
  !> outputs:
  real*8,dimension(lens_inout%n_x,n_samples),intent(out) :: x_pos
  !>
  integer :: ii
  !$omp parallel do default(shared) firstprivate(n_samples) private(ii)
  do ii=1,n_samples
    x_pos(:,ii) = lens_inout%center
  enddo
  !$omp end parallel do
end subroutine pinhole_sampling

!> return return the pdf for a set of sampled positions
!> Due to the fact that the pdf of a pinhole camera is a 
!> delta-Dirac distribution, 1.d0 is returned by default
!> inputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   n_samples:  (integer) number of samples to be generated
!>   x_pos:      (n_x,n_samples) arrays of the sampled lens points
!> outputs:
!>   lens_inout: (pinhole_lens) the pinhole lens class
!>   pdf:        (n_samples) arrays contaning the sample pdfs
subroutine pinhole_pdf(lens_inout,n_samples,x_pos,pdf)
  implicit none
  !> inputs-outputs
  class(pinhole_lens),intent(inout) :: lens_inout
  !> inputs
  integer,intent(in) :: n_samples
  real*8,dimension(lens_inout%n_x,n_samples),intent(in)  :: x_pos
  !> outputs:
  real*8,dimension(n_samples),intent(out) :: pdf
  !>
  integer :: ii
  pdf = 1.d0
end subroutine pinhole_pdf 

!>---------------------------------------------------------------------
end module mod_pinhole_lens
