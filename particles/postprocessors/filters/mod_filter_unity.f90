!> mod_filter_unity is a dummy filter which returns 1
!> whatever is the input
module mod_filter_unity
use mod_filter, only: filter
implicit none

private 
public :: filter_unity

!> Variables and basic type definitions -------------------
!> define the filter unity type
type,extends(filter) :: filter_unity
  contains
  procedure,pass(filter_inout) :: init_filter => init_filter_unity
  procedure,pass(filter_inout) :: compute_filter_from_position => compute_filter_unity_position
  procedure,pass(filter_inout) :: compute_filter_from_position_vectorial => &
  compute_filter_unity_position_vectorial
  procedure,pass(filter_inout) :: compute_filter_from_position_matrix => &
  compute_filter_unity_position_matrix
end type filter_unity

contains
!> Procedures ---------------------------------------------
!> initialise the unity filter
!> inouts:
!>   filter_unity:     (filter_unity) unity filter to be initialised
!>   n_dimensions:     (integer) dimension of the filter 1d,2d,3d,...
!>   stencil_shape_in: (integer)(n_dimensions)(optional) filter stencil shape (not used)
!> outputs:
!>   filter_unity: (filter_unity) initialised unity filter
subroutine init_filter_unity(filter_inout,n_dimensions,stencil_shape_in)
  implicit none
  class(filter_unity),intent(inout) :: filter_inout
  integer,intent(in) :: n_dimensions
  integer,dimension(n_dimensions),intent(in),optional :: stencil_shape_in
  filter_inout%n_dimensions = n_dimensions
end subroutine init_filter_unity

!> return the unity filter value or 1
!> inputs:
!>   filter_inout: (filter_unity) unity filter
!>   pos:          (real8)(n_dimensions) position where the filter weight is computed
!> outputs:
!>   filter_inout: (filter_unity) unity filter
!>   weight:       (real8) filter weight = 1
subroutine compute_filter_unity_position(filter_inout,pos,weight)
  implicit none
  !> inputs-outputs:
  class(filter_unity),intent(inout) :: filter_inout
  !> inputs:
  real*8,dimension(filter_inout%n_dimensions),intent(in) :: pos
  !> outputs:
  real*8,intent(out) :: weight
  weight = 1.d0
end subroutine compute_filter_unity_position

!> return the unity filter value or 1: vectorial version
!> inputs:
!>   filter_inout: (filter_unity) unity filter
!>   n_positions:  (integer) number of positions
!>   pos:          (real8)(n_dimensions,n_positions) positions where the 
!>                        filter weights are computed
!> outputs:
!>   filter_inout: (filter_unity) unity filter
!>   weights:      (real8)(n_positions) filter weight = 1
subroutine compute_filter_unity_position_vectorial(filter_inout,n_positions,pos,weights)
  implicit none
  !> inputs-outputs:
  class(filter_unity),intent(inout) :: filter_inout
  !> inputs:
  integer,intent(in) :: n_positions
  real*8,dimension(filter_inout%n_dimensions,n_positions),intent(in) :: pos
  !> outputs:
  real*8,dimension(n_positions),intent(out) :: weights
  weights = 1.d0
end subroutine compute_filter_unity_position_vectorial

!> return the unity filter value or 1: matrix version
!> inputs:
!>   filter_inout: (filter_unity) unity filter
!>   n_positions:  (integer) number of positions
!>   pos:          (real8)(n_dimensions,n_positions,n_intervals) 
!>                        positions where the filter weights are computed
!> outputs:
!>   filter_inout: (filter_unity) unity filter
!>   weights:      (real8)(n_positions,n_intervals) filter weight = 1
subroutine compute_filter_unity_position_matrix(filter_inout,n_intervals,n_positions,pos,weights)
  implicit none
  !> inputs-outputs:
  class(filter_unity),intent(inout) :: filter_inout
  !> inputs:
  integer,intent(in) :: n_positions,n_intervals
  real*8,dimension(filter_inout%n_dimensions,n_positions,n_intervals),intent(in) :: pos
  !> outputs:
  real*8,dimension(n_positions,n_intervals),intent(out) :: weights
  weights = 1.d0
end subroutine compute_filter_unity_position_matrix

!>---------------------------------------------------------
end module mod_filter_unity

