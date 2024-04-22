!> mod_filter_test contains variables an procedures used for
!> testing the filter class. The filter class is an abstract
!> class hence the filter_unity class is used instead.
module mod_filter_test
use fruit
implicit none

private
public :: run_fruit_filter

!> Variables and and parameters -------------------------------
character(len=13),parameter                                        :: input_file='filter_inputs'
integer,parameter                                                  :: read_unit=43
integer,parameter                                                  :: n_dimensions_sol=2
integer,parameter                                                  :: n_intervals_sol=3
integer,parameter                                                  :: n_positions_sol=13
real*8,parameter                                                   :: tol_real8=5.d-16
integer,dimension(n_dimensions_sol),parameter                      :: filter_shape_sol=(/5,10/)
real*8,dimension(n_dimensions_sol),parameter                       :: pos_lowbnd=(/-3.d1,1.d-2/)
real*8,dimension(n_dimensions_sol),parameter                       :: pos_uppbnd=(/3.d1,5.d0/)
real*8,dimension(n_dimensions_sol,n_positions_sol)                 :: positions_sol
real*8,dimension(n_dimensions_sol,n_positions_sol,n_intervals_sol) :: positions_2d_sol
real*8,dimension(n_positions_sol)                                  :: weights_sol
real*8,dimension(n_positions_sol,n_intervals_sol)                  :: weights_2d_sol

!> Interfaces -------------------------------------------------
contains
!> Fruit basket -----------------------------------------------
!> fruit basket containing all set-up, test and tear-down 
!> procedure to be executed
subroutine run_fruit_filter()
  implicit none
  write(*,*) "  ... setting-up filter tests"
  call setup
  call write_filter_input_file
  write(*,*) "  ... running: filter tests"
  call run_test_case(test_read_filter_inputs,'test_read_filter_inputs')
  call run_test_case(test_filter_allocation_deallocation,&
  'test_filter_allocation_deallocation')
  call run_test_case(test_initialisation_filter_unity,&
  'test_initialisation_filter_unity')
  call run_test_case(test_compute_weights_filter_unity,&
  'test_compute_weights_filter_unity')
  call run_test_case(test_compute_weights_filter_unity_vectorial,&
  'test_compute_weights_filter_unity_vectorial')
  call run_test_case(test_compute_weights_filter_unity_matrix,&
  'test_compute_weights_filter_unity_matrix')
  write(*,*) "  ... tearing-down: filter tests"
  call teardown
end subroutine run_fruit_filter

!> Set-up and tear-down ---------------------------------------
!> initiliase test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> variables
  integer :: ii,jj
  !> initilaise positions and weights
  do ii=1,n_positions_sol
    call gnu_rng_interval(n_dimensions_sol,pos_lowbnd,pos_uppbnd,positions_sol(:,ii))
  enddo 
  do jj=1,n_intervals_sol
    do ii=1,n_positions_sol
      call gnu_rng_interval(n_dimensions_sol,pos_lowbnd,pos_uppbnd,positions_2d_sol(:,ii,jj))
    enddo
  enddo
  weights_sol = 1.d0;  weights_2d_sol = 1.d0;
end subroutine setup

!> method used for generating an input file
subroutine write_filter_input_file()
  implicit none
  !> variables
  integer :: ifail
  !> write file
  open(read_unit,file=input_file,status='unknown',action='write',iostat=ifail)
  write(read_unit,'(/A)') '&filter_in'
  write(read_unit,'(/A,I10)') 'n_dimensions = ',n_dimensions_sol
  write(read_unit,'(/A)') '/'
  write(read_unit,'(/A)') '&filter_stencil_in'
  write(read_unit,'(/A,I10)') 'stencil_shape(1) = ',filter_shape_sol(1)
  write(read_unit,'(/A,I10)') 'stencil_shape(2) = ',filter_shape_sol(2)
  write(read_unit,'(/A)') '/'
  close(read_unit)
end subroutine write_filter_input_file

!> tear down test features
subroutine teardown()
  implicit none
  call system('rm '//input_file)
end subroutine teardown

!> Tests ------------------------------------------------------
!> test read filter from output procedures
subroutine test_read_filter_inputs()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  type(filter_unity)                    :: filter_test
  integer                               :: rank,ifail
  integer,dimension(2)                  :: n_inputs
  integer,dimension(1+n_dimensions_sol) :: int_param_sol
  integer,dimension(:),allocatable      :: int_param
  real*8,dimension(:),allocatable       :: real_param
  !> initialisationsi
  rank = 0; int_param_sol(1) = n_dimensions_sol;
  int_param_sol(2:n_dimensions_sol+1) = filter_shape_sol
  n_inputs = filter_test%return_n_filter_inputs()
  !> read data
  open(read_unit,file=input_file,status='old',action='read',iostat=ifail)
  call filter_test%read_filter_inputs(rank,read_unit,n_inputs,&
  int_param,real_param)
  close(read_unit)
  !> checks
  call assert_equals((/n_dimensions_sol+1,0/),n_inputs,2,&
  "Error filter unity read inputs: N# inputs mismatch!")
  call assert_equals(int_param_sol,int_param,1+n_dimensions_sol,&
  "Error filter unity read inputs: integer inputs mismatch!")
  call assert_false(allocated(real_param),&
  "Error filter unity read inputs: real input array allocated!")
  !> cleanup
  if(allocated(int_param)) deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine test_read_filter_inputs

!> test the filter allocation and deallocation methods
subroutine test_filter_allocation_deallocation()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_filter_unity,        only: filter_unity
  implicit none
  !> variables
  type(filter_unity) :: filter_test
  !> tests allocation
  call filter_test%allocate_filter(n_dimensions_sol)
  call assert_equals(n_dimensions_sol,filter_test%n_dimensions,&
  "Error allocate filter: n_dimensions mistmatech!")
  call assert_equals_allocatable_arrays(n_dimensions_sol,&
  filter_test%stencil_shape,0,"Error allocate filter: stencil_shape")
  !> test deallocation
  call filter_test%deallocate_filter
  call assert_equals(0,filter_test%n_dimensions,&
  "Error deallocate filter: n_dimensions not reset to 0!")
  call assert_false(allocated(filter_test%stencil_shape),&
  "Error deallocate filter: stencil shape not deallocated!")
end subroutine test_filter_allocation_deallocation

!> test initialisation of filter unity
subroutine test_initialisation_filter_unity()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  type(filter_unity) :: filter_test
  !> test initialisation without stencil shape
  call filter_test%init_filter(n_dimensions_sol)
  call assert_equals(n_dimensions_sol,filter_test%n_dimensions,&
  "Error filter unity initialisation: n dimensions mismatch")
  call filter_test%deallocate_filter
  !> test initialisation with stencil shape
  call filter_test%init_filter(n_dimensions_sol,filter_shape_sol)
  call assert_equals(n_dimensions_sol,filter_test%n_dimensions,&
  "Error filter unity initialisation stencil shape: n dimensions mismatch")
end subroutine test_initialisation_filter_unity 

!> test compute weights filter unity
subroutine test_compute_weights_filter_unity()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(n_positions_sol) :: weights
  type(filter_unity) :: filter_test
  !> initialisation
  call filter_test%init_filter(n_dimensions_sol)
  !> test
  do ii=1,n_positions_sol
    call filter_test%compute_filter_from_position(positions_sol(:,ii),weights(ii))
  enddo
  call assert_equals(weights_sol,weights,n_positions_sol,tol_real8,&
  "Error filter unity compute filter form positions: filter weights mistmatch!")
  !> deallocation
  call filter_test%deallocate_filter
end subroutine test_compute_weights_filter_unity

!> test compute weights filter unity: vectorial version
subroutine test_compute_weights_filter_unity_vectorial()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  real*8,dimension(n_positions_sol) :: weights
  type(filter_unity) :: filter_test
  !> initialisation
  call filter_test%init_filter(n_dimensions_sol)
  !> test
  call filter_test%compute_filter_from_position_vectorial(n_positions_sol,positions_sol,weights)
  call assert_equals(weights_sol,weights,n_positions_sol,tol_real8,&
  "Error filter unity compute filter form positions (vectorial): filter weights mistmatch!")
  !> deallocation
  call filter_test%deallocate_filter
end subroutine test_compute_weights_filter_unity_vectorial

!> test compute weights filter unity: matrix version
subroutine test_compute_weights_filter_unity_matrix()
  use mod_filter_unity, only: filter_unity
  implicit none
  !> variables
  real*8,dimension(n_positions_sol,n_intervals_sol) :: weights_2d
  type(filter_unity) :: filter_test
  !> initialisation
  call filter_test%init_filter(n_dimensions_sol)
  !> test
  call filter_test%compute_filter_from_position_matrix(n_positions_sol,n_intervals_sol,&
  positions_2d_sol,weights_2d)
  call assert_equals(weights_2d_sol,weights_2d,n_positions_sol,n_intervals_sol,tol_real8,&
  "Error filter unity compute filter form positions (matrix): filter weights mistmatch!")
  !> deallocation
  call filter_test%deallocate_filter
end subroutine test_compute_weights_filter_unity_matrix

!> Tools ------------------------------------------------------
!>-------------------------------------------------------------
end module mod_filter_test
