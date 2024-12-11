module mod_neighbours_mpi_test
use fruit
use fruit_mpi
use data_structure
implicit none
private
public :: run_fruit_neighbours_mpi
!> Fruit basket -----------------------------------
integer,parameter :: message_len=100
integer,parameter :: n_sides=4
integer,parameter :: ny=2
real*8,parameter  :: R_begin=0d0
real*8,parameter  :: Z_begin=0d0
real*8,parameter  :: Z_end=1.d0
real*8,parameter  :: tol=1d-12
logical,parameter :: fill_boundary=.false.
type(type_node_list)    :: test_nodes
type(type_element_list) :: test_elements
integer :: rank_loc,n_tasks_loc,ifail_loc
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_neighbours_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: neighbours mpi"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: neighbours mpi"
  call run_test_case(test_nearby_elements,'test_nearby_elements')
  call run_test_case(test_update_neighbours,'test_update_neighbours')
  call run_test_case(test_coordinate_in_neighbour_square,&
  'test_coordinate_in_neighbour_square')
  write(*,'(/A)') "  ... tearing-down: neighbours mpi"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_neighbours_mpi

subroutine setup(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=rank; n_tasks_loc=n_tasks; ifail_loc=ifail;
  call init_node_list(test_nodes, n_nodes_max, test_nodes%n_dof, n_var)
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=-1; n_tasks_loc=-1; ifail=ifail_loc;
  call dealloc_node_list(test_nodes)
end subroutine teardown

!> Tests ------------------------------------------
subroutine test_nearby_elements()
  use mod_element_rtree,                 only: nearby_elements
  use mod_projection_helpers_test_tools, only: default_square_grid
  implicit none
  integer,parameter :: nx=4
  real*8,parameter  :: R_end=3d0
  integer,dimension(nx-1),parameter   :: n_nearby=(/2,3,2/)
  integer,dimension(3,nx-1),parameter :: expected_nearby_ids=&
                                         reshape((/1,2,0,1,2,3,&
                                         2,3,0/),(/3,nx-1/))
  integer                          :: ii
  integer,dimension(:),allocatable :: nearby_ids
  character(len=message_len)       :: message
  ! Construct a few specific set of elements
  !
  ! 5-------6-------7-------8 x(2) = 1
  ! |       |       |       |
  ! |   1   |   2   |   3   |
  ! |       |       |       |
  ! 1-------2-------3-------4 x(2) = 0
  ! 0       1       2       3
  !x(1)
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,test_nodes,&
  test_elements,ifail_loc,R_begin,R_end,Z_begin,Z_end)
  !> Returns the insertion order
  do ii=1,nx-1
    call nearby_elements(test_nodes,test_elements,ii,nearby_ids)
    write(*,*) 'ii: ',ii,' nearby_ids: ',nearby_ids
    write(message,'(A,I0,A)') 'Error nearby element test: N# neighbours element ',&
    ii,' mismatch!'
    call assert_equals(n_nearby(ii),size(nearby_ids),trim(message))
    write(message,'(A,I0,A)') 'Error nearby element test: neighbours array element ',&
    ii,' incorrect!'
    call assert_equals(expected_nearby_ids(1:n_nearby(ii),ii),nearby_ids,&
    n_nearby(ii),trim(message))   
    if(allocated(nearby_ids)) deallocate(nearby_ids)
  enddo
end subroutine test_nearby_elements

!> Test a simple case for the update_neighbours routine
!> TODO:
!>  * Test degenerate elements (like axis)
!>  * Test multi-node points (like axis, xpoint)
!>  * Test elements with different orientation (co/counter)
!>
!> When manually creating grids make sure to set sizes 
!> correctly, they are used to determing the element bounding
!> boxes and this is nearby_elements.
subroutine test_update_neighbours
  use mod_element_rtree,                 only: rtree_initialized
  use mod_neighbours,                    only: neighbours
  use mod_neighbours,                    only: update_neighbours
  use mod_projection_helpers_test_tools, only: default_square_grid
  implicit none
  integer,parameter :: side_elm1=2
  integer,parameter :: side_elm2=4
  integer,parameter :: nx=3
  real*8,parameter  :: R_end=2d0
  integer,dimension(n_sides),parameter :: nb_elm1=(/0,2,0,0/)
  integer,dimension(n_sides),parameter :: nb_elm2=(/0,0,0,1/)
  integer           :: ii,inb_i,inb_j
  ! Construct a few specific set of elements
  !
  ! 4-------5-------6 x(2) = 1
  ! |       |       |
  ! |   1   |   2   |
  ! |       |       |
  ! 1-------2-------3 x(2) = 0
  ! 0       1       2
  !x(1)
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,test_nodes,&
  test_elements,ifail_loc,R_begin,R_end,Z_begin,Z_end)
  !> test neighbour side from element 1 to 2
  call assert_true(neighbours(test_nodes,test_elements%element(1),&
  test_elements%element(2),inb_i,inb_j),&
  'Error update neighbours: element 2 not neighbour element 1!')
  call assert_equals(side_elm1,inb_i,&
  'Error update neighbours: wrong side element i (1->2)!')
  call assert_equals(side_elm2,inb_j,&
  'Error update neighbours: wrong side element j (1->2)!')
  !> test neighbour side from element 2 to 1
  call assert_true(neighbours(test_nodes,test_elements%element(2),&
  test_elements%element(1),inb_i,inb_j),&
  'Error update neighbours: element 1 not neighbour element 2!')
  call assert_equals(side_elm2,inb_i,&
  'Error update neighbours: wrong side element i (2->1)!')
  call assert_equals(side_elm1,inb_j,&
  'Error update neighbours: wrong side element j (2->1)!')
  !> re-initialise tree and update neighbours
  rtree_initialized=.false.; !< force tree construction
  call update_neighbours(test_nodes,test_elements) 
  do ii=1,n_sides
    call assert_equals(nb_elm1(ii),test_elements%element(1)%neighbours(ii),&
    'Error update neighbours: element 1 not next to 2 after update!')
    call assert_equals(nb_elm2(ii),test_elements%element(2)%neighbours(ii),&
    'Error update neighbours: element 2 not next to 1 after update!')
  enddo
end subroutine test_update_neighbours

subroutine test_coordinate_in_neighbour_square()
  use mod_settings,                      only: n_vertex_max
  use mod_element_rtree,                 only: rtree_initialized
  use mod_neighbours,                    only: neighbours
  use mod_neighbours,                    only: update_neighbours
  use mod_neighbours,                    only: coord_in_neighbour
  use mod_projection_helpers_test_tools, only: default_square_grid
  implicit none
  integer,parameter :: nx=3
  integer,parameter :: neighbour_sol=2
  integer,dimension(n_vertex_max,2),parameter :: vertices_elm1=reshape(&
                                                 (/1,2,5,4,2,5,4,1/),&
                                                 (/n_vertex_max,2/))
  integer,dimension(n_vertex_max,3),parameter :: vertices_elm2=reshape(&
                                                 (/2,3,6,5,3,6,5,2,3,2,5,6/),&
                                                 (/n_vertex_max,3/))
  real*8,parameter  :: R_end=2d0
  real*8,dimension(3),parameter  :: s_sol=(/0d0,1d0,3d-1/)
  real*8,dimension(3),parameter  :: t_sol=(/0d0,1d0,3d-1/)
  real*8,dimension(2),parameter  :: s_test=(/1d0,3d-1/)
  real*8,dimension(2),parameter  :: t_test=(/3d-1,0d0/)
  integer             :: neighbour_id
  real*8,dimension(2) :: st
  
  ! Construct a specific set of elements
  ! 4-------5-------6 x(2) = 1
  ! |       |       |
  ! |   1   |   2   |
  ! |       |       |
  ! 1-------2-------3 x(2) = 0
  ! 0       1       2
  !x(1)
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,test_nodes,&
  test_elements,ifail_loc,R_begin,R_end,Z_begin,Z_end)
  ! test default element orientation
  ! 4-------5-------6 (node number)
  ! |4     3|4     3| (vertex number)
  ! |   1   |   2   |
  ! |1     2x1     2|
  ! 1-------2-------3
  test_elements%element(1)%vertex = vertices_elm1(:,1)
  test_elements%element(2)%vertex = vertices_elm2(:,1)
  call update_neighbours(test_nodes,test_elements)
  !> test 1 position on the boundary and check if 
  !> the transformation is correct
  st = (/s_test(1),t_test(1)/)
  call coord_in_neighbour(test_nodes,test_elements,1,neighbour_id,st)
  call assert_equals(neighbour_sol,neighbour_id,&
  'Error coordinate in neighbour square test: right element id is not 1 (default)!')
  call assert_equals([s_sol(1),t_sol(3)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 1 (default)!')
  call coord_in_neighbour(test_nodes,test_elements,2,neighbour_id,st)
  call assert_equals([s_sol(2),t_sol(3)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 2 (default)!')
  ! Rotate the right element and retest
  ! 4-------5-------6 (node number)
  ! |4     3|3     2| (vertex number)
  ! |   1   |   2   |
  ! |1     2x4     1|
  ! 1-------2-------3
  test_elements%element(2)%vertex = vertices_elm2(:,2)
  rtree_initialized = .false.
  call update_neighbours(test_nodes,test_elements)
  st = (/s_test(1),t_test(1)/)
  call coord_in_neighbour(test_nodes,test_elements,1,neighbour_id,st)
  call assert_equals(neighbour_sol,neighbour_id,&
  'Error coordinate in neighbour square test: right element id is not 2 (right rotation)!')
  call assert_equals([s_sol(3),t_sol(2)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 1 (right rotation)!')
  call coord_in_neighbour(test_nodes,test_elements,2,neighbour_id,st)
  call assert_equals([s_sol(2),t_sol(3)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 2 (right rotation)!')
  ! Rotate the left element and retest
  ! 4-------5-------6 (node number)
  ! |3     2|3     2| (vertex number)
  ! |   1   |   2   |
  ! |4     1x4     1|
  ! 1-------2-------3
  test_elements%element(1)%vertex = vertices_elm1(:,2)
  rtree_initialized = .false.
  call update_neighbours(test_nodes,test_elements)
  st = (/s_test(2),t_test(2)/)
  call coord_in_neighbour(test_nodes,test_elements,1,neighbour_id,st)
  call assert_equals(neighbour_sol,neighbour_id,&
  'Error coordinate in neighbour square test: right element id is not 3 (left rotation)!')
  call assert_equals([s_sol(3),t_sol(2)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 1 (left rotation)!')
  call coord_in_neighbour(test_nodes,test_elements,2,neighbour_id,st)
  call assert_equals([s_sol(3),t_sol(1)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 2 (left rotation)!')
  ! Reverse orientation of right element and retest
  ! 4-------5-------6 (node number)
  ! |3     2|3     4| (vertex number)
  ! |   1   |   2   |
  ! |4     1x2     1|
  ! 1-------2-------3
  test_elements%element(2)%vertex = vertices_elm2(:,3)
  rtree_initialized = .false.
  call update_neighbours(test_nodes,test_elements)
  st = (/s_test(2),t_test(2)/)
  call coord_in_neighbour(test_nodes,test_elements,1,neighbour_id,st)
  call assert_equals(neighbour_sol,neighbour_id,&
  'Error coordinate in neighbour square test: right element id is not 4 (right reverse)!')
  call assert_equals([s_sol(2),t_sol(3)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 1 (right reverse)!')
  call coord_in_neighbour(test_nodes,test_elements,2,neighbour_id,st)
  call assert_equals([s_sol(3),t_sol(1)],st,2,tol,&
  'Error coordinate in neighbour square test: st mismatch elm 2 (right reverse)!')  
end subroutine test_coordinate_in_neighbour_square

!> ------------------------------------------------
end module mod_neighbours_mpi_test
