!> Verify the generation of few simple grids
module mod_grid_generation_mpi_test
use fruit
use fruit_mpi
use data_structure

implicit none
private
public :: run_fruit_grid_generation_mpi
!> Variables --------------------------------------
integer,parameter       :: string_len=40
integer,parameter       :: nx=10
integer,parameter       :: ny=10
integer,parameter       :: n_corners=4
real*8,parameter        :: R_begin=5.d-1
real*8,parameter        :: R_end=1.5d0
real*8,parameter        :: Z_begin=-5.d-1
real*8,parameter        :: Z_end=5.d-1
real*8,parameter        :: tol_node_pos=1d-20
type(type_node_list)    :: test_nodes
type(type_element_list) :: test_elements
integer :: rank_loc,n_tasks_loc,ifail_loc
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_grid_generation_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: grid generation mpi"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: grid generation mpi"
  call run_test_case(test_sqaure_grid,'test_sqaure_grid')
  write(*,'(/A)') "  ... tearing-down: grid generation mpi"
  call teardown(rank,n_tasks,ifail) 
end subroutine run_fruit_grid_generation_mpi

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
!> Test some reasonable properties of our simple grid.
!> Note that we should not use knowledge of the grid
!> construction algorithm and we should not use the
!> node ordering or element ordering in any way.
subroutine test_sqaure_grid
  use mod_projection_helpers_test_tools, only: default_square_grid
  implicit none
  character(len=22),parameter :: message='Error square grid test'
  character(len=string_len)   :: string
  !> generate the square grid
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,test_nodes,test_elements,&
  ifail_loc,R_begin,R_end,Z_begin,Z_end)
  !> Verify all nodes in box R_begin,R_end,Z_begin,Z_end
  call assert_true(minval(test_nodes%node(1:test_nodes%n_nodes)%x(1,1,1)).ge.R_begin,&
  message//': found node wth R < R_begin!')
  call assert_true(maxval(test_nodes%node(1:test_nodes%n_nodes)%x(1,1,1)).le.R_end,&
  message//': found node wth R > R_end!')
  call assert_true(minval(test_nodes%node(1:test_nodes%n_nodes)%x(1,1,2)).ge.Z_begin,&
  message//': found node wth Z < Z_begin!')
  call assert_true(maxval(test_nodes%node(1:test_nodes%n_nodes)%x(1,1,2)).le.Z_end,&
  message//': found node wth Z > Z_end!') 
  !> Test if every element has valid nodes
  call check_all_valid_nodes(message,test_nodes,test_elements)
  !> Verify there are no duplicate node positions
  call check_no_duplicate_node_positions(tol_node_pos,message,test_nodes)
  !> Test every node is used by 1-4 elments
  call check_all_nodes_used(message,test_nodes,test_elements)
  !> Test there are no nodes inside the quadrilateral defined by element vertices 
  !> https://stackoverflow.com/questions/5922027/how-to-determine-if-a-point-is-within-a-quadrilateral
  call check_node_inside_element(message,test_nodes,test_elements)
end subroutine test_sqaure_grid

!> Tools ------------------------------------------
!> check that all nodes of all elements are valid
subroutine check_all_valid_nodes(message,node_list,element_list)
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list
  character(len=*),intent(in)           :: message
  integer                   :: ii,jj
  character(len=string_len) :: string
  do ii=1,element_list%n_elements
    do jj=1,n_vertex_max
      write(string,'(A,I3,A,I3)') 'element ',ii,', vertex ',jj
      call assert_true(element_list%element(ii)%vertex(jj).le.node_list%n_nodes,&
      trim(adjustl(message))//': element has not vertex '//trim(adjustl(string)))
      call assert_true(element_list%element(ii)%vertex(jj).gt.0,&
      trim(adjustl(message))//': element has not vertex '//trim(adjustl(string)))
    enddo
  enddo
end subroutine check_all_valid_nodes

!> Verify there are no duplicate node positions
subroutine check_no_duplicate_node_positions(pos_tol,message,node_list)
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  real*8,intent(in)                     :: pos_tol
  character(len=*),intent(in)           :: message
  integer                   :: ii,jj
  character(len=string_len) :: string
  do ii=1,node_list%n_nodes
    do jj=1,node_list%n_nodes
      if(ii.ne.jj) then
        write(string,'(I3,A,I3)') ii,',',jj
        call assert_false(norm2(node_list%node(ii)%x(1,1,:)-node_list%node(jj)%x(1,1,:)).le.&
        pos_tol,trim(adjustl(message))//': nodes '//trim(adjustl(string))//&
        ' have the same position!')
      endif
    enddo
  enddo
end subroutine check_no_duplicate_node_positions

!> Test every node is used by 1-4 elments
subroutine check_all_nodes_used(message,node_list,element_list)
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list
  character(len=*),intent(in)           :: message
  integer                   :: ii,jj,cc
  character(len=string_len) :: string
  do ii=1,node_list%n_nodes
    cc=0
    do jj=1,n_vertex_max
      cc = cc + count(element_list%element(1:element_list%n_elements)%vertex(jj).eq.ii)
      write(string,'(A,I3,A,I3)') 'node ',ii,', vertex ',jj
    enddo
      call assert_true(cc.gt.0,trim(adjustl(message))//&
      ': orphan node found for '//trim(adjustl(string)))
      call assert_true(cc.le.n_vertex_max,trim(adjustl(message))//&
      ': overused node for '//trim(adjustl(string)))
  enddo
end subroutine check_all_nodes_used

!> Test there are no nodes inside the quadrilateral defined by element vertices 
!> https://stackoverflow.com/questions/5922027/how-to-determine-if-a-point-is-within-a-quadrilateral
subroutine check_node_inside_element(message,node_list,element_list)
  implicit none
  type(type_node_list),intent(inout)    :: node_list
  type(type_element_list),intent(inout) :: element_list
  character(len=*),intent(in)           :: message
  integer                     :: ii,jj
  real*8                      :: area,quad_area
  real*8,dimension(2)         :: xt
  real*8,dimension(n_corners) :: x,z
  character(len=string_len)   :: string
  do ii=1,element_list%n_elements
    x = node_list%node(element_list%element(ii)%vertex(:))%x(1,1,1)
    z = node_list%node(element_list%element(ii)%vertex(:))%x(1,1,2)
    !> Calculate area of quad
    xt = [sum(x)/4, sum(z)/4] ! test point in quad
    quad_area = triangle_area([x(1),z(1)], [x(2),z(2)], xt) + &
                triangle_area([x(2),z(2)], [x(3),z(3)], xt) + &
                triangle_area([x(3),z(3)], [x(4),z(4)], xt) + &
                triangle_area([x(4),z(4)], [x(1),z(1)], xt)
    do jj=1,node_list%n_nodes
      if (.not. any(element_list%element(ii)%vertex .eq. jj)) then
        xt = node_list%node(jj)%x(1,1,:)
        !> Calculate sum of areas of triangles between point and consecutive pairs of vertices
        area = triangle_area([x(1),z(1)], [x(2),z(2)], xt) + &
               triangle_area([x(2),z(2)], [x(3),z(3)], xt) + &
               triangle_area([x(3),z(3)], [x(4),z(4)], xt) + &
               triangle_area([x(4),z(4)], [x(1),z(1)], xt)
        write(string,'(A,i3,A,i3)') 'element', ii, ', node ', jj
        call assert_true(area .gt. quad_area, trim(adjustl(message))//&
        ': node outside of quad of other element '//string)
      end if
    end do
  end do
end subroutine check_node_inside_element

!> Use Heron's formula to calculate the area
function triangle_area(x1,x2,x3) result(area)
  real*8,dimension(2),intent(in) :: x1,x2,x3 !< corners
  real*8 :: s,a,b,c
  real*8 :: area
  a=norm2(x1-x2); b=norm2(x2-x3); c=norm2(x3-x1);
  s = (a+b+c)*5.d-1
  area = sqrt(max(s*(s-a)*(s-b)*(s-c), 0.d0))
end function triangle_area
!> ------------------------------------------------
end module mod_grid_generation_mpi_test
