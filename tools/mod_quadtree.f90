module mod_quadtree
!< this module is derived from the octree module (from wall_collisions), simplified to an quad_tree in 2D
!< The routine quadtree_find gives the node which contains a given point. 
!< The node contains the list of elements that are present within the bounding box of this node.
!<
use data_structure
use mod_interp

implicit none

integer, parameter :: max_depth = 7
real*8, parameter  :: tolerance = 1.d-6

!> Quad representing one element.
type quadtree_quad
  integer      :: quad_id     !< Unique identifier
  real(kind=8) :: vertex(4,2) !< R,Z coordinates of boundingbox
end type quadtree_quad

!> Node / Area belonging to the quadtree
type quadtree_node
integer      :: depth = 0                        !< The depth of this node, for root depth = 1
real(kind=8) :: boundary(2,2)                    !< [x_min, x_max; y_min, y_max, z_min, z_max] of the volume belonging this node
type(quadtree_quad), allocatable :: contained(:) !< Quads contained within this node.
type(quadtree_node), pointer     :: children(:)  !< Child nodes this element has.
end type quadtree_node

type(quadtree_node)              :: quadtree    !< root node of the quadtree

contains

!< Initialize quadtree
subroutine quadtree_init(node_list, element_list)
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
    
type(quadtree_quad), allocatable :: quads(:)
real*8                           :: boundary(2,2) !< Global bounding box
integer                          :: i, n, n_elements
real*8                           :: rmin, rmax, zmin, zmax 
logical  :: do_plot = .false.

  write(*,'(A)')      '*************************************'
  write(*,'(A,i3,A)') '* Quadtree_init, max depth : ',max_depth,'    *'
  write(*,'(A)')      '*************************************'

  call quadtree_free(quadtree) 

  quadtree%depth    = 1
  
  n_elements = element_list%n_elements

  allocate(quads(n_elements))
  boundary(1,:) =  1.d33
  boundary(2,:) = -1.d33

  do i=1,n_elements
    call RZ_minmax(node_list, element_list, i, rmin, rmax, zmin, zmax)
    quads(i)%quad_id = i
    quads(i)%vertex(1,:) = (/ rmin, zmin /)
    quads(i)%vertex(2,:) = (/ rmax, zmin /)
    quads(i)%vertex(3,:) = (/ rmax, zmax /)
    quads(i)%vertex(4,:) = (/ rmin, zmax /)
    boundary(1,1) = min(boundary(1,1), rmin)
    boundary(1,2) = min(boundary(1,2), zmin)
    boundary(2,1) = max(boundary(2,1), rmax)
    boundary(2,2) = max(boundary(2,2), zmax)
  end do
  boundary(1,:) = boundary(1,:) - tolerance
  boundary(2,:) = boundary(2,:) + tolerance

  quadtree%boundary = boundary

  call quadtree_construction(quads, max_depth, quadtree)

  if (do_plot) then
    call plot_node(quadtree, .false.)

    do i = 1, size(quads)
      call plot_quad(quads(i)%vertex)
    enddo
  endif

  deallocate(quads)
    
end subroutine quadtree_init
    
!< Routine for initializing nodes recursively
recursive subroutine quadtree_construction(quads, max_depth, node)
implicit none
    
type(quadtree_quad), intent(in)    :: quads(:)
integer, intent(in)                :: max_depth
type(quadtree_node), intent(inout) :: node
integer :: i, j
integer :: contained_quad

! Find how many quads are contained within this node
contained_quad = 0

do i = 1, size(quads)
    if ( .not. contains_quad(node%boundary, quads(i)) ) cycle
    contained_quad = contained_quad + 1
end do

allocate(node%contained(contained_quad))
    
! No quads in this node -> no need to make children
if( contained_quad .eq. 0 ) then
  return
endif    
! Repeat the loop, this time assigning all the quads that are within this node
j=1

do i = 1, size(quads)
  if ( .not. contains_quad(node%boundary, quads(i)) ) cycle
  node%contained(j)%quad_id = quads(i)%quad_id
  node%contained(j)%vertex  = quads(i)%vertex
  j = j + 1
end do

! Return if we have reached the max depth
if(node%depth==max_depth) return
    
! Make children and free the memory allocated for quads as those will be contained within child nodes
call node_division(node)
call quadtree_construction(node%contained, max_depth, node%children(1))
call quadtree_construction(node%contained, max_depth, node%children(2))
call quadtree_construction(node%contained, max_depth, node%children(3))
call quadtree_construction(node%contained, max_depth, node%children(4))

deallocate(node%contained)
    
end subroutine quadtree_construction
    
!< Function for checking if a given axis aligned box contains the given quad
!<
!< This function is conservative. In reality it doesn't check if the quad is
!< inside the box but simply checks if the bounding box of the quad and the given
!< box overlap. For now this results in some extra collision checks so feel free to improve
!< this algorithm if you want.
logical function contains_quad(nbox, quad)
real*8, intent(in)              :: nbox(2,2)
type(quadtree_quad), intent(in) :: quad
real*8                          :: bbox(2,2) ! Bounding box for the quad

bbox(1,1) = minval(quad%vertex(:,1))  
bbox(2,1) = maxval(quad%vertex(:,1))
bbox(1,2) = minval(quad%vertex(:,2))
bbox(2,2) = maxval(quad%vertex(:,2))

contains_quad = ( &
         ( ( bbox(1,1) .le. nbox(1,1) .and. nbox(1,1) .le. bbox(2,1) ) .or.    &
           ( bbox(1,1) .le. nbox(2,1) .and. nbox(2,1) .le. bbox(2,1) ) .or.    &
           ( nbox(1,1) .le. bbox(1,1) .and. bbox(1,1) .le. nbox(2,1) ) .or.    &
           ( nbox(1,1) .le. bbox(2,1) .and. bbox(2,1) .le. nbox(2,1) ) ) .and. &
         ( ( bbox(1,2) .le. nbox(1,2) .and. nbox(1,2) .le. bbox(2,2) ) .or.    &
           ( bbox(1,2) .le. nbox(2,2) .and. nbox(2,2) .le. bbox(2,2) ) .or.    &
           ( nbox(1,2) .le. bbox(1,2) .and. bbox(1,2) .le. nbox(2,2) ) .or.    &
           ( nbox(1,2) .le. bbox(2,2) .and. bbox(2,2) .le. nbox(2,2) ) )       &
           )

end function contains_quad
    
!< Divide node volume and allocated and initialize the children
subroutine node_division(node)
implicit none
    
type(quadtree_node), intent(inout), target :: node
integer :: i, j, k
real(kind=8) :: offset
    
! Add just a little bit of offset so that the boxes overlap
! and we don't end up with (numerical) gaps
offset = tolerance
    
allocate(node%children(4))

k = 1
do j = 1,2
  do i = 1,2
    node%children(k)%depth = node%depth + 1
               
    node%children(k)%boundary(1,1) = node%boundary(1,1) + (i-1) * (node%boundary(2,1) - node%boundary(1,1)) * 0.5d0 - offset
    node%children(k)%boundary(2,1) = node%boundary(2,1) - (2-i) * (node%boundary(2,1) - node%boundary(1,1)) * 0.5d0 + offset
    node%children(k)%boundary(1,2) = node%boundary(1,2) + (j-1) * (node%boundary(2,2) - node%boundary(1,2)) * 0.5d0 - offset
    node%children(k)%boundary(2,2) = node%boundary(2,2) - (2-j) * (node%boundary(2,2) - node%boundary(1,2)) * 0.5d0 + offset
    
    k = k + 1
  end do
end do
end subroutine node_division
   

!< Find the node containing the point
subroutine quadtree_find_point(node, x, nodeout, i_err)
implicit none
type(quadtree_node), target, intent(in)   :: node      !< The quadtree root node
real(kind=8), intent(in)                  :: x(2)      !< x,y,z coordinates of the query point
type(quadtree_node), pointer, intent(inout) :: nodeout !< Node containg the point
integer, intent(out)                      :: i_err     !< Non-zero if the point is outside the quadtree

if ( (x(1)) > node%boundary(1,1) .and. (x(1)) < node%boundary(2,1) .and. &
     (x(2)) > node%boundary(1,2) .and. (x(2)) < node%boundary(2,2) ) then
  i_err = 0
  call quadtree_traverse(node, x, nodeout)
else
  i_err = 1
endif

end subroutine quadtree_find_point
    

!< Find recursively the node where given point belongs to
recursive subroutine quadtree_traverse(node, x, nodeout)
implicit none
type(quadtree_node), target, intent(in)     :: node
real(kind=8), intent(in)                    :: x(2)
type(quadtree_node), pointer, intent(inout) :: nodeout
    
integer :: i

! This node has no children, return the node
if ( .not. associated(node%children) ) then
  nodeout => node
else
! Find the child to which the given point belongs to
  do i = 1,4
    if((x(1)) > node%children(i)%boundary(1,1) .and. &
       (x(1)) < node%children(i)%boundary(2,1) .and. &
       (x(2)) > node%children(i)%boundary(1,2) .and. &
       (x(2)) < node%children(i)%boundary(2,2) ) then
    
       ! Child found (no need to loop further), repeat the process for the child node
      call quadtree_traverse(node%children(i), x, nodeout)
      exit
    endif
  enddo
endif
    
end subroutine quadtree_traverse

    
!< Free resources used by the quadtree
!<
!< This routine should be called only for the root node.
subroutine quadtree_free(node)
implicit none
type(quadtree_node), intent(inout) :: node

call clean_node(node)
node%depth    = 0
node%boundary = 0

if (allocated(node%contained)) then
  deallocate(node%contained)
endif

end subroutine quadtree_free
    
!< Free node resources (children and quads) recursively
recursive subroutine clean_node(node)
implicit none
    
type(quadtree_node), intent(inout) :: node
integer :: i
    
if (associated(node%children)) then
  if (size(node%children) .gt. 0) then
    call clean_node(node%children(1))
    if ( allocated(node%children(1)%contained) ) deallocate(node%children(1)%contained)
    call clean_node(node%children(2))
    if ( allocated(node%children(2)%contained) ) deallocate(node%children(2)%contained)
    call clean_node(node%children(3))
    if ( allocated(node%children(3)%contained) ) deallocate(node%children(3)%contained)
    call clean_node(node%children(4))
    if ( allocated(node%children(4)%contained) ) deallocate(node%children(4)%contained)    
  endif
  deallocate(node%children)
  if (allocated(node%contained)) deallocate(node%contained)
endif
end subroutine clean_node


subroutine plot_quad(vertex)
real*8 :: vertex(4,2)

call lplot6(1,1,(/vertex(1,1),vertex(2,1),vertex(3,1),vertex(4,1),vertex(1,1) /), &
                (/vertex(1,2),vertex(2,2),vertex(3,2),vertex(4,2),vertex(1,2) /), -5, ' ')
return
end


recursive subroutine plot_node(node,initialised)
type(quadtree_node), intent(inout) :: node
logical :: initialised
integer, save :: i_color, i

if (.not. initialised) then
  call lincol(0)
  call lplot6(1,1,(/node%boundary(1,1),node%boundary(2,1),node%boundary(2,1),node%boundary(1,1),node%boundary(1,1)/), &
                  (/node%boundary(1,2),node%boundary(1,2),node%boundary(2,2),node%boundary(2,2),node%boundary(1,2)/),+5,' ')  
else
  i_color = node%depth-1 
  call lincol(i_color)
  call lplot6(1,1,(/node%boundary(1,1),node%boundary(2,1),node%boundary(2,1),node%boundary(1,1),node%boundary(1,1)/), &
                  (/node%boundary(1,2),node%boundary(1,2),node%boundary(2,2),node%boundary(2,2),node%boundary(1,2)/),-5,' ')  
endif

if (associated(node%children)) then      ! replacing the 4 statements by a do loop does not seem to work
   if (size(node%children) .ge. 1) call plot_node(node%children(1),.true.)
   if (size(node%children) .ge. 2) call plot_node(node%children(2),.true.)
   if (size(node%children) .ge. 3) call plot_node(node%children(3),.true.)
   if (size(node%children) .ge. 4) call plot_node(node%children(4),.true.)
else
  call lincol(0)
  do i=1, size(node%contained)
!    call plot_quad(node%contained(i)%vertex)
  enddo
endif

end

recursive subroutine print_node(node)
implicit none
type(quadtree_node), intent(inout) :: node        !< Root node

write(11,'(A,i3,4f8.4,L)') 'print node ', node%depth,  node%boundary, associated(node%children)

if (associated(node%children)) then
  write(11,*) node%depth,' number of children : ',size(node%children)
  if (size(node%children) .ge. 1) call print_node(node%children(1))
  if (size(node%children) .ge. 2) call print_node(node%children(2))
  if (size(node%children) .ge. 3) call print_node(node%children(3))
  if (size(node%children) .ge. 4) call print_node(node%children(4))
else
  write(11,*) ' no more children'
  write(11,*) node%depth,' contained : ',node%contained(:)%quad_id
endif
return
end

! try to find all the elements in the box continaing i_elm
subroutine nearby_elements_quadtree2(node_list, element_list, i_elm, i_nearby)
  implicit none
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  integer                             :: i_elm
  integer, allocatable, intent(out)   :: i_nearby(:)
  
  integer, dimension(n_vertex_max)    :: vertices
  real*8                              :: x1(2), x2(2), x3(2), x4(2), rmin, rmax, zmin, zmax
  type(quadtree_node), pointer        :: nodeout_1, nodeout_2, nodeout_3, nodeout_4    
  integer, allocatable                :: i_nearby_1(:), i_nearby_2(:), i_nearby_3(:), i_nearby_4(:)
  integer                             :: i, i_err_1, i_err_2, i_err_3, i_err_4, n, n1, n2, n3, n4
  integer, dimension(:), allocatable  :: i_elms_1, i_elms_2, i_elms_3, i_elms_4

  nullify(nodeout_1,nodeout_2,nodeout_3,nodeout_4)

  call RZ_minmax(node_list, element_list, i_elm, rmin, rmax, zmin, zmax)

  x1 = (/ rmin, zmin /)
  x2 = (/ rmax, zmin /)
  x3 = (/ rmax, zmax /)
  x4 = (/ rmin, zmax /)

  call quadtree_find_point(quadtree, x1, nodeout_1, i_err_1)
  call quadtree_find_point(quadtree, x2, nodeout_2, i_err_2)
  call quadtree_find_point(quadtree, x3, nodeout_3, i_err_3)
  call quadtree_find_point(quadtree, x4, nodeout_4, i_err_4)

  if (allocated(i_nearby)) deallocate(i_nearby)

  if (i_err_1 .eq. 0) then
    if (allocated(nodeout_1%contained)) then
       i_nearby_1 = pack( nodeout_1%contained(:)%quad_id,                  &
                         (nodeout_1%contained(:)%quad_id .gt. 0)           &
                      .and. (rmin .le. nodeout_1%contained(:)%vertex(3,1)) &
                      .and. (rmax .ge. nodeout_1%contained(:)%vertex(1,1)) &
                      .and. (zmin .le. nodeout_1%contained(:)%vertex(3,2)) &
                      .and. (zmax .ge. nodeout_1%contained(:)%vertex(1,2)) )
    endif
  endif
  if (i_err_2 .eq. 0) then
    if (allocated(nodeout_2%contained)) then
       i_nearby_2 = pack( nodeout_2%contained(:)%quad_id,                  &
                         (nodeout_2%contained(:)%quad_id .gt. 0)           & 
                      .and. (rmin .le. nodeout_2%contained(:)%vertex(3,1)) &
                      .and. (rmax .ge. nodeout_2%contained(:)%vertex(1,1)) &
                      .and. (zmin .le. nodeout_2%contained(:)%vertex(3,2)) &
                      .and. (zmax .ge. nodeout_2%contained(:)%vertex(1,2)) )
    endif
  endif
  if (i_err_3 .eq. 0) then
    if (allocated(nodeout_3%contained)) then
        i_nearby_3 = pack( nodeout_3%contained(:)%quad_id,                  &
                          (nodeout_3%contained(:)%quad_id .gt. 0)           &
                       .and. (rmin .le. nodeout_3%contained(:)%vertex(3,1)) &
                       .and. (rmax .ge. nodeout_3%contained(:)%vertex(1,1)) &
                       .and. (zmin .le. nodeout_3%contained(:)%vertex(3,2)) &
                       .and. (zmax .ge. nodeout_3%contained(:)%vertex(1,2)) )
    endif
  endif
  if (i_err_4 .eq. 0) then
    if (allocated(nodeout_4%contained)) then
      i_nearby_4 = pack( nodeout_4%contained(:)%quad_id,                  &
                        (nodeout_4%contained(:)%quad_id .gt. 0)           &
                     .and. (rmin .le. nodeout_4%contained(:)%vertex(3,1)) &
                     .and. (rmax .ge. nodeout_4%contained(:)%vertex(1,1)) &
                     .and. (zmin .le. nodeout_4%contained(:)%vertex(3,2)) &
                     .and. (zmax .ge. nodeout_4%contained(:)%vertex(1,2)) )
    endif
  endif
  
  allocate(i_nearby(size(i_nearby_1)+size(i_nearby_2)+size(i_nearby_3)+size(i_nearby_4)))

  i_nearby = [i_nearby_1, i_nearby_2, i_nearby_3, i_nearby_4]
  
end

!try to find all the elements in the box continaing i_elm
subroutine nearby_elements_quadtree(node_list, element_list, i_elm, i_nearby)
  implicit none
  type(type_node_list), intent(in)    :: node_list
  type(type_element_list), intent(in) :: element_list
  integer                             :: i_elm
  integer, allocatable, intent(out)   :: i_nearby(:)
  
  integer, dimension(n_vertex_max)    :: vertices
  real*8                              :: x1(2), x2(2), x3(2), x4(2)
  integer                             :: i, j, k, i_err
  integer, dimension(:), allocatable  :: i_elms_1, i_elms_2, i_elms_3, i_elms_4

  call interp_RZ(node_list, element_list, i_elm, 0.5d0, 0.0d0, x1(1), x1(2))
  call interp_RZ(node_list, element_list, i_elm, 1.0d0, 0.5d0, x2(1), x2(2))
  call interp_RZ(node_list, element_list, i_elm, 0.5d0, 1.0d0, x3(1), x3(2))
  call interp_RZ(node_list, element_list, i_elm, 0.0d0, 0.5d0, x4(1), x4(2))

  call elements_containing_point_quadtree(x1(1), x1(2), i_elms_1, .true.) 
  call elements_containing_point_quadtree(x2(1), x2(2), i_elms_2, .true.) 
  call elements_containing_point_quadtree(x3(1), x3(2), i_elms_3, .true.) 
  call elements_containing_point_quadtree(x4(1), x4(2), i_elms_4, .true.) 

  allocate(i_nearby(size(i_elms_1)+size(i_elms_2)+size(i_elms_3)+size(i_elms_4)))

  if (i_elm .eq. 1) then
    i_nearby = [[i_elm+1], i_elms_1, i_elms_2, i_elms_3, i_elms_4]
  else
    i_nearby = [[i_elm-1, i_elm+1], i_elms_1, i_elms_2, i_elms_3, i_elms_4]
  endif

end

recursive subroutine n_contained(node,n)
implicit none
type(quadtree_node), intent(inout) :: node        !< Root node
integer :: n, n1, n2, n3, n4

n1 = 0; n2 = 0; n3 = 0; n4 = 0
if (associated(node%children)) then
  call n_contained(node%children(1),n1)
  call n_contained(node%children(2),n2)
  call n_contained(node%children(3),n3)
  call n_contained(node%children(4),n4)

  n = n + n1 + n2 + n3 + n4
  write(*,'(A,5i4)') node%depth,' sum : ', n, n1, n2, n3, n4
else
  if (allocated(node%contained)) n = size(node%contained)
  write(*,*) node%depth,' n_contained : ',n
endif

end


subroutine elements_containing_point_quadtree(R, Z, i_elms, filter)
implicit none
real*8, intent(in)                              :: R, Z
integer, dimension(:), allocatable, intent(out) :: i_elms
logical, intent(in) , optional                  :: filter

type(quadtree_node), pointer :: nodeout     
integer                      :: i_err, i
logical                      :: do_filter = .true.
real*8                       :: tol = 1.d-6

if (present(filter)) then
  do_filter = filter
endif

call quadtree_find_point(quadtree, (/R, Z/), nodeout, i_err)

if (i_err .eq.0) then
  if (allocated(nodeout%contained)) then

    if (do_filter)  then

      i_elms = pack(nodeout%contained(:)%quad_id,  &
                     ((R .ge. nodeout%contained(:)%vertex(1,1) - tol) &
                .and. (R .le. nodeout%contained(:)%vertex(3,1) + tol) & 
                .and. (Z .ge. nodeout%contained(:)%vertex(1,2) - tol) &
                .and. (Z .le. nodeout%contained(:)%vertex(3,2) + tol)))

    else
      allocate(i_elms(size(nodeout%contained)))
      i_elms = nodeout%contained(:)%quad_id
    endif
    
  endif
endif
end

subroutine plot_element(node_list,element_list,i_elm)
implicit none
type(type_node_list), intent(in)    :: node_list
type(type_element_list), intent(in) :: element_list
integer                             :: i_elm

integer :: i, j, ip, np, iuv, inode_0, inode_p
real*8  :: xx_0(2), uv_0(2), huv_0, xx_p(2), uv_p(2), huv_p, s
real*8, allocatable :: xb(:,:), xp(:,:)

np = 5

allocate(xb(4,np), xp(np,2))

do i=1,4                         ! over the 4 edges

    iuv = mod(i+1,2)+1             ! the direction vector corresponding to this edge (i)
 
    inode_0 = element_list%element(i_elm)%vertex(i)
    xx_0    = node_list%node(inode_0)%x(1,1,:)
    uv_0    = node_list%node(inode_0)%x(1,iuv+1,:)
    huv_0   = element_list%element(i_elm)%size(i,iuv+1)
 
    ip      = mod(i,4)+1
    inode_p = element_list%element(i_elm)%vertex(ip)
    xx_p    = node_list%node(inode_p)%x(1,1,:)
    uv_p    = node_list%node(inode_p)%x(1,iuv+1,:)
    huv_p   = element_list%element(i_elm)%size(ip,iuv+1)
  
    xb      = 0.d0
 
    do j=1,np
      s = (real(j-1,8)/real(np-1,8))
      xb(1,:) = xx_0
      xb(2,:) = xx_0+uv_0*huv_0
      xb(3,:) = xx_p+uv_p*huv_p
      xb(4,:) = xx_p
 
      call bezier_1d(n_dim, s, xb, xp(j,:))
    enddo
 
    call lplot6(1,1,xp(1:np,1),xp(1:np,2),-np,' ')
  
  enddo
end

end
