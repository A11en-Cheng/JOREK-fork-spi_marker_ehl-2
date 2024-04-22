!> mod_geometry_test contains variables and procedures
!> used for testing the routines contained in mod_geometry
module mod_geometry_test
use fruit
use constants, only: PI,TWOPI
implicit none

private
public :: run_fruit_geometry

!> Variables ----------------------------------------------
integer,parameter              :: n_planes=11
integer,parameter              :: n_nodes_plane=3
integer,parameter              :: n_lines_per_test=23
integer,parameter              :: n_tests=5
integer,parameter              :: n_lines=n_tests*n_lines_per_test
integer,dimension(2),parameter :: mirror_xy_interval=(/0,1/) 
real*8,parameter               :: tol_real8=7.5d-10
real*8,parameter               :: radius=2.d2
real*8,dimension(2),parameter  :: len_interval=(/2.d0,5.d2/)
real*8,dimension(2),parameter  :: cos_theta_interval=(/-1.d0,1.d0/)
real*8,dimension(2),parameter  :: phi_interval=(/0.d0,TWOPI/)
real*8,dimension(3),parameter  :: half_angle_lowbnd=(/PI/21.d0,PI/11.d0,0d0/)
real*8,dimension(3),parameter  :: half_angle_uppbnd=(/PI/2.d0,TWOPI/5.d0,TWOPI/)
real*8,dimension(3),parameter  :: ppxyz_lowbound=(/-1.d0,2.d0,-4.d1/)
real*8,dimension(3),parameter  :: ppxyz_uppbound=(/5.d0,2.5d1,-3.2d1/)
real*8,dimension(3),parameter  :: plxyz_lowbound=(/-3.d0,1.d1,4.d1/)
real*8,dimension(3),parameter  :: plxyz_uppbound=(/7.d0,5.d0,3.2d1/)
real*8,dimension(3),parameter  :: pl0_lowbound=(/-2.1d1,4.5d0,-9.3d1/)
real*8,dimension(3),parameter  :: pl0_uppbound=(/8.d1,1.1d1,5.2d1/)
real*8,dimension(5,n_planes),parameter :: zeros_5Xn_planes_r8=0.d0
real*8,dimension(4,n_planes),parameter :: ones_4Xn_planes_r8=1.d0
real*8,dimension(2,n_tests),parameter :: st_lowbound=&
       reshape((/0.d0,0.d0,-3.d1,-4.d0,-5.1d1,1.1d0,&
       1.2d0,-8.d1,1.1d0,1.05d0/),shape(st_lowbound))
real*8,dimension(2,n_tests),parameter :: st_uppbound=&
       reshape((/1.d0,1.d0,-5.d-2,-4.d-3,-1.d-1,5.d1,&
       9.2d1,-1.d-3,1.1d1,3.05d0/),shape(st_uppbound))
logical,dimension(n_lines,n_planes)        :: intersect_sol
integer,dimension(2,n_planes)              :: mirror_xy_sol
real*8,dimension(3,n_planes)               :: half_angles_sol
real*8,dimension(3,n_planes)               :: origins_sol,rthetaphi_sol
real*8,dimension(3,n_nodes_plane,n_planes) :: pp_sol
real*8,dimension(3,n_lines,n_planes)       :: stq_sol,pos_intersect_sol
real*8,dimension(3,2,n_lines,n_planes)     :: pl_sol
!> Interfaces ---------------------------------------------

contains
!> Fruit basket -------------------------------------------
!> basket having all set-up, tests and tearing-down routines
subroutine run_fruit_geometry()
  implicit none
  write(*,'(/A)') "  ... setting-up: geometry tests"
  call setup_plane_line_intersection_test
  call setup_plane_definition_from_half_angle
  write(*,'(/A)') "  ... running: geometry tests"
  call run_test_case(test_compute_test_line_intersect_cart_points,&
  'test_compute_test_line_intersect_cart_points')
  call run_test_case(test_define_standard_plane_from_half_angle,&
  'test_define_standard_plane_from_half_angle')
  call run_test_case(test_define_direction_plane_from_half_angle,&
  'test_define_direction_plane_from_half_angle')
  call run_test_case(test_define_direction_plane_from_half_angle_origin,&
  'test_define_direction_plane_from_half_angle_origin')
  call run_test_case(test_vertex_spherical_coord_standard,&
  'test_vertex_spherical_coord_standard')
  call run_test_case(test_vertex_spherical_coord_origin,&
  'test_vertex_spherical_coord_origin')
  call run_test_case(test_compute_edge_length,&
  'test_compute_edge_length')
  write(*,'(/A)') "  ... tearing-down: geometry tests"
end subroutine run_fruit_geometry

!> Set-up and tear-down -----------------------------------
!> set-up all the variables for testing the plane-line intersections
subroutine setup_plane_line_intersection_test()
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_gnu_rng,  only: gnu_rng_interval
  implicit none
  !> variables
  integer :: ii,jj,kk,id
  real*8              :: length
  real*8,dimension(2) :: st
  real*8,dimension(3) :: pos

  do ii=1,n_planes
    !> generate plane nodes
    do jj=1,n_nodes_plane
      call gnu_rng_interval(3,ppxyz_lowbound,&
      ppxyz_uppbound,pp_sol(:,jj,ii))
    enddo
    !> generate lines
    do jj=1,n_tests
      do kk=1,n_lines_per_test
        !> generate first line node
        id = (jj-1)*n_lines_per_test+kk
        call gnu_rng_interval(3,plxyz_lowbound,&
        plxyz_uppbound,pl_sol(:,1,id,ii))
        !> generate second line node
        call gnu_rng_interval(2,st_lowbound(:,jj),&
        st_uppbound(:,jj),st)
        call gnu_rng_interval(len_interval,length)
        pos = compute_global_cart_coord_plane_points(pp_sol(:,:,ii),st)
        pl_sol(:,2,id,ii) = pl_sol(:,1,id,ii) + &
        (pos-pl_sol(:,1,id,ii))*length
        !> check if intersection
        intersect_sol(id,ii)=(all(st.ge.0.d0).and.all(st.le.1.d0))
        !> store intersection position values
        pos_intersect_sol(:,id,ii) = pos
        stq_sol(:,id,ii) = (/st(1),st(2),norm2(pos-pl_sol(:,1,id,ii))/&
        norm2(pl_sol(:,2,id,ii)-pl_sol(:,1,id,ii))/)
      enddo
    enddo
  enddo
end subroutine setup_plane_line_intersection_test

!> set-up all the variables for testing the routines
!> generating planes in cartesian coordinates
subroutine setup_plane_definition_from_half_angle()
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(3) :: rand

  !> generate plane coordinates, widths and heights
  call gnu_rng_interval(2,n_planes,mirror_xy_interval,mirror_xy_sol)
  do ii=1,n_planes
    call random_number(rand)
    rthetaphi_sol(:,ii) = sample_uniform_sphere(radius,&
    cos_theta_interval,phi_interval,rand)
    call gnu_rng_interval(3,pl0_lowbound,pl0_uppbound,origins_sol(:,ii))
    call gnu_rng_interval(3,half_angle_lowbnd,half_angle_uppbnd,half_angles_sol(:,ii))
  enddo
end subroutine setup_plane_definition_from_half_angle

!> Tests --------------------------------------------------
!> test the procedure for computing vertex given some spherical
!> coordinates
subroutine test_vertex_spherical_coord_standard()
  use constants,    only: TWOPI
  use mod_geometry, only: define_vertex_spherical_coord
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(3)          :: vertex_test
  real*8,dimension(3,n_planes) :: rthetaphi_test
  !> execute tests
  do ii=1,n_planes
    call define_vertex_spherical_coord(rthetaphi_sol(:,ii),vertex_test)
    rthetaphi_test(1,ii) = norm2(vertex_test)
    rthetaphi_test(2:3,ii) = (/acos(vertex_test(3)/rthetaphi_test(1,ii)),&
    atan2(vertex_test(2),vertex_test(1))/)
  enddo
  where(rthetaphi_test(3,:).lt.0.d0) rthetaphi_test(3,:) = TWOPI + rthetaphi_test(3,:)
  call assert_equals(rthetaphi_sol,rthetaphi_test,3,n_planes,tol_real8,&
  "Error vertex spherical coordinates standard: spherical coord. mismatch!")
end subroutine test_vertex_spherical_coord_standard

!> test the procedure for computing the vertex give some spherical
!> coordinates and origin vertices
subroutine test_vertex_spherical_coord_origin()
  use mod_geometry, only: define_vertex_spherical_coord
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(3)          :: vertex_test
  real*8,dimension(3,n_planes) :: rthetaphi_test
  !> execute tests
  do ii=1,n_planes
    call define_vertex_spherical_coord(rthetaphi_sol(:,ii),origins_sol(:,ii),vertex_test)
    vertex_test = vertex_test - origins_sol(:,ii)
    rthetaphi_test(1,ii) = norm2(vertex_test)
    rthetaphi_test(2:3,ii) = (/acos(vertex_test(3)/rthetaphi_test(1,ii)),&
    atan2(vertex_test(2),vertex_test(1))/)
  enddo
  where(rthetaphi_test(3,:).lt.0.d0) rthetaphi_test(3,:) = TWOPI + rthetaphi_test(3,:)
  call assert_equals(rthetaphi_sol,rthetaphi_test,3,n_planes,tol_real8,&
  "Error vertex spherical coordinates origins: spherical coord. mismatch!")
end subroutine test_vertex_spherical_coord_origin

!> test the procedure for finding the intersection between
!> a line and a plane
subroutine test_compute_test_line_intersect_cart_points()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_geometry, only: compute_plane_line_intersection_cart
  implicit none
  integer :: ii,jj
  logical :: intersect
  real*8,dimension(3) :: stq,pos

  do ii=1,n_planes
    do jj=1,n_lines
      !> compute intersection point in local and global coordinates
      call compute_plane_line_intersection_cart(&
      pp_sol(:,:,ii),pl_sol(:,:,jj,ii),intersect,stq)
      pos = compute_global_cart_coord_plane_points(&
      pp_sol(:,:,ii),stq(1:2))
      !> check solutions
      call assert_equals_rel_error(3,stq_sol(:,jj,ii),stq,tol_real8,&
      "Error compute plane line intersection cart: stq mismatch!")
      call assert_equals_rel_error(3,pos_intersect_sol(:,jj,ii),pos,tol_real8,&
      "Error compute plane line intersection cart: pos intersect mismatch!")
      call assert_equals(intersect_sol(jj,ii),intersect,&
      "Error compute plane line intersection cart: intersect mismatch!")
    enddo
  enddo
end subroutine test_compute_test_line_intersect_cart_points

!> Test define standard plane from half angles 
subroutine test_define_standard_plane_from_half_angle()
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(3)            :: normal,origin
  real*8,dimension(n_planes)     :: plane_angle
  real*8,dimension(4,n_planes)   :: distance_loc
  real*8,dimension(2,n_planes)   :: half_width_angle_loc,half_height_angle_loc
  real*8,dimension(5,n_planes)   :: orthogonal
  real*8,dimension(3,3)          :: plane
  !> initialise quantities
  normal = (/0.d0,0.d0,1.d0/); origin = (/0.d0,0.d0,0.d0/);
  !> compute plane and quantities to test
  do ii=1,n_planes
    call define_plane_from_half_angles(mirror_xy_sol(:,ii),half_angles_sol(:,ii),plane)
    call extract_values_for_plane_definition_test(normal,origin,mirror_xy_sol(:,ii),&
    plane,half_height_angle_loc(:,ii),half_width_angle_loc(:,ii),&
    distance_loc(:,ii),orthogonal(:,ii),plane_angle(ii))
  enddo
  !> check values 
  call assert_equals(half_angles_sol(1,:),half_width_angle_loc(1,:),n_planes,&
  tol_real8,"Error define standard plane from half angles: width angle 1 mismatch!")
  call assert_equals(half_angles_sol(1,:),half_width_angle_loc(2,:),n_planes,&
  tol_real8,"Error define standard plane from half angles: width angle 2 mismatch!")
  call assert_equals(half_angles_sol(2,:),half_height_angle_loc(1,:),n_planes,&
  tol_real8,"Error define standard plane from half angles: height angle 1 mismatch!")
  call assert_equals(half_angles_sol(2,:),half_height_angle_loc(2,:),n_planes,&
  tol_real8,"Error define standard plane from half angles: height angle 2 mismatch!")
  call assert_equals(ones_4Xn_planes_r8,distance_loc,4,n_planes,&
  tol_real8,"Error define standard plane from half angles: distance is not unity!")
  call assert_equals(zeros_5Xn_planes_r8,orthogonal,5,n_planes,&
  tol_real8,"Error define standard plane from half angles: axis not orthogonal!")
  call assert_equals(half_angles_sol(3,:),plane_angle,n_planes,&
  tol_real8,"Error define standard plane from half angles: plane angle mismatch!")
end subroutine test_define_standard_plane_from_half_angle

!> Test define plane from half angles with direction and distance
subroutine test_define_direction_plane_from_half_angle()
  use mod_coordinate_transforms, only: spherical_colatitude_to_cartesian
  use mod_geometry,              only: define_plane_from_half_angles
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(3)            :: normal,origin
  real*8,dimension(n_planes)     :: plane_angle
  real*8,dimension(4,n_planes)   :: distance_loc
  real*8,dimension(2,n_planes)   :: half_width_angle_loc,half_height_angle_loc
  real*8,dimension(5,n_planes)   :: orthogonal
  real*8,dimension(3,3)          :: plane
  !> initialise quantities
  origin = (/0.d0,0.d0,0.d0/);
  !> compute plane and quantities to test
  do ii=1,n_planes
    normal = spherical_colatitude_to_cartesian(rthetaphi_sol(:,ii),origin)
    distance_loc(:,ii) = rthetaphi_sol(1,ii)
    call define_plane_from_half_angles(mirror_xy_sol(:,ii),&
    half_angles_sol(:,ii),rthetaphi_sol(:,ii),plane)
    call extract_values_for_plane_definition_test(normal,origin,mirror_xy_sol(:,ii),&
    plane,half_height_angle_loc(:,ii),half_width_angle_loc(:,ii),&
    distance_loc(:,ii),orthogonal(:,ii),plane_angle(ii))
  enddo
  !> check values 
  call assert_equals(half_angles_sol(1,:),half_width_angle_loc(1,:),n_planes,&
  tol_real8,"Error define direction plane from half angles: width angle 1 mismatch!")
  call assert_equals(half_angles_sol(1,:),half_width_angle_loc(2,:),n_planes,&
  tol_real8,"Error define direction plane from half angles: width angle 2 mismatch!")
  call assert_equals(half_angles_sol(2,:),half_height_angle_loc(1,:),n_planes,&
  tol_real8,"Error define direction plane from half angles: height angle 1 mismatch!")
  call assert_equals(half_angles_sol(2,:),half_height_angle_loc(2,:),n_planes,&
  tol_real8,"Error define direction plane from half angles: height angle 2 mismatch!")
  call assert_equals(distance_loc,distance_loc,4,n_planes,&
  tol_real8,"Error define direction plane from half angles: distance is not unity!")
  call assert_equals(zeros_5Xn_planes_r8,orthogonal,5,n_planes,&
  tol_real8,"Error define direction plane from half angles: axis not orthogonal!")
  call assert_equals(half_angles_sol(3,:),plane_angle,n_planes,&
  tol_real8,"Error define direction plane from half angles: plane angle mismatch!")
end subroutine test_define_direction_plane_from_half_angle

!> Test define plane from half angles with direction, distance and origin
subroutine test_define_direction_plane_from_half_angle_origin()
  use mod_coordinate_transforms, only: spherical_colatitude_to_cartesian
  use mod_geometry,              only: define_plane_from_half_angles
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(3)            :: normal,origin
  real*8,dimension(n_planes)     :: plane_angle
  real*8,dimension(4,n_planes)   :: distance_loc
  real*8,dimension(2,n_planes)   :: half_width_angle_loc,half_height_angle_loc
  real*8,dimension(5,n_planes)   :: orthogonal
  real*8,dimension(3,3)          :: plane
  !> compute plane and quantities to test
  do ii=1,n_planes
    origin = origins_sol(:,ii)
    normal = spherical_colatitude_to_cartesian(rthetaphi_sol(:,ii),origin)
    distance_loc(:,ii) = rthetaphi_sol(1,ii)
    call define_plane_from_half_angles(mirror_xy_sol(:,ii),half_angles_sol(:,ii),&
    rthetaphi_sol(:,ii),origin,plane)
    call extract_values_for_plane_definition_test(normal,origin,mirror_xy_sol(:,ii),&
    plane,half_height_angle_loc(:,ii),half_width_angle_loc(:,ii),&
    distance_loc(:,ii),orthogonal(:,ii),plane_angle(ii))
  enddo
  !> check values 
  call assert_equals(half_angles_sol(1,:),half_width_angle_loc(1,:),n_planes,&
  tol_real8,"Error define direction plane from half angles origin: width angle 1 mismatch!")
  call assert_equals(half_angles_sol(1,:),half_width_angle_loc(2,:),n_planes,&
  tol_real8,"Error define direction plane from half angles origin: width angle 2 mismatch!")
  call assert_equals(half_angles_sol(2,:),half_height_angle_loc(1,:),n_planes,&
  tol_real8,"Error define direction plane from half angles origin: height angle 1 mismatch!")
  call assert_equals(half_angles_sol(2,:),half_height_angle_loc(2,:),n_planes,&
  tol_real8,"Error define direction plane from half angles origin: height angle 2 mismatch!")
  call assert_equals(distance_loc,distance_loc,4,n_planes,&
  tol_real8,"Error define direction plane from half angles origin: distance is not unity!")
  call assert_equals(zeros_5Xn_planes_r8,orthogonal,5,n_planes,&
  tol_real8,"Error define direction plane from half angles origin: axis not orthogonal!")
  call assert_equals(half_angles_sol(3,:),plane_angle,n_planes,&
  tol_real8,"Error define direction plane from half angles origin: plane angle mismatch!")
end subroutine test_define_direction_plane_from_half_angle_origin

!> test if the lenght of the plane sides are computed correctly 
subroutine test_compute_edge_length()
  use mod_geometry, only: compute_plane_edge_length
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(2,n_planes) :: solution_length,test_length
  !> compute test and solution lengths
  do ii=1,n_planes
    call compute_plane_edge_length(pp_sol(:,:,ii),test_length(:,ii))
    solution_length(:,ii) = [norm2(pp_sol(:,2,ii)-pp_sol(:,1,ii)),&
    norm2(pp_sol(:,3,ii)-pp_sol(:,1,ii))]
  enddo
  !> test results 
  call assert_equals(solution_length,test_length,2,n_planes,&
  tol_real8,"Error computing length of plane edges: edge length mismatch!")
end subroutine test_compute_edge_length

!> Tools --------------------------------------------------
!> extract values for testing the definition of plane
subroutine extract_values_for_plane_definition_test(n_vector,origin,mirror_xy,&
plane,half_height_angle,half_width_angle,distance,orthogonal,plane_angle)
  use constants,                 only: TWOPI
  use mod_coordinate_transforms, only: cartesian_to_spherical_colatitude
  use mod_coordinate_transforms, only: vectors_cartesian_to_spherical
  use mod_coordinate_transforms, only: mirror_around_cart_y
  use mod_coordinate_transforms, only: mirror_around_cart_x
  use mod_geometry, only: compute_global_cart_coord_plane_points
  implicit none
  !> inputs
  integer,dimension(2),intent(in)  :: mirror_xy
  real*8,dimension(3),intent(in)   :: origin,n_vector
  real*8,dimension(3,3),intent(in) :: plane 
  !> outputs:
  real*8,intent(out)              :: plane_angle
  real*8,dimension(2),intent(out) :: half_height_angle
  real*8,dimension(2),intent(out) :: half_width_angle
  real*8,dimension(4),intent(out) :: distance
  real*8,dimension(5),intent(out) :: orthogonal
  !> variables
  integer :: ii
  real*8  :: norm_n_origin_vector
  real*8,dimension(3)   :: A1,A2,B1,B2,C,P4,normal,rthetachi
  real*8,dimension(3,4) :: plane_test

  !> initialise data
  normal = n_vector - origin
  norm_n_origin_vector = norm2(normal)
  normal = normal/norm_n_origin_vector
  P4 = compute_global_cart_coord_plane_points(plane,(/1.d0,1.d0/))
  A1 = 5.d-1*(plane(:,2)+plane(:,1));
  A2 = 5.d-1*(P4+plane(:,3))
  B1 = 5.d-1*(plane(:,1)+plane(:,3))
  B2 = 5.d-1*(plane(:,2)+P4)
  !> extract unsigned half height and half width
  half_height_angle(1) = asin(norm2(A1-n_vector)/norm2(A1-origin))
  half_height_angle(2) = asin(norm2(A2-n_vector)/norm2(A2-origin))
  half_width_angle(1)  = asin(norm2(B1-n_vector)/norm2(B1-origin))
  half_width_angle(2)  = asin(norm2(B2-n_vector)/norm2(B2-origin))
  !> extract distances
  distance(1) = dot_product(A1-origin,normal)
  distance(2) = dot_product(A2-origin,normal)
  distance(3) = dot_product(B1-origin,normal)
  distance(4) = dot_product(B2-origin,normal) 
  !> extract orthogonality
  orthogonal(1) = dot_product((A1-n_vector)/norm2(A1-n_vector),normal)
  orthogonal(2) = dot_product((A2-n_vector)/norm2(A2-n_vector),normal)
  orthogonal(3) = dot_product((B1-n_vector)/norm2(B1-n_vector),normal)
  orthogonal(4) = dot_product((B2-n_vector)/norm2(B2-n_vector),normal)
  orthogonal(5) = dot_product((B2-B1)/norm2(B2-B1),(A2-A1)/norm2(A2-A1))
  !> compute the normal spherical coordinates
  rthetachi = cartesian_to_spherical_colatitude(n_vector,origin)
  !> anti-transform the plane nodes
  do ii=1,size(plane,dim=2)
    plane_test(:,ii) = plane(:,ii)-origin;
  enddo
  plane_test(:,4) = P4-origin
  call vectors_cartesian_to_spherical(size(plane_test,dim=2),rthetachi,plane_test)
  !> anti-mirror the plane position if required
  if(mirror_xy(2).eq.1) call mirror_around_cart_y(size(plane_test,dim=2),plane_test)
  if(mirror_xy(1).eq.1) call mirror_around_cart_x(size(plane_test,dim=2),plane_test)
  !> compute the rotation angle
  B1 = 5d-1*(plane_test(:,2)+plane_test(:,4)); 
  C = 2.5d-1*sum(plane_test,dim=2); B1 = B1-C;
  plane_angle = atan2(-B1(2),B1(1))
  if(plane_angle.lt.0d0) plane_angle = TWOPI+plane_angle
end subroutine extract_values_for_plane_definition_test

!>---------------------------------------------------------
end module mod_geometry_test

