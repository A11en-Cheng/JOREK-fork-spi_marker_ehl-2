!> mod_geometry contains variables and procedures for solving
!> basic / common geometrical problems such as intersections
module mod_geometry
implicit none

private
public :: compute_plane_edge_length
public :: compute_global_cart_coord_plane_points
public :: compute_plane_line_intersection_cart
public :: define_plane_from_half_angles
public :: define_vertex_spherical_coord

!> Interfaces -------------------------------------------------
interface compute_plane_edge_length
  module procedure compute_plane_edge_length_r8
end interface compute_plane_edge_length
interface compute_global_cart_coord_plane_points
  module procedure compute_global_cart_coord_plane_points_r8
end interface compute_global_cart_coord_plane_points

interface compute_plane_line_intersection_cart
  module procedure compute_plane_line_intersect_cart_points_r8
end interface compute_plane_line_intersection_cart

interface define_plane_from_half_angles
  module procedure define_direction_origin_plane_from_half_angles
  module procedure define_direction_plane_from_half_angles
  module procedure define_standard_plane_from_half_angles
end interface define_plane_from_half_angles

interface define_vertex_spherical_coord
  module procedure define_vertex_spherical_coord_standard
  module procedure define_vertex_spherical_coord_origin
end interface define_vertex_spherical_coord

contains
!> Procedures -------------------------------------------------
!> compute the width (P1->P2) and the height (P1->P3) of a plane
!> defined by the points
!> P1 -> s -> P2
!> |           |
!> v           v
!> t           t
!> |           |
!> v           v
!> P3 -> s -> P4
!> inputs:
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> origin
!>                                 pp(:,2) -> s=1,t=0 node
!>                                 pp(:,3) -> s=0,t=1 node
!> outputs:
!>   lengths: (real8)(2) 1) width (P1->P2) and 2) height
!>                       (P3->P4) of a plane
subroutine compute_plane_edge_length_r8(pp,lengths)
  implicit none
  real*8,dimension(3,3),intent(in) :: pp
  real*8,dimension(2),intent(out)  :: lengths
  lengths = [norm2(pp(:,2)-pp(:,1)),norm2(pp(:,3)-pp(:,1))]
end subroutine compute_plane_edge_length_r8

!> compute the global cartesian coordinates of a point on 
!> a plane given the nodes defining the plane and the point
!> local coordinates (double precision)
!> P1 -> s -> P2
!> |           |
!> v           v
!> t           t
!> |           |
!> v           v
!> P3 -> s -> P4
!> inputs:
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> origin
!>                                 pp(:,2) -> s=1,t=0 node
!>                                 pp(:,3) -> s=0,t=1 node
!>   st: (real8)(2) position of a point on a plane in the plane
!>                  local coordinates
!> outpus:
!>   x: (real8)(3) position of a point on a plane in the global
!>                 cartesian coordinates
pure function compute_global_cart_coord_plane_points_r8(pp,st) result(x)
  implicit none
  real*8,dimension(2),intent(in)   :: st
  real*8,dimension(3,3),intent(in) :: pp
  real*8,dimension(3)              :: x
  x = pp(:,1)*(1.d0-st(1)-st(2))+pp(:,2)*st(1)+pp(:,3)*st(2)
end function compute_global_cart_coord_plane_points_r8

!> compute_plane_line_intersect_cart_points_r8 computes the 
!> intersection between a line and a plane in double precision
!> the plane must be defined by three points while the line
!> must be defined by two points. The function returns the
!> intersection coordinates in the local plane (s,t) and lines
!> (q) coordinate system. The intersection is found if the
!> local coordinates of the intersection points are in [0,1].
!> The north-east node is the origin of the plane
!> P1 -> s -> P2
!> |           |
!> v           v
!> t           t
!> |           |
!> v           v
!> P3 -> s -> P4 
!> inputs:
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> origin
!>                                 pp(:,2) -> s=1,t=0 node
!>                                 pp(:,3) -> s=0,t=1 node
!>   pl: (real8)(3,2) points defining a line in cartesian
!>                    coordinates: pl(:,1) -> origin
!>                                 pl(:,2) -> q=1 node
!> outputs:
!>   intersect: (bool) if true intersection found
!>   stq:       (real8)(3) plane(s,t) and line (q) local
!>              coordinates of the intersection
subroutine compute_plane_line_intersect_cart_points_r8(pp,pl,&
intersect,stq)
  use mod_math_operators, only: solve_3x3_linear_problem
  implicit none
  real*8,dimension(3,3),intent(in) :: pp
  real*8,dimension(3,2),intent(in) :: pl
  logical,intent(out)              :: intersect
  real*8,dimension(3),intent(out)  :: stq
  real*8,dimension(3,3) :: matrix
  
  matrix(:,1) = pp(:,2)-pp(:,1); matrix(:,2) = pp(:,3)-pp(:,1);
  matrix(:,3) = pl(:,1)-pl(:,2);
  call  solve_3x3_linear_problem(matrix,pl(:,1)-pp(:,1),stq)
  intersect = (all(stq.ge.0.d0)).and.(all(stq.le.1.d0))
end subroutine compute_plane_line_intersect_cart_points_r8

!> define a plane vertices given the half width, half height,
!> the spherical coordinates of the plane mid point in
!> the origin coordinate system and the origin coordinates
!> inputs:
!>   mirror_xy:   (integer)(2) mirror plane vertices w.r.t.
!>                the x and y axes.
!>                1) if 1 mirror w.r.t. the x-axis
!>                2) if 1 mirror w.r.t. the y-axis
!>
!>   half_angles: (real8)(2) 1:half_width 2:half_height
!>                3: orientation - alignment of the 
!>                P1,P2,P3 points w.r.t. the Z axis
!> 
!>   rthetaphi:   (real8)(3) spherical coordinates of the
!>                plane midpoin: rthetaphi(1): distance
!>                               rthetaphi(2): colatitude
!>                               rthetaphi(3): azimuth
!>   origin:      (real8)(3) position of the sphere origin
!> outputs
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> P1
!>                                 pp(:,2) -> P2
!>                                 pp(:,3) -> P3
subroutine define_direction_origin_plane_from_half_angles(&
mirror_xy,half_angles,rthetaphi,origin,pp)
  implicit none
  integer,dimension(2),intent(in) :: mirror_xy
  real*8,dimension(3),intent(in)  :: half_angles
  real*8,dimension(3),intent(in)  :: rthetaphi,origin
  !> outputs:
  real*8,dimension(3,3),intent(out) :: pp
  !> variables
  integer :: ii
  !> compute directional plane vertices
  call define_direction_plane_from_half_angles(&
  mirror_xy,half_angles,rthetaphi,pp)
  !> translate vertices
  do ii=1,3
    pp(:,ii) = origin + pp(:,ii)
  enddo
end subroutine define_direction_origin_plane_from_half_angles


!> define a plane vertices given the half width, half height
!> and the spherical coordinates of the plane mid point in
!> the origin coordinate system
!> inputs:
!>   mirror_xy:   (integer)(2) mirror plane vertices w.r.t.
!>                the x and y axes.
!>                1) if 1 mirror w.r.t. the x-axis
!>                2) if 1 mirror w.r.t. the y-axis
!>
!>   half_angles: (real8)(2) 1:half_width 2:half_height
!>                3: orientation - alignment of the 
!>                P1,P2,P3 points w.r.t. the Z axis
!> 
!>   rthetaphi:   (real8)(3) spherical coordinates of the
!>                plane midpoin: rthetaphi(1): distance
!>                               rthetaphi(2): colatitude
!>                               rthetaphi(3): azimuth
!> outputs
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> P1
!>                                 pp(:,2) -> P2
!>                                 pp(:,3) -> P3
subroutine define_direction_plane_from_half_angles(mirror_xy,&
half_angles,rthetaphi,pp)
  use mod_coordinate_transforms, only: vectors_spherical_to_cartesian
  implicit none
  !> inputs:
  integer,dimension(2),intent(in) :: mirror_xy
  real*8,dimension(3),intent(in)  :: half_angles
  real*8,dimension(3),intent(in)  :: rthetaphi
  !> outputs:
  real*8,dimension(3,3),intent(out) :: pp  !> compute rotation transform
  !> compute plane vertices
  call define_standard_plane_from_half_angles(mirror_xy,half_angles,pp)
  !> rotate the points defining the plane in the cartesian reference
  call vectors_spherical_to_cartesian(3,rthetaphi,pp)
end subroutine define_direction_plane_from_half_angles

!> define a plane vertices given the half width and half height
!> angles in the standard reference systme: origin = (/0,0,0/)
!> and plane normal z=(/0,0,1/)
!> P1 -> width -> P2
!>   |           |
!>   v           v
!> height      height
!>   |           |
!>   v           v
!> P3 -> width -> P4 
!> inputs:
!>   mirror_xy:   (integer)(2) mirror plane vertices w.r.t.
!>                the x and y axes.
!>                1) if 1 mirror w.r.t. the x-axis
!>                2) if 1 mirror w.r.t. the y-axis
!>
!>   half_angles: (real8)(2) 1:half_width 2:half_height
!>                3: orientation - alignment of the 
!>                P1,P2,P3 points w.r.t. the Z axis
!>                 
!> outputs
!>   pp: (real8)(3,3) points defining a plane in cartesian 
!>                    coordinates: pp(:,1) -> P1
!>                                 pp(:,2) -> P2
!>                                 pp(:,3) -> P3
subroutine define_standard_plane_from_half_angles(mirror_xy,half_angles,pp)
  use mod_coordinate_transforms, only: mirror_around_cart_x
  use mod_coordinate_transforms, only: mirror_around_cart_y
  use mod_coordinate_transforms, only: rotate_vectors_cart_z
  implicit none
  !> inputs:
  integer,dimension(2),intent(in) :: mirror_xy
  real*8,dimension(3),intent(in)  :: half_angles
  !> outputs:
  real*8,dimension(3,3),intent(out) :: pp
  !> variables
  real*8,dimension(2)   :: tan_angles
  !> compute plane vertex coordinates
  tan_angles = tan(half_angles(1:2))
  pp(:,1) = (/-tan_angles(1),tan_angles(2),1.d0/)
  pp(:,2) = (/tan_angles(1),tan_angles(2),1.d0/)
  pp(:,3) = (/-tan_angles(1),-tan_angles(2),1.d0/)
  !> change the orientation of the plane points
  call rotate_vectors_cart_z(3,half_angles(3),pp)
  !> mirror w.r.t. x and y axis if required
  if(mirror_xy(1).eq.1) call mirror_around_cart_x(3,pp)
  if(mirror_xy(2).eq.1) call mirror_around_cart_y(3,pp)
end subroutine define_standard_plane_from_half_angles

!> define vertex given its spherical coordinates and an origin
!> inputs:
!>   rthetaphi:   (real8)(3) spherical coordinates of the
!>                plane midpoin: rthetaphi(1): distance
!>                               rthetaphi(2): colatitude
!>                               rthetaphi(3): azimuth
!>   origin:      (real8)(3) position of the sphere origin
!> outputs: 
!>   vertex:      (real8)(3) vertex
subroutine define_vertex_spherical_coord_origin(rthetaphi,origin,vertex)
  implicit none
  !> inputs:
  real*8,dimension(3),intent(in)  :: rthetaphi,origin
  !> outputs:
  real*8,dimension(3),intent(out) :: vertex
  !> compute vertex
  call define_vertex_spherical_coord_standard(rthetaphi,vertex)
  vertex = origin + vertex
end subroutine define_vertex_spherical_coord_origin

!> define vertex its the spherical coordinates
!> inputs:
!>   rthetaphi:   (real8)(3) spherical coordinates of the
!>                plane midpoin: rthetaphi(1): distance
!>                               rthetaphi(2): colatitude
!>                               rthetaphi(3): azimuth
!> outputs: 
!>   vertex:      (real8)(3) vertex
subroutine define_vertex_spherical_coord_standard(rthetaphi,vertex)
  implicit none
  !> inputs:
  real*8,dimension(3),intent(in)  :: rthetaphi
  !> outputs:
  real*8,dimension(3),intent(out) :: vertex
  !> compute vertex
  vertex = rthetaphi(1)*(/cos(rthetaphi(3)),sin(rthetaphi(3)),cos(rthetaphi(2))/)
  vertex(1:2) = vertex(1:2)*sin(rthetaphi(2))
end subroutine define_vertex_spherical_coord_standard

!>-------------------------------------------------------------
end module mod_geometry
