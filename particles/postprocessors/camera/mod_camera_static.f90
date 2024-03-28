!> mod_camera_static contains all variables and 
!> procedures defining a camera with static properties
!> (do not vary with time)
module mod_camera_static
use mod_camera, only: camera
implicit none

private
public :: camera_static

!> Variables and type definitions ----------------------------
type,extends(camera),abstract :: camera_static
 contains
  procedure,pass(camera_inout) :: generate_points_on_lens_pdf => &
  generate_points_on_lens_static
  procedure,pass(camera_inout) :: find_ray_image_plane_intersection => &
  find_ray_image_plane_intersection_static
  procedure,pass(camera_inout) :: define_image_plane_pixel_size
end type camera_static

!> Interfaces ------------------------------------------------
contains
!> Procedures ------------------------------------------------
!> generate points on lens and retrive their pdf
!> inputs:
!>   camera_inout: (camera_static) camera with initialised points on lens
!>   lens_inout:   (lens) lens model for generating points
!>   n_points_in:  (integer)(optional) number of points to sample
!>                                     default: 1000, pinhole: 1
!> outputs:
!>   camera_out:   (camera_static) camera with initialised points on lens
subroutine generate_points_on_lens_static(camera_inout,&
lens_inout,n_points_in)
  use mod_lens,         only: lens
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> inputs-outputs:
  class(camera_static),intent(inout) :: camera_inout
  class(lens),intent(inout)          :: lens_inout
  !> inputs:
  integer,intent(in),optional        :: n_points_in
  !> initialisations
  camera_inout%n_vertices = 1000
  if(present(n_points_in)) camera_inout%n_vertices = n_points_in
  !> check if the lens is a pinhole, overwrite n_points_lens
  select type(ln=>lens_inout)
    type is(pinhole_lens)
    camera_inout%n_vertices = 1
  end select
  camera_inout%n_active_vertices = camera_inout%n_vertices
  !> sample the lens
  if(allocated(camera_inout%x)) then
    if(camera_inout%n_vertices.ne.size(camera_inout%x,dim=2)) &
    deallocate(camera_inout%x)
  endif
  if(allocated(camera_inout%properties)) then
    if(camera_inout%n_vertices.ne.size(camera_inout%properties,dim=2)) &
    deallocate(camera_inout%properties)
  endif
  if(.not.allocated(camera_inout%x)) &
  allocate(camera_inout%x(camera_inout%n_x,camera_inout%n_vertices,1))
  if(.not.allocated(camera_inout%properties)) &
  allocate(camera_inout%properties(camera_inout%n_property_vertex,camera_inout%n_vertices,1))
  call lens_inout%sampling(camera_inout%n_vertices,camera_inout%x(:,:,1))
  call lens_inout%pdf(camera_inout%n_vertices,camera_inout%x(:,:,1),&
  camera_inout%properties(1,:,1))
end subroutine generate_points_on_lens_static

!> generate the image plane direction, vertices and compute the pixel width and height
!> inputs:
!>  camera_inout: (camera_static) camera with unallocated image plane
!>  n_int_param:  (integer) number of integer parameters
!>  n_real_param: (integer) number of real parameters
!>  int_param:    (integer)(n_int_param) integer parameters
!>                4) mirror w.r.t. the x axis
!>                5) mirror w.r.t. the y axis
!>  real_param:   (real8)(n_real_param) real parameters, order:
!>                1:3) plane width, height half angles and plane orientation
!>                4:6) plane position w.r.t. the pupil in spherical coordinates
!>                7:9) position of the pupil in cartesian coordinates
!> outputs:
!>  camera_inout: (camera_static) camera with defined image plane
subroutine define_image_plane_pixel_size(camera_inout,n_int_param,n_real_param,&
int_param,real_param)
  use mod_geometry, only: compute_plane_edge_length
  use mod_geometry, only: define_vertex_spherical_coord
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> inputs-outputs:
  class(camera_static),intent(inout)        :: camera_inout
  !> inputs:
  integer,intent(in)                        :: n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in) :: int_param
  real*8,dimension(n_real_param),intent(in) :: real_param
  !> define the image plane direction and store it
  call define_vertex_spherical_coord((/1.d0,real_param(5),real_param(6)/),&
  camera_inout%image_plane_direction(:,1))
  !> define the plane from the width/height half angles, the distance
  !> from the pupil and the pupil position
  call define_plane_from_half_angles(int_param(4:5),real_param(1:3),&
  real_param(4:6),real_param(7:9),camera_inout%image_plane(:,:,1))
  !> compute the image plane lenght of the edges
  call compute_plane_edge_length(camera_inout%image_plane(:,:,1),camera_inout%plane_edge_length(:,1))
  !> compute the pixel width and heigh in the pixel reference system
  camera_inout%pixel_size = (/1.d0,1.d0/)/real(camera_inout%n_pixels_spectra(2:3),kind=8)
end subroutine define_image_plane_pixel_size

!> find intersection between a camera ray and the image plane
!> inputs:
!>   camera_inout: (camera_perspective_static) initialised camera object
!>   x_pos:        (real8)(3) coordinated of the point defining a ray
!>   x_lens_id:    (integer) index of the point on lens to be treated
!>   time_id_in:   (integer)(optional) time index (not used)
!> outputs:
!>   camera_inout: (camera_perspective_static) initialised camera object
!>   intersect:    (bool) if true an intersection is found
!>   local_coords: (real8)(3) plane (s,t) and ray (q) local 
!>                            coordinates of the intersection
subroutine find_ray_image_plane_intersection_static(camera_inout,&
x_pos,x_lens_id,intersect,local_coords,time_id_in)
  use mod_geometry, only: compute_plane_line_intersection_cart
  implicit none
  !> inputs-outpus:
  class(camera_static),intent(inout)             :: camera_inout
  !> inputs:
  integer,intent(in)                             :: x_lens_id
  real*8,dimension(camera_inout%n_x),intent(in)  :: x_pos
  integer,intent(in),optional                    :: time_id_in
  !> outputs:
  logical,intent(out)                            :: intersect
  real*8,dimension(3),intent(out)                :: local_coords
  !> variables
  real*8,dimension(camera_inout%n_x,3)           :: plane
  real*8,dimension(camera_inout%n_x,2)           :: ray
  !> compute intersection between the ray and the image plane
  ray(:,1) = camera_inout%x(:,x_lens_id,1); ray(:,2) = x_pos;
  call compute_plane_line_intersection_cart(camera_inout%image_plane(:,:,1),&
  ray,intersect,local_coords)
end subroutine find_ray_image_plane_intersection_static

!>------------------------------------------------------------
end module mod_camera_static
