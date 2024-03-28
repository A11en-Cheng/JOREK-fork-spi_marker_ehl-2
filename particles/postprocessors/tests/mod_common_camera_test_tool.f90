!> mod_common_camera_test_tool contains variables and
!> procedures used by multiple camera unit test modules
module mod_common_camera_test_tool
implicit none
private
public :: generate_image_plane_variables
public :: generate_ray_variables_from_origin_plane

!> Variable and data types -------------------------
!> Interfaces --------------------------------------
contains
!> Procedures --------------------------------------
!> sample the image plane variables
!> inputs:
!>   n_x:                (integer) size of the position vector
!>   n_planes:           (integer) number of image plane to test
!>   mirror_xy_interval: (integer)(2) interval (0,1) for selecting if a plane
!>                       should be mirrored or not
!>   plane_distance:     (real8) distance of the image plane from the pupil
!>   costheta_interval:  (real8)(2) selection interval of the plane center colatitude
!>   phi_interval:       (real8)(2) selection interval of the plane center azimuth
!>   half_angle_lowbnd:  (real8)(3) image plane dimension lowerbounds
!>   half_angle_uppbnd:  (real8)(3) image plane dimension upperbounds
!>   center_pos_lowbnd:  (real8)(3) pupil position lowerbounds
!>   center_pos_uppbnd:  (real8)(3) pupil position upperbounds
!> outputs:
!>   mirror_xy:         (integer)(2,n_planes) if a plane should be mirrored or not
!>   half_angle:        (real8)(3,n_panes) image plane angular coordinates
!>   image_plane_coord: (real8)(n_x,n_planes) positions of the image plane
!>                       centers in spherical coordinates w.r.t. the pupil
!>   pupil_positions:   (real8)(n_x,n_planes) position in xyz coordinates
!>                      of the camera pupil for each plane
subroutine generate_image_plane_variables(n_x,n_planes,mirror_xy_interval,&
plane_distance,costheta_interval,phi_interval,half_angle_lowbnd,&
half_angle_uppbnd,center_pos_lowbnd,center_pos_uppbnd,mirror_xy,&
half_angle,image_plane_coords,pupil_positions)
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  !> inputs:
  integer,intent(in)               :: n_x,n_planes
  integer,dimension(2),intent(in)  :: mirror_xy_interval
  real*8,intent(in)                :: plane_distance
  real*8,dimension(2),intent(in)   :: costheta_interval,phi_interval
  real*8,dimension(3),intent(in)   :: half_angle_lowbnd,half_angle_uppbnd 
  real*8,dimension(n_x),intent(in) :: center_pos_lowbnd,center_pos_uppbnd
  !> outputs:
  integer,dimension(2,n_planes),intent(out)  :: mirror_xy
  real*8,dimension(3,n_planes),intent(out)   :: half_angle
  real*8,dimension(n_x,n_planes),intent(out) :: image_plane_coords
  real*8,dimension(n_x,n_planes),intent(out) :: pupil_positions
  !> variables:
  integer             :: ii
  real*8,dimension(3) :: rand
  !> generation routine
  call gnu_rng_interval(2,n_planes,mirror_xy_interval,mirror_xy)
  do ii=1,n_planes
    call random_number(rand)
    image_plane_coords(:,ii) = sample_uniform_sphere(&
    plane_distance,costheta_interval,phi_interval,rand)
    call gnu_rng_interval(3,half_angle_lowbnd,&
    half_angle_uppbnd,half_angle(:,ii))
    call gnu_rng_interval(n_x,center_pos_lowbnd,&
    center_pos_uppbnd,pupil_positions(:,ii))
  enddo
end subroutine generate_image_plane_variables

!> sample the vertex of a ray given an origin and a plane
!> inputs:
!>   n_x:             (integer) size of the position vector
!>   n_st:            (integer) size of the plane local coordinates
!>   n_stq:           (integer) number of local plane/ray coordinates
!>   n_rays:          (integer) number of ray to generate
!>   mirror_xy:       (integer)(2) check if the plane should be mirrored 
!>   st_false_lowbnd: (real8)(2) lower bound of the local plane coordinates
!>   st_false_uppbnd: (real8)(2) upper bound of the local plane coordinates
!>   ray_q_interval:  (real8)(2) interval of the local ray coordinate
!>   half_width:      (real8)(3) plane dimensions
!>   origin:          (real8)(n_x) position of the pupil
!>   plane_coords:    (real8)(n_x) position of the plane in the pupil reference
!>   x_lens:          (real8)(n_x) position on the camera lens
!>                                 intersected by the ray
!> outputs:
!>   accept_ray:      (logical)(n_rays) if true plane and ray intersect
!>   ray_stq:         (real8)(n_stq,n_rays) local plane-ray coordinates
!>   ray_vertices:    (n_x,n_rays) vertices of the rays
subroutine generate_ray_variables_from_origin_plane(n_x,n_st,n_stq,n_rays,&
st_false_lowbnd,st_false_uppbnd,ray_q_interval,mirror_xy,half_width,&
origin,plane_coords,x_lens,accept_ray,ray_stq,ray_vertices)
  use mod_geometry, only: define_plane_from_half_angles
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_gnu_rng,  only: gnu_rng_interval
  implicit none
  !> parameters:
  real*8,parameter  :: accept_threshold=5.d-1
  !> inputs:
  integer,intent(in)                 :: n_x,n_st,n_stq,n_rays
  integer,dimension(2),intent(in)    :: mirror_xy
  real*8,dimension(3),intent(in)     :: half_width
  real*8,dimension(2),intent(in)     :: ray_q_interval
  real*8,dimension(2),intent(in)     :: st_false_lowbnd,st_false_uppbnd
  real*8,dimension(n_x),intent(in)   :: origin,x_lens
  real*8,dimension(n_x),intent(in)   :: plane_coords
  !> outputs:
  logical,dimension(n_rays),intent(out)      :: accept_ray
  real*8,dimension(n_stq,n_rays),intent(out) :: ray_stq
  real*8,dimension(n_x,n_rays),intent(out)   :: ray_vertices
  !> variables:
  integer :: ii
  real*8 ::  rand
  real*8,dimension(n_st)  :: st_value
  real*8,dimension(n_x)   :: plane_pos
  real*8,dimension(n_x,3) :: plane
  !> intiialisation
  call define_plane_from_half_angles(mirror_xy,half_width,plane_coords,origin,plane)
  ray_stq(1:2,:) = 5.d-1
  !> loop on the number of rays
  do ii=1,n_rays
    !> check if a ray with or without intersection must be generated
    call random_number(rand)
    if(rand.gt.accept_threshold) then
      call random_number(ray_stq(1:2,ii))
      accept_ray(ii) = .true.
    else
      do while((all(ray_stq(1:2,ii).ge.0.d0).and.all(ray_stq(1:2,ii).le.1.d0)))
        call gnu_rng_interval(2,st_false_lowbnd,st_false_uppbnd,ray_stq(1:2,ii))
      enddo
      accept_ray(ii) = .false.
    endif
    !> sample a ray lenght
    call gnu_rng_interval(ray_q_interval,ray_stq(3,ii))
    !> compute the position on the plane
    plane_pos = compute_global_cart_coord_plane_points(plane,ray_stq(1:2,ii))
    !> compute and store the ray vertex and its local coordinate
    ray_vertices(:,ii) = x_lens+(plane_pos-x_lens)/ray_stq(3,ii)
  enddo
end subroutine generate_ray_variables_from_origin_plane
!> -------------------------------------------------
end module mod_common_camera_test_tool
