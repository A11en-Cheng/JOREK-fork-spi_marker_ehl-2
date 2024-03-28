!> mod_camera_static_test contains all variables and 
!> procedures used for unit testing the routines
!> shared by all static camera
module mod_camera_static_test
use fruit
use mod_assert_equals_tools,       only: assert_equals_allocatable_arrays
use constants,                     only: PI,TWOPI
use mod_pinhole_lens,              only: pinhole_lens
use mod_spectra_deterministic,     only: spectrum_integrator_2nd
use mod_camera_perspective_static, only: camera_perspective_static
implicit none

private
public :: run_fruit_camera_static

!> Variable and data types -------------------------
integer,parameter :: n_int_param=5
integer,parameter :: n_real_param=9
integer,parameter :: n_planes=11
integer,parameter :: n_properties=1
integer,parameter :: n_st_sol=2
integer,parameter :: n_stq_sol=3
integer,parameter :: n_x_sol=3
integer,parameter :: n_times_sol=1
integer,parameter :: n_plane_vertices=3
integer,parameter :: n_pixels_x=512
integer,parameter :: n_pixels_y=512
integer,parameter :: n_points_on_lens_sol=2345
integer,parameter :: n_lines_per_spectrum=13
integer,parameter :: n_spectra=2
integer,parameter :: n_rays_sol=123
integer,parameter :: test_plane_pupil_id=1
integer,dimension(2),parameter :: mirror_xy_interval=(/0,1/)
real*8,parameter  :: tol_real8=5.d-15
real*8,parameter  :: tol_real8_2=5.d-9
real*8,parameter  :: plane_distance=2.5d1
real*8,dimension(2),parameter :: costheta_interval=(/-1.d0,1.d0/)
real*8,dimension(2),parameter :: phi_interval=(/0.d0,TWOPI/)
real*8,dimension(2),parameter :: ray_q_interval=(/3.1d-2,9.5d-1/)
real*8,dimension(3),parameter :: half_angle_lowbnd=(/PI/1.d1,PI/4.d0,0d0/)
real*8,dimension(3),parameter :: half_angle_uppbnd=(/PI/1.9d0,PI/3.d0,TWOPI/)
real*8,dimension(3),parameter :: center_pos_lowbnd=(/-2.d-1,4.d1,-7.d0/)
real*8,dimension(3),parameter :: center_pos_uppbnd=(/3.4d1,3.d2,5.d0/)
real*8,dimension(2),parameter :: st_false_lowbnd=(/-5.4d1,-9.d0/)
real*8,dimension(2),parameter :: st_false_uppbnd=(/1.4d1,4.5d2/)
real*8,dimension(3),parameter :: test_points_lowbnd=(/-7.4d-1,2.9d1,-5.d0/)
real*8,dimension(3),parameter :: test_points_uppbnd=(/1.4d1,4.5d2,9.d0/)
real*8,dimension(n_spectra),parameter :: min_wlen=(/3.d0-6,2.5d-7/)
real*8,dimension(n_spectra),parameter :: max_wlen=(/3.5d0-6,4.2d-7/)
logical,dimension(n_rays_sol)          :: accept_ray_sol
integer,dimension(2,n_planes)          :: mirror_xy_sol
real*8,dimension(n_st_sol)             :: pixel_size_sol
real*8,dimension(3,n_planes)           :: half_angle_sol
real*8,dimension(n_x_sol,n_planes)     :: image_plane_coords
real*8,dimension(n_x_sol,n_planes)     :: pupil_positions
real*8,dimension(n_x_sol,n_rays_sol)   :: test_ray_vertices
real*8,dimension(n_stq_sol,n_rays_sol) :: test_ray_stq_sol
real*8,dimension(n_x_sol,n_planes)     :: test_points
real*8,dimension(:,:,:),allocatable    :: points_on_lens
real*8,dimension(:,:,:),allocatable    :: pdf_points_on_lens
type(camera_perspective_static) :: camera_sol
type(pinhole_lens)              :: pinhole_sol
type(spectrum_integrator_2nd)   :: spectrum_sol
!> Interfaces --------------------------------------
contains
!> Fruit basket ------------------------------------
!> fruit basket executes all set-up, test and 
!> tear-down procedures
subroutine run_fruit_camera_static()
  implicit none
  write(*,'(/A)') "  ... setting-up: camera static tests"
  call setup
  write(*,'(/A)') "  ... running: camera static tests"
  call run_test_case(test_points_on_lens_pdf_pinhole,&
  'test_points_on_lens_pdf_pinhole')
  call run_test_case(test_image_plane_pixel_size_definitions,&
  'test_image_plane_pixel_size_definitions')
  call run_test_case(test_find_ray_image_plane_intersection,&
  'test_find_ray_image_plane_intersection')
  write(*,'(/A)') "  ... tearing-up: camera static tests"
  call teardown
end subroutine run_fruit_camera_static

!> Set-up and tear-down ----------------------------
!> set-up unit test features
!> set-up unit test features
subroutine setup()
  use mod_gnu_rng,                 only: gnu_rng_interval
  use mod_common_camera_test_tool, only: generate_image_plane_variables
  use mod_common_camera_test_tool, only: generate_ray_variables_from_origin_plane
  implicit none
  integer :: ii
  real*8,dimension(3) :: center
  !> compute the solution pixel size
  pixel_size_sol = (/1.d0,1.d0/)/real((/n_pixels_x,n_pixels_y/),kind=8)
  !> initialise the pinhole lens
  call gnu_rng_interval(3,center_pos_lowbnd,center_pos_uppbnd,center)
  call pinhole_sol%init_pinhole(n_x_sol,center)
  !> initialise camera spectrum
  spectrum_sol = spectrum_integrator_2nd(n_lines_per_spectrum,&
  n_spectra,min_wlen,max_wlen)
  !> initialise image plane variables
  call generate_image_plane_variables(n_x_sol,n_planes,mirror_xy_interval,&
  plane_distance,costheta_interval,phi_interval,half_angle_lowbnd,&
  half_angle_uppbnd,center_pos_lowbnd,center_pos_uppbnd,mirror_xy_sol,&
  half_angle_sol,image_plane_coords,pupil_positions)
  !> initialise camera
  camera_sol%n_plane_points = n_plane_vertices
  !> initialise the test points
  do ii=1,n_planes
    call gnu_rng_interval(n_x_sol,test_points_lowbnd,&
    test_points_uppbnd,test_points(:,ii))
  enddo
  !> generate the ray variables for one image plane
  call generate_ray_variables_from_origin_plane(n_x_sol,n_st_sol,n_stq_sol,n_rays_sol,&
  st_false_lowbnd,st_false_uppbnd,ray_q_interval,mirror_xy_sol(:,test_plane_pupil_id),&
  half_angle_sol(:,test_plane_pupil_id),pupil_positions(:,test_plane_pupil_id),&
  image_plane_coords(:,test_plane_pupil_id),center,accept_ray_sol,&
  test_ray_stq_sol,test_ray_vertices)
end subroutine setup

!> tearing-down unit test features
subroutine teardown()
  implicit none
  pixel_size_sol = 0.d0
  call pinhole_sol%deallocate_lens; call spectrum_sol%deallocate_spectrum;
  if(allocated(points_on_lens))     deallocate(points_on_lens)
  if(allocated(pdf_points_on_lens)) deallocate(pdf_points_on_lens)
end subroutine teardown

!> Tests -------------------------------------------
!> test generation of points and pdf on lens using
!> a pinhole lens
subroutine test_points_on_lens_pdf_pinhole()
  implicit none
  !> allocate and initialise variables
  camera_sol%n_property_vertex = n_properties
  call camera_sol%allocate_vertices(1,1)
  allocate(points_on_lens(n_x_sol,1,1))
  allocate(pdf_points_on_lens(1,1,1))
  call pinhole_sol%sampling(1,points_on_lens)
  call pinhole_sol%pdf(1,points_on_lens(:,:,1),pdf_points_on_lens(:,1,1))
  !> test generation without input number of points
  call camera_sol%generate_points_on_lens_pdf(pinhole_sol)
  call assert_equals(1,camera_sol%n_vertices,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: x")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens")
  !> test generation with input number points
  call camera_sol%generate_points_on_lens_pdf(pinhole_sol,n_points_on_lens_sol)
  call assert_equals(1,camera_sol%n_vertices,&
  "Error camera perspective static generate points on pinhole: n vertices not 1!")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: points on lens")
  call assert_equals_allocatable_arrays(n_properties,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error camera perspective static generate points on pinhole: pdf points on lens")
  !> deallocate variable
  call camera_sol%deallocate_vertices
  deallocate(points_on_lens); deallocate(pdf_points_on_lens);
end subroutine test_points_on_lens_pdf_pinhole

!> test the generation of image planes and the calculation of the pixel size
subroutine test_image_plane_pixel_size_definitions()
  use constants,    only: TWOPI
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> variables
  integer                                    :: ii
  integer,dimension(n_int_param)             :: int_param
  real*8,dimension(n_x_sol)                  :: direction_sol,direction_test
  real*8,dimension(n_x_sol,n_plane_vertices) :: image_plane_sol
  real*8,dimension(n_x_sol,n_plane_vertices) :: image_plane_std_sol
  real*8,dimension(n_real_param)             :: real_param
  !> initialisation
  camera_sol%n_property_vertex = n_properties;
  call camera_sol%allocate_camera_perspective_static(spectrum_sol,3,&
  (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/))
  !> test the definition of the image plane and pixel size
  do ii=1,n_planes
    !> store plane value in parameters
    int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y,&
    mirror_xy_sol(1,ii),mirror_xy_sol(2,ii)/)
    real_param(1:3) = half_angle_sol(:,ii)
    real_param(4:6) = image_plane_coords(:,ii)
    real_param(7:9) = pupil_positions(:,ii)
    !> generate camera and solution planes
    call camera_sol%define_image_plane_pixel_size(n_int_param,n_real_param,int_param,real_param)
    call define_plane_from_half_angles(mirror_xy_sol(:,ii),half_angle_sol(:,ii),&
    image_plane_coords(:,ii),pupil_positions(:,ii),image_plane_sol)
    call define_plane_from_half_angles(mirror_xy_sol(:,ii),half_angle_sol(:,ii),&
    image_plane_coords(:,ii),image_plane_std_sol)
    direction_test = (/1.d0,acos(camera_sol%image_plane_direction(3,1)),&
    atan2(camera_sol%image_plane_direction(2,1),camera_sol%image_plane_direction(1,1))/)
    if(direction_test(3).lt.0.d0) direction_test(3) = TWOPI + direction_test(3)
    direction_sol = (/1.d0,real_param(5),real_param(6)/)
    !> test plane and pixel size
    call assert_equals(direction_sol,direction_test,n_x_sol,tol_real8,&
    "Error camera perspective static define image plane: image plane direction mismatch!")
    call assert_equals(image_plane_sol,camera_sol%image_plane(:,:,1),n_x_sol,&
    n_plane_vertices,tol_real8,&
    "Error camera perspective static define image plane: image plane mismatch!")
    call assert_equals(pixel_size_sol,camera_sol%pixel_size,2,tol_real8,&
    "Error camera perspective static define image plane: pixel size mismatch!")
  enddo
  !> deallocate camera
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_image_plane_pixel_size_definitions

!> test the method for finding image plane - ray intersections
!> for simplicity, only one image plane is tested
subroutine test_find_ray_image_plane_intersection()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_int_param)         :: int_param
  logical,dimension(n_rays_sol)          :: test_intersection
  real*8,dimension(n_real_param)         :: real_param
  real*8,dimension(n_stq_sol,n_rays_sol) :: test_ray_stq
  !> initialisation
  camera_sol%n_property_vertex = n_properties;
  int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y,&
  mirror_xy_sol(1,test_plane_pupil_id),mirror_xy_sol(2,test_plane_pupil_id)/)
  real_param(1:3) = half_angle_sol(:,test_plane_pupil_id)
  real_param(4:6) = image_plane_coords(:,test_plane_pupil_id)
  real_param(7:9) = pupil_positions(:,test_plane_pupil_id)
  call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
  n_int_param,n_real_param,int_param,real_param)
  !> compute rays and intersections
  do ii=1,n_rays_sol
    call camera_sol%find_ray_image_plane_intersection(&
    test_ray_vertices(:,ii),test_plane_pupil_id,&
    test_intersection(ii),test_ray_stq(:,ii))
  enddo
  !> check solutions
  call assert_equals(accept_ray_sol,test_intersection,n_rays_sol,&
  "Error find ray image plane intersection: intersections mismatch!")
  call assert_equals(test_ray_stq_sol,test_ray_stq,n_stq_sol,n_rays_sol,&
  tol_real8_2,"Error find ray image plane intersection: local coordinates mismatch!")
  !> cleanup
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_find_ray_image_plane_intersection

!> Tools -------------------------------------------
!>--------------------------------------------------
end module mod_camera_static_test
