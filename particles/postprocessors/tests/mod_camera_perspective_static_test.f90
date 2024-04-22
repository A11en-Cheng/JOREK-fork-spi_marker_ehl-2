!> mod_camera_perspective_static_test contains all 
!> variables and procedures used for unit testing 
!> the camera perspective static model
module mod_camera_perspective_static_test
use fruit
use mod_assert_equals_tools,       only: assert_equals_allocatable_arrays
use constants,                     only: PI,TWOPI
use mod_pinhole_lens,              only: pinhole_lens
use mod_spectra_deterministic,     only: spectrum_integrator_2nd
use mod_camera_perspective_static, only: camera_perspective_static
implicit none

private
public :: run_fruit_camera_perspective_static

!> Variable and data types -------------------------
integer,parameter :: n_int_param=5
integer,parameter :: n_real_param=9
integer,parameter :: n_planes=11
integer,parameter :: n_properties=1
integer,parameter :: n_st_sol=2
integer,parameter :: n_x_sol=3
integer,parameter :: n_times_sol=1
integer,parameter :: n_plane_vertices=3
integer,parameter :: n_pixels_x=512
integer,parameter :: n_pixels_y=512
integer,parameter :: n_points_on_lens_sol=2345
integer,parameter :: n_lines_per_spectrum=13
integer,parameter :: n_spectra=2
integer,parameter :: test_plane_pupil_id=1
integer,dimension(2),parameter :: mirror_xy_interval=(/0,1/)
real*8,parameter  :: tol_real8=5.d-15 
real*8,parameter  :: plane_distance=2.5d1
real*8,dimension(2),parameter :: costheta_interval=(/-1.d0,1.d0/)
real*8,dimension(2),parameter :: phi_interval=(/0.d0,TWOPI/)
real*8,dimension(3),parameter :: half_angle_lowbnd=(/PI/1.d1,PI/4.d0,0d0/)
real*8,dimension(3),parameter :: half_angle_uppbnd=(/PI/1.9d0,PI/3.d0,TWOPI/)
real*8,dimension(3),parameter :: center_pos_lowbnd=(/-2.d-1,4.d1,-7.d0/)
real*8,dimension(3),parameter :: center_pos_uppbnd=(/3.4d1,3.d2,5.d0/)
real*8,dimension(3),parameter :: test_points_lowbnd=(/-7.4d-1,2.9d1,-5.d0/)
real*8,dimension(3),parameter :: test_points_uppbnd=(/1.4d1,4.5d2,9.d0/)
real*8,dimension(n_spectra),parameter :: min_wlen=(/3.d0-6,2.5d-7/)
real*8,dimension(n_spectra),parameter :: max_wlen=(/3.5d0-6,4.2d-7/)
integer,dimension(2,n_planes)         :: mirror_xy_sol
real*8,dimension(n_st_sol)            :: pixel_size_sol
real*8,dimension(3,n_planes)          :: half_angle_sol
real*8,dimension(n_x_sol,n_planes)    :: image_plane_coords
real*8,dimension(n_x_sol,n_planes)    :: pupil_positions
real*8,dimension(n_x_sol,n_planes)    :: test_points
real*8,dimension(:,:,:),allocatable   :: points_on_lens
real*8,dimension(:,:,:),allocatable   :: pdf_points_on_lens
type(camera_perspective_static) :: camera_sol
type(pinhole_lens)              :: pinhole_sol
type(spectrum_integrator_2nd)   :: spectrum_sol

!> Interfaces --------------------------------------
interface compute_cos_angle_two_vectors
  module procedure compute_cos_angle_two_vectors_origin_points
end interface compute_cos_angle_two_vectors

contains
!> Fruit basket ------------------------------------
!> fruit basket executes all set-up, test and 
!> tear-down procedures
subroutine run_fruit_camera_perspective_static
  implicit none
  write(*,'(/A)') "  ... setting-up: camera perspective static tests"
  call setup
  write(*,'(/A)') "  ... running: camera perspective static tests"
  call run_test_case(test_de_allocation_camera_perspective_static,&
  'test_de_allocation_camera_perspective_static')
  call run_test_case(test_init_camera_perspective_static_pinhole,&
  'test_init_camera_perspective_static_pinhole')
  call run_test_case(test_cosine_view_angle_static,&
  'test_cosine_view_angle_static')
  call run_test_case(test_material_funct_perspective_static,&
  'test_material_funct_perspective_static') 
  write(*,'(/A)') "  ... tearing-up: camera perspective static tests"
  call teardown
end subroutine run_fruit_camera_perspective_static

!> Set-up and tear-down ----------------------------
!> set-up unit test features
subroutine setup()
  use mod_gnu_rng,                 only: gnu_rng_interval
  use mod_common_camera_test_tool, only: generate_image_plane_variables
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
!> test the initialisation camera perspective static
subroutine test_init_camera_perspective_static_pinhole()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_geometry, only: define_plane_from_half_angles
  implicit none
  !> variables
  integer,dimension(n_int_param)           :: int_param
  real*8,dimension(n_x_sol)                :: direction_sol,direction_test
  real*8,dimension(n_real_param)           :: real_param
  real*8,dimension(n_x_sol,n_plane_vertices) :: image_plane_sol
  !> initialisation, only one image plane is considered
  int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y,&
  mirror_xy_sol(1,test_plane_pupil_id),mirror_xy_sol(2,test_plane_pupil_id)/)
  real_param(1:3) = half_angle_sol(:,test_plane_pupil_id)
  real_param(4:6) = image_plane_coords(:,test_plane_pupil_id)
  real_param(7:9) = pupil_positions(:,test_plane_pupil_id)
  allocate(points_on_lens(n_x_sol,1,1))
  allocate(pdf_points_on_lens(1,1,1))
  call pinhole_sol%sampling(1,points_on_lens)
  call pinhole_sol%pdf(1,points_on_lens(:,:,1),pdf_points_on_lens(:,1,1))
  !> initialise the camera perspective static
  call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
  n_int_param,n_real_param,int_param,real_param)
  !> perform tests
  call assert_equals(1,camera_sol%n_times,&
  "Error init camera perspective static pinhole: n_times is not 1")
  call assert_equals(1,camera_sol%n_vertices,&
  "Error init camera perspective static pinhole: n vertices is not 1")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,&
  "Error init camera perspective static pinhole: n active vertices")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%times,&
  "Error init camera perspective static pinhole: times")
  call assert_equals_allocatable_arrays(n_x_sol,1,&
  n_times_sol,camera_sol%x,"Error init camera perspective static pinhole: x")
  call assert_equals_allocatable_arrays(1,n_properties,n_times_sol,&
  camera_sol%properties,"Error init camera perspective static pinhole: properties")
  call assert_equals((/n_spectra,n_pixels_x,n_pixels_y/),camera_sol%n_pixels_spectra,&
  3,"Error init camera perspective static pinhole: n pixels / spectra")
  call assert_equals_allocatable_arrays(n_x_sol,camera_sol%n_vertices,&
  camera_sol%n_times,camera_sol%x,points_on_lens,tol_real8,&
  ":Error init camera perspective static pinhole: points on lens")
  call assert_equals_allocatable_arrays(camera_sol%n_vertices,n_properties,&
  camera_sol%n_times,camera_sol%properties,pdf_points_on_lens,tol_real8,&
  "Error init camera perspective static pinhole: pdf points on lens") 
  call assert_equals_allocatable_arrays(n_x_sol,n_times_sol,camera_sol%image_plane_direction,&
  "Error init camera perspective static pinhole: image plane direction")
  if(allocated(camera_sol%image_plane_direction)) then
    direction_test = (/1.d0,acos(camera_sol%image_plane_direction(3,1)),&
    atan2(camera_sol%image_plane_direction(2,1),camera_sol%image_plane_direction(1,1))/)
    if(direction_test(3).lt.0.d0) direction_test(3) = TWOPI + direction_test(3)
    direction_sol = (/1.d0,image_plane_coords(2,1),image_plane_coords(3,1)/)
    call assert_equals(direction_sol,direction_test,n_x_sol,tol_real8,&
    "Error init camera perspective static pinhole: image plane mismatch!")
  endif
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_vertices,n_times_sol,&
  camera_sol%image_plane,"Error init camera perspective static pinhole: image plane")
  if(allocated(camera_sol%image_plane)) then
    call define_plane_from_half_angles(mirror_xy_sol(:,1),half_angle_sol(:,1),&
    image_plane_coords(:,1),pupil_positions(:,1),image_plane_sol)
    call assert_equals(image_plane_sol,camera_sol%image_plane(:,:,1),n_x_sol,&
    n_plane_vertices,tol_real8,&
    "Error init camera perspective static pinhole: image plane mismatch!")
  endif
  call assert_equals(camera_sol%pixel_size,pixel_size_sol,n_st_sol,tol_real8,&
  "Error init camera perspective static pinhole: pixel size mismatch!")
  !> clean up
  deallocate(points_on_lens); deallocate(pdf_points_on_lens);
  call camera_sol%deallocate_camera_perspective_static
end subroutine test_init_camera_perspective_static_pinhole

!> test allocation and deallocation of camera_perspective_static
!> attributs (only)
subroutine test_de_allocation_camera_perspective_static()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> initialise property size
  camera_sol%n_property_vertex = n_properties
  !> allocate camera perspective static
  call camera_sol%allocate_camera_perspective_static(spectrum_sol,3,&
  (/n_points_on_lens_sol,n_pixels_x,n_pixels_y/))
  !> test allocation
  call assert_equals(camera_sol%n_vertices,n_points_on_lens_sol,&
  "Error allocate camera perspective static: n vertices mismatch")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,&
  "Error allocate camera perspective static: n active vertices")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%times,&
  "Error allocate camera perspective static: times")
  call assert_equals_allocatable_arrays(n_x_sol,n_points_on_lens_sol,&
  n_times_sol,camera_sol%x,"Error allocate camera perspective static: x")
  call assert_equals_allocatable_arrays(n_properties,n_points_on_lens_sol,n_times_sol,&
  camera_sol%properties,"Error allocate camera perspective static: properties")
  call assert_equals((/n_spectra,n_pixels_x,n_pixels_y/),camera_sol%n_pixels_spectra,&
  3,"Error allocate camera perspective static: n pixels / spectra")
  call assert_equals_allocatable_arrays(n_x_sol,n_times_sol,camera_sol%image_plane_direction,&
  "Error allocate camera perspective static: image plane direction")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_vertices,n_times_sol,&
  camera_sol%image_plane,"Error allocate_camera perspective static: image plane") 
  !> deallocate camera perspective static
  call camera_sol%deallocate_camera_perspective_static
  call assert_equals(camera_sol%n_vertices,0,&
  "Error deallocate camera perspective static: n vertices not 0!")
  call assert_equals((/0,0,0/),camera_sol%n_pixels_spectra,3,&
  "Error deallocate camera perspective static: n pixels spectra not 0!")
  call assert_false(allocated(camera_sol%n_active_vertices),&
  "Error deallocate camera perspective static: n active vertices not deallocated!")
  call assert_false(allocated(camera_sol%times),&
  "Error deallocate camera perspective static: times not deallocated!")
  call assert_false(allocated(camera_sol%x),&
  "Error deallocate camera perspective static: x not deallocated!")
  call assert_false(allocated(camera_sol%properties),&
  "Error deallocate camera perspective static: properties not deallocated!")
  call assert_false(allocated(camera_sol%image_plane_direction),&
  "Error deallocate camera perspective static: image plane direction not deallocated!")
  call assert_false(allocated(camera_sol%image_plane),&
  "Error deallocate camera perspective static: image plane not deallocated!")
  !> nullify property size
  camera_sol%n_property_vertex = 0
end subroutine test_de_allocation_camera_perspective_static

!> test the calculation of the cosinus between the image plane direction and a ray
subroutine test_cosine_view_angle_static()
  use mod_geometry,     only: define_vertex_spherical_coord
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> variables
  type(pinhole_lens)             :: pinhole
  integer :: ii
  integer,dimension(n_int_param) :: int_param
  real*8,dimension(n_x_sol)      :: vertex_1
  real*8,dimension(n_real_param) :: real_param
  real*8,dimension(n_planes)     :: cos_view_angle_sol
  real*8,dimension(n_planes)     :: cos_view_angle_test 
  !> test the calculation of the view angle cosinus
  do ii=1,n_planes
    !> store plane value in parameters
    int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y,&
    mirror_xy_sol(1,ii),mirror_xy_sol(2,ii)/)
    real_param(1:3) = half_angle_sol(:,ii)
    real_param(4:6) = image_plane_coords(:,ii)
    real_param(7:9) = pupil_positions(:,ii)
    call pinhole%init_pinhole(n_x_sol,real_param(7:9))
    call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
    n_int_param,n_real_param,int_param,real_param)
    !> compute solution
    call define_vertex_spherical_coord(real_param(4:6),real_param(7:9),vertex_1)
    call compute_cos_angle_two_vectors(real_param(7:9),vertex_1,test_points(:,ii),cos_view_angle_sol(ii))
    !> compute test value
    call camera_sol%generate_points_on_lens_pdf(pinhole)
    call camera_sol%define_image_plane_pixel_size(n_int_param,n_real_param,int_param,real_param)
    call camera_sol%cos_view_angle_static(test_points(:,ii),1,cos_view_angle_test(ii))
    !> deallocate camera perspective static
    call camera_sol%deallocate_camera_perspective_static
  enddo
  !> check solutions
  call assert_equals(cos_view_angle_sol,cos_view_angle_test,n_planes,tol_real8,&
  "Error computation cosine view angle perspective static: cosine view angles mismatch!")
  !> deallocate lens
  call pinhole%deallocate_lens
end subroutine test_cosine_view_angle_static 

!> test the calculation of the physical material function for the perspective static camera
subroutine test_material_funct_perspective_static()
  use mod_geometry,     only: define_vertex_spherical_coord
  use mod_pinhole_lens, only: pinhole_lens
  implicit none
  !> variables
  type(pinhole_lens)             :: pinhole
  integer :: ii
  integer,dimension(n_int_param) :: int_param
  real*8,dimension(n_x_sol)      :: vertex_1
  real*8,dimension(n_real_param) :: real_param
  real*8,dimension(n_planes)     :: material_sol
  real*8,dimension(n_planes)     :: material_test
  !> test the calculation of the view angle cosinus
  do ii=1,n_planes
    !> store plane value in parameters
    int_param = (/n_points_on_lens_sol,n_pixels_x,n_pixels_y,&
    mirror_xy_sol(1,ii),mirror_xy_sol(2,ii)/)
    real_param(1:3) = half_angle_sol(:,ii)
    real_param(4:6) = image_plane_coords(:,ii)
    real_param(7:9) = pupil_positions(:,ii)
    !> compute solution
    call camera_sol%init_camera(pinhole_sol,spectrum_sol,&
    n_int_param,n_real_param,int_param,real_param)
    call define_vertex_spherical_coord(real_param(4:6),real_param(7:9),vertex_1)
    call compute_cos_angle_two_vectors(real_param(7:9),vertex_1,test_points(:,ii),material_sol(ii))
    !> compute test value
    call pinhole%init_pinhole(n_x_sol,real_param(7:9))
    call camera_sol%generate_points_on_lens_pdf(pinhole)
    call camera_sol%define_image_plane_pixel_size(n_int_param,n_real_param,int_param,real_param)
    call camera_sol%physical_material_funct(test_points(:,ii),1,material_test(ii),1)
    !> deallocate camera perspective static
    call camera_sol%deallocate_camera_perspective_static
  enddo
  !> check solutions
  call assert_equals(material_sol,material_test,n_planes,tol_real8,&
  "Error computation physical material funct perspective static: cosine view angles mismatch!")
  !> deallocate lens
  call pinhole%deallocate_lens
end subroutine test_material_funct_perspective_static

!> Tools -------------------------------------------
!> compute angle between two vectors with same origin
subroutine compute_cos_angle_two_vectors_origin_points(&
origin,vertex_1,vertex_2,cos_angle)
  implicit none
  !> inputs:
  real*8,dimension(n_x_sol) :: origin,vertex_1,vertex_2
  !> outputs:
  real*8 :: cos_angle
  !> comput angle
  cos_angle = dot_product(vertex_2-origin,vertex_1-origin)
  cos_angle = cos_angle/(norm2(vertex_2-origin)*norm2(vertex_1-origin))
end subroutine compute_cos_angle_two_vectors_origin_points

!>--------------------------------------------------
end module mod_camera_perspective_static_test

