!> mod_camera_test implements all variables and procedures
!> used for testing the camera model. Camera is an abstract
!> type => the type camera_perspective_static is used instead.
module mod_camera_test
use fruit
use mod_camera_perspective_static, only: camera_perspective_static
implicit none

private
public :: run_fruit_camera

!> Variables --------------------------------------------------------
integer,parameter :: n_x_sol=3
integer,parameter :: n_st_sol=2
integer,parameter :: n_properties_sol=11
integer,parameter :: n_times_sol=12
integer,parameter :: n_plane_points_sol=3
integer,parameter :: n_pixels_x_sol=256
integer,parameter :: n_pixels_y_sol=512
integer,parameter :: n_points_per_pixel_sol=11
integer,parameter :: n_spectra_sol=3
integer,parameter :: n_vertices_sol=1245
integer,parameter :: n_times_2_sol=33
integer,parameter :: n_plane_points_2_sol=3
integer,parameter :: n_pixels_x_2_sol=1024
integer,parameter :: n_pixels_y_2_sol=128
integer,parameter :: n_spectra_2_sol=12
integer,parameter :: n_vertices_2_sol=234
integer,parameter :: n_pixels_x_3_sol=31
integer,parameter :: n_pixels_y_3_sol=12
real*8,parameter  :: tol_real8=7.5e-15 
real*8,parameter  :: tol_real8_rel=3.d-6
real*8,dimension(2),parameter :: plane_edge_length_lowbnd=(/2.4d-3,4.3d-1/)
real*8,dimension(2),parameter :: plane_edge_length_uppbnd=(/3.25d2,4.74d1/) 
real*8,dimension(n_st_sol)             :: pixel_size_sol
real*8,dimension(2,n_times_sol)        :: plane_edge_length_sol
integer,dimension(:,:,:,:),allocatable :: pixel_ids
integer,dimension(:,:,:,:),allocatable :: pixel_ids_sol
real*8,dimension(:,:,:,:),allocatable  :: st_point_on_pixels
real*8,dimension(:,:,:,:),allocatable  :: st_point_on_pixels_sol
type(camera_perspective_static) :: camera_sol
!> Interfaces -------------------------------------------------------
contains
!> Fruit basket -----------------------------------------------------
!> fruit basket runs all the set-up, test and tear-down procedures
subroutine run_fruit_camera
  implicit none
  write(*,'(/A)') "  ... setting-up: camera tests"
  call setup()
  write(*,'(/A)') "  ... running: camera tests"
  call run_test_case(test_de_allocate_camera,'test_de_allocate_camera')
  call run_test_case(test_computation_pixel_ids_st_plane_point,&
  'test_computation_pixel_ids_st_plane_point')
  call run_test_case(test_computation_pixel_to_scaled_plane_point,&
  'test_computation_pixel_to_scaled_plane_point')
  call run_test_case(test_compute_scaled_pixel_position_on_plane,&
  'test_compute_scaled_pixel_position_on_plane')
  write(*,'(/A)') "  ... tearing-down: camera tests"
  call teardown()
end subroutine run_fruit_camera

!> Set-up and tear-down ---------------------------------------------
!> set-up test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> variables
  integer :: ii
  !> set n_properties
  camera_sol%n_property_vertex = n_properties_sol
  !> set the number of points of a plane
  camera_sol%n_plane_points = n_plane_points_sol
  pixel_size_sol = (/1.d0,1.d0/)/real((/n_pixels_x_sol,n_pixels_y_sol/),kind=8)
  !> generate random plane_edge_length for all times
  do ii=1,n_times_sol
    call gnu_rng_interval(2,plane_edge_length_lowbnd,plane_edge_length_uppbnd,&
    plane_edge_length_sol(:,ii))
  enddo
end subroutine setup

!> tear-down test features
subroutine teardown()
  implicit none
  camera_sol%n_property_vertex=0; camera_sol%n_plane_points=0;
  pixel_size_sol = 0d0;
  if(allocated(pixel_ids))          deallocate(pixel_ids)
  if(allocated(pixel_ids_sol))      deallocate(pixel_ids_sol)
  if(allocated(st_point_on_pixels)) deallocate(st_point_on_pixels)
  if(allocated(st_point_on_pixels)) deallocate(st_point_on_pixels_sol)
end subroutine teardown 

!> Tests ------------------------------------------------------------
!> test allocation and deallocation of the camera classes
subroutine test_de_allocate_camera()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> test allocation from allocated
  call camera_sol%allocate_camera(n_times_sol,n_vertices_sol,&
  n_pixels_x_sol,n_pixels_y_sol,n_spectra_sol)
  call assert_equals(n_times_sol,camera_sol%n_times,&
  "Error allocate camera from unallocated: n_times mismatch!")
  call assert_equals(n_vertices_sol,camera_sol%n_vertices,&
  "Error allocate camera from unallocated: n_vertices mismatch!")
  call assert_equals((/n_spectra_sol,n_pixels_x_sol,n_pixels_y_sol/),&
  camera_sol%n_pixels_spectra,3,&
  "Error allocate camera from unallocated: n_pixels_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_times_sol,camera_sol%n_active_vertices,0,&
  "Error allocate camera from unallocated: n_active_vertices")  
  call assert_equals_allocatable_arrays(n_times_sol,0.d0,camera_sol%times,&
  "Error allocate camera from unallocated: times") 
  call assert_equals_allocatable_arrays(n_x_sol,n_vertices_sol,n_times_sol,&
  0.d0,camera_sol%x,"Error allocate camera from unallocated: x")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_sol,n_times_sol,&
  0.d0,camera_sol%properties,"Error allocate camera from unallocated: properties")
  call assert_true(camera_sol%exposure_time.eq.0.d0,&
  "Error allocate camera from unallocated: exposure_time not zero!")
  call assert_equals_allocatable_arrays(2,n_times_sol,0.d0,camera_sol%plane_edge_length,&
  "Error allocate camera from unallocated: plane_edge_length")
  call assert_equals_allocatable_arrays(n_x_sol,n_times_sol,0.d0,camera_sol%image_plane_direction,&
  "Error allocate camera from unallocated: image_plane_direction")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_points_sol,n_times_sol,&
  0.d0,camera_sol%image_plane,"Error allocate camera from unallocated: image_plane_direction")
  !> test allocation from allocated
  camera_sol%n_plane_points = n_plane_points_2_sol
  call camera_sol%allocate_camera(n_times_2_sol,n_vertices_2_sol,&
  n_pixels_x_2_sol,n_pixels_y_2_sol,n_spectra_2_sol)
  call assert_equals(n_times_2_sol,camera_sol%n_times,&
  "Error allocate camera from allocated: n_times mismatch!")
  call assert_equals(n_vertices_2_sol,camera_sol%n_vertices,&
  "Error allocate camera from allocated: n_vertices mismatch!")
  call assert_equals((/n_spectra_2_sol,n_pixels_x_2_sol,n_pixels_y_2_sol/),&
  camera_sol%n_pixels_spectra,3,&
  "Error allocate camera from allocated: n_pixels_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_times_2_sol,camera_sol%n_active_vertices,0,&
  "Error allocate camera from allocated: n_active_vertices")  
  call assert_equals_allocatable_arrays(n_times_2_sol,0.d0,camera_sol%times,&
  "Error allocate camera from allocated: times") 
  call assert_equals_allocatable_arrays(n_x_sol,n_vertices_2_sol,n_times_2_sol,&
  0.d0,camera_sol%x,"Error allocate camera from allocated: x")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_2_sol,n_times_2_sol,&
  0.d0,camera_sol%properties,"Error allocate camera from allocated: properties")
  call assert_equals_allocatable_arrays(2,n_times_2_sol,0.d0,camera_sol%plane_edge_length,&
  "Error allocate camera from allocated: plane_edge_length")
  call assert_equals_allocatable_arrays(n_x_sol,n_times_2_sol,0.d0,camera_sol%image_plane_direction,&
  "Error allocate camera from allocated: image_plane_direction")
  call assert_equals_allocatable_arrays(n_x_sol,n_plane_points_2_sol,n_times_2_sol,&
  0.d0,camera_sol%image_plane,"Error allocate camera from allocated: image_plane_direction")
  !> test deallocation
  call camera_sol%deallocate_camera
  call assert_equals(0,camera_sol%n_times,&
  "Error deallocate camera: n_times not zero!")
  call assert_equals(0,camera_sol%n_vertices,&
  "Error deallocate camera: n_vertices not zero!")
  call assert_equals(0,camera_sol%n_plane_points,&
  "Error deallocate camera: n_plane_points not zero!")
  call assert_equals((/0,0,0/),camera_sol%n_pixels_spectra,3,&
  "Error deallocate camera: n_pixels_spectra not zero!")
  call assert_equals((/0.d0,0.d0/),camera_sol%pixel_size,2,&
  "Error deallocate camera: pixels_size not zero!")
  call assert_false(allocated(camera_sol%n_active_vertices),&
  "Error deallocate camera: n_active_vertices allocated!")
  call assert_false(allocated(camera_sol%times),&
  "Error deallocate camera: times size allocated!")
  call assert_false(allocated(camera_sol%x),&
  "Error deallocate camera: x allocated!")
  call assert_false(allocated(camera_sol%properties),&
  "Error deallocate camera: properties allocated!")
  call assert_true(camera_sol%exposure_time.eq.0.d0,&
  "Error deallocate camera: exposure time not zero!")
  call assert_false(allocated(camera_sol%plane_edge_length),&
  "Error deallocate camera: plane_edge_length allocated!")
  call assert_false(allocated(camera_sol%image_plane_direction),&
  "Error deallocate camera: image_plane_direction allocated!")
  call assert_false(allocated(camera_sol%image_plane),&
  "Error deallocate camera: image_plane allocated!")
end subroutine test_de_allocate_camera

!> test the calculation of the calculation of the position of a point on
!> a plane in the pixel coordinate system
subroutine test_computation_pixel_ids_st_plane_point()
  use mod_assert_equals_tools, only: assert_equals_extended
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer :: ii,jj,kk
  real*8,dimension(n_st_sol) :: st_plane
  !> initialisations
  camera_sol%pixel_size = pixel_size_sol
  allocate(st_point_on_pixels_sol(n_st_sol,n_points_per_pixel_sol,n_pixels_x_sol,n_pixels_y_sol))
  allocate(st_point_on_pixels(n_st_sol,n_points_per_pixel_sol,n_pixels_x_sol,n_pixels_y_sol))
  allocate(pixel_ids(n_st_sol,n_points_per_pixel_sol,n_pixels_x_sol,n_pixels_y_sol))
  allocate(pixel_ids_sol(n_st_sol,n_points_per_pixel_sol,n_pixels_x_sol,n_pixels_y_sol))
  !> compute local and global coordinates of the points
  do ii=1,n_pixels_y_sol
    do jj=1,n_pixels_x_sol
      do kk=1,n_points_per_pixel_sol
        pixel_ids_sol(:,kk,jj,ii) = (/jj,ii/)
        call random_number(st_point_on_pixels_sol(:,kk,jj,ii))
        st_plane = (st_point_on_pixels_sol(:,kk,jj,ii) + &
        real(pixel_ids_sol(:,kk,jj,ii)-1,kind=8))*pixel_size_sol
        call camera_sol%plane_to_pixel_local_coord(st_plane,&
        pixel_ids(:,kk,jj,ii),st_point_on_pixels(:,kk,jj,ii))
      enddo
    enddo
  enddo
  !> test results
  call assert_equals_extended(n_st_sol,n_points_per_pixel_sol,n_pixels_x_sol,&
  n_pixels_y_sol,pixel_ids_sol,pixel_ids,&
  "Error computation position in pixel coordinates: pixel ids mismatch!")
  call assert_equals_rel_error(n_st_sol,n_points_per_pixel_sol,n_pixels_x_sol,&
  n_pixels_y_sol,st_point_on_pixels_sol,st_point_on_pixels,tol_real8_rel,&
  "Error computation position in pixel coordinates: pixel coords. mismatch!")
  !> cleanup
  camera_sol%pixel_size = 0.d0
  deallocate(st_point_on_pixels); deallocate(pixel_ids_sol);
  deallocate(st_point_on_pixels_sol); deallocate(pixel_ids);
end subroutine test_computation_pixel_ids_st_plane_point

!> test the computation of the scaled plane coordinates from
!> the pixel coordinates
subroutine test_computation_pixel_to_scaled_plane_point()
  use mod_assert_equals_tools, only: assert_equals_extended
  implicit none
  !> variables
  integer :: ii,jj,kk
  real*8,dimension(2,n_pixels_x_3_sol,n_pixels_y_3_sol,n_times_sol) :: st_out,st_in
  !> initialise camera and test data
  call camera_sol%allocate_camera(n_times_sol,n_vertices_sol,&
  n_pixels_x_3_sol,n_pixels_y_3_sol,n_spectra_sol) 
  camera_sol%pixel_size = real(1d0/(/n_pixels_x_3_sol,n_pixels_y_3_sol/),kind=8)
  camera_sol%plane_edge_length = plane_edge_length_sol
  !> loop for computing the particle coordinates
  do kk=1,n_times_sol
    do jj=1,n_pixels_y_3_sol
      do ii=1,n_pixels_x_3_sol
        call random_number(st_in(:,ii,jj,kk))
        call camera_sol%pixel_local_to_scaled_plane_coord(kk,[ii,jj],st_in(:,ii,jj,kk),st_out(:,ii,jj,kk))
        st_out(:,ii,jj,kk) = (st_out(:,ii,jj,kk)/(camera_sol%pixel_size*camera_sol%plane_edge_length(:,kk))) - &
        real([ii-1,jj-1],kind=8)
      enddo
    enddo
  enddo
  !> check test results
  call assert_equals_extended(2,n_pixels_x_3_sol,n_pixels_y_3_sol,&
  n_times_sol,st_in,st_out,tol_real8,&
  "Error computation scaled position: local pixel coord. mismatch!")
  !> cleanup
  call camera_sol%deallocate_camera
end subroutine test_computation_pixel_to_scaled_plane_point

!> test the generation of the scaled pixel positions for all times
subroutine test_compute_scaled_pixel_position_on_plane()
  implicit none
  !> variables
  integer :: ii
  real*8,dimension(n_pixels_x_3_sol)             :: x_pixel_ids
  real*8,dimension(n_pixels_y_3_sol)             :: y_pixel_ids
  real*8,dimension(n_pixels_x_3_sol,n_times_sol) :: x_pixels,dot5_x
  real*8,dimension(n_pixels_y_3_sol,n_times_sol) :: y_pixels,dot5_y
  !> initialisations
  dot5_x = 5d-1; dot5_y = 5d-1;
  do ii=1,n_pixels_x_3_sol
    x_pixel_ids(ii) = real(ii-1,kind=8)
  enddo
  do ii=1,n_pixels_y_3_sol
    y_pixel_ids(ii) = real(ii-1,kind=8)
  enddo
  call camera_sol%allocate_camera(n_times_sol,n_vertices_sol,&
  n_pixels_x_3_sol,n_pixels_y_3_sol,n_spectra_sol) 
  camera_sol%pixel_size = real(1d0/(/n_pixels_x_3_sol,n_pixels_y_3_sol/),kind=8)
  camera_sol%plane_edge_length = plane_edge_length_sol
  !> loop for computing the mesh in local pixel coordinates
  call camera_sol%compute_scaled_pixel_positions_on_plane(x_pixels,y_pixels)
  do ii=1,n_times_sol
    x_pixels(:,ii) = (x_pixels(:,ii)/(camera_sol%pixel_size(1)*&
    camera_sol%plane_edge_length(1,ii))) - x_pixel_ids 
    y_pixels(:,ii) = (y_pixels(:,ii)/(camera_sol%pixel_size(2)*&
    camera_sol%plane_edge_length(2,ii))) - y_pixel_ids
  enddo
  !> test results
  call assert_equals(dot5_x,x_pixels,n_pixels_x_3_sol,n_times_sol,tol_real8,&
  "Error computation scaled pixel positions: x pixel positions mismatch!")
  call assert_equals(dot5_y,y_pixels,n_pixels_y_3_sol,n_times_sol,tol_real8,&
  "Error computation scaled pixel positions: y pixel positions mismatch!")
  !> cleanup
  call camera_sol%deallocate_camera
end subroutine test_compute_scaled_pixel_position_on_plane

!> Tools ------------------------------------------------------------
!>-------------------------------------------------------------------
end module mod_camera_test
