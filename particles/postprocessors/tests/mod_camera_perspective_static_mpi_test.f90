!> mod_camera_perspective_static_mpi_test cantains all variables
!> and procedures used for testing mpi enabled procedures of the
!> camera_perspective_static class.
module mod_camera_perspective_static_mpi_test
use fruit
use fruit_mpi
use constants,                           only: PI,TWOPI
use mod_spectra_deterministic,           only: spectrum_integrator_2nd
use mod_filter_unity,                    only: filter_unity
use mod_pinhole_lens,                    only: pinhole_lens
use mod_camera_perspective_static,       only: camera_perspective_static
use mod_omnidirectional_gaussian_lights, only: omnidirectional_gaussian_lights
implicit none
private
public :: run_fruit_camera_perspective_static_mpi

!> Variable and data types ---------------------------------------------------
character(len=32),parameter                 :: input_file='camera_perspective_static_inputs'
integer,parameter                           :: read_unit=43
integer,parameter                           :: n_st_sol=2
integer,parameter                           :: n_x_sol=3
integer,parameter                           :: n_property_light=2
integer,parameter                           :: n_spectra_sol=2
integer,parameter                           :: n_lines_per_spectrum_sol=53
integer,parameter                           :: n_times=3
integer,parameter                           :: n_pixels_x=64
integer,parameter                           :: n_pixels_y=32
integer,parameter                           :: n_active_light_time_rank=345
integer,dimension(2),parameter              :: rng_seed_interval=(/642,12456324/)
integer,dimension(2),parameter              :: n_inputs_sol=(/5,9/)
integer,dimension(2),parameter              :: mirror_xy_interval=(/0,1/)
real*8,parameter                            :: accept_threshold=9d-1
real*8,parameter                            :: plane_distance_sol=4.23d0
real*8,parameter                            :: dt_sol=1.d-2
real*8,parameter                            :: tol_real8=5.d-16
real*8,parameter                            :: tol_real8_rel=5.d-15
real*8,parameter                            :: light_intensity=1.d13
real*8,dimension(2),parameter               :: st_outbnd=(/-2.d2,5.d3/)
real*8,dimension(2),parameter               :: q_interval=(/-3.1d-1,2.5d-1/)
real*8,dimension(2),parameter               :: weight_interval=(/2d-1,5d-2/)
real*8,dimension(2),parameter               :: property_interval=(/5.d-10,5.d-8/)
real*8,dimension(2),parameter               :: width_gaussian=(/5.d-5,5.d-1/)
real*8,dimension(2),parameter               :: costheta_interval=(/-1.d0,1.d0/)
real*8,dimension(2),parameter               :: phi_interval=(/0.d0,TWOPI/)
real*8,dimension(3),parameter               :: half_angle_lowbnd=(/PI/1.d1,PI/4.d0,0d0/)
real*8,dimension(3),parameter               :: half_angle_uppbnd=(/PI/1.9d0,PI/3.d0,TWOPI/)
real*8,dimension(n_spectra_sol),parameter   :: min_wlen=(/3.d-6,2.5d-7/)
real*8,dimension(n_spectra_sol),parameter   :: max_wlen=(/3.5d-6,4.2d-7/)
real*8,dimension(n_x_sol),parameter         :: center_pos_lowbnd=(/-9.d-1,6.d1,-3.2d0/)
real*8,dimension(n_x_sol),parameter         :: center_pos_uppbnd=(/3.5d0,9.5d1,1.2d0/)
type(filter_unity)                          :: filter_pixel_sol
type(filter_unity)                          :: filter_time_sol
type(spectrum_integrator_2nd)               :: spectra_sol
type(pinhole_lens)                          :: pinhole_sol
type(camera_perspective_static)             :: camera_sol
type(omnidirectional_gaussian_lights)       :: lights_rank
type(filter_unity),dimension(n_spectra_sol) :: filter_spectra_sol
integer                                     :: n_particles_time
integer                                     :: rank_loc,n_tasks_loc,ifail_loc
integer,dimension(5)                        :: int_param_sol
real*8                                      :: pixel_area_sol
real*8,dimension(9)                         :: real_param_sol
real*8,dimension(n_spectra_sol,2,n_pixels_x,n_pixels_y,1) :: image_filter_rank
real*8,dimension(:,:,:,:),allocatable       :: image_sol
real*8,dimension(:,:,:,:,:),allocatable     :: image_filter_sol

!> Interfaces ----------------------------------------------------------------
contains
!> Fruit basket --------------------------------------------------------------
!> fruit basket executes all set-up, test and tear-down procedures
subroutine run_fruit_camera_perspective_static_mpi(rank,n_tasks,ifail)
  use mpi
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  if(rank.eq.0) write(*,*) "  ... setting-up: camera perspective static mpi tests"
  call setup_baseline(rank,n_tasks,ifail)
  call setup_spectra_filters
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  call setup_camera
  call write_camera_input_file
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  call setup_lights
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  call setup_image
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  if(rank.eq.0) write(*,*) "  ... running: camera perspective static mpi tests"
  call run_test_case(test_camera_perspective_static_inputs,&
  'test_camera_perspective_static_inputs')
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  call run_test_case(test_reduce_particle_light_image_static,&
  'test_reduce_particle_light_image_static')
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  call run_test_case(test_compute_images,'test_compute_images')
  call MPI_Barrier(MPI_COMM_WORLD,ifail_loc)
  if(rank.eq.0) write(*,*) "  ... tearing-down: camera perspective static mpi tests"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_camera_perspective_static_mpi

!> Set-up and tear-down ------------------------------------------------------
subroutine setup_baseline(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs
  integer,intent(inout) :: ifail
  !> inputs
  integer,intent(in)    :: rank,n_tasks
  rank_loc = rank; n_tasks_loc = n_tasks; ifail_loc = ifail;
end subroutine setup_baseline

!> set-up spectra and filter types
subroutine setup_spectra_filters()
  implicit none
  !> variables
  integer :: ii,n_1d,n_2d
  !> initialisations
  n_1d = 1; n_2d = 2;
  spectra_sol = spectrum_integrator_2nd(n_lines_per_spectrum_sol,&
  n_spectra_sol,min_wlen,max_wlen)
  call spectra_sol%generate_spectrum
  call filter_time_sol%init_filter(n_1d)
  call filter_pixel_sol%init_filter(n_2d)
  do ii=1,n_spectra_sol
    call filter_spectra_sol(ii)%init_filter(n_1d)
  enddo
end subroutine setup_spectra_filters

!> set-up the pinhole and camera features
subroutine setup_camera
  use mpi
  use mod_gnu_rng,  only: gnu_rng_interval
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  !> variables:
  integer :: ii,n_int_param,n_real_param
  integer,dimension(:),allocatable :: int_param
  real*8,dimension(n_x_sol)        :: center
  real*8,dimension(:),allocatable  :: real_param
  !> initialisation
  n_int_param = 5; n_real_param = 9;
  if(.not.allocated(int_param))  allocate(int_param(n_int_param))
  if(.not.allocated(real_param)) allocate(real_param(n_real_param))
  !> initialise pinhole object 
  if(rank_loc.eq.0) then 
    int_param(1:3) = (/n_lines_per_spectrum_sol,n_pixels_x,n_pixels_y/)
    call gnu_rng_interval(2,mirror_xy_interval,int_param(4:5))
    call gnu_rng_interval(n_x_sol,center_pos_lowbnd,center_pos_uppbnd,center)
  endif 
  call MPI_Bcast(int_param,n_int_param,MPI_INTEGER,0,MPI_COMM_WORLD,ifail_loc)
  call MPI_Bcast(center,n_x_sol,MPI_DOUBLE,0,MPI_COMM_WORLD,ifail_loc)
  call pinhole_sol%init_pinhole(n_x_sol,center)
  !> initialise camera object
  int_param_sol = (/n_lines_per_spectrum_sol,n_pixels_x,n_pixels_y,int_param(4),int_param(5)/)
  if(rank_loc.eq.0) then
    call gnu_rng_interval(3,half_angle_lowbnd,half_angle_uppbnd,real_param_sol(1:3))
    call random_number(real_param_sol(4:6))
    real_param_sol(4:6) = sample_uniform_sphere(plane_distance_sol,&
    costheta_interval,phi_interval,real_param_sol(4:6))
    call gnu_rng_interval(n_x_sol,center_pos_lowbnd,center_pos_uppbnd,real_param_sol(7:9))
    real_param = real_param_sol
  endif
  call MPI_Bcast(real_param,n_real_param,MPI_DOUBLE,0,MPI_COMM_WORLD,ifail_loc)
  call camera_sol%init_camera(pinhole_sol,spectra_sol,n_int_param,&
  n_real_param,int_param,real_param)
  !> compute the pixel area for squared pixels
  pixel_area_sol = (norm2(camera_sol%image_plane(:,2,1)-camera_sol%image_plane(:,1,1))*&
  norm2(camera_sol%image_plane(:,3,1)-camera_sol%image_plane(:,1,1)))/(n_pixels_x*n_pixels_y)
  !> cleanup
  if(allocated(int_param))  deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine setup_camera

!> write the camera input file
subroutine write_camera_input_file
  implicit none
  if(rank_loc.eq.0) then
    open(read_unit,file=input_file,status='unknown',action='write',iostat=ifail_loc)
    write(read_unit,'(/A)') '&camera_in'
    write(read_unit,'(/A,I10)') 'n_lens_samples = ',                     int_param_sol(1)
    write(read_unit,'(/A,I10)') 'n_pixels_x = ',                         int_param_sol(2)
    write(read_unit,'(/A,I10)') 'n_pixels_y = ',                         int_param_sol(3)
    write(read_unit,'(/A,I10)') 'mirror_x = ',                           int_param_sol(4)
    write(read_unit,'(/A,I10)') 'mirror_y = ',                           int_param_sol(5)
    write(read_unit,'(/A,F20.16)') 'image_plane_half_width = ',          real_param_sol(1)
    write(read_unit,'(/A,F20.16)') 'image_plane_half_height = ',         real_param_sol(2)
    write(read_unit,'(/A,F20.16)') 'image_plane_orientation = ',         real_param_sol(3)
    write(read_unit,'(/A,F20.16)') 'image_plane_focal_point_distance = ',real_param_sol(4) 
    write(read_unit,'(/A,F20.16)') 'image_plane_colatitude = ',          real_param_sol(5)
    write(read_unit,'(/A,F20.16)') 'image_plane_azimuth = ',             real_param_sol(6)
    write(read_unit,'(/A,F20.16)') 'focal_point_x_pos = ',               real_param_sol(7)
    write(read_unit,'(/A,F20.16)') 'focal_point_y_pos = ',               real_param_sol(8)
    write(read_unit,'(/A,F20.16)') 'focal_point_z_pos = ',               real_param_sol(9)
    write(read_unit,'(/A)') '/'
    close(read_unit)
  endif
end subroutine write_camera_input_file

!> set-up the light features
subroutine setup_lights
  use mpi
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_particle_common_test_tools, only: fill_particle_simulations_no_init
  use mod_particle_types,             only: particle_kinetic_relativistic_id
  use mod_particle_sim,               only: particle_sim
  implicit none
  !> variables:
  integer :: ii

  !> allocate light properties
  lights_rank%n_property_vertex = n_property_light
  lights_rank%light_intensity = light_intensity
  call lights_rank%allocate_vertices(n_times,n_active_light_time_rank)
  lights_rank%n_active_vertices = n_active_light_time_rank
  if(rank_loc.eq.0) then
    do ii=1,n_times
      call gnu_rng_interval((/(ii-1)*dt_sol,ii*dt_sol/),lights_rank%times(ii))
    enddo
  endif
  call MPI_Bcast(lights_rank%times,n_times,MPI_DOUBLE,0,MPI_COMM_WORLD,ifail_loc)
  !> compute and store the exposure time
  camera_sol%exposure_time = lights_rank%times(n_times)-lights_rank%times(1)
  call fill_omnidirectional_light(rank_loc,n_tasks_loc,ifail_loc) !< prepare particle lights
end subroutine setup_lights

!> set up the solution image
subroutine setup_image
  implicit none
  !> compute the image-filter for each rank
  image_filter_rank = 0.d0
  call compute_contribution_image_plane(rank_loc,n_tasks_loc,ifail_loc)
  !> compute the final image in the rank 0 tasks
  if(rank_loc.eq.0) then
    if(.not.allocated(image_sol)) allocate(image_sol(n_spectra_sol,n_pixels_x,n_pixels_y,1))
    if(.not.allocated(image_filter_sol)) allocate(image_filter_sol(&
    n_spectra_sol,2,n_pixels_x,n_pixels_y,1))
    image_sol = 0.d0; image_filter_sol = 0.d0;
  endif
  call compute_solution_image(rank_loc,n_tasks_loc,ifail_loc)
end subroutine setup_image
 
!> tear-down all test features
subroutine teardown(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  integer :: ii
  call spectra_sol%deallocate_spectrum !< clean spectra
  !> clean filter
  call filter_time_sol%deallocate_filter
  call filter_pixel_sol%deallocate_filter
  do ii=1,n_spectra_sol
    call filter_spectra_sol(ii)%deallocate_filter
  enddo
  call pinhole_sol%deallocate_lens
  call camera_sol%deallocate_camera_perspective_static
  if(allocated(image_sol))        deallocate(image_sol)
  if(allocated(image_filter_sol)) deallocate(image_filter_sol)
  if(rank.eq.0) call system("rm "//input_file)
  rank_loc = -1; n_tasks_loc = -1; ifail = ifail_loc;
end subroutine teardown

!> Tests ---------------------------------------------------------------------
!> test the procedure used for reading lght input files
subroutine test_camera_perspective_static_inputs
  implicit none
  !> variables
  integer,dimension(2) :: n_inputs
  integer,dimension(:),allocatable :: int_param
  real*8,dimension(:),allocatable  :: real_param
  !> read values
  if(rank_loc.eq.0) then
    n_inputs = camera_sol%return_n_camera_inputs()
    open(read_unit,file=input_file,status='old',action='read',iostat=ifail_loc)
    call camera_sol%read_camera_inputs(rank_loc,read_unit,int_param,real_param)
    close(read_unit)
  else
    n_inputs = n_inputs_sol; int_param = int_param_sol;
    real_param = real_param_sol
  endif
  !> checks
  call assert_equals(n_inputs_sol,n_inputs,2,&
  "Error read input camera perspective: N# of inputs mismatch!")
  call assert_equals(int_param_sol,int_param,n_inputs_sol(1),&
  "Error read input camera perspective: integer inputs mismatch!")
  call assert_equals(real_param_sol,real_param,n_inputs_sol(2),tol_real8,&
  "Error read input camera perspective: real inputs mismatch!")
  !> cleanup
  if(allocated(int_param))  deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine test_camera_perspective_static_inputs 

!> test the reduction of all light sources and filters on image plane
subroutine test_reduce_particle_light_image_static
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> inputs:
  real*8,dimension(n_spectra_sol,2,n_pixels_x,n_pixels_y,1) :: error
  real*8,dimension(n_spectra_sol,2,n_pixels_x,n_pixels_y,1) :: pixel_intensities_test
  error = 0.d0
  !> execute reduce particle using omnidirectional light sources
  call camera_sol%reduce_light_image(lights_rank,spectra_sol,filter_pixel_sol,&
  filter_spectra_sol,filter_time_sol,rank_loc,pixel_intensities_test,ifail_loc)
  !> test reduction
  if(rank_loc.eq.0) then
  call assert_equals_rel_error(n_spectra_sol,2,n_pixels_x,n_pixels_y,1,&
  image_filter_sol,pixel_intensities_test,tol_real8_rel,&
  "Error camera reduction particle light image static: pixel intensity mismatch!")
  endif
end subroutine test_reduce_particle_light_image_static

!> test image generation
subroutine test_compute_images
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  real*8,dimension(n_spectra_sol,n_pixels_x,n_pixels_y,1) :: test_image,error
  error = 0.d0
  !> generate image from distribution of omnidirectional light sources
  call camera_sol%compute_images(lights_rank,spectra_sol,filter_pixel_sol,&
  filter_spectra_sol,filter_time_sol,rank_loc,test_image,ifail_loc)
  if(rank_loc.eq.0) call assert_equals_rel_error(n_spectra_sol,n_pixels_x,n_pixels_y,1,&
  image_sol,test_image,tol_real8_rel,"Error camera compute images: images mismatch!")
end subroutine test_compute_images

!> Tools ---------------------------------------------------------------------
!> fill omnidirectional light
subroutine fill_omnidirectional_light(rank,n_tasks,ifail)
  !$ use omp_lib
  use mod_geometry, only: compute_global_cart_coord_plane_points
  use mod_gnu_rng,  only: set_seed_sys_time
  use mod_gnu_rng,  only: gnu_rng_interval
  implicit none
  !> inputs:
  integer,intent(in) :: rank,n_tasks,ifail
  !> variables
  integer                                 :: ii,jj,thread_id
  real*8                                  :: rand,weight
  real*8,dimension(n_x_sol)               :: stq_sol
  real*8,dimension(n_property_light)      :: properties_sol
  real*8,dimension(n_x_sol)               :: x_sol,plane_pos
  real*8,dimension(spectra_sol%n_spectra) :: light_integral
  !> initialise omnidirectional light
  !$omp parallel default(shared) private(ii,jj,thread_id,rand,&
  !$omp x_sol,stq_sol,properties_sol,plane_pos)
  stq_sol=5.d-1; thread_id = 0;
  !$ thread_id = omp_get_thread_num()
  call set_seed_sys_time(rng_seed_interval,rank,thread_id)
  !$omp do collapse(2)
  do ii=1,n_times
    do jj=1,n_active_light_time_rank
      !> random number for deciding shadowed lights
      call random_number(rand)
      if(rand.lt.accept_threshold) then
        call random_number(stq_sol)
      else
        do while(all(stq_sol.ge.0.d0).and.all(stq_sol.le.1.d0))
          call gnu_rng_interval(n_st_sol,st_outbnd,stq_sol(1:n_st_sol))
        end do
      endif
      !> sample a ray lenght
      call gnu_rng_interval(q_interval,stq_sol(n_x_sol))
      !> compute the position on plane
      plane_pos = compute_global_cart_coord_plane_points(&
      camera_sol%image_plane(:,:,ii),stq_sol(1:n_st_sol))
      !> compute and store the ray vertex and its local coordinates
      lights_rank%x(:,jj,ii) = pinhole_sol%center + (plane_pos-pinhole_sol%center)/stq_sol(n_x_sol)
      !> sample the particle weight
      call gnu_rng_interval(weight_interval,weight)
      !> compute light properties
      call gnu_rng_interval(property_interval,lights_rank%properties(1,jj,ii))
      lights_rank%properties(2,jj,ii) = weight/sqrt(2.d0*PI*lights_rank%properties(1,jj,ii))
    enddo
  enddo
  !$omp end do
  !$omp end parallel
end subroutine fill_omnidirectional_light

!> compute light contribution to the image plane
subroutine compute_contribution_image_plane(rank,n_tasks,ifail)
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in)    :: rank,n_tasks
  !> variables:
  logical                                                  :: found
  integer                                                  :: ii,jj,kk
  integer,dimension(n_st_sol)                              :: i_pixel
  real*8 :: material_value,pixel_filter,geometry
  real*8,dimension(n_st_sol)                               :: st_coords
  real*8,dimension(n_x_sol)                                :: stq_coords
  real*8,dimension(n_times)                                :: time_filter
  real*8,dimension(n_lines_per_spectrum_sol,n_spectra_sol) :: spectral_filter
  real*8,dimension(n_lines_per_spectrum_sol,n_spectra_sol) :: irradiance
  real*8,dimension(n_spectra_sol)                          :: light_integral
  !> initialise time filter
  call filter_time_sol%compute_filter_from_position_vectorial(&
  n_times,lights_rank%times,time_filter)
  !> initialise spectral filter
  do kk=1,n_spectra_sol
    call filter_spectra_sol(kk)%compute_filter_from_position_vectorial(&
    n_lines_per_spectrum_sol,spectra_sol%points(:,kk),spectral_filter(:,kk))
  enddo
  !> compute and integrate light
  do kk=1,camera_sol%n_vertices
    do jj=1,n_times
      do ii=1,n_active_light_time_rank
        !> compute the material function
        call camera_sol%physical_material_funct(lights_rank%x(:,ii,jj),kk,material_value)
        if(material_value.le.0.d0) cycle
        !> intersection with the image plane
        call camera_sol%find_ray_image_plane_intersection(lights_rank%x(:,ii,jj),kk,&
        found,stq_coords)
        if(.not.found) cycle
        call camera_sol%plane_to_pixel_local_coord(stq_coords(1:2),i_pixel,st_coords)
        !> compute pixel filter
        call filter_pixel_sol%compute_filter_from_position(st_coords,pixel_filter)
        !> compute distance
        geometry = dot_product(lights_rank%x(:,ii,jj)-camera_sol%x(:,kk,1),&
        lights_rank%x(:,ii,jj)-camera_sol%x(:,kk,1))
        !> compute the spectral irradiance
        call lights_rank%spectral_irradiance(spectra_sol,jj,ii,camera_sol%x(:,kk,1),irradiance)
        !> integrate the spectral irradiance
        call spectra_sol%integrate_data(irradiance*spectral_filter,light_integral)
        !> reduction of the irradiance and the filter function
        image_filter_rank(:,1,i_pixel(1),i_pixel(2),1) = &
        image_filter_rank(:,1,i_pixel(1),i_pixel(2),1) + (light_integral*material_value*&
        pixel_area_sol*pixel_filter*time_filter(jj)*camera_sol%exposure_time)/geometry
        image_filter_rank(:,2,i_pixel(1),i_pixel(2),1) = &
        image_filter_rank(:,2,i_pixel(1),i_pixel(2),1) + time_filter(jj)*pixel_filter
      enddo
    enddo
  enddo
end subroutine compute_contribution_image_plane 

subroutine compute_solution_image(rank,n_tasks,ifail)
  use mpi
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  !> variables:
  integer :: ii,n_count,stat
  real*8,dimension(:,:,:,:,:,:),allocatable   :: image_buffer
  !> receive all contributions from all tasks
  n_count = n_spectra_sol*2*n_pixels_x*n_pixels_y*1
  allocate(image_buffer(n_spectra_sol,2,n_pixels_x,n_pixels_y,1,n_tasks))
  image_buffer = 0.d0
  !> retrive all the image planes
  call MPI_Gather(image_filter_rank,n_count,MPI_REAL8,image_buffer,n_count,&
  MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  if(rank.eq.0) then
    do ii=1,n_tasks
      image_filter_sol = image_filter_sol + image_buffer(:,:,:,:,:,ii)
    enddo
    where(image_filter_sol(:,2,:,:,:).ne.0.d0) &
    image_sol = image_filter_sol(:,1,:,:,:)/image_filter_sol(:,2,:,:,:)
  endif
end subroutine compute_solution_image

!>----------------------------------------------------------------------------
end module mod_camera_perspective_static_mpi_test

