!> mod_camera_perspective_static contains all variables and
!> procedures defining a perspective camera which does not
!> vary with time (static)
module mod_camera_perspective_static
use mod_camera_static, only: camera_static
implicit none
  
private
public :: camera_perspective_static

!> Variables and type definitions ----------------------------
type,extends(camera_static) :: camera_perspective_static
  contains
  procedure,pass(camera_in)    :: return_n_camera_inputs => &
  return_n_camera_perspective_static_inputs
  procedure,pass(camera_inout) :: read_camera_inputs => &
  read_camera_perspective_static_inputs
  procedure,pass(camera_inout) :: init_camera => init_camera_perspective_static
  procedure,pass(camera_inout) :: reduce_light_image => &
  reduce_particle_light_image_static
  procedure,pass(camera_inout) :: physical_material_funct => &
  physical_material_funct_perspective_static
  procedure,pass(camera_inout) :: allocate_camera_perspective_static
  procedure,pass(camera_inout) :: deallocate_camera_perspective_static
  procedure,pass(camera_inout) :: cos_view_angle_static
end type camera_perspective_static

!> Interfaces ------------------------------------------------
contains
!> Procedures ------------------------------------------------
!> reduce the direct illumination, geometric and structural term 
!> from a set of particles in physical units.
!> inputs:
!>   camera_inout:         (camera) initialised camera perspective static
!>   light_inout:          (light_vertices) initialised particle light vertices
!>   spectra_inout:        (spectrum_inout) camera spectrum object
!>   filter_image_inout:   (filter) 2d image plane filter
!>   filter_spectra_inout: (filter)(n_spectra) set of 1d spectral filters
!>   filter_time_inout:    (filter) filter in time
!>   rank:                 (integer) mpi rank
!>   ierr:                 (integer) error for the mpi procedures
!> outputs:
!>   camera_inout:         (camera) camera perspective static with reduced light
!>   light_inout:          (light_vertices) initialised particle light vertices
!>   spectra_inout:        (spectrum_inout) camera spectrum object
!>   filter_image_inout:   (filter) 2d image plane filter
!>   filter_spectra_inout: (filter)(n_spectra) set of 1d spectral filters
!>   filter_time_inout:    (filter) filter in time
!>   pixel_intensities:    (real8)(n_spectra,2,n_pixels_x,n_pixels_y,n_times) array
!>                         containing the reduction of the pixel intensities (:,1,:,:,:)
!>                         and the reduction of the filter values (:,2,:,:,:)
!>   ierr:                 (integer) error for the mpi procedures
subroutine reduce_particle_light_image_static(camera_inout,light_inout,&
spectra_inout,filter_image_inout,filter_spectra_inout,filter_time_inout,rank,&
pixel_intensities,ierr)
  use mpi
  use mod_light_vertices, only: light_vertices
  use mod_spectra,        only: spectrum_base
  use mod_filter,         only: filter
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout)                 :: camera_inout
  class(light_vertices),intent(inout)                            :: light_inout
  class(spectrum_base),intent(inout)                             :: spectra_inout
  class(filter),intent(inout)                                    :: filter_image_inout
  class(filter),dimension(spectra_inout%n_spectra),intent(inout) :: filter_spectra_inout
  class(filter),intent(inout)                                    :: filter_time_inout
  integer,intent(inout)                                          :: ierr
  !> inputs:
  integer,intent(in) :: rank
  !> outputs:
  real*8,dimension(camera_inout%n_pixels_spectra(1),2,camera_inout%n_pixels_spectra(2),&
  camera_inout%n_pixels_spectra(3),camera_inout%n_times),intent(out) :: pixel_intensities
  !> variables
  logical :: intersect
  integer :: ii,jj,kk
  integer,dimension(2) :: i_pixel
  real*8  :: material_value,pixel_filter_weight,visibility_geometry,pixel_area
  real*8,dimension(2) :: pixel_coords
  real*8,dimension(3) :: plane_line_coords
  real*8,dimension(light_inout%n_times) :: light_time_weights !< time filter
  real*8,dimension(spectra_inout%n_spectra) :: integrated_irradiance
  real*8,dimension(spectra_inout%n_points,spectra_inout%n_spectra) :: spectral_weights !< spectral filter
  real*8,dimension(spectra_inout%n_points,spectra_inout%n_spectra) :: spectral_irradiance

  !> compute and store the exposure time and the pixel area
  pixel_intensities = 0.d0; spectral_irradiance = 0d0; integrated_irradiance = 0d0;
  camera_inout%exposure_time = light_inout%times(light_inout%n_times) - light_inout%times(1)
  if(camera_inout%exposure_time.le.0.d0) camera_inout%exposure_time = 1.d0
  pixel_area = (norm2(camera_inout%image_plane(:,2,1)-camera_inout%image_plane(:,1,1))*&
               norm2(camera_inout%image_plane(:,3,1)-camera_inout%image_plane(:,1,1)))/&
               (real(camera_inout%n_pixels_spectra(2)*camera_inout%n_pixels_spectra(3),kind=8))
  !> compute the time filter weights
  call filter_time_inout%compute_filter_from_position_vectorial(light_inout%n_times,&
  light_inout%times,light_time_weights)
  !> compute the spectral filter weights
  do ii=1,spectra_inout%n_spectra
    call filter_spectra_inout(ii)%compute_filter_from_position_vectorial(&
    spectra_inout%n_points,spectra_inout%points(:,ii),spectral_weights(:,ii))
  enddo

  !> loop on the light times
  do jj=1,light_inout%n_times
    !$omp parallel do default(none) firstprivate(jj,pixel_area,&
    !$omp light_time_weights,spectral_weights,spectral_irradiance,&
    !$omp integrated_irradiance) &
    !$omp shared(light_inout,camera_inout,spectra_inout,&
    !$omp filter_image_inout,filter_spectra_inout,filter_time_inout) &
    !$omp private(ii,kk,material_value,intersect,plane_line_coords,&
    !$omp i_pixel,pixel_coords,pixel_filter_weight,visibility_geometry) &
    !$omp reduction(+:pixel_intensities) collapse(2)
    !> loop on the points on lens
    do ii=1,camera_inout%n_vertices
      !> loop on the active lights per time
      do kk=1,light_inout%n_active_vertices(jj)
        !> compute the physical material function (skip if negative because light behind camera)
        call camera_inout%physical_material_funct(light_inout%x(:,kk,jj),ii,material_value)
        if(material_value.le.0.d0) cycle
        !> compute the intersection of the camera-light ray with the image plane 
        call camera_inout%find_ray_image_plane_intersection(light_inout%x(:,kk,jj),ii,&
        intersect,plane_line_coords)
        if(.not.intersect) cycle
        !> compute the intersection point pixel coordinates
        call camera_inout%plane_to_pixel_local_coord(plane_line_coords(1:2),&
        i_pixel,pixel_coords)
        !> compute the pixel filter weight
        call filter_image_inout%compute_filter_from_position(pixel_coords,pixel_filter_weight)
        !> compute the geometry and visibility functions
        call camera_inout%visibility_geometry_funct(light_inout,1,jj,ii,kk,visibility_geometry)
        !> compute the spectral irradiance
        call light_inout%spectral_irradiance(spectra_inout,jj,kk,camera_inout%x(:,ii,1),&
        spectral_irradiance)
        !> integrate the weighted spectral irradiance for each spectral interval
        call spectra_inout%integrate_data(spectral_weights*spectral_irradiance,integrated_irradiance)
        !> accumulate the irradiance of the pixel
        pixel_intensities(:,1,i_pixel(1),i_pixel(2),1) = &
        pixel_intensities(:,1,i_pixel(1),i_pixel(2),1) + &
        (light_time_weights(jj)*pixel_filter_weight*integrated_irradiance*&
        visibility_geometry*material_value*pixel_area*camera_inout%exposure_time)
        !> accumulate the overall filter functions of the pixel (E.Veach, PhD thesis, 1997)
        pixel_intensities(:,2,i_pixel(1),i_pixel(2),1) = &
        pixel_intensities(:,2,i_pixel(1),i_pixel(2),1) + &
        light_time_weights(jj)*pixel_filter_weight
      enddo
    enddo
    !$omp end parallel do
  enddo
  !> reduce all mpi images into the root image
  if(rank.eq.0) then
    call MPI_Reduce(MPI_IN_PLACE,pixel_intensities,&
    size(pixel_intensities),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  else
    call MPI_Reduce(pixel_intensities,pixel_intensities,&
    size(pixel_intensities),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  endif
end subroutine reduce_particle_light_image_static

!> procedure use for reading the camera perspective static inputs
!> inputs:
!>   camera_inout: (camera_perspective_static) camera object
!>   my_id:        (integer) mpi rank
!>   r_unit:       (integer) file unit id
!> outputs:
!>   camera_inout: (camera_perspective_static) camera object
!>   int_param:    (integer)(:) integer input parameters
!>   real_param:    (real8)(:) real input parameters
subroutine read_camera_perspective_static_inputs(camera_inout,&
my_id,r_unit,int_param,real_param)
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> inputs:
  integer,intent(in) :: my_id,r_unit
  !> outputs:
  integer,dimension(:),allocatable,intent(out) :: int_param
  real*8,dimension(:),allocatable,intent(out)  :: real_param
  !> variables:
  integer :: n_lens_samples,n_pixels_x,n_pixels_y,mirror_x,mirror_y
  integer,dimension(2) :: n_inputs
  real*8  :: image_plane_half_width,image_plane_half_height
  real*8  :: image_plane_orientation,image_plane_focal_point_distance
  real*8  :: image_plane_colatitude,image_plane_azimuth
  real*8  :: focal_point_x_pos,focal_point_y_pos,focal_point_z_pos
  !> definitions and initialisations
  namelist /camera_in/ n_lens_samples,n_pixels_x,n_pixels_y,mirror_x,mirror_y,&
                       image_plane_half_width,image_plane_half_height,&
                       image_plane_orientation,image_plane_focal_point_distance,&
                       image_plane_colatitude,image_plane_azimuth,&
                       focal_point_x_pos,focal_point_y_pos,focal_point_z_pos

  !> if master, read input file
  if(my_id.eq.0) then
    n_inputs = camera_inout%return_n_camera_inputs()
    if(allocated(int_param)) deallocate(int_param)
    allocate(int_param(n_inputs(1)))
    if(allocated(real_param)) deallocate(real_param)
    allocate(real_param(n_inputs(2)))
    read(r_unit,camera_in)
    int_param = (/n_lens_samples,n_pixels_x,n_pixels_y,mirror_x,mirror_y/)
    real_param = (/image_plane_half_width,image_plane_half_height,&
    image_plane_orientation,image_plane_focal_point_distance,&
    image_plane_colatitude,image_plane_azimuth,focal_point_x_pos,&
    focal_point_y_pos,focal_point_z_pos/)
  endif
end subroutine read_camera_perspective_static_inputs

!> return the input array size
!> inputs:
!>   camera_in: (camera_perspective_static) camera object 
!> outputs:
!>   n_inputs: (integer)(2) size of the integer 
!>             and real input arrays
function return_n_camera_perspective_static_inputs(&
camera_in) result(n_inputs)
  implicit none
  !> inputs:
  class(camera_perspective_static),intent(in) :: camera_in
  !> outputs:
  integer,dimension(2) :: n_inputs
  n_inputs = (/5,9/)
end function return_n_camera_perspective_static_inputs

!> procedure used for initialising a static perspective camera
!> inputs:
!>   camera_inout:   (camera) unallocated camera
!>   lens_inout:     (lens) camera lens object
!>   spectrum_inout: (spectrum_inout) camera spectrum object
!>   n_int_param:    (integer) number of integer parameter: 3
!>   n_real_param:   (integer) number of real parameter: 8
!>   int_praram:     (integer)(n_int_param) integer parameters:
!>                            1) number of lens samples
!>                            2) number of pixels in the x-direction
!>                            3) number of pixels in the y-direction
!>                            4) mirror image plane along x if 1
!>                            5) mirror image plane along y if 1
!>  real_param:      (real8)(n_int_param) real parameters:
!>                            1:3) image plane half widht, half
!>                                 height and orientation  angles 
!>                                 in the focal reference
!>                            4:6) distance betweem the image plane
!>                                 and the camera focal point
!>                            7:9) camera focal point position
!> outputs:
!>   camera_inout:   (camera) allocated camera
!>   lens_inout:     (lens) camera lens object
!>   spectrum_inout: (spectrum_inout) camera spectrum object
subroutine init_camera_perspective_static(camera_inout,lens_inout,&
spectrum_inout,n_int_param,n_real_param,int_param,real_param)
  use mod_spectra,  only: spectrum_base
  use mod_geometry, only: define_plane_from_half_angles
  use mod_lens,     only: lens
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout) :: camera_inout
  class(lens),intent(inout)                      :: lens_inout
  class(spectrum_base),intent(inout)             :: spectrum_inout
  !> inputs:
  integer,intent(in)                             :: n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in)      :: int_param
  real*8,dimension(n_real_param),intent(in)      :: real_param

  !> set variables
  camera_inout%n_times=1;
  camera_inout%n_property_vertex=1; camera_inout%n_plane_points=3;
  !> initialise attributes
  !> lens samples are stored as x positions of the
  !> vertices while their pdfs as property
  call camera_inout%allocate_camera_perspective_static(spectrum_inout,&
  n_int_param,int_param)
  !> define the image plane characteristics
  call camera_inout%define_image_plane_pixel_size(n_int_param,n_real_param,&
  int_param,real_param)
  !> sample the lens
  call camera_inout%generate_points_on_lens_pdf(lens_inout,int_param(1))
end subroutine init_camera_perspective_static

!> procedure used for allocating all attributes of camera perspective static
!> inputs:
!>   camera_inout:   (camera) unallocated camera
!>   spectrum_inout: (spectrum_inout) camera spectrum object
!>   n_int_param:    (integer) number of integer parameter: 3
!>   int_praram:     (integer)(n_int_param) integer parameters:
!>                            1) number of lens samples
!>                            2) number of pixels in the x-direction
!>                            3) number of pixels in the y-direction
!> outputs:
!>   camera_inout:   (camera) allocated camera
!>   spectrum_inout: (spectrum_inout) camera spectrum object
subroutine allocate_camera_perspective_static(camera_inout,&
spectrum_inout,n_int_param,int_param)
  use mod_spectra, only: spectrum_base
  implicit none
  !> inputs-outputs
  class(camera_perspective_static),intent(inout) :: camera_inout
  class(spectrum_base),intent(inout)             :: spectrum_inout
  !> inputs
  integer,intent(in)                        :: n_int_param
  integer,dimension(n_int_param),intent(in) :: int_param
  !> allocate camera base type
  call camera_inout%allocate_camera(1,int_param(1),&
  int_param(2),int_param(3),spectrum_inout%n_spectra)
end subroutine allocate_camera_perspective_static

!> procedure used for deallocating all attributes of camera perspective static
!> inputs:
!>   camera_inout: (camera) allocated camera
!> outputs:
!>   camera_inout: (camera) deallocated camera
subroutine deallocate_camera_perspective_static(camera_inout)
  implicit none
  !> inputs-outputs
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> deallocate variables
  call camera_inout%deallocate_camera
end subroutine deallocate_camera_perspective_static

!> compute the physical importance function for the camera perspective static
!> inputs:
!>   camera_inout:   (camera_perspective_static) initialised camera object
!>   x_pos:          (real8)(3) coordinated of the point defining a ray
!>   x_lens_id:      (integer) index of the point on lens to be treated
!>   time_id_in:     (integer)(optional) time index (not used)
!> outputs:
!>   camera_inout:   (camera_perspective_static) camera object
!>   material_value: (real8) camera physical importance
subroutine physical_material_funct_perspective_static(camera_inout,x_pos,x_lens_id,&
material_value,time_id_in)
  implicit none
  !> inputs-outputs:
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> inputs:
  integer,intent(in)                            :: x_lens_id
  integer,intent(in),optional                   :: time_id_in
  real*8,dimension(camera_inout%n_x),intent(in) :: x_pos
  !> outputs
  real*8,intent(out)                            :: material_value
  !> compute the physical material function which at the moment it is 
  !> only the cosinus of the view angle
  call camera_inout%cos_view_angle_static(x_pos,x_lens_id,material_value)
end subroutine physical_material_funct_perspective_static

!> compute the cosinus of the view angle for static cameras.
!> the cosinus of the view angle is defined as the scalar product between
!> a ray and the lens-image plane vector. The ray is defined by a position
!> on the lens and a given position in the space
!> inputs:
!>   camera_inout: (camera_perspective_static) initialised camera object
!>   x_pos:     (real8)(3) coordinated of the point defining a ray
!>   x_lens_id: (integer) index of the point on lens to be treated
!> outputs:
!>   camera_inout:   (camera_perspective_static) camera object
!>   cos_view_angle: (real8) cosinus of the view angle. Negative values mean
!>                           that the point is behind the lens
subroutine cos_view_angle_static(camera_inout,x_pos,x_lens_id,cos_view_angle)
  implicit none
  !> inputs-outputs
  class(camera_perspective_static),intent(inout) :: camera_inout
  !> inputs
  integer,intent(in)                             :: x_lens_id
  real*8,dimension(camera_inout%n_x),intent(in)  :: x_pos
  !> outputs
  real*8,intent(out)                             :: cos_view_angle
  !> variables
  real*8,dimension(camera_inout%n_x)             :: ray
  !> compute the cosinus of the view angle
  ray = x_pos - camera_inout%x(:,x_lens_id,1)
  ray = ray/sqrt(ray(1)*ray(1)+ray(2)*ray(2)+ray(3)*ray(3))
  cos_view_angle = ray(1)*camera_inout%image_plane_direction(1,1) + &
                   ray(2)*camera_inout%image_plane_direction(2,1) + &
                   ray(3)*camera_inout%image_plane_direction(3,1)
end subroutine cos_view_angle_static

!>------------------------------------------------------------
end module mod_camera_perspective_static
