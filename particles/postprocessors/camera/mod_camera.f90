!> mod_camera extends mod_vertices containing the basic
!> variables, datatypes and procedures common to the
!> different camera models. Each camera is supposed to have
!> one pupil and one visual plane hence, camera systems having
!> multiple pupils with multiple visual planes must be
!> simulated as different synthetic cameras
module mod_camera
use mod_vertices, only: vertices
implicit none

private
public :: camera

!> Variable and type definitions -----------------------------------
type,abstract,extends(vertices) :: camera
   !> the x variables contains the global position of the 
   !> visual points on the lens / film for each frame
   !> number of spectra and pixels
   integer              :: n_plane_points !< number of points of a plane
   integer,dimension(3) :: n_pixels_spectra
   real*8,dimension(2)  :: pixel_size !< width and height of a pixel 
   real*8,dimension(:,:),allocatable   :: plane_edge_length !< length of the image plane edges
   real*8,dimension(:,:,:),allocatable :: image_plane !< vertices of the image plane
   real*8,dimension(:,:),allocatable   :: image_plane_direction !< view direction of the image plane
   real*8               :: exposure_time !< exposure time for each camera frame
  contains
   procedure(int_ret_n_cam_inps),pass(camera_in),deferred       :: return_n_camera_inputs
   procedure(int_read_cam_inps),pass(camera_inout),deferred     :: read_camera_inputs
   procedure(int_init_camera),pass(camera_inout),deferred       :: init_camera
   procedure(int_gen_points_lens),pass(camera_inout),deferred   :: generate_points_on_lens_pdf
   procedure(int_reduce_light_img),pass(camera_inout),deferred  :: reduce_light_image
   procedure(int_phys_mat_funct),pass(camera_inout),deferred    :: physical_material_funct
   procedure(int_ray_plane_inter),pass(camera_inout),deferred   :: find_ray_image_plane_intersection
   procedure,pass(camera_inout)                                 :: allocate_camera
   procedure,pass(camera_inout)                                 :: deallocate_camera
   procedure,pass(camera_inout)                                 :: compute_images
   procedure,pass(camera_inout)                                 :: plane_to_pixel_local_coord
   procedure,pass(camera_inout)                                 :: pixel_local_to_scaled_plane_coord
   procedure,pass(camera_inout)                                 :: compute_scaled_pixel_positions_on_plane
end type camera

!> Interfaces ------------------------------------------------------
interface
  !> interfarce for procedure used for reading the camera inputs
  !> inputs:
  !>   camera_inout: (camera) camera object
  !>   my_id:        (integer) mpi rank
  !>   r_unit:       (integer) file unit id
  !> outputs:
  !>   camera_inout: (camera) camera object
  !>   int_param:    (integer)(:) integer input parameters
  !>   real_param:   (real8)(:) real input parameters
  subroutine int_read_cam_inps(camera_inout,my_id,r_unit,&
  int_param,real_param)
  IMPORT :: camera
  implicit none
  !> inputs-outputs:
  class(camera),intent(inout) :: camera_inout
  !> inputs:
  integer,intent(in) :: my_id,r_unit
  !> outputs:
  integer,dimension(:),allocatable,intent(out) :: int_param
  real*8,dimension(:),allocatable,intent(out)  :: real_param
  end subroutine int_read_cam_inps

  !> interface of the function used for returning the input
  !> array size
  !> inputs:
  !>   camera_in: (camera) input camera object
  !> outputs:
  !>   n_inputs:  (integer)(2) size of the integer and real
  !>              input arrays
  function int_ret_n_cam_inps(camera_in) result(n_inputs)
    IMPORT :: camera
    !> inputs:
    class(camera),intent(in) :: camera_in
    !> outputs:
    integer,dimension(2) :: n_inputs
  end function int_ret_n_cam_inps

  !> interface of the deferred procedure init_camera which initialises a 
  !> camera type.
  !> inputs:
  !>   camera_inout:  (camera) camera to be initialised
  !>   lens_inout:    (lens) the lens to be used for the camera
  !>   spectrum_inout (spectrum base) camera sensitivity spectra
  !>   n_int_param:   (integer) number of integer parameters to be passed
  !>   n_real_param:  (integer) number of real parameters to be passed
  !>   int_param:     (integer)(n_int_param) integer parameter array
  !>   real_param:    (real8)(n_real_param) real parameter array
  !> outputs:
  !>   camera_inout:  (camera) initialised camera
  subroutine int_init_camera(camera_inout,lens_inout,spectrum_inout,&
  n_int_param,n_real_param,int_param,real_param)
    use mod_spectra, only: spectrum_base
    use mod_lens,    only: lens
    IMPORT :: camera
    implicit none
    !> inputs-outputs
    class(camera),intent(inout)        :: camera_inout
    class(lens),intent(inout)          :: lens_inout
    class(spectrum_base),intent(inout) :: spectrum_inout
    !> inputs
    integer,intent(in)                        :: n_int_param,n_real_param
    integer,dimension(n_int_param),intent(in) :: int_param
    real*8,dimension(n_real_param),intent(in) :: real_param
  end subroutine int_init_camera

  !> generate a set of points on a lens
  !> inputs:
  !>   camera_inout: (camera) camera model used for sampling
  !>   lens_inout:   (lens) lens from which points are sampled
  !>   n_points_in:  (integer)(optional) number of points to sample
  !>                                     default: 1000, pinhole: 1
  !> outputs:
  !>   camera_inout: (camera) camera with sampled points on lens
  !>   lens_inout:   (lens) lens from which points are sampled
  subroutine int_gen_points_lens(camera_inout,lens_inout,n_points_in)
    use mod_lens, only: lens
    IMPORT :: camera
    implicit none
    !> inputs-outputs:
    class(camera),intent(inout) :: camera_inout
    class(lens),intent(inout)   :: lens_inout
    !> inputs:
    integer,intent(in),optional :: n_points_in
  end subroutine int_gen_points_lens

  !> reduce the particle lights on the image plane
  !> inputs:
  !>   camera_inout:         (camera) initialised camera
  !>   light_inout:          (light_vertices) initialised particle light vertices
  !>   spectra_inout:        (spectrum_base) camera spectrum object
  !>   filter_image_inout:   (filter) 2d image plane filter
  !>   filter_spectra_inout: (filter)(n_spectra) set of 1d spectral filter
  !>   filter_time_inout:    (filter) filter in time
  !>   rank:                 (integer) mpi rank
  !>   ierr:                 (integer) error for the mpi procedures
  !> outputs:
  !>   camera_inout:         (camera) camera with reduced light
  !>   light_inout:          (light_vertices) initialised particle light vertices
  !>   spectra_inout:        (spectrum_base) camera spectrum object
  !>   filter_image_inout:   (filter) 2d image plane filter
  !>   filter_spectra_inout: (filter)(n_spectra) set of 1d spectral filter
  !>   filter_time_inout:    (filter) filter in time
  !>   pixel_intensities:    (real8)(n_spectra,2,n_pixels_x,n_pixels_y,n_times) array
  !>                         containing the reduction of the pixel intensities (:,1,:,:,:)
  !>                         and the reduction of the filter values (:,2,:,:,:)
  !>   ierr:                 (integer) error for the mpi procedures
  subroutine int_reduce_light_img(camera_inout,light_inout,spectra_inout,&
  filter_image_inout,filter_spectra_inout,filter_time_inout,rank,pixel_intensities,ierr)
  use mod_light_vertices, only: light_vertices
  use mod_spectra,        only: spectrum_base
  use mod_filter,         only: filter
  IMPORT :: camera
  implicit none
  !> imputs-outputs
  class(camera),intent(inout)                                    :: camera_inout
  class(light_vertices),intent(inout)                            :: light_inout
  class(spectrum_base),intent(inout)                             :: spectra_inout
  class(filter),intent(inout)                                    :: filter_image_inout
  class(filter),dimension(spectra_inout%n_spectra),intent(inout) :: filter_spectra_inout
  class(filter),intent(inout)                                    :: filter_time_inout
  integer,intent(inout)                                          :: ierr
  !> inputs:
  integer,intent(in) :: rank
  !> outputs:
  real*8,dimension(camera_inout%n_pixels_spectra(1),2,&
  camera_inout%n_pixels_spectra(2),camera_inout%n_pixels_spectra(3),&
  camera_inout%n_times),intent(out) :: pixel_intensities
  end subroutine int_reduce_light_img

  !> structure function: it returns the physical camera importance normalised to 1
  !> inputs:
  !>   camera_inout:   (camera) camera with reduced light
  !>   x_pos:          (real8)(3) position used for defining the importance ray
  !>   x_lesn_id:      (integer) index of points on the lens
  !>   time_id_in:     (integer)(optional) index of the treated time
  !> outputs:
  !>   camera_inout:   (camera) camera with reduced light
  !>   material_value: (real8) physical camera physical importance
  subroutine int_phys_mat_funct(camera_inout,x_pos,x_lens_id,&
  material_value,time_id_in)
    IMPORT :: camera
    implicit none
    !> inputs-outputs:
    class(camera),intent(inout)                   :: camera_inout
    !> inputs:
    integer,intent(in)                            :: x_lens_id
    integer,intent(in),optional                   :: time_id_in
    real*8,dimension(camera_inout%n_x),intent(in) :: x_pos
    !> outputs:
    real*8,intent(out)                            :: material_value
  end subroutine int_phys_mat_funct

  !> find the intersection between a camera ray and the image plane
  !> inputs:
  !>   camera_inout: (camera) camera object
  !> outputs:
  !>   camera_inout: (camera) camera object
  !>   intersect:    (bool) if true an intersection is found
  !>   local_coords: (real8)(3) plane (s,t) and ray (q) local
  !>                            coordinates of the intersection
  subroutine int_ray_plane_inter(camera_inout,x_pos,&
  x_lens_id,intersect,local_coords,time_id_in)
    use mod_geometry, only: compute_plane_line_intersection_cart
    IMPORT :: camera
    implicit none
    !> inputs-outputs:
    class(camera),intent(inout)                   :: camera_inout
    !> inputs:
    integer,intent(in)                            :: x_lens_id
    real*8,dimension(camera_inout%n_x),intent(in) :: x_pos
    integer,intent(in),optional                   :: time_id_in
    !> outputs:
    logical,intent(out)                           :: intersect
    real*8,dimension(3),intent(out)               :: local_coords
  end subroutine int_ray_plane_inter
end interface

contains
!> Procedures ------------------------------------------------------
!> allocate camera arrays and initialise them to zero
!> inputs:
!>   camera_inout: (camera) camera to be allocated
!>   n_times:      (integer) number of times
!>   n_vertices:   (integer) number of vertices
!>   n_pixels_x:   (integer) number of pixels in the x-direction
!>   n_pixels_y:   (integer) number of pixels in the y-direction
!>   n_spectra:    (integer) number of spectrum intervals
!> outputs:
!>   camera_inout: (camera) allocated camera
subroutine allocate_camera(camera_inout,n_times,n_vertices,&
n_pixels_x,n_pixels_y,n_spectra)
  implicit none
  !> inputs-outputs:
  class(camera),intent(inout) :: camera_inout
  !> inputs:
  integer,intent(in) :: n_times,n_vertices
  integer,intent(in) :: n_spectra,n_pixels_x,n_pixels_y
  
  !> allocate vertices
  call camera_inout%allocate_vertices(n_times,n_vertices)
  if(allocated(camera_inout%plane_edge_length)) then
    if(size(camera_inout%plane_edge_length,2).ne.n_times) deallocate(camera_inout%plane_edge_length)
  endif
  if(.not.allocated(camera_inout%plane_edge_length)) allocate(camera_inout%plane_edge_length(2,n_times))
  if(allocated(camera_inout%image_plane_direction)) then
    if(size(camera_inout%image_plane_direction,2).ne.n_times) &
    deallocate(camera_inout%image_plane_direction)
  endif
  if(.not.allocated(camera_inout%image_plane_direction)) &
  allocate(camera_inout%image_plane_direction(camera_inout%n_x,n_times))
  if(allocated(camera_inout%image_plane)) then
    if(.not.((size(camera_inout%image_plane,2).eq.camera_inout%n_plane_points).and.&
    (size(camera_inout%image_plane,3).eq.n_times))) deallocate(camera_inout%image_plane)
  endif
  if(.not.allocated(camera_inout%image_plane)) &
  allocate(camera_inout%image_plane(camera_inout%n_x,camera_inout%n_plane_points,n_times))
  !> allocate camera
  camera_inout%exposure_time = 0d0; camera_inout%plane_edge_length = 0d0;
  camera_inout%image_plane_direction = 0d0; camera_inout%image_plane = 0d0;
  camera_inout%pixel_size = 0d0;
  camera_inout%n_pixels_spectra = (/n_spectra,n_pixels_x,n_pixels_y/)
end subroutine allocate_camera

!> deallocate camera arrays and reset counters to 0
!> inputs:
!>   camera_inout: (camera) camera to deallocated
!> outputs:
!>   camera_inout: (camera) deallocated camera
subroutine deallocate_camera(camera_inout)
  implicit none
  !> inputs-outputs
  class(camera),intent(inout) :: camera_inout
  !> deallocate everything and reset counters
  call camera_inout%deallocate_vertices
  if(allocated(camera_inout%plane_edge_length))     deallocate(camera_inout%plane_edge_length)
  if(allocated(camera_inout%image_plane_direction)) deallocate(camera_inout%image_plane_direction)
  if(allocated(camera_inout%image_plane))           deallocate(camera_inout%image_plane)
  camera_inout%n_pixels_spectra = 0; camera_inout%exposure_time = 0.d0;
  camera_inout%n_plane_points = 0;  camera_inout%pixel_size = 0d0;
end subroutine deallocate_camera

!> compute the pixel number and the position in the pixel local coordinates
!> of a point on the image plane (in the plane local coordinates)
!> inputs:
!>  camera_inout: (camera) camera vertices
!>  st_plane:     (real8)(2) position in the plane local coordinates
!> outputs:
!>  camera_inout: (camera) camera vertices
!>  i_pixel:      (integer)(2) pixel indices of the point on the plane (s,t)
!>  st_pixel:     (real8)(2) position in the pixel local coordinates
subroutine plane_to_pixel_local_coord(camera_inout,st_plane,i_pixel,st_pixel)
  implicit none
  !> inputs-outputs
  class(camera),intent(inout)      :: camera_inout
  !> inputs
  real*8,dimension(2),intent(in)   :: st_plane
  !> outputs
  integer,dimension(2),intent(out) :: i_pixel
  real*8,dimension(2),intent(out)  :: st_pixel
  !> find the local pixel coordinates and find the position in the pixel local coordinates
  st_pixel = st_plane/camera_inout%pixel_size
  i_pixel = floor(st_pixel)
  st_pixel = st_pixel - real(i_pixel,kind=8)
  i_pixel = i_pixel + 1
end subroutine plane_to_pixel_local_coord

!> compute the coordinate in the scaled plane given the local pixel coordinates
!> inputs:
!>  camera_inout: (camera) camera vertices
!>  id_time:      (integer) index of the considered time
!>  i_pixel:      (integer)(2) pixel indices of the point on the plane (s,t)
!>  st_pixel:     (real8)(2) position in the pixel local coordinates
!> outputs: 
!>  camera_inout: (camera) camera vertices
!>  st_scale:     (real8)(2) position rescaled by the plane size
subroutine pixel_local_to_scaled_plane_coord(camera_inout,id_time,i_pixel,st_pixel,st_scale)
  implicit none
  !> inputs-outputs
  class(camera),intent(inout)     :: camera_inout
  !> inputs:
  integer,intent(in)              :: id_time
  integer,dimension(2),intent(in) :: i_pixel
  real*8,dimension(2),intent(in)  :: st_pixel
  !> outputs:
  real*8,dimension(2),intent(out) :: st_scale
  !> compute the plane local coordinate rescaled by the plane size
  st_scale = camera_inout%plane_edge_length(:,id_time)*camera_inout%pixel_size*&
            (real(i_pixel-1,kind=8)+st_pixel)
end subroutine pixel_local_to_scaled_plane_coord

!> generate the pixel mesh grid scaled as the image plane
!> inputs:
!>  camera_inout: (camera) camera vertices
!> outputs:
!>  x_pos:        (real8)(n_pixels_x,n_times) x pixel positions for all times
!>  y_pos:        (real8)(n_pixels_y,n_times) y pixel positions for all times
!>  camera_inout: (camera) camera vertices
subroutine compute_scaled_pixel_positions_on_plane(camera_inout,x_pos,y_pos)
  implicit none
  !> inputs-outputs
  class(camera),intent(inout)  :: camera_inout
  !> outputs: 
  real*8,dimension(camera_inout%n_pixels_spectra(2),camera_inout%n_times),intent(out) :: x_pos
  real*8,dimension(camera_inout%n_pixels_spectra(3),camera_inout%n_times),intent(out) :: y_pos
  !> variables:
  integer :: ii,jj
  real*8,dimension(2) :: st
  do jj=1,camera_inout%n_times
    !> compute x mesh
    do ii=1,camera_inout%n_pixels_spectra(2)
      call camera_inout%pixel_local_to_scaled_plane_coord(jj,[ii,1],[5d-1,5d-1],st)
      x_pos(ii,jj) = st(1)
    enddo
    !> compute y mesh
    do ii=1,camera_inout%n_pixels_spectra(3)
      call camera_inout%pixel_local_to_scaled_plane_coord(jj,[1,ii],[5d-1,5d-1],st)
      y_pos(ii,jj) = st(2)
    enddo
  enddo
end subroutine compute_scaled_pixel_positions_on_plane

!> procedure used for the rendering of images
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
!>   images:
!>   ierr:                 (integer) error for the mpi procedures
subroutine compute_images(camera_inout,light_inout,spectra_inout,&
  filter_image_inout,filter_spectra_inout,filter_time_inout,rank,&
  images,ierr)
  use mod_light_vertices, only: light_vertices
  use mod_spectra,        only: spectrum_base
  use mod_filter,         only: filter
  implicit none
  !> inputs-outputs:
  class(camera),intent(inout)                                    :: camera_inout
  class(light_vertices),intent(inout)                            :: light_inout
  class(spectrum_base),intent(inout)                             :: spectra_inout
  class(filter),intent(inout)                                    :: filter_image_inout
  class(filter),dimension(spectra_inout%n_spectra),intent(inout) :: filter_spectra_inout
  class(filter),intent(inout)                                    :: filter_time_inout
  integer,intent(inout)                                          :: ierr
  !> inputs:
  integer,intent(in) :: rank
  !> outputs:
  real*8,dimension(camera_inout%n_pixels_spectra(1),&
  camera_inout%n_pixels_spectra(2),camera_inout%n_pixels_spectra(3),&
  camera_inout%n_times),intent(out) :: images
  !> variables:
  real*8,dimension(camera_inout%n_pixels_spectra(1),2,&
  camera_inout%n_pixels_spectra(2),camera_inout%n_pixels_spectra(3),&
  camera_inout%n_times) :: pixel_intensities

  !> reduce all lights contribution
  call camera_inout%reduce_light_image(light_inout,spectra_inout,&
  filter_image_inout,filter_spectra_inout,filter_time_inout,rank,&
  pixel_intensities,ierr)
  !> generate images
  if(rank.eq.0) then
    images = 0.d0
    where(pixel_intensities(:,2,:,:,:).ne.0.d0) images = &
    pixel_intensities(:,1,:,:,:)/pixel_intensities(:,2,:,:,:)
  endif
end subroutine compute_images
!>------------------------------------------------------------------
end module mod_camera
