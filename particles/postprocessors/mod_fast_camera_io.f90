!> The mod_fast_camera_io module contains procedures
!> for writing and reading fast camera data 
!> in from HDF5 file
module mod_fast_camera_io
#ifdef USE_HDF5
implicit none

private 
public :: write_pixel_intensity_hdf5

contains

!> Procedures first order integrator--------------------------------------
!> Write pixel intensity / filter structure in HDF5 file
!> inputs:
!>   filename:  (char) string contaning the hdf5 filename
!>   n_spectra: (integer) number of spectral lines
!>   n_values:  (integer) number of values stored per spectral
!>              line, per pixel and per time
!>   n_pixel_x: (integer) number of pixels in the x direction
!>   n_pixel_y: (integer) number of pixels in the y direction
!>   x_pixel_coord: (real8)(n_pixel_x,n_times) pixel x-coordinates
!>   y_pixel_coord: (real8)(n_pixel_y,n_times) pixel y-coordinates
!>   pixel_filter_array: (real8)(n_spectra,n_values,n_pixel_x,
!>                       n_pixel_y,n_times) camera image and filter data
!> outputs:
!>   ierr: (integer) error resulting from the HDF5 opening procedure
subroutine write_pixel_intensity_hdf5(filename,n_spectra,&
n_values,n_pixel_x,n_pixel_y,n_times,x_pixel_coord,y_pixel_coord,&
pixel_filter_array,ierr)
  use hdf5
  use hdf5_io_module, only: HDF5_open_or_create,HDF5_close
  use hdf5_io_module, only: HDF5_array2D_saving
  use hdf5_io_module, only: HDF5_array5D_saving
  implicit none
  !> Inputs:
  character(len=*),intent(in) :: filename
  integer,intent(in)          :: n_spectra,n_values,n_pixel_x
  integer,intent(in)          :: n_pixel_y,n_times
  real*8,dimension(n_pixel_x,n_times),intent(in) :: x_pixel_coord
  real*8,dimension(n_pixel_y,n_times),intent(in) :: y_pixel_coord
  real*8,dimension(n_spectra,n_values,n_pixel_x,n_pixel_y,n_times),intent(in) :: pixel_filter_array
  !> Outputs:
  integer,intent(out)         :: ierr
  !> Variables:
  integer(HID_T)              :: file_id
  !> open / store / close image
  call HDF5_open_or_create((trim(filename)//'.h5'),file_id,ierr,file_access=H5F_ACC_TRUNC_F)
  call HDF5_array2D_saving(file_id,x_pixel_coord,n_pixel_x,n_times,'x_pixel_coordinates')
  call HDF5_array2D_saving(file_id,y_pixel_coord,n_pixel_y,n_times,'y_pixel_coordinates')
  call HDF5_array5D_saving(file_id,pixel_filter_array,n_spectra,&
  n_values,n_pixel_x,n_pixel_y,n_times,'pixel_filter_intensities')
  call HDF5_close(file_id)
end subroutine write_pixel_intensity_hdf5

!>---------------------------------------------------------------
#endif
end module mod_fast_camera_io

