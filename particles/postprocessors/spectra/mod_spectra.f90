!> mod_spectra contains base variables and procedures for
!> the spectra modules
module mod_spectra
implicit none

private
public :: spectrum_base

!> Variable and type definitions ---------------------
!> spectrum: abstract class containing the basic types
!> and procedures for generating radiation spectra
type,abstract :: spectrum_base
  integer :: n_points  !< number of wave lengths or colors
  integer :: n_spectra !< number of wave length or color intervals
  !> wave lengths or colors
  real*8,dimension(:,:),allocatable :: points
  contains
  procedure(int_integral),pass(spectrum),deferred :: integrate_data
  procedure,pass(spectrum)                        :: read_spectrum_inputs
  procedure,pass(spectrum)                        :: return_n_spectrum_inputs
  procedure,pass(spectrum)                        :: allocate_spectrum_base
  procedure,pass(spectrum)                        :: deallocate_spectrum_base
end type spectrum_base

!> Interfaces ---------------------------------------
interface
  !> interface of the procedure used for integrating spectral data
  !> inputs:
  !>   spectrum: (spectrum_base) spectrum object
  !>   dat:      (n_points,n_intervals) spectral data to be integrated
  !> outputs:
  !>   spectrum:  (spectrum_base) spectrum object
  !>   integrals: (n_spectra) integral value
  subroutine int_integral(spectrum,dat,integrals)
    IMPORT :: spectrum_base
    implicit none
    !> inputs-outputs
    class(spectrum_base),intent(inout) :: spectrum
    !> inputs:
    real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(in) :: dat
    !> outputs:
    real*8,dimension(spectrum%n_spectra),intent(out) :: integrals
  end subroutine int_integral
end interface

contains
!> Procedures spectrum base -------------------------
!> read the spectrum base inputs
!> inputs:
!>   spectrum: (spectrum_base) generate and integrated
!>   my_id:    (integer) mpi rank id
!>   r_unit:   (integer) read unit id
!>   n_inputs: (integer)(2) size of the integer and real input arrays
!> outputs:
!>   spectrum:   (spectrum_base) generate and integrated
!>   n_inputs:   (integer)(2) size of the integer and real input arrays
!>   int_param:  (integer)(:) integer input parameters
!>   real_param: (real8)(:) real8 input parameters
subroutine read_spectrum_inputs(spectrum,my_id,r_unit,n_inputs,&
int_param,real_param)
  implicit none
  !> inputs-outputs:
  class(spectrum_base),intent(inout) :: spectrum
  integer,dimension(2),intent(inout) :: n_inputs
  !> inputs:
  integer,intent(in) :: my_id,r_unit
  !> outputs:
  integer,dimension(:),allocatable,intent(out) :: int_param
  real*8,dimension(:),allocatable,intent(out)  :: real_param
  !> variables:
  integer :: n_points,n_spectra
  real*8,dimension(:),allocatable :: min_max_wavelength
  !> initialisation and definitions
  namelist /spectrum_in/ n_points,n_spectra
  namelist /wavelength_in/ min_max_wavelength
  !> read spectrum inputs
  if(my_id.eq.0) then
    n_inputs = spectrum%return_n_spectrum_inputs()
    if(allocated(int_param))  deallocate(int_param)
    allocate(int_param(n_inputs(1)))
    read(r_unit,spectrum_in)
    n_inputs = n_inputs + (/0,2*n_spectra+1/)
    int_param = (/n_points,n_spectra/)
    allocate(min_max_wavelength(2*n_spectra))
    read(r_unit,wavelength_in)
    if(allocated(real_param)) deallocate(real_param)
    allocate(real_param(2*n_spectra))
    real_param = min_max_wavelength
  endif
  !> cleanup
  if(allocated(min_max_wavelength)) deallocate(min_max_wavelength)
end subroutine read_spectrum_inputs 

!> return the number of integer and real inputs
!> inputs:
!>   spectrum: (spectrum_integrator_2nd) spectrum object
!> outputs:
!>   n_inputs: (integer)(2) size of the input parameters
!>             if -1, size to be determined at reading
function return_n_spectrum_inputs(spectrum) result(n_inputs)
  implicit none
  !> inputs:
  class(spectrum_base),intent(in) :: spectrum
  !> outputs:
  integer,dimension(2) :: n_inputs
  n_inputs = (/2,-1/)
end function return_n_spectrum_inputs

!> allocate spectrum base and set counters
subroutine allocate_spectrum_base(spectrum,n_points,n_spectra)
  implicit none
  !> input-outputs
  class(spectrum_base),intent(inout) :: spectrum
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  if(.not.allocated(spectrum%points).and.((n_points.gt.0).and.&
  (n_spectra.gt.0))) &
  allocate(spectrum%points(n_points,n_spectra))
  spectrum%n_points  = n_points
  spectrum%n_spectra = n_spectra
end subroutine allocate_spectrum_base

!> deallocate spectrum base and set counters to 0
subroutine deallocate_spectrum_base(spectrum)
  implicit none
  !> inputs-outputs
  class(spectrum_base),intent(inout) :: spectrum
  !> cleanup
  if(allocated(spectrum%points)) deallocate(spectrum%points)
  spectrum%n_points  = -1
  spectrum%n_spectra = -1
end subroutine deallocate_spectrum_base

!>---------------------------------------------------

end module mod_spectra
