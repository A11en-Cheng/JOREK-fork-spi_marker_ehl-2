!> mod_spectra_deterministic contains variables and procedures
!> for generating deterministic spectra and integatre values 
!> on spectra using deterministic methods
module mod_spectra_deterministic
use mod_spectra, only: spectrum_base
implicit none

private
public :: spectrum_integrator_2nd

!> Variables and type definitions -------------------------------
!> type generating and integrating data on spectra using
!> the central point method
type,extends(spectrum_base) :: spectrum_integrator_2nd
  real*8,dimension(:),allocatable :: min_wlen  !< minimum wavelength
  real*8,dimension(:),allocatable :: wbin_size !< wlenght interval size
  contains
  procedure,pass(spectrum) :: allocate_spectrum     => allocate_spectrum_integrator_2nd
  procedure,pass(spectrum) :: set_spectrum_interval => set_spectrum_int_2nd_properties
  procedure,pass(spectrum) :: generate_spectrum     => generate_midpoint_spectra 
  procedure,pass(spectrum) :: integrate_data        => integrate_spectrum_rectangle
  procedure,pass(spectrum) :: deallocate_spectrum   => deallocate_spectrum_integrator_2nd
end type spectrum_integrator_2nd

!> Interfaces----------------------------------------------------
interface spectrum_integrator_2nd
  module procedure construct_spectrum_integarator_2nd
end interface

contains

!> Constructors--------------------------------------------------
!> constructor of the spectrum_integrator_2nd class
!> inputs:
!>   n_points:   (integer) number of spectral points
!>   n_spectrum: (integer) number of spectra
!>   min_wlen:   (real8)(n_spectra) spectrum intervals lower bound
!>   max_wlen:   (real8)(n_spectra) spectrum intervals upper bound
!> outputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
function construct_spectrum_integarator_2nd(n_points,n_spectra,&
min_wlen,max_wlen) result(spectrum)
  implicit none
  !> inputs:
  integer,intent(in) :: n_points,n_spectra
  real*8,dimension(n_spectra),intent(in),optional :: min_wlen,max_wlen
  !> outputs
  type(spectrum_integrator_2nd),target :: spectrum
  !> variables
  real*8,dimension(2*n_spectra) :: real8_param
  if(present(min_wlen).and.present(max_wlen)) then
    real8_param(1:n_spectra) = min_wlen
    real8_param(n_spectra+1:2*n_spectra) = max_wlen
    call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  else
    call spectrum%allocate_spectrum(n_points,n_spectra)
  endif
end function construct_spectrum_integarator_2nd

!> Procedures first order integrator-----------------------------
!> allocate spectrum integrator 2nd order datatype
!> inputs:
!>   spectrum:    (spectrum_integrator_2nd) generate and integrated
!>                variables using the central point rule
!>   n_points:    (integer) number of spectral points
!>   n_spectra:   (integer) number of spectra
!>   real8_param: (real8)(2*n_spectra)(optional) double parameters:
!>                1:n_spectra             -> minimum wavelength
!>                n_spectra+1:2*n_spectra -> maximum wavelength
!>   int_param:   (integer)(0)(optional) integer parameters
!> outputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
subroutine allocate_spectrum_integrator_2nd(spectrum,n_points,&
n_spectra,real8_param,int_param)
  implicit none
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  real*8,dimension(2*n_spectra),intent(in),optional :: real8_param
  integer,dimension(0),intent(in),optional          :: int_param
  !> inputs-outputs
  class(spectrum_integrator_2nd),intent(inout) :: spectrum

  !> allocate arrays
  call spectrum%allocate_spectrum_base(n_points,n_spectra)
  if(present(real8_param)) call spectrum%set_spectrum_interval(&
  n_spectra,real8_param(1:n_spectra),&
  real8_param(n_spectra+1:2*n_spectra))
end subroutine allocate_spectrum_integrator_2nd

!> set the spectrum integrator 2nd order properties
!> inputs:
!>   spectrum:   (spectrum_integrator_2nd) generate and integrated
!>               variables using the central point rule
!>   n_spectrum: (integer) number of spectra
!>   min_wlen:   (real8)(n_spectra) spectrum intervals lower bound
!>   max_wlen:   (real8)(n_spectra) spectrum intervals upper bound
!> outputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
subroutine set_spectrum_int_2nd_properties(spectrum,n_spectra,&
min_wlen,max_wlen)
  implicit none
  !> inputs
  integer,intent(in) :: n_spectra
  real*8,dimension(n_spectra) :: min_wlen,max_wlen
  !> inputs-outputs
  class(spectrum_integrator_2nd),intent(inout) :: spectrum
  !> set values
  if(spectrum%n_spectra.ne.n_spectra) then
    if(allocated(spectrum%min_wlen))  deallocate(spectrum%min_wlen)
    if(allocated(spectrum%wbin_size)) deallocate(spectrum%wbin_size)
    spectrum%n_spectra = n_spectra
    write(*,*) 'WARNING: n_spectra is changed -> regenerate spectrum points!'
  endif
  if(.not.allocated(spectrum%min_wlen))  allocate(spectrum%min_wlen(n_spectra))
  if(.not.allocated(spectrum%wbin_size)) allocate(spectrum%wbin_size(n_spectra))
  spectrum%min_wlen = min_wlen
  spectrum%wbin_size = (max_wlen-min_wlen)/spectrum%n_points
end subroutine set_spectrum_int_2nd_properties

!> generate spectral points using the mid-point rule 
!> (mesh are the mid points of the intervals)
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
!> outputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
subroutine generate_midpoint_spectra(spectrum)
  implicit none
  !> inputs-outputs
  class(spectrum_integrator_2nd),intent(inout) :: spectrum
  !> variables
  integer :: ii,jj
  !$omp parallel do default(private) shared(spectrum) collapse(2)
  do jj=1,spectrum%n_spectra
    do ii=1,spectrum%n_points
      spectrum%points(ii,jj) = spectrum%min_wlen(jj) + &
      (real(ii,kind=8)-5.d-1)*spectrum%wbin_size(jj)
    enddo
  enddo
  !$omp end parallel do
end subroutine generate_midpoint_spectra

!> integrate data on the spectrum interval using the rectangle rule
!> inputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
!>   dat:      (real8)(n_points,n_spectra) integrand values evaluated
!>             at the mid-point of each interval
!> outputs:
!>   spectrum:  (spectrum_integrator_2nd) generate and integrated
!>              variables using the central point rule
!>   integrals: (real8)(n_spectra) integrals of midpoint data for each spectrum
subroutine integrate_spectrum_rectangle(spectrum,dat,integrals)
  !$ use omp_lib
  implicit none
  !> input-outputs
  class(spectrum_integrator_2nd),intent(inout) :: spectrum
  !> inputs
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(in) :: dat
  !> outputs
  real*8,dimension(spectrum%n_spectra),intent(out) :: integrals
  !> variables
  integer :: ii,thread_id,n_threads,n_points_per_thread,residual_id
  logical :: in_parallel
  !> integration
  in_parallel = .true.
  thread_id = 0
  n_threads = 1
  integrals = 0.d0
  !$ in_parallel = omp_in_parallel()
  if(in_parallel) then
    !> reduction clause in taskloop seems not to be supported by all compilers
    !> so a single loop is implemented with no parallelization on the summatory
#ifdef USE_TASKLOOP
    !$omp taskloop default(shared) private(ii) !nogroup
#endif
    do ii=1,spectrum%n_spectra
      integrals(ii) = sum(dat(:,ii))
    enddo
#ifdef USE_TASKLOOP
    !$omp end taskloop
#endif
  else
    !$omp parallel default(private) shared(spectrum,dat,residual_id) &
    !$omp firstprivate(n_threads) reduction(+:integrals)
    !$ n_threads = omp_get_num_threads()
    !$ thread_id = omp_get_thread_num()
    n_points_per_thread = spectrum%n_points/n_threads
    residual_id = spectrum%n_points - n_points_per_thread*n_threads
    do ii=1,spectrum%n_spectra
      integrals(ii) = integrals(ii) + sum(dat(&
      thread_id*n_points_per_thread+1:(thread_id+1)*n_points_per_thread,ii))
    enddo
    !$omp end parallel
    if(residual_id.gt.0) then
      do ii=1,spectrum%n_spectra
        integrals(ii) = integrals(ii) + sum(dat(&
        spectrum%n_points-residual_id+1:spectrum%n_points,ii))
      enddo
    endif
  endif
  integrals = integrals*spectrum%wbin_size
end subroutine integrate_spectrum_rectangle

!> deallocate spectrum intgrator 2nd order datatype
!> inputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
!> outputs:
!>   spectrum: (spectrum_integrator_2nd) generate and integrated
!>             variables using the central point rule
subroutine deallocate_spectrum_integrator_2nd(spectrum)
  implicit none
  !> inputs-outputs
  class(spectrum_integrator_2nd),intent(inout) :: spectrum

  !> deallocate everything and set variables to default
  call spectrum%deallocate_spectrum_base
  if(allocated(spectrum%min_wlen))  deallocate(spectrum%min_wlen)
  if(allocated(spectrum%wbin_size)) deallocate(spectrum%wbin_size)

end subroutine deallocate_spectrum_integrator_2nd

!>---------------------------------------------------------------

end module mod_spectra_deterministic
