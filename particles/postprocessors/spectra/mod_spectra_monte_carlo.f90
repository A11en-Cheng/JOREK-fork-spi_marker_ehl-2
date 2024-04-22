!> mod_spectra_monte_carlo contains all variables and procedures
!> for generating spectra and integrating values using Monte-Carlo
!> integrators
module mod_spectra_monte_carlo
use mod_spectra, only: spectrum_base
implicit none

private
public :: spectrum_rng_uniform

!> generate a set of spectral points for integration
!> using a uniform distribution
type,extends(spectrum_base) :: spectrum_rng_uniform
  real*8,dimension(:),allocatable :: min_wlen !< lower wavelenght of the interval
  real*8,dimension(:),allocatable :: i_pdf !< 1 over probability density function
  contains
  procedure,pass(spectrum) :: allocate_spectrum      => allocate_spectrum_rng_uniform
  procedure,pass(spectrum) :: set_spectrum_interval  => set_uniform_spectrum_interval
  procedure,pass(spectrum) :: generate_spectrum      => generate_uniform_rng_spectrum
  procedure,pass(spectrum) :: integrate_data         => integrate_rng_uniform
  procedure,pass(spectrum) :: deallocate_spectrum    => deallocate_spectrum_rng_uniform
end type spectrum_rng_uniform

!> Interfaces ---------------------------------------
interface spectrum_rng_uniform
  module procedure construct_spectrum_rng_uniform
end interface

contains

!> Constructors -------------------------------------
!> construct a uniform random spectrum
!> inputs:
!>   n_points:  (integer) number of discrete random variables
!>   n_spectra: (integer) number of spectra
!>   min_wlen:  (real8) minimum wavelength (interval lowerbound)
!>   max_wlen:  (real8) maximum wavelength (interval upperbound)
!> outputs:
!>   spectrum: (spectrum_rng_uniform) uniform random spectrum
function construct_spectrum_rng_uniform(n_points,n_spectra,&
min_wlen,max_wlen) result(spectrum)
  implicit none
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  real*8,dimension(n_spectra),intent(in),optional :: min_wlen,max_wlen
  !> outputs
  type(spectrum_rng_uniform),target :: spectrum
  !> variables
  real*8,dimension(2*n_spectra) :: real8_param 
  if(present(min_wlen).and.present(max_wlen)) then
    real8_param(1:n_spectra) = min_wlen
    real8_param(n_spectra+1:2*n_spectra) = max_wlen
    call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  else
    call spectrum%allocate_spectrum(n_points,n_spectra)
  endif
end function construct_spectrum_rng_uniform

!> Procedures spectrum rng uniform ------------------
!> allocate the spectrum_rng_uniform datatype
!> inputs:
!>   spectrum:    (spectrum_rng_uniform) generates and integrates
!>                variables along a uniform spectral distribution
!>   n_points:    (integer) number of spectral points
!>   n_spectra:   (integer) number of spectral intervals
!>   real8_param: (real8)(2*n_spectra)(optional) double parameters:
!>                1:n_spectral            -> minimum wavelengths
!>                n_spectra+1:2*n_spectra ->maximum wavelengths
!>   int_param:   (integer)(0)(optional) integer parameters
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates
!>             variables along a uniform spectral distribution
subroutine allocate_spectrum_rng_uniform(spectrum,n_points,&
n_spectra,real8_param,int_param)
  implicit none
  !> inputs
  integer,intent(in) :: n_points,n_spectra
  real*8,dimension(2*n_spectra),intent(in),optional  :: real8_param
  integer,dimension(0),intent(in),optional :: int_param
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> allocated all variables
  call spectrum%allocate_spectrum_base(n_points,n_spectra)
  if(present(real8_param)) &
  call spectrum%set_spectrum_interval(n_spectra,real8_param(1:n_spectra),&
  real8_param(n_spectra+1:2*n_spectra))
end subroutine allocate_spectrum_rng_uniform

!> set the spectrum properties
!> inputs:
!>   spectrum:  (spectrum_rng_uniform) generates and integrates
!>              variables along a uniform spectral distribution
!>   n_spectra: (integer) number of spectral intervals
!>   min_wlen:  (real8)(n_spectra) minimum wavelength
!>   max_wlen:  (real8)(n_spectra) maximum wavelength
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates
!>             variables along a uniform spectral distribution
subroutine set_uniform_spectrum_interval(spectrum,n_spectra,min_wlen,max_wlen)
  implicit none
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> inputs
  integer,intent(in) :: n_spectra
  real*8,dimension(n_spectra),intent(in) :: min_wlen,max_wlen
  !> set values
  if(.not.allocated(spectrum%min_wlen)) allocate(spectrum%min_wlen(n_spectra))
  if(.not.allocated(spectrum%i_pdf))    allocate(spectrum%i_pdf(n_spectra))
  if(spectrum%n_spectra.ne.n_spectra) then
    deallocate(spectrum%min_wlen); deallocate(spectrum%i_pdf);
    allocate(spectrum%min_wlen(n_spectra))
    allocate(spectrum%i_pdf(n_spectra))
    spectrum%n_spectra = n_spectra
    write(*,*) 'WARNING: n_spectra is changed -> regenerate spectrum points!'
  endif
  spectrum%min_wlen = min_wlen
  spectrum%i_pdf = max_wlen - min_wlen
end subroutine set_uniform_spectrum_interval

!> generate uniform random spectrum
!> inputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates 
!>             variables along a uniform spectral distribution
!>   rngs:     (type_rng)(n_omp_threads,n_spectra) random number generators
!> outputs:
!>   spectrum: (spectrum_rng_uniform) generates and integrates 
!>             variables along a uniform spectral distribution
subroutine generate_uniform_rng_spectrum(spectrum,rngs)
  use mod_rng
  !$ use omp_lib
  implicit none
  !> inputs-outpus
  class(spectrum_rng_uniform),intent(inout)   :: spectrum
  class(type_rng),dimension(:),allocatable,intent(inout) :: rngs
  !> variables
  integer :: n_points_per_tile=25
  integer :: ii,jj,n_tiles,n_residual,thread_id
  real*8,dimension(:),allocatable :: rands
  !> generate spectrum from uniform random number distribution
  thread_id = 1
  n_tiles = int(spectrum%n_points/n_points_per_tile)
  allocate(rands(n_points_per_tile)); rands=0.d0;
  !$omp parallel default(private) shared(spectrum,rngs) &
  !$omp firstprivate(n_tiles,n_points_per_tile)
  !$ thread_id = omp_get_thread_num()+1
  !$omp do collapse(2)
  do ii=1,spectrum%n_spectra
    do jj=1,n_tiles
      call rngs(thread_id)%next(rands)
      spectrum%points((jj-1)*n_points_per_tile+1:&
      jj*n_points_per_tile,ii) = spectrum%min_wlen(ii) + &
      spectrum%i_pdf(ii)*rands
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  n_residual = spectrum%n_points-n_tiles*n_points_per_tile
  if(n_residual.gt.0.d0) then
    do ii=1,spectrum%n_spectra
      call rngs(ii)%next(rands(1:n_residual))
      spectrum%points(n_tiles*n_points_per_tile+1:spectrum%n_points,ii) = &
      spectrum%min_wlen(ii) + spectrum%i_pdf(ii)*rands(1:n_residual)
    enddo
  endif
  deallocate(rands)
  end subroutine generate_uniform_rng_spectrum

!> integrate data computed using uniformly distributed spectra
!> inputs:
!>   spectrum:     (spectrum_rng_uniform) generates and integrates 
!>                 variables along a uniform spectral distribution
!>   dat:          (real8)(n_points,n_spectra) data obtained from
!>                 a uniform spectral distribution
!> outputs:
!>   spectrum:  (spectrum_rng_uniform) generates and integrates 
!>              variables along a uniform spectral distribution
!>   integrals: (n_spectra) integrals of uniform_data for each spectrum
subroutine integrate_rng_uniform(spectrum,dat,integrals)
  !$ use omp_lib
  implicit none
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  !> inputs
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(in) :: dat
  !> outputs
  real*8,dimension(spectrum%n_spectra),intent(out) :: integrals
  !> variables
  integer :: ii,thread_id,n_threads,n_points_per_thread,residual_id
  logical :: in_parallel
  !> integrate
  in_parallel = .true.; thread_id = 0; 
  n_threads = 1; integrals = 0d0;
  !$ in_parallel = omp_in_parallel()
  if(in_parallel) then
    !> Not all compilers support reduction clauses with taskloop
    !> hence only the outer loop is parallelised
#ifdef USE_TASKLOOP
    !$omp taskloop default(shared) private(ii)
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
    !$ thread_id = omp_get_thread_num()
    !$ n_threads = omp_get_num_threads()
    n_points_per_thread = spectrum%n_points/n_threads
    residual_id = n_points_per_thread*n_threads
    do ii=1,spectrum%n_spectra
      integrals(ii) = integrals(ii) + sum(dat(&
      thread_id*n_points_per_thread+1:(thread_id+1)*n_points_per_thread,ii))
    enddo 
    !$omp end parallel
    if(residual_id.lt.spectrum%n_points) then
      do ii=1,spectrum%n_spectra
        integrals(ii) = integrals(ii) + sum(dat(residual_id+1:spectrum%n_points,ii))
      enddo
    endif  
  endif
  integrals = integrals*spectrum%i_pdf/spectrum%n_points
end subroutine integrate_rng_uniform

!> deallocate spectrum rng uniform and cleanup
subroutine deallocate_spectrum_rng_uniform(spectrum)
  implicit none
  !> inputs-outputs
  class(spectrum_rng_uniform),intent(inout) :: spectrum
  call spectrum%deallocate_spectrum_base
  if(allocated(spectrum%min_wlen)) deallocate(spectrum%min_wlen)
  if(allocated(spectrum%i_pdf))    deallocate(spectrum%i_pdf)
end subroutine deallocate_spectrum_rng_uniform

!>---------------------------------------------------

end module mod_spectra_monte_carlo 
