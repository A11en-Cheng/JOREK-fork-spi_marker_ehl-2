!> mod_spectra_monte_carlo_test contains all setup, teadown and test function
!> for verifying the correctness of the procedure defining spectra based on 
!> Monte-Carlo method
module mod_spectra_monte_carlo_test
use fruit
use mod_rng, only: type_rng
implicit none

private
public :: run_fruit_spectra_monte_carlo

!> Variables -----------------------------------------------------------------
logical,parameter :: use_xor_time_pid=.true.
character(len=33),parameter :: input_file='spectrum_deterministic_2nd_inputs'
integer,parameter :: read_unit=43     !< file read unit
integer,parameter :: n_trials=50      !< number of trials for computing std deviation
integer,parameter :: n_convergence=5  !< number of points for convergence
integer,parameter :: n_points=5512437 !< number of points
integer,parameter :: n_spectra=2      !< number of spectra
integer,parameter :: n_bins=500       !< number of bins for histograms
real*8,parameter  :: tol_real8=1.d-16 !< tolerance for assert
real*8,parameter  :: tol_updf=5.d-4   !< tolerance on the uniform probability
real*8,parameter  :: expected_std_conv_coeff=-5.d-1 !< expected std convergence rate
real*8,parameter  :: tol_std_conv_coeff=5.d-2       !< tolerance on the convergence rate
real*8,parameter  :: tol_min_int_error=7.5d-4       !< tolerance on the minimum int error
!> n_points for convergence study
integer,dimension(n_convergence),parameter :: n_points_conv=(/997,100725,100000,1003757,10023947/)
real*8,dimension(2),parameter :: min_wlen=(/3.d-6,3.0d-7/) !< minimum wavelength
real*8,dimension(2),parameter :: max_wlen=(/3.5d-6,4.d-7/) !< maximum wavelength
real*8,dimension(2),parameter :: min_angle=(/6.d-1,2.3d0/) !< minimum angle for integration
real*8,dimension(2),parameter :: max_angle=(/3.6d0,3.4d0/) !< maximum angle for integration
class(type_rng),dimension(:),allocatable :: rngs !< random number generators
integer                                  :: n_threads,thread_id !< N# and id omp threads
real*8,dimension(n_spectra)              :: i_pdf !< 1/pdf=min_wlen,max_wlen

contains
!> Test basket ---------------------------------------------------------------
!> test basket for executing the simulation set-up, tests and tear-down
subroutine run_fruit_spectra_monte_carlo()
implicit none
  write(*,'(/A)') "  ... setting-up: spectra Monte-Carlo tests"
  call setup
  call write_spectrum_inputs
  write(*,'(/A)') "  ... running: spectra Monte-Carlo tests"
  call run_test_case(test_spectrum_input_reader,'test_spectrum_input_reader')
  call run_test_case(test_spectrum_rng_uniform_construction_noinit,&
  'test_spectrum_rng_uniform_construction_noinit')
  call run_test_case(test_set_uniform_spectrum_interval,&
  'test_set_uniform_spectrum_interval')
  call run_test_case(test_spectrum_rng_uniform_construction_init,&
  'test_spectrum_rng_uniform_construction_init')
  call run_test_case(test_spectrum_generation_rng_uniform,&
  'test_spectrum_generation_rng_uniform')
  call run_test_case(test_spectrum_integration_rng_uniform,&
  'test_spectrum_integration_rng_uniform')
  write(*,'(/A)') "  ... tearing-down: spectra Monte-Carlo tests"
  call teardown
end subroutine run_fruit_spectra_monte_carlo

!> Set-up and tear-down ------------------------------------------------------
!> Set-up test variables
subroutine setup()
  use mod_common_test_tools, only: omp_initialize_rngs
  use mod_pcg32_rng,  only: pcg32_rng
  !$ use omp_lib
  implicit none
  !> variables
  integer :: n_threads
  !> set the inverse of the pdf
  i_pdf = max_wlen-min_wlen
  !> initialise the rngs using the pcg32
  n_threads = 1
  !$ n_threads = omp_get_max_threads()
  allocate(pcg32_rng::rngs(n_threads))
  call omp_initialize_rngs(n_points,n_threads,rngs,&
  use_xor_time_pid_in=use_xor_time_pid)
end subroutine setup

!> write an input for testing the reading routines
subroutine write_spectrum_inputs()
  implicit none
  !> variables
  integer :: ifail
  open(read_unit,file=input_file,status='unknown',action='write',iostat=ifail)
  write(read_unit,'(/A)') '&spectrum_in'
  write(read_unit,'(/A,I10)') 'n_points = ',n_points
  write(read_unit,'(/A,I10)') 'n_spectra = ',n_spectra
  write(read_unit,'(/A)') '/'
  write(read_unit,'(/A)') '&wavelength_in'
  write(read_unit,'(/A,F20.16)') 'min_max_wavelength(1) = ',min_wlen(1)
  write(read_unit,'(/A,F20.16)') 'min_max_wavelength(2) = ',min_wlen(2)
  write(read_unit,'(/A,F20.16)') 'min_max_wavelength(3) = ',max_wlen(1)
  write(read_unit,'(/A,F20.16)') 'min_max_wavelength(4) = ',max_wlen(2)
  write(read_unit,'(/A)') '/'
  close(read_unit)
end subroutine write_spectrum_inputs

!> clean up all test variables
subroutine teardown()
  implicit none
  i_pdf = 0.d0
  deallocate(rngs)
  call system("rm "//input_file)
end subroutine teardown

!> Tests ---------------------------------------------------------------------
!> test input reading parameters
subroutine test_spectrum_input_reader()
  use mod_spectra_monte_carlo, only: spectrum_rng_uniform
  implicit none 
  !> variables
  type(spectrum_rng_uniform) :: spectrum
  integer :: rank,ifail
  integer,dimension(2)             :: n_inputs
  integer,dimension(:),allocatable :: int_param
  real*8,dimension(:),allocatable  :: real_param
  !> initialisation
  rank = 0
  !> read parameters from file
  open(read_unit,file=input_file,status='old',action='read',iostat=ifail)
  call spectrum%read_spectrum_inputs(rank,read_unit,n_inputs,&
  int_param,real_param)
  close(read_unit)
  !> checks
  call assert_equals((/2,2*n_spectra/),n_inputs,2,&
  "Error spectrum monte carlo uniform input reader: N# inputs mismatch!")
  call assert_equals((/n_points,n_spectra/),int_param,2,&
  "Error spectrum monte carlo uniform input reader: N# integer parameters mismatch!")
  call assert_equals((/min_wlen(1),min_wlen(2),max_wlen(1),max_wlen(2)/),real_param,4,&
  "Error spectrum monte carlo uniform input reader: N# real parameters mismatch!")
  !> cleanup
  if(allocated(int_param)) deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine test_spectrum_input_reader

!> test allocation, deallocation and construction of spectrum_base class
subroutine test_spectrum_rng_uniform_construction_noinit()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_spectra_monte_carlo, only: spectrum_rng_uniform
  implicit none
  !> variables 
  type(spectrum_rng_uniform) :: spectrum

  !> try allocation
  call spectrum%allocate_spectrum(n_points,n_spectra)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum base allocation: n_points mismatch!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum base allocation: n_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_points,n_spectra,&
  spectrum%points,"Error spectrum base allocation: points")

  !> try deallocate
  call spectrum%deallocate_spectrum
  call assert_equals(-1,spectrum%n_points,&
  "Error spectrum base deallocation: failed cleaning n_points!")
  call assert_equals(-1,spectrum%n_spectra,&
  "Error spectrum base deallocation: failed cleaning n_spectra!")
  call assert_false(allocated(spectrum%points),&
  "Error spectrum base deallocation: points still allocated!")

  !> try construction
  spectrum = spectrum_rng_uniform(n_points,n_spectra) 
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum base construction: n_points mismatch!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum base construction: n_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_points,n_spectra,&
  spectrum%points,"Error spectrum base construction: points")
  call spectrum%deallocate_spectrum
  
end subroutine test_spectrum_rng_uniform_construction_noinit

!> test set and change of the spectrum interval
subroutine test_set_uniform_spectrum_interval()
  use mod_spectra_monte_carlo, only: spectrum_rng_uniform
  implicit none
  !> variables
  type(spectrum_rng_uniform)            :: spectrum
  real*8,dimension(n_spectra),parameter :: min_wlen_2=(/1.d-5,5.d-9/)
  real*8,dimension(n_spectra),parameter :: max_wlen_2=(7.4d-5,1.25d-8)
  real*8,dimension(2*n_spectra)         :: min_wlen_3,max_wlen_3

  !> test setting of min_wlen,max_wlen
  spectrum = spectrum_rng_uniform(n_points,n_spectra)
  call spectrum%set_spectrum_interval(n_spectra,min_wlen,max_wlen)
  call assert_equals(n_spectra,spectrum%n_spectra,"Error set uniform interval: n_spectra mismatch!")
  call assert_equals(min_wlen,spectrum%min_wlen,n_spectra,tol_real8,&
  "Error set uniform interval: min_wlen mismatch")
  call assert_equals(i_pdf,spectrum%i_pdf,n_spectra,tol_real8,&
  "Error set uniform interval: i_pdf mismatch")

  !> test change with equal number of spectral intervals
  call spectrum%set_spectrum_interval(n_spectra,min_wlen_2,max_wlen_2)
  call assert_equals(min_wlen_2,spectrum%min_wlen,n_spectra,tol_real8,&
  "Error set uniform interval: min_wlen_2 mismatch")
  call assert_equals(max_wlen_2-min_wlen_2,spectrum%i_pdf,n_spectra,tol_real8,&
  "Error set uniform interval: i_pdf_2 mismatch")

  !> test change with different number of intervals
  min_wlen_3(1:n_spectra) = min_wlen; min_wlen_3(n_spectra+1:2*n_spectra) = min_wlen_2;
  max_wlen_3(1:n_spectra) = max_wlen; max_wlen_3(n_spectra+1:2*n_spectra) = max_wlen_2;
  call spectrum%set_spectrum_interval(2*n_spectra,min_wlen_3,max_wlen_3)
  call assert_equals(2*n_spectra,spectrum%n_spectra,&
  "Error set uniform interval: n_spectra_2 mismatch!")
  call assert_equals(min_wlen_3,spectrum%min_wlen,2*n_spectra,tol_real8,&
  "Error set uniform interval: min_wlen_3 mismatch")
  call assert_equals(max_wlen_3-min_wlen_3,spectrum%i_pdf,2*n_spectra,tol_real8,&
  "Error set uniform interval: i_pdf_3 mismatch")  

  !> deallocate variables
  call spectrum%deallocate_spectrum()
end subroutine test_set_uniform_spectrum_interval

!> test allocation and construction with initialisation
subroutine test_spectrum_rng_uniform_construction_init()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_spectra_monte_carlo, only: spectrum_rng_uniform
  implicit none
  !> variables
  type(spectrum_rng_uniform)    :: spectrum
  real*8,dimension(2*n_spectra) :: real8_param

  !> test allocation with initialisation
  real8_param(1:n_spectra) = min_wlen
  real8_param(n_spectra+1:2*n_spectra) = max_wlen
  spectrum = spectrum_rng_uniform(n_points,n_spectra)
  call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum rng uniform allocation: n_points mismatch!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum rng uniform allocation: n_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_points,n_spectra,&
  spectrum%points,"Error spectrum rng uniform allocation: points")
  call assert_equals_allocatable_arrays(n_spectra,min_wlen,spectrum%min_wlen,&
  tol_real8,"Error spectrum rng uniform allocation: min_wlen")
  call assert_equals_allocatable_arrays(n_spectra,i_pdf,spectrum%i_pdf,&
  tol_real8,"Error spectrum rng uniform allocation: i_pdf")
  call spectrum%deallocate_spectrum

  !> test constructor with initialisation
  spectrum = spectrum_rng_uniform(n_points,n_spectra,min_wlen,max_wlen)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum rng uniform constructor: n_points mismatch!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum rng uniform constructor: n_spectra mismatch!")
  call assert_equals_allocatable_arrays(n_points,n_spectra,&
  spectrum%points,"Error spectrum rng uniform constructor: points")
  call assert_equals_allocatable_arrays(n_spectra,min_wlen,spectrum%min_wlen,&
  tol_real8,"Error spectrum rng uniform constructor: min_wlen")
  call assert_equals_allocatable_arrays(n_spectra,i_pdf,spectrum%i_pdf,&
  tol_real8,"Error spectrum rng uniform constructor: i_pdf")
  call spectrum%deallocate_spectrum

end subroutine test_spectrum_rng_uniform_construction_init

!> test the random generation of uniform spectra within interval
subroutine test_spectrum_generation_rng_uniform()
  use mod_spectra_monte_carlo, only: spectrum_rng_uniform
  implicit none
  !> variables
  type(spectrum_rng_uniform)          :: spectrum
  integer                             :: ii,jj
  integer,dimension(n_bins,n_spectra) :: histogram
  real*8                              :: min_w,delta_bin_w
  real*8,dimension(n_spectra)         :: max_variation


  !> initialise and construct the spectrum class
  spectrum = spectrum_rng_uniform(n_points,n_spectra,min_wlen,max_wlen)
  !> generate new random spectra, check intervals bounds
  call spectrum%generate_spectrum(rngs)
  do ii=1,n_spectra
    call assert_true(minval(spectrum%points(:,ii)).ge.min_wlen(ii),&
    "Error generate uniform random spectrum: minmum point value out-of-bound")
    call assert_true(maxval(spectrum%points(:,ii)).le.max_wlen(ii),&
    "Error generate uniform random spectrum: maximum point value out-of-bound")
  enddo
  !> generating histograms for checking pdf, if the probability is uniform
  !> then for a very large number of realization the histograms of all bins
  !> are almost equal. We admit a difference between the number of points
  !> in each bin < tol_var for success
  do jj=1,n_spectra
    delta_bin_w = spectrum%i_pdf(jj)/n_bins
    min_w = spectrum%min_wlen(jj)
    do ii=1,n_bins
      histogram(ii,jj) = count((spectrum%points(:,jj).ge.min_w).and.&
      (spectrum%points(:,jj).lt.(min_w+delta_bin_w)))
      min_w = min_w + delta_bin_w  
    enddo
      max_variation(jj) = real((maxval(histogram(:,jj))-minval(histogram(:,jj))),kind=8)/&
      real(spectrum%n_points,kind=8)
  enddo
  call assert_true(all(max_variation.le.tol_updf),&
  "Error generate uniform random spectrum: probability density not uniform!")
  !> clean-up the spectrum class
  call spectrum%deallocate_spectrum

end subroutine test_spectrum_generation_rng_uniform

!> test the integration via Monte-Carlo method (uniform distribution)
subroutine test_spectrum_integration_rng_uniform()
  use constants,               only: PI
  use mod_test_functions,      only: sin2x,int_sin2x
  use mod_common_test_tools,   only: omp_initialize_rngs
  use mod_linear_reg,          only: linear_regression
  use mod_spectra_monte_carlo, only: spectrum_rng_uniform
  !$ use omp_lib
  implicit none
  !> variables
  type(spectrum_rng_uniform)  :: spectrum
  integer                     :: ii,kk,pp,n_threads
  !$ integer                  :: jj
  real*8,dimension(2)         :: std_conv_coeff    !< linear regression coeff. of the std deviation
  real*8,dimension(n_spectra) :: integral,solution !< required for avoiding wrong memory accesses
  real*8,dimension(:,:),allocatable         :: integrands
  real*8,dimension(n_spectra,n_trials)      :: integrals,integrals_taskloop
  real*8,dimension(n_spectra,n_convergence) :: avg_integrals,&
  avg_integrals_taskloop !< average value of the integrals
  real*8,dimension(n_spectra,n_convergence) :: std_dev,std_dev_taskloop !< integral std deviation
  real*8,dimension(n_spectra,n_convergence) :: int_error,int_error_taskloop !< integral error

  !> initialisation
  std_dev = 0.d0; n_threads = 1;
  !$ n_threads = omp_get_max_threads()
  !> loop on the number of convergence points
  do kk=1,n_convergence
    !> initialise structures and grids
    allocate(integrands(n_points_conv(kk),n_spectra))
    call omp_initialize_rngs(n_points_conv(kk),n_threads,rngs)!< change the size of the rng output
    spectrum = spectrum_rng_uniform(n_points_conv(kk),n_spectra,PI*min_angle,PI*max_angle)
    do pp=1,n_trials
      call spectrum%generate_spectrum(rngs)
#ifdef _OPENMP
      !$omp parallel do default(private) shared(spectrum,integrands) collapse(2)
      do jj=1,spectrum%n_spectra
        do ii=1,spectrum%n_points
          integrands(ii,jj) = sin2x(spectrum%points(ii,jj))
        enddo
      enddo
      !$omp end parallel do
#else
      do ii=1,spectrum%n_spectra
        integrands(:,ii) = sin2x(spectrum%n_points,spectrum%points(:,ii))
      enddo
#endif
      !> integrate  and compute error
      call spectrum%integrate_data(integrands,integral)
      integrals(:,pp) = integral
      !> integrate using the taskloop parallelism
      !$omp parallel default(private) shared(spectrum,integrands,integral)
      !$omp single
      call spectrum%integrate_data(integrands,integral)
      !$omp end single
      !$omp end parallel
      integrals_taskloop(:,pp) = integral
    enddo

    !> compute integrand standard deviation and error
    solution = int_sin2x(n_spectra,PI*max_angle)-int_sin2x(n_spectra,PI*min_angle)
    call compute_average_std_dev_error(n_spectra,n_trials,integrals,solution,&
    avg_integrals(:,kk),std_dev(:,kk),int_error(:,kk))
    call compute_average_std_dev_error(n_spectra,n_trials,integrals_taskloop,solution,&
    avg_integrals_taskloop(:,kk),std_dev_taskloop(:,kk),int_error_taskloop(:,kk))

    !> clean up everything
    deallocate(integrands)
    call spectrum%deallocate_spectrum
  enddo

  !> compute the comvergence rate of the standard deviation
  std_dev = sqrt(std_dev/real(n_trials,kind=8))
  std_dev_taskloop = sqrt(std_dev_taskloop/real(n_trials,kind=8))
  do ii=1,n_spectra
    call linear_regression(n_convergence,log10(real(n_points_conv,kind=8)),&
    log10(std_dev(ii,:)),std_conv_coeff)
    !> check that the convergence coefficient of the standard deviation ~0.5
    call assert_equals(expected_std_conv_coeff,std_conv_coeff(1),tol_std_conv_coeff,&
    "Error integrate uniform random spectrum: std convergence rate is not 1/2!")
    call assert_true(int_error(ii,n_convergence).lt.tol_min_int_error,&
    "Error integrate uniform random spectrum: minimum error larger than expected!") 
    call linear_regression(n_convergence,log10(real(n_points_conv,kind=8)),&
    log10(std_dev_taskloop(ii,:)),std_conv_coeff)
    !> check that the convergence coefficient of the standard deviation ~0.5
    call assert_equals(expected_std_conv_coeff,std_conv_coeff(1),tol_std_conv_coeff,&
    "Error integrate uniform random spectrum taskloop: std convergence rate is not 1/2!")
    call assert_true(int_error_taskloop(ii,n_convergence).lt.tol_min_int_error,&
    "Error integrate uniform random spectrum taskloop: minimum error larger than expected!") 
  enddo

end subroutine test_spectrum_integration_rng_uniform

!> Tools ---------------------------------------------------------------------
!> method for computing averages, standard deviatons and error
subroutine compute_average_std_dev_error(n_values,n_trials,integrals,solution,&
avg_integrals,std_dev,error)
  implicit none
  !> inputs:
  integer,intent(in)                             :: n_values,n_trials
  real*8,dimension(n_values),intent(in)          :: solution
  real*8,dimension(n_values,n_trials),intent(in) :: integrals
  !> outputs:
  real*8,dimension(n_values),intent(out) :: avg_integrals,std_dev,error
  !> variables
  integer :: jj
  !> initialisation
  std_dev = 0d0; error = 0d0; avg_integrals = 0d0;
  !> compute the average solution
  avg_integrals = sum(integrals,dim=2)/real(n_trials,kind=8)
  !> compute standard deviation
  do jj=1,n_trials
    std_dev = std_dev + ((integrals(:,jj)-avg_integrals)*&
    (integrals(:,jj)-avg_integrals))
  enddo
  !> compute error
  error = abs(avg_integrals-solution)
end subroutine compute_average_std_dev_error
!>----------------------------------------------------------------------------
end module mod_spectra_monte_carlo_test
