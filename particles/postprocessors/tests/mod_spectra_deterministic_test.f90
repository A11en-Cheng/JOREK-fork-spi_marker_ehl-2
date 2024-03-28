!> mod_spectra_deterministic_test contains variables and procedures
!> for testing the deterministic spectrum integrators
module mod_spectra_deterministic_test
use fruit
implicit none

private
public :: run_fruit_spectra_deterministic

!> Variables --------------------------------------------------------
character(len=33),parameter :: input_file='spectrum_deterministic_2nd_inputs'
integer,parameter :: read_unit=43         !< file read unit
integer,parameter :: n_convergence=5      !< number of points for convergence
integer,parameter :: n_points=512437      !< number of points
integer,parameter :: n_spectra=2          !< number of spectra
real*8,parameter  :: tol_grid=3.d-16      !< tolerance for grid check
real*8,parameter  :: accuracy_order=-2.d0 !< accuracy order
real*8,parameter  :: tol_accuracy=5.d-2   !< tolerance on the accuray order
real*8,parameter  :: tol_int_error=5.d-12 !< tolerance on the minim integratio error
!> n_points for convergence study
integer,dimension(n_convergence) :: n_points_conv=(/997,10725,100000,1003757,10023947/)
real*8,dimension(2),parameter :: min_wlen=(/3.0d-6,3.0d-7/) !< minimum wavelength
real*8,dimension(2),parameter :: max_wlen=(/3.5d-6,4.0d-7/) !< maximum wavelength
real*8,dimension(2),parameter :: min_angle=(/1.75d-1,8.4d1/) !< minimum angle for integration
real*8,dimension(2),parameter :: max_angle=(/3.25d1,1.75d2/) !< maximum angle for integration
real*8,dimension(n_spectra)   :: wbin_size !> size of the wavelength integration interval

!> Interfaces -------------------------------------------------------

contains

!> Fruit basket -----------------------------------------------------
!> run_fruit_spectra_deterministic executes the set-up, runs the tests
!> and clean-up all test features
subroutine run_fruit_spectra_deterministic()
  implicit none
  write(*,'(/A)') "  ... setting-up: spectra deterministic integrator tests"
  call setup
  call write_spectrum_inputs
  write(*,'(/A)') "  ... running: spectra deterministic integrator tests"
  call run_test_case(test_spectrum_input_reader,'test_spectrum_input_reader')
  call run_test_case(test_deterministic_allocation_noinit,&
  'test_deterministic_allocation_noinit')
  call run_test_case(test_deterministic_allocation_init,&
  'test_deterministic_allocation_init')
  call run_test_case(test_set_spectrum_int_2nd_properties,&
  'test_set_spectrum_int_2nd_properties')
  call run_test_case(test_generate_midpoint_spectra,&
  'test_generate_midpoint_spectra')
  call run_test_case(test_2nd_order_rectangle_integrator,&
  'test_2nd_order_rectangle_integrator')
  write(*,'(/A)') "  ... tearing-down: spectra deterministic integrator tests"
  call teardown
end subroutine run_fruit_spectra_deterministic

!> Set-up and tear-down ---------------------------------------------
!> Set-up the test variables
subroutine setup()
  implicit none
  !> set spectrum integration interval
  wbin_size = (max_wlen-min_wlen)/n_points
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

!> Tear-down the test variables
subroutine teardown()
  implicit none
  wbin_size = 0.d0
  call system("rm "//input_file)
end subroutine teardown

!> Tests ------------------------------------------------------------
!> test input reading parameters
subroutine test_spectrum_input_reader()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none 
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
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
  "Error spectrum deterministic 2nd input reader: N# inputs mismatch!")
  call assert_equals((/n_points,n_spectra/),int_param,2,&
  "Error spectrum deterministic 2nd input reader: N# integer parameters mismatch!")
  call assert_equals((/min_wlen(1),min_wlen(2),max_wlen(1),max_wlen(2)/),real_param,4,&
  "Error spectrum deterministic 2nd input reader: N# real parameters mismatch!")
  !> cleanup
  if(allocated(int_param)) deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine test_spectrum_input_reader

!> test the allocation, deallocation and construction of the
!> spectrum_integrator_2nd without initialisation
subroutine test_deterministic_allocation_noinit()
  use mod_assert_equals_tools,   only: assert_equals_allocatable_arrays
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum

  !> test allocation and deallocation
  call spectrum%allocate_spectrum(n_points,n_spectra)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum integration 2nd allocation: n_points do not match!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum integration 2nd allocation: n_spectra do not match!")
  call assert_equals_allocatable_arrays(n_points,n_spectra,spectrum%points,&
  "Error spectrum integration 2nd allocation: points")
  call spectrum%deallocate_spectrum
  call assert_equals(-1,spectrum%n_points,&
  "Error spectrum integration 2nd deallocation: n_points not set to default!")
  call assert_equals(-1,spectrum%n_spectra,&
  "Error spectrum integration 2nd deallocation: n_spectra not set to default!")
  call assert_false(allocated(spectrum%points),&
  "Error spectrum integration 2nd deallocation: points not deallocated!")

  !> test constructor
  spectrum = spectrum_integrator_2nd(n_points,n_spectra)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum integration 2nd construction: n_points do not match!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum integration 2nd construction: n_spectra do not match!")
  call assert_equals_allocatable_arrays(n_points,n_spectra,spectrum%points,&
  "Error spectrum integration 2nd construction: points")
  call spectrum%deallocate_spectrum

end subroutine test_deterministic_allocation_noinit

!> test allocation, deallocation and construction of the
!> spectrum_integrator_2nd with initialisation
subroutine test_deterministic_allocation_init()
  use mod_assert_equals_tools,   only: assert_equals_allocatable_arrays
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  real*8,dimension(2*n_spectra) :: real8_param

  !> initialisation
  real8_param(1:n_spectra) = min_wlen
  real8_param(n_spectra+1:2*n_spectra) = max_wlen
  !> test allocation and deallocation with initialisation
  call spectrum%allocate_spectrum(n_points,n_spectra,real8_param)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum integration 2nd allocation init: n_points do not match!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum integration 2nd allocation init: n_spectra do not match!")
  call assert_equals_allocatable_arrays(n_spectra,min_wlen,spectrum%min_wlen,&
  tol_grid,"Error spectrum integration 2nd allocation init: min wavelengths")
  call assert_equals_allocatable_arrays(n_spectra,spectrum%wbin_size,&
  "Error spectrum integration 2nd allocation init: wavelengths bin size")
  call assert_equals_allocatable_arrays(n_points,n_spectra,spectrum%points,&
  "Error spectrum integration 2nd allocation init: points")
  call spectrum%deallocate_spectrum
  call assert_equals(-1,spectrum%n_points,&
  "Error spectrum integration 2nd deallocation init: n_points not set to default!")
  call assert_equals(-1,spectrum%n_spectra,&
  "Error spectrum integration 2nd deallocation init: n_spectra not set to default!")
  call assert_false(allocated(spectrum%points),&
  "Error spectrum integration 2nd deallocation init: points not deallocated!")
  call assert_false(allocated(spectrum%min_wlen),&
  "Error spectrum integration 2nd deallocation init: min_wlen not deallocated!")
  call assert_false(allocated(spectrum%wbin_size),&
  "Error spectrum integration 2nd deallocation init: wbin_size not deallocated!")

  !> test construction with initialisation
  spectrum = spectrum_integrator_2nd(n_points,n_spectra,min_wlen,max_wlen)
  call assert_equals(n_points,spectrum%n_points,&
  "Error spectrum integration 2nd construction init: n_points do not match!")
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum integration 2nd construction init: n_spectra do not match!")
  call assert_equals_allocatable_arrays(n_spectra,min_wlen,spectrum%min_wlen,&
  tol_grid,"Error spectrum integration 2nd construction init: min wavelengths")
  call assert_equals_allocatable_arrays(n_spectra,spectrum%wbin_size,&
  "Error spectrum integration 2nd construction init: wavelengths bin size")
  call assert_equals_allocatable_arrays(n_points,n_spectra,spectrum%points,&
  "Error spectrum integration 2nd construction init: points")
  call spectrum%deallocate_spectrum
end subroutine test_deterministic_allocation_init

!> test set properties deterministic integrator 2nd order 
subroutine test_set_spectrum_int_2nd_properties()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  real*8,dimension(n_spectra),parameter :: min_wlen_2=(/1.d-5,5.d-9/)
  real*8,dimension(n_spectra),parameter :: max_wlen_2=(/7.4d-5,1.25d-8/)
  real*8,dimension(n_spectra)           :: wbin_size_2
  real*8,dimension(2*n_spectra)         :: min_wlen_3,max_wlen_3,wbin_size_3

  !> initialisation
  wbin_size_2 = (max_wlen_2-min_wlen_2)/n_points
  min_wlen_3(1:n_spectra) = min_wlen; min_wlen_3(n_spectra+1:2*n_spectra) = min_wlen_2;
  max_wlen_3(1:n_spectra) = max_wlen; max_wlen_3(n_spectra+1:2*n_spectra) = max_wlen_2;
  wbin_size_3 = (max_wlen_3-min_wlen_3)/n_points

  !> test settings on unallocated properties
  spectrum = spectrum_integrator_2nd(n_points,n_spectra)
  call spectrum%set_spectrum_interval(n_spectra,min_wlen,max_wlen)
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum integration set interval not allocated: n_spectra do not match!")
  call assert_equals(min_wlen,spectrum%min_wlen,n_spectra,&
  "Error spectrum integration set interval not allocated: min_wlen mismatch!")
  call assert_equals(wbin_size,spectrum%wbin_size,n_spectra,&
  "Error spectrum integration set interval not allocated: wbin_size mismatch!")

  !> Test change of interval alrady allocated
  call spectrum%set_spectrum_interval(n_spectra,min_wlen_2,max_wlen_2)
  call assert_equals(n_spectra,spectrum%n_spectra,&
  "Error spectrum integration set interval allocated: n_spectra do not match!")
  call assert_equals(min_wlen_2,spectrum%min_wlen,n_spectra,&
  "Error spectrum integration set interval allocated: min_wlen mismatch!")
  call assert_equals(wbin_size_2,spectrum%wbin_size,n_spectra,&
  "Error spectrum integration set interval allocated: wbin_size mismatch!")

  !> Test change of interval alrady with re-allocated
  call spectrum%set_spectrum_interval(2*n_spectra,min_wlen_3,max_wlen_3)
  call assert_equals(spectrum%n_spectra,2*n_spectra,&
  "Error spectrum integration set interval reallocated: n_spectra do not match!")
  call assert_equals(min_wlen_3,spectrum%min_wlen,2*n_spectra,&
  "Error spectrum integration set interval reallocated: min_wlen mismatch!")
  call assert_equals(wbin_size_3,spectrum%wbin_size,2*n_spectra,&
  "Error spectrum integration set interval reallocated: wbin_size mismatch!")

  !> cleanup
  call spectrum%deallocate_spectrum
end subroutine test_set_spectrum_int_2nd_properties

!> test the generation of midpoint grids
subroutine test_generate_midpoint_spectra()
  use mod_spectra_deterministic, only: spectrum_integrator_2nd
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  integer :: ii,jj
  real*8,dimension(n_points+1) :: interval_nodes,grid_nodes

  !> initialise variables
  spectrum = spectrum_integrator_2nd(n_points,n_spectra,min_wlen,max_wlen)
  !> generate grids and check correctness
  call spectrum%generate_spectrum
  do jj=1,spectrum%n_spectra
    do ii=1,spectrum%n_points+1
      interval_nodes(ii) = min_wlen(jj) +  wbin_size(jj)*real(ii-1,kind=8)
    enddo
    grid_nodes(1:n_points) = (spectrum%points(:,jj)-5.d-1*spectrum%wbin_size(jj))
    grid_nodes(n_points+1) = (spectrum%points(n_points,jj)+5.d-1*spectrum%wbin_size(jj))
    grid_nodes = grid_nodes/interval_nodes
    interval_nodes = 1.d0
    call assert_equals(interval_nodes,grid_nodes,n_points,tol_grid,&
    "Error spectrum integration generate spectrum: spectral grid mismatch!")
  enddo
  !> cleanaup
  call spectrum%deallocate_spectrum

end subroutine test_generate_midpoint_spectra

!> test the 2nd order rectangle method integrator. The reltive error 
!> is used due to the large value of the integral.
subroutine test_2nd_order_rectangle_integrator()
  use omp_lib
  use constants,                 only: PI
  use mod_test_functions,        only: expxsin2x,int_expxsin2x
  use mod_linear_reg,            only: linear_regression
  use mod_spectra_deterministic, only: spectrum_integrator_2nd 
  implicit none
  !> variables
  type(spectrum_integrator_2nd) :: spectrum
  integer                       :: ii,kk
  !$ integer                    :: jj
  real*8,dimension(2)           :: conv_coeff !< convergence coeff. from linear regression
  real*8,dimension(n_spectra)   :: integral,integral_taskloop,solution
  real*8,dimension(:,:),allocatable         :: integrands
  real*8,dimension(n_spectra,n_convergence) :: rel_int_error !< integretion error
  real*8,dimension(n_spectra,n_convergence) :: rel_int_error_taskloop !< integretion error
  !> convergence loop
  do kk=1,n_convergence
    !> initialise
    allocate(integrands(n_points_conv(kk),n_spectra))
    spectrum = spectrum_integrator_2nd(n_points_conv(kk),n_spectra,min_angle,max_angle)
    call spectrum%generate_spectrum
#ifdef _OPENMP
   !$omp parallel do default(private) shared(spectrum,integrands) collapse(2)
   do jj=1,spectrum%n_spectra
     do ii=1,spectrum%n_points
       integrands(ii,jj) = expxsin2x(spectrum%points(ii,jj))
     enddo
   enddo
   !$omp end parallel do
#else
    do ii=1,spectrum%n_spectra
      integrands(:,ii) = expxsin2x(spectrum%n_points,spectrum%points(:,ii))
    enddo
#endif
    call spectrum%integrate_data(integrands,integral) !< integrate values
    !> integral values in taskloop
    !$omp parallel default(private) shared(spectrum,integrands,integral_taskloop)
    !$omp single
    call spectrum%integrate_data(integrands,integral_taskloop)
    !$omp end single
    !$omp end parallel
    solution = int_expxsin2x(n_spectra,max_angle)-int_expxsin2x(n_spectra,min_angle)
    rel_int_error(:,kk) = abs((integral-solution)/solution)!< compute the error  
    rel_int_error_taskloop(:,kk) = abs((integral_taskloop-solution)/solution)!< compute the error  
    deallocate(integrands)
    call spectrum%deallocate_spectrum
  enddo

  !> compute convergence rate and check for error
  do ii=1,n_spectra
    call linear_regression(n_convergence,log10(real(n_points_conv,kind=8)),&
    log10(rel_int_error(ii,:)),conv_coeff)
    call assert_equals(accuracy_order,conv_coeff(1),tol_accuracy,&
    "Error spectrum integration: expected accuracy order not matched!")
    call assert_true(rel_int_error(ii,n_convergence).lt.tol_int_error,&
    "Error spectrum integration: expected minimum error not achieved!")
    call linear_regression(n_convergence,log10(real(n_points_conv,kind=8)),&
    log10(rel_int_error_taskloop(ii,:)),conv_coeff)
    call assert_equals(accuracy_order,conv_coeff(1),tol_accuracy,&
    "Error spectrum integration taskloop: expected accuracy order not matched!")
    call assert_true(rel_int_error_taskloop(ii,n_convergence).lt.tol_int_error,&
    "Error spectrum integration taskloop: expected minimum error not achieved!")
  enddo
end subroutine test_2nd_order_rectangle_integrator


!>-------------------------------------------------------------------
end module mod_spectra_deterministic_test
