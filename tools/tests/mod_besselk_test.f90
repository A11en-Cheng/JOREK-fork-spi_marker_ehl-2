! mod_besselk_test contains procedures for testing procedured
! computing the modified bessel functions of fractional order
module mod_besselk_test
use fruit 
use mod_dynamic_array_tools, only: allocate_check
use mod_dynamic_array_tools, only: deallocate_check
implicit none

private
public :: run_fruit_besselk

!> variables common to all tests
integer,parameter               :: max_it=100000   !< maximum number of iterations
real*8,parameter                :: test_tol=1.d-13 !< tolerance for success
! order of the bessel function for testing
integer :: Nnu=3
real*8,dimension(3),parameter :: nu=(/1.d0/3.d0,2.d0/3.d0,5.d0/3.d0/)
! values of x used for testing the bessel function for x>=x_val=2
integer                             :: Nx_ge_2=10
real*8,dimension(10),parameter :: x_ge_2=(/7.12d0,8.17d0,2.22d0,19.2d0,&
                                       4.12d0,3.7d0,6.05d0,21.18d0,14.4d0,2.d0/)
! values of x used for testing the bessel function for x<x_val=2
integer                             :: Nx_lt_2=10
real*8,dimension(10),parameter :: x_lt_2=(/1.d-4,1.82d0,0.36d0,5.28d-3,&
                                       2.9d-6,0.27d0,1.74d0,1.d0,1.d-8,1.99d0/)
! values of the bessel function 2nd kind and its derivatives for x_ge_2 and 
! nu=1/3,2/3,5/2 as computed by matlab (mat) and python (py)
real*8,dimension(30),parameter :: bknu_ge_2_mat = (/&
     (/3.764170850828786d-4,1.231043600701636d-4,8.899156804745512d-2,&
     1.307444784232888d-9,9.879203508537923d-3,1.584159892509702d-2,&
     1.188697476014876d-3,1.719285416494302d-10,1.832364180921524d-7,&
     1.165449612961646d-1/),(/3.847728624956395d-4,1.255008926911841d-4,&
     9.476954484002706d-2,1.318559733980535d-9,1.024573759899036d-2,&
     1.649057589808570d-2,1.219489768014882d-3,1.732560154686813d-10,&
     1.852997068247446d-7,1.248389274881278d-1/),(/4.484719282468561d-4,&
     1.435859419047102d-4,1.459102135970209d-1,1.399011432425981d-9,&
     1.319497619429532d-2,2.178414879828106d-2,1.457455551610443d-3,&
     1.828354359187838d-10,2.003937983537028d-7,1.997709129549173d-1/)/)
real*8,dimension(30),parameter :: bknu_ge_2_py = (/&
     (/3.7641708508287861d-4,1.2310436007016362d-4,8.8991568047455116d-2,&
       1.3074447842328883d-9,9.8792035085379232d-3,1.5841598925097022d-2,&
       1.1886974760148762d-3,1.7192854164943015d-10,1.8323641809215239d-7,&
       1.1654496129616169d-1,3.8477286249563949d-4,1.2550089269118409d-4,&
       9.4769544840027059d-2,1.3185597339805353d-9,1.0245737598990364d-2,&
       1.6490575898085696d-2,1.2194897680148818d-3,1.7325601546868132d-10,&
       1.8529970682474454d-7,1.2483892748813481d-1/),(/4.4847192824685610d-4,&
       1.4358594190471018d-4,1.4591021359702089d-1,1.3990114324259808d-9,&
       1.3194976194295318d-2,2.1784148798281056d-2,1.4574555516104425d-3,&
       1.8283543591878372d-10,2.0039379835370282d-7,1.9977091295491711d-1/)/)
! values of the bessel function 2nd kind and its derivatives for x_lt_2
! and nu=1/3,2/3,5/2 as computed by matlab (mat) and python (py)
real*8,dimension(30),parameter :: bknu_lt_2_mat = (/&
     (/3.628396070100760d1,1.459197599310183d-1,1.313842449151210d0,&
     9.411032704078833d0,1.183202403880544d2,1.626582338799409d0,&
     1.614737084812082d-1,4.384306334415338d-1,7.833229062430851d2,&
     1.179972353239948d-1/),(/4.988585910043111d2,1.571933403240230d-1,&
     1.679684673610099d0,3.540754405283199d1,5.285034257170110d3,&
     2.179872025239006d0,1.744438023694717d-1,4.944750621042081d-1,&
     2.315509105323807d5,1.264314529006558d-1/),(/6.651484164018186d6,&
     2.610797162123548d-1,7.534896795855283d0,8.950710035944485d3,&
     2.429900926214546d9,1.239138246343648d1,2.951471202585817d-1,&
     1.097730716247145d0,3.087345473843414d13,2.027084265472155d-1/)/)
real*8,dimension(30),parameter :: bknu_lt_2_py = (/&
     (/3.6283960701007608d1,1.4591975993101533d-1,1.3138424491512091d0,&
     9.4110327040788331d0,1.1832024038805436d2,1.6265823387994081d0,&
     1.6147370848120510d-1,4.3843063344153255d-1,7.8332290624308507d2,&
     1.1799723532399106d-1/),(/4.9885859100431156d2,1.5719334032402899d-1,&
     1.6796846736101023d0,3.5407544052832023d1,5.2850342571701158d3,&
     2.1798720252390091d0,1.7444380236947751d-1,4.9447506210421149d-1,&
     2.3155091053238098d5,1.2643145290066365d-1/),(/6.6514841640181877d6,&
     2.6107971621235415d-1,7.5348967958552855d0,8.9507100359444867d3,&
     2.4299009262145467d9,1.2391382463436482d1,2.9514712025858114d-1,&
     1.0977307162471459d0,3.0873454738434148d13,2.0270842654721480d-1/)/)
integer                           :: Nx ! size of the x array
real*8,dimension(:,:),allocatable :: x_all            !< input array x
real*8,dimension(:,:),allocatable :: bknu_mat,bknu_py !< shuffled matlab/python solutions

contains

! Tests baskets -----------------------------------------------------

! run_fruit_besselk performs the setup,
! run the tests and performs the teardown of
! the besselk functions
subroutine run_fruit_besselk()
  implicit none

  write(*,'(/A)') "  ... setting-up: besselk tests" 
  call setup !< setup test variables
  write(*,'(/A)') "  ... running: besselk tests"
  call run_test_case(test_f_besselk,'test_f_besselk')
  call run_test_case(test_besselk,'test_besselk')
  call run_test_case(test_besselk_x_array,'test_besselk_x_array')
  call run_test_case(test_besselk_nu_array,'test_besselk_nu_array')
  call run_test_case(test_besselk_x_nu_array,'test_besselk_x_nu_array')
#ifdef UNIT_TESTS
  call run_test_case(test_handle_float_exceptions,&
  'test_handle_float_exceptions')
#endif
  write(*,'(/A)') "  ... tearing-down: besselk tests" 
  call teardown !< cleanup test variables  
  
end subroutine run_fruit_besselk

! Set-up and tear-down ---------------------------------------------

! setup initialises the module variables
! inputs:
!   Nx_val:    (integer)(optional) size of the x array
!   Nx_lt_val: (integer)(optional) number of x lower than threshold
!   Nx_ge_val: (integer)(optional) number of x higher than threshold
!   x_val:     (doiuble)(optional) threshold value
!   x_min:     (double)(optional) lower bound of x
!   x_max:     (double)(optional) upper bound of x
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_array_norep
  implicit none

  ! variables
  integer :: ii,ierr
  integer,dimension(:),allocatable :: ids
  real*8,dimension(:),allocatable  :: rnd

  ! set default inputs
  ierr = 0
  Nx = Nx_ge_2 + Nx_lt_2

  ! allocate arrays for tests
  call allocate_check(Nx,rnd)
  call allocate_check(Nx,ids)
  call allocate_check(Nx,Nnu,x_all)
  call allocate_check(Nx,Nnu,bknu_mat)
  call allocate_check(Nx,Nnu,bknu_py)
  ! initialise arrays
  rnd=0.d0;

  ! randomize matlab and python solutions for testing
  do ii=1,Nnu
    call gnu_rng_array_norep(Nx,(/1,Nx/),ids,ierr,max_it)
    call assert_equals(0,ierr,&
    "Error, setup: generate random indices failed reached max_it")
    x_all(ids(1:Nx_lt_2),ii)    = x_lt_2
    x_all(ids(Nx_lt_2+1:Nx),ii) = x_ge_2
    bknu_mat(ids(1:Nx_lt_2),ii)    = bknu_lt_2_mat((ii-1)*Nx_lt_2+1:ii*Nx_lt_2)
    bknu_mat(ids(Nx_lt_2+1:Nx),ii) = bknu_ge_2_mat((ii-1)*Nx_ge_2+1:ii*Nx_ge_2)
    bknu_py(ids(1:Nx_lt_2),ii)     = bknu_lt_2_py((ii-1)*Nx_lt_2+1:ii*Nx_lt_2)
    bknu_py(ids(Nx_lt_2+1:Nx),ii)  = bknu_ge_2_py((ii-1)*Nx_ge_2+1:ii*Nx_ge_2)
  enddo


  ! deallocate arrays
  call deallocate_check(ids)
  call deallocate_check(rnd)

end subroutine setup

! teardown cleans-up the unit tests
subroutine teardown()
  implicit none
  call deallocate_check(x_all)
  call deallocate_check(bknu_mat)
  call deallocate_check(bknu_py)
end subroutine teardown

! Tests ------------------------------------------------------------
! test_besselk tests the computation of the modified bessel
! function of the second kind and fractional order using
! randomized matlab and python solutions. Test the function.
subroutine test_f_besselk()
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_besselk, only: f_besselk
  implicit none

  ! variables
  integer,parameter :: idnu(2)=(/1,2/)
  integer,dimension(2) :: ids
  real*8  :: bknu_1,bknu_2
  real*8  :: nu_1,nu_2
  real*8  :: x_1,x_2
  real*8,dimension(2) :: bknu_ref_mat
  real*8,dimension(2) :: bknu_ref_py

  ! initialisation
  call gnu_rng_interval(2,(/1,Nx/),ids)
  x_1 = x_all(ids(1),idnu(1)); x_2 = x_all(ids(2),idnu(2));
  nu_1 = nu(idnu(1)); nu_2 = nu(idnu(2));
  bknu_ref_mat = (/bknu_mat(ids(1),idnu(1)),bknu_mat(ids(2),idnu(2))/);
  bknu_ref_py = (/bknu_py(ids(1),idnu(1)),bknu_py(ids(2),idnu(2))/);

  ! compute modified bessel function for multiple nu
  bknu_1 = f_besselk(nu_1,x_1)
  bknu_2 = f_besselk(nu_2,x_2)

  ! check solution
  call assert_equals(1.d0,bknu_1/bknu_ref_mat(1),test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (function)")
  call assert_equals(1.d0,bknu_1/bknu_ref_py(1),test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (function)")
  call assert_equals(1.d0,bknu_2/bknu_ref_mat(2),test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (function)")
  call assert_equals(1.d0,bknu_2/bknu_ref_py(2),test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (function)")

end subroutine test_f_besselk

! test_besselk tests the computation of the modified bessel
! function of the second kind and fractional order using
! randomized matlab and python solutions. Test the subroutine.
subroutine test_besselk()
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_besselk, only: besselk
  implicit none

  ! variables
  integer,parameter :: idnu(2)=(/1,2/)
  integer,dimension(2) :: ids
  real*8  :: bknu_1,bknu_2
  real*8  :: nu_1,nu_2
  real*8  :: x_1,x_2
  real*8,dimension(2) :: bknu_ref_mat
  real*8,dimension(2) :: bknu_ref_py

  ! initialisation
  call gnu_rng_interval(2,(/1,Nx/),ids)
  x_1 = x_all(ids(1),idnu(1)); x_2 = x_all(ids(2),idnu(2));
  nu_1 = nu(idnu(1)); nu_2 = nu(idnu(2));
  bknu_ref_mat = (/bknu_mat(ids(1),idnu(1)),bknu_mat(ids(2),idnu(2))/);
  bknu_ref_py = (/bknu_py(ids(1),idnu(1)),bknu_py(ids(2),idnu(2))/);

  ! compute modified bessel function for multiple nu
  call besselk(nu_1,x_1,bknu_1)
  call besselk(nu_2,x_2,bknu_2)

  ! check solution
  call assert_equals(1.d0,bknu_1/bknu_ref_mat(1),test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (subroutine)")
  call assert_equals(1.d0,bknu_1/bknu_ref_py(1),test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (subroutine)")
  call assert_equals(1.d0,bknu_2/bknu_ref_mat(2),test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (subroutine)")
  call assert_equals(1.d0,bknu_2/bknu_ref_py(2),test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (subroutine)")

end subroutine test_besselk

! test_besselk_x_array tests the computation of the modified
! bessel function of the second kind and fractional order when 
! applied to arrays of x
subroutine test_besselk_x_array()
  use mod_besselk, only: besselk
  implicit none
  
  ! variables
  integer,parameter :: Nnu_loc=2
  integer,parameter :: idnu(2)=(/1,2/)
  integer :: ii
  real*8  :: nu_loc
  real*8,dimension(Nx)         :: x_loc
  real*8,dimension(Nx*Nnu_loc) :: bknu
  real*8,dimension(Nx*Nnu_loc) :: bknu_mat_loc,bknu_py_loc

  ! compute modified bessel function for multiple nu
  do ii=1,Nnu_loc
    nu_loc = nu(idnu(ii))
    x_loc = x_all(:,idnu(ii))
    bknu_mat_loc((ii-1)*Nx+1:ii*Nx) = bknu_mat(:,idnu(ii))
    bknu_py_loc((ii-1)*Nx+1:ii*Nx) = bknu_py(:,idnu(ii))
    call besselk(Nx,nu_loc,x_loc,bknu((ii-1)*Nx+1:ii*Nx))
  enddo

  ! check solution
  call assert_equals(bknu_mat_loc/bknu_mat_loc,bknu/bknu_mat_loc,Nx*Nnu_loc,test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (x-array)")
  call assert_equals(bknu_py_loc/bknu_py_loc,bknu/bknu_py_loc,Nx*Nnu_loc,test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (x-array)")

end subroutine test_besselk_x_array

! test_besselk_nu_array tests the computation of the modified
! bessel function of the second kind and fractional order when 
! applied to arrays of nu
subroutine test_besselk_nu_array()
  use mod_gnu_rng, only: gnu_rng_interval
  use mod_besselk, only: besselk
  implicit none

  ! variables
  integer :: ii,id
  real*8 :: x_loc
  real*8,dimension(Nnu) :: bknu,bknu_loc

  ! initialization
  call gnu_rng_interval((/1,Nx/),id)
  
  ! compute solution
  do ii=1,Nnu
     x_loc = x_all(id,ii)
     call besselk(Nnu,nu,x_loc,bknu_loc)
     bknu(ii) = bknu_loc(ii)
  enddo

  ! check solutions
  call assert_equals(bknu_mat(id,:)/bknu_mat(id,:),bknu/bknu_mat(id,:),Nnu,test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (nu-array)")
  call assert_equals(bknu_py(id,:)/bknu_py(id,:),bknu/bknu_py(id,:),Nnu,test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (nu-array)")

end subroutine test_besselk_nu_array

! test_besselk_x_nu_array tests the computation of the modified
! bessel function of the second kind and fractional order when 
! applied to arrays of x and nu
subroutine test_besselk_x_nu_array()
  use mod_besselk, only: besselk
  implicit none

  ! variables
  integer :: ii
  real*8,dimension(Nx,Nnu) :: bknu,bknu_loc
  
  ! compute solution
  do ii=1,Nnu
    call besselk(Nx,Nnu,nu,x_all(:,ii),bknu_loc)
    bknu(:,ii) = bknu_loc(:,ii)
  enddo

  ! check solution
  call assert_equals(bknu_mat/bknu_mat,bknu/bknu_mat,Nx,Nnu,test_tol,&
  "Error: no match between Matlab and JOREK modified bessel function 2nd kind (x-nu-arrays)")
  call assert_equals(bknu_mat/bknu_mat,bknu/bknu_mat,Nx,Nnu,test_tol,&
  "Error: no match between Python and JOREK modified bessel function 2nd kind (x-nu-arrays)")

end subroutine test_besselk_x_nu_array

#ifdef UNIT_TESTS
! test the floating point exception handlers of the besselk functions
subroutine test_handle_float_exceptions()
  use, intrinsic :: ieee_arithmetic, only: IEEE_Value,IEEE_QUIET_NAN
  use, intrinsic :: ieee_arithmetic, only: IEEE_POSITIVE_INF,IEEE_NEGATIVE_INF
  use mod_besselk, only: handle_float_exceptions
  implicit none
  ! tests
  call assert_equals(0d0,handle_float_exceptions(IEEE_Value(0d0,IEEE_QUIET_NAN)),0d0,&
  "Error: handle of nan by besselk float exception handler failed!")
  call assert_equals(0d0,handle_float_exceptions(IEEE_Value(0d0,IEEE_POSITIVE_INF)),0d0,&
  "Error: handle of positive overflows by besselk float exception handler failed!") 
  call assert_equals(0d0,handle_float_exceptions(IEEE_Value(0d0,IEEE_NEGATIVE_INF)),0d0,&
  "Error: handle of negative overflows by besselk float exception handler failed!")
  call assert_equals(0d0,handle_float_exceptions(tiny(0d0)*tiny(0d0)),0d0,&
  "Error: handle of positive underflows by besselk float exception handler failed!")
  call assert_equals(0d0,handle_float_exceptions(-tiny(0d0)*tiny(0d0)),0d0,&
  "Error: handle of negative underflows by besselk float exception handler failed!")
  call assert_equals(2.35d0,handle_float_exceptions(2.35d0),0d0,&
  "Error: return of variables by besselk float exception handler failed!")
end subroutine test_handle_float_exceptions 
#endif

end module mod_besselk_test
