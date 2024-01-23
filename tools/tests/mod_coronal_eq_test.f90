!> Module testing Coronal equilibrium for selected elements
!> IMPORTANT NOTE: the original set of adas data have been 
!> shortened specifically, 96_li, 96_n 50_w, 96_c adas datafile
!> have been removed since the respective xcd, ccd files
!> and/or prc datafiles downloaded from the ADAS presents errors
module mod_coronal_eq_test
use fruit
use mod_openadas, only: ADF11_all
implicit none
private
public :: run_fruit_coronal_eq
!> Variables --------------------------------------
integer,parameter                       :: set_len=200
integer,parameter                       :: n_sets=4
real*8,parameter                        :: sleep_time=0.5
logical,parameter                       :: write_coronal=.false.
character(len=5),dimension(n_sets),parameter :: sets=(/'89_b ','96_n ',&
                                                       '89_ar','96_he'/)
character(len=2),dimension(n_sets),parameter :: type_imp=(/'B ','N ','Ar','He'/)
type(ADF11_all),dimension(n_sets) :: adas
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_coronal_eq
  implicit none
  write(*,'(/A)') "  ... setting-up: coronal eq"
  call setup
  write(*,'(/A)') "  ... running: coronal eq"
  call run_test_case(test_coronal_Z_increasing,'test_coronal_Z_increasing')
  write(*,'(/A)') "  ... tearing-down: coronal eq"
  call teardown
end subroutine run_fruit_coronal_eq

!> Donwload, read and store necessary files from the open-adas website
subroutine setup()
  use phys_module,  only: imp_type
  use mod_openadas, only: read_adf11
  implicit none
  integer :: ii
  character(len=set_len) :: sleep_command
  !> copy adas data from regression test folder
  call system('cp reg_tests/unit_tests/adas_files/*.dat .')
  !> wait until system finishes
  write(sleep_command,'(A,F0.3)') 'sleep ',sleep_time
  call system(sleep_command)
  do ii=1,n_sets
    adas(ii) = read_adf11(0,trim(sets(ii)))
  enddo
  !> set the impurity type
  imp_type(1:n_sets) = type_imp
end subroutine setup

subroutine teardown()
  use phys_module,  only: imp_type
  implicit none
  integer :: ii
  imp_type(1:n_sets) = ''
  do ii=1,n_sets
    call system('rm *'//trim(sets(ii))//'.dat')
  enddo
end subroutine teardown

!> Tests ------------------------------------------

!> For a single density test that <Z> is increasing with T
subroutine test_coronal_Z_increasing
  use mod_coronal,  only: coronal
  use mod_impurity, only: output_coronal
  implicit none
  integer,parameter           :: T_string_len=10
  type(coronal)               :: cor 
  integer                     :: ii,jj,set_id
  real*8                      :: Z_eff,Z_eff_old
  character(len=T_string_len) :: T_string
  do set_id=1,n_sets
    cor = coronal(adas(set_id))
    if(write_coronal) call output_coronal(cor,set_id)
    Z_eff_old = 0.d0
    do ii=1,size(cor%temperature,1)
      Z_eff = sum(cor%Z(1,ii,:)*(/(jj,jj=0,cor%n_Z)/))
      write(T_string,'(G10.3)') 10d0**(cor%temperature(ii))
      call assert_true(Z_eff.ge.Z_eff_old,&
      'Error coronal Z increasing test: <Z> must increase with T ('//&
      trim(sets(set_id))//') T_e='//trim(T_string)//'K')
      Z_eff_old = Z_eff
    enddo
  enddo
end subroutine test_coronal_Z_increasing

!> ------------------------------------------------
end module mod_coronal_eq_test
