module mod_moivre_series_test
use fruit
use fruit_mpi
#ifndef UNIT_TESTS
use mod_parameters, only: n_tor,n_period
#endif
implicit none
private
public :: run_fruit_moivre_series
integer,parameter :: id_string_len=3
integer,parameter :: n_tor_test=7
integer,parameter :: n_period_test=16
integer,parameter :: n_max=40
integer,parameter :: n_phi=4
real*8,parameter  :: tol=1d-14
real*8,dimension(n_phi),parameter :: phi0=(/1d-1,5d-1,1d0,1.2d1/)
#ifdef UNIT_TESTS
real*8,dimension(n_tor_test) :: HZ,dHZ
#else
real*8,dimension(n_tor)      :: HZ,dHZ
#endif
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_moivre_series
  implicit none
  write(*,'(/A)') "  ... setting-up: moivre series"
  write(*,'(/A)') "  ... running: moivre series"
  call run_test_case(test_moivre_sincos,'test_moivre_sincos')
  call run_test_case(test_mode_moivre,'test_mode_moivre')
  !call run_test_case(,'')
  write(*,'(/A)') "  ... tearing-down: moivre series"
end subroutine run_fruit_moivre_series

!> Tests ------------------------------------------
!> Not clear which function this test is testing
subroutine test_moivre_sincos()
  implicit none
  integer :: ii,jj
  real*8  :: e1r,e1i,enr,eni,enr_tmp
  character(len=id_string_len) :: id_string
  do ii=1,n_phi
    e1r=cos(phi0(ii)); e1i=sin(phi0(ii));
    enr=1.d0; eni=0.d0;
    do jj=1,n_max
      enr_tmp=enr*e1r-eni*e1i; eni=eni*e1r+enr*e1i;
      enr=enr_tmp; write(id_string,'(I3)') jj
      call assert_equals(cos(real(jj,kind=8)*phi0(ii)),enr,tol,&
      'Error moivre sincos test: cosine mismatch id: '//id_string//'!')
      call assert_equals(sin(real(jj,kind=8)*phi0(ii)),eni,tol,&
      'Error moivre sincos test: sine mismatch id: '//id_string//'!')
    enddo
  enddo
end subroutine test_moivre_sincos

#ifdef UNIT_TESTS
subroutine test_mode_moivre
  use mod_interp, only: mode_moivre
  implicit none
  integer                      :: ii,jj
  real*8                       :: angle
  character(len=id_string_len) :: id_string
  do jj=1,n_phi
    call mode_moivre(phi0(jj),HZ,n_tor_test,n_period_test)
    call assert_equals(1.d0,HZ(1),tol,&
    'Error mode moivre: cos n=0 mismatch!')
    if(size(HZ)>1) then
      do ii=1,(n_tor_test-1)/2
        write(id_string,'(I3)') ii
        angle = real(n_period_test*ii,kind=8)*phi0(jj)
        call assert_equals(cos(angle),HZ(2*ii),tol,&
       'Error mode moivre: cos mismatch id: '//id_string//'!')
        call assert_equals(sin(angle),HZ(2*ii+1),tol,&
       'Error mode moivre: sin mismatch id: '//id_string//'!')
      enddo
    endif
  enddo 
end subroutine test_mode_moivre
#else
subroutine test_mode_moivre
  use mod_interp, only: mode_moivre
  implicit none
  integer                      :: ii,jj
  real*8                       :: angle
  character(len=id_string_len) :: id_string
  do jj=1,n_phi
    call mode_moivre(phi0(jj),HZ)
    call assert_equals(1.d0,HZ(1),tol,&
    'Error mode moivre: cos n=0 mismatch!')
    if(size(HZ)>1) then
      do ii=1,(n_tor-1)/2
        write(id_string,'(I3)') ii
        angle = real(n_period*ii,kind=8)*phi0(jj)
        call assert_equals(cos(angle),HZ(2*ii),tol,&
       'Error mode moivre: cos mismatch id: '//id_string//'!')
        call assert_equals(sin(angle),HZ(2*ii+1),tol,&
       'Error mode moivre: sin mismatch id: '//id_string//'!')
      enddo
    endif
  enddo 
end subroutine test_mode_moivre
#endif

#ifdef UNIT_TESTS
subroutine test_sincosperiod_moivre
  use mod_interp, only: sincosperiod_moivre
  implicit none
  integer                      :: ii,jj
  real*8                       :: cosphi,sinphi
  character(len=id_string_len) :: id_string 
  do jj=1,n_phi
    call sincosperiod_moivre(phi0(jj),HZ,dHZ,n_tor_test,n_period_test)
    call assert_equals(1.d0,HZ(1),tol,&
    'Error mode sincosperiod moivre: cos n=0 mismatch!')
    call assert_equals(0.d0,dHZ(1),tol,&
    'Error mode sincosperiod moivre: dcos n=0 mismatch!')
    if(size(HZ)>1) then
      do ii=1,(n_tor_test-1)/2
        write(id_string,'(I3)') ii
        cosphi = cos(real(n_period_test*ii,kind=8)*phi0(jj))
        sinphi = sin(real(n_period_test*ii,kind=8)*phi0(jj))
        call assert_equals(cosphi,HZ(2*ii),tol,&
        'Error mode sincosperiod moivre: cos mismatch id: '//id_string//'!')
        call assert_equals(sinphi,HZ(2*ii+1),tol,&
        'Error mode sincosperiod moivre: sin mismatch id: '//id_string//'!')
        call assert_equals(n_period_test*sinphi,dHZ(2*ii),tol,&
        'Error mode sincosperiod moivre: dcos mismatch id: '//id_string//'!')
        call assert_equals(-n_period_test*cosphi,dHZ(2*ii+1),tol,&
        'Error mode sincosperiod moivre: dsin mismatch id: '//id_string//'!')
      enddo
    endif
  enddo
end subroutine test_sincosperiod_moivre
#else
subroutine test_sincosperiod_moivre
  use mod_interp, only: sincosperiod_moivre
  implicit none
  integer                      :: ii,jj
  real*8                       :: cosphi,sinphi
  character(len=id_string_len) :: id_string 
  do jj=1,n_phi
    call sincosperiod_moivre(phi0(jj),HZ,dHZ)
    call assert_equals(1.d0,HZ(1),tol,&
    'Error mode sincosperiod moivre: cos n=0 mismatch!')
    call assert_equals(0.d0,dHZ(1),tol,&
    'Error mode sincosperiod moivre: dcos n=0 mismatch!')
    if(size(HZ)>1) then
      do ii=1,(n_tor-1)/2
        write(id_string,'(I3)') ii
        cosphi = cos(real(n_period*ii,kind=8)*phi0(jj))
        sinphi = sin(real(n_period*ii,kind=8)*phi0(jj))
        call assert_equals(cosphi,HZ(2*ii),tol,&
        'Error mode sincosperiod moivre: cos mismatch id: '//id_string//'!')
        call assert_equals(sinphi,HZ(2*ii+1),tol,&
        'Error mode sincosperiod moivre: sin mismatch id: '//id_string//'!')
        call assert_equals(n_period*sinphi,dHZ(2*ii),tol,&
        'Error mode sincosperiod moivre: dcos mismatch id: '//id_string//'!')
        call assert_equals(-n_period*cosphi,dHZ(2*ii+1),tol,&
        'Error mode sincosperiod moivre: dsin mismatch id: '//id_string//'!')
      enddo
    endif
  enddo
end subroutine test_sincosperiod_moivre
#endif
!> ------------------------------------------------
end module mod_moivre_series_test
