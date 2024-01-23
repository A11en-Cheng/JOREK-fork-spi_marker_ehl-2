module mod_gauss_test
use fruit
use gauss, only: n_gauss
implicit none
private
public :: run_fruit_gauss

!> Variables --------------------------------------
integer,parameter                      :: n_len=2
integer,parameter                      :: s_len=20
integer,parameter                      :: n_schemes=1
integer,dimension(n_schemes),parameter :: n_gauss_p=(/n_gauss/)
integer,dimension(n_schemes),parameter :: n_tests=(/14/)
real*8,dimension(n_schemes),parameter  :: tol_test=(/2d-3/)
real*8,dimension(n_schemes),parameter  :: tol_exact=(/n_gauss*1d-16/)
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_gauss
  implicit none
  write(*,'(/A)') "  ... setting-up: gauss "
  write(*,'(/A)') "  ... running: gauss "
  call run_test_case(test_gauss_p,'test_gauss_p')
  write(*,'(/A)') "  ... tearing-down: gauss "
end subroutine run_fruit_gauss

!> Tests ------------------------------------------
subroutine test_gauss_p
  use gauss, only: Xgauss,Wgauss
  implicit none
  integer :: ii,jj,kk
  real*8  :: a,tol
  character(len=n_len) :: n
  character(len=s_len) :: s
  do kk=1,n_schemes
    do jj=1,n_tests(kk)
      a = 0d0
      do ii=1,n_gauss_p(kk)
        a = a + (Xgauss(ii)**jj)*Wgauss(ii)
      enddo
      if(jj.lt.2*n_gauss) then
        tol =  tol_exact(kk)
      else
        tol = tol_test(kk)
      endif
      write(n,'(i2)') jj
      write(s,'(g20.12)') 1d0/real(jj+1,kind=8)
      call assert_equals(1d0/real(jj+1,kind=8),a,tol,&
      'Error gauss p test: integral x**'//trim(adjustl(n))//' must be'//s//'!')
    enddo
  enddo
end subroutine test_gauss_p
!> ------------------------------------------------
end module mod_gauss_test
