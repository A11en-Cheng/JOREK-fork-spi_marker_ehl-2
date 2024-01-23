module mod_basisfunctions2_test
use fruit
implicit none
private
public :: run_fruit_basisfunctions2
!> Variables --------------------------------------
integer,parameter :: message_len=200
real*8,parameter  :: tol=1d-14
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_basisfunctions2
  implicit none
  write(*,'(/A)') "  ... setting-up: basisfunctions2"
  write(*,'(/A)') "  ... running: basisfunctions2"
  call run_test_case(test_basisfunctions_2D_1_basisfunctions_2D_2,&
  'test_basisfunctions_2D_1_basisfunctions_2D_2')
  call run_test_case(test_basisfunctions_2D_vs_transpose,&
  'test_basisfunctions_2D_vs_transpose')
  call run_test_case(test_basisfunction_properties,&
  'test_basisfunction_properties')
  write(*,'(/A)') "  ... tearing-down: basisfunctions2"
end subroutine run_fruit_basisfunctions2

!> Tests ------------------------------------------

!> Check that the first-derivative and second-derivative versions of the basisfunctions match at the gaussian points
subroutine test_basisfunctions_2D_1_basisfunctions_2D_2
  use gauss
  use mod_basisfunctions
  use mod_parameters, only: n_vertex_max, n_degrees

  integer :: k, l, m, n
  real*8  :: s, t
  real*8, dimension(n_vertex_max,n_degrees) :: H, H_s, H_t, H_st, H_ss, H_tt, G, G_s, G_t, G_st, F, F_s, F_t

  ! Verify against other expression
  do k=1,n_gauss
    s = xgauss(k)
    do l=1,n_gauss
      t = xgauss(l)
      call basisfunctions(s,t,H, H_s, H_t, H_st, H_ss, H_tt)
      call basisfunctions(s,t,G, G_s, G_t, G_st)
      call basisfunctions(s,t,F, F_s, F_t)

      do m=1,n_vertex_max
        do n=1,n_degrees
          call assert_equals(H(m,n), G(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H,G values mismatch!')
          call assert_equals(H(m,n), F(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H,F values mismatch!')
          call assert_equals(H_s(m,n), G_s(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H_s,G_s values mismatch!')
          call assert_equals(H_s(m,n), F_s(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H_s,F_s values mismatch!')
          call assert_equals(H_t(m,n), G_t(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H_t,G_t values mismatch!')
          call assert_equals(H_t(m,n), F_t(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H_t,F_t values misatch!')
          call assert_equals(H_st(m,n), G_st(m,n), tol, 'Error basisfunctions 2D-1 and 2D-2 test: H_st,G_st values mismatch!')
        enddo
      enddo
    enddo
  enddo
end subroutine test_basisfunctions_2D_1_basisfunctions_2D_2

!> Check that the transposed version of the basisfunctions is equal to the transpose of the normal versino
subroutine test_basisfunctions_2D_vs_transpose
  use gauss
  use mod_basisfunctions
  use mod_parameters, only: n_vertex_max, n_degrees

  integer :: k, l, m, n
  real*8  :: s, t
  real*8, dimension(n_vertex_max,n_degrees) :: H, H_s, H_t, G, G_s, G_t

  do k=1,n_gauss
    s = xgauss(k)
    do l=1,n_gauss
      t = xgauss(l)
      call basisfunctions(s,t,H, H_s, H_t)
      call basisfunctions_T(s,t,G, G_s, G_t)

      do m=1,n_vertex_max
        do n=1,n_degrees
          call assert_equals(H(m,n), G(n,m), tol, &
          'Error basisfunctions 2D vs transpose test: H,G  mismatch!')
          call assert_equals(H_s(m,n), G_s(n,m), tol, &
          'Error basisfunctions 2D vs transpose test: H_s,G_s  mismatch!')
          call assert_equals(H_t(m,n), G_t(n,m), tol, &
          'Error basisfunctions 2D vs transpose test: H_t,G_t  mismatch!')
        enddo
      enddo
    enddo
  enddo
end subroutine test_basisfunctions_2D_vs_transpose

!> Test if the basisfunctions are 1 on their node
subroutine test_basisfunction_properties
  use mod_basisfunctions
  use mod_parameters, only: n_vertex_max, n_degrees

  integer :: j, k
  real*8  :: s, t
  real*8, dimension(n_vertex_max,n_degrees) :: H
  character(len=message_len)                :: message

  do k=1,n_vertex_max

    select case (k)
    case (1)
      s = 0
      t = 0
    case (2)
      s = 1
      t = 0
    case (3)
      s = 1
      t = 1
    case (4)
      s = 0
      t = 1
    end select
    call basisfunctions(s,t,H)
    write(message,'(A,F0.3,A,F0.3)') 'Error basisfunction properties test for s: ',s,' t: ',t
    call assert_equals(1.d0,H(k,1), tol, &
    trim(adjustl(message))//': right node value not 1!')
    do j=0,2
      call assert_equals(0.d0,H(mod(k+j,4)+1,1), tol, &
      trim(adjustl(message))//': other node value not 0!')
    enddo
  enddo
end subroutine test_basisfunction_properties

!> ------------------------------------------------
end module mod_basisfunctions2_test
