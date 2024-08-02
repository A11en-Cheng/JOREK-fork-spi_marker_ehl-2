module mod_SolveMN
implicit none
contains
subroutine SolveMN(A, info)
!-------------------------------------------------------------------------------
! Solves a matrix system by invering an NxN matrix using the Gauss-Jordan method
! Useful only for small matrices, where calling PASTIX makes no sense
!-------------------------------------------------------------------------------
  implicit none
  real*8, dimension(:,:), allocatable, intent(inout) :: A
  integer,                             intent(out)   :: info
  
  real*8, dimension(:,:), allocatable :: Ainv
  real*8, dimension(:),   allocatable :: tmp
  integer                             :: N, i, j, ml
  
  N = size(A,1)
  allocate(Ainv(N,N)); allocate(tmp(N))
  
  if (size(A,2) .ne. N) then
    info = 3
    return
  end if
  
  Ainv = 0.d0
  do i=1,N
    Ainv(i,i) = 1.d0
  end do
  
  do j=1,N
    ml = maxloc(abs(A(j:,j)),1) + j - 1
    if (A(ml,j) .eq. 0.d0) then
      info = 1
      return
    end if
    if (ml .ne. j) then
      tmp = Ainv(j,:)
      Ainv(j,:) = Ainv(ml,:)
      Ainv(ml,:) = tmp
      
      tmp = A(j,:)
      A(j,:) = A(ml,:)
      A(ml,:) = tmp
    end if
    
    Ainv(j,:) = Ainv(j,:)/A(j,j)
    A(j,:) = A(j,:)/A(j,j)
    
    do i=1,N
      if (i .ne. j) then
        Ainv(i,:) = Ainv(i,:) - A(i,j)*Ainv(j,:)
        A(i,:) = A(i,:) - A(i,j)*A(j,:)
      end if
    end do
  end do
  
  info = 0
  A = Ainv
end subroutine SolveMN
end module mod_SolveMN
