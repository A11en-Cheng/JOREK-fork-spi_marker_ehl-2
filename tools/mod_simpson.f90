!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
!> A module for evaluating one-dimensional definite integrals with the adaptive Simpson rule
!<
module mod_simpson

  implicit none
 
  integer, private, parameter :: simpson_default_max_depth = 10


  !> Interface for passing double precision 1D functions as arguments.
  !> To use these, declare
  !>
  !> procedure(func_real8_1d), pointer :: my_function_pointer => null()
  !>
  !> and then set the my_function_pointer to point to the function
  !> by declaring
  !>
  !> my_function_pointer => my_function
  !>
  !<
  abstract interface
     function func_real8_1D (x)
       import
       real*8 :: func_real8_1d
       real*8, intent (in) :: x
     end function func_real8_1D
  end interface
  
  public :: simpson_adaptive, func_real8_1D

  private
  
contains
    
  !> Helper routine for "simpson_adaptive"
  recursive function simpson_helper(f, a, b, eps, S, fa, fb, fc, bottom) result (val)
    implicit none
    real*8, intent(in) :: a,b,eps,s,fa,fb,fc
    integer, intent(in) :: bottom
    real*8 :: val
    real*8 :: c, h, d, e,fd,fe, Sleft, Sright, S2
    procedure(func_real8_1d) :: f

    c  = ( a + b ) / 2
    h  = b - a
    d  = ( a + c ) / 2
    e  = ( c + b ) / 2
    fd = f(d) 
    fe = f(e);  
    Sleft  = ( h / 12 ) * ( fa + 4 * fd + fc )
    Sright = ( h / 12 ) * ( fc + 4 * fe + fb )
    S2     = Sleft + Sright
        
    if (bottom .lt. 0 .or. abs(S2 - S) .le. eps*abs(S)) then
       val = S2 + ( S2 - S ) / 15;      
    else
       val = simpson_helper(f, a, c, eps, Sleft,  fa, fc, fd, bottom-1) &
            &+ simpson_helper(f, c, b, eps, Sright, fc, fb, fe, bottom-1)  
    end if
  end function simpson_helper
  
  !> Adaptive Simpsons rule for integral
  !<
  !< This function uses recursive splitting of the interval
  !< until either maximum number of intervals is reached or
  !< given relative error tolerance is obtained
  function simpson_adaptive(f, a, b, eps, maxDepth) result (val)
    implicit none
    procedure(func_real8_1d) :: f          !< (1D) integrand function
    real*8, intent(in) :: a                !< lower limit for the interval
    real*8, intent(in) :: b                !< upper limit for the interval
    real*8, intent(in) :: eps              !< absolute tolerance
    integer, optional, intent(in) :: maxDepth !< maximum number of splits for the interval (default 10)
    
    real*8  :: val
    integer :: depth
    
    ! Helper variables
    real*8 :: c,h,fa,fb,fc,S
    
    if(.not. present(maxDepth)) then 
       depth = simpson_default_max_depth
    else
       depth = maxDepth
    end if

    c = ( a + b ) / 2
    h = b - a                                      
    fa = f(a) 
    fb = f(b) 
    fc = f(c)                                      
    S = ( h / 6 ) * ( fa + 4 * fc + fb );
    val = simpson_helper(f, a, b, eps, S, fa, fb, fc, depth);        
  end function simpson_adaptive

end module mod_simpson
