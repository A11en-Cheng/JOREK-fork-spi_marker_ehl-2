!> This module is used in some models for automatic code generation of the RHS and AMAT in the 
!! matrix construction routines. It contains algebraic types and operators that can be used to
!! construct the matrix contributions in a compact, legible format that is then converted into 
!! code
module mod_semianalytical
  implicit none

  ! A simple algebraic expression with two operands and an arithmetical operator (+, -, * or \)
  ! The operands can either be basics variables or algebraic expressions themselves
  type algexpr
    type(algexpr), pointer    :: operand1 => NULL(), operand2 => NULL()
    character                 :: oprtr = ''
    integer                   :: dx = 0, dy = 0, dp = 0    ! Orders of differential operators acting on this expression
    logical                   :: basic = .false.           ! .true. if there are no more sub-operands in the expression
    integer                   :: var = 0                   ! The index of the basic variable. Only initialized if basic .eq. .true.
    real*8                    :: factor = 1.d0, add = 0.d0 ! A numerical multiplicative factor and an additive constant
    character(:), allocatable :: factcode, addcode         ! Fortran code for the multiplicative and additive constants
    
  end type algexpr
  
  ! A constant derived from a Fortran variable that does not vary in space or time
  type const
    real*8                    :: value ! The numerical value of the constant
    character(:), allocatable :: token ! The name of the corresponding Fortran variable
  end type const

  
  type(algexpr), parameter :: one = algexpr(basic = .true.,var=1,factor=0,add=1)

  ! Operators for making algebraic expressions
  interface operator (+)
    procedure addexpr  ! algexpr + algexpr
    procedure addexprn ! algexpr + real
    procedure addnexpr ! real + algexpr
    procedure addexprc ! algexpr + const
    procedure addcexpr ! const + algexpr
    procedure addcc    ! const + const
    procedure addcn    ! const + real
    procedure addnc    ! real + const
  end interface operator (+)
  
  interface operator (-)
    procedure subexpr  ! algexpr - algexpr
    procedure subexprn ! algexpr - real
    procedure subnexpr ! real - algexpr
    procedure subexprc ! algexpr - const
    procedure subcexpr ! const - algexpr
    procedure subcc    ! const - const
    procedure subcn    ! const - real
    procedure subnc    ! real - const
    procedure negate   ! -algexpr
    procedure negatec  ! -const
  end interface operator (-)
  
  interface operator (*)
    procedure multexpr  ! algexpr*algexpr
    procedure multexprn ! algexpr*real
    procedure multnexpr ! real*algexpr
    procedure multexprc ! algexpr*const
    procedure multcexpr ! const*algexpr
    procedure multcc    ! const*const
    procedure multcn    ! const*real
    procedure multnc    ! real*const
  end interface operator (*)
  
  interface operator (/)
    procedure divexpr  ! algexpr/algexpr
    procedure divexprn ! algexpr/real
    procedure divnexpr ! real/algexpr
    procedure divexprc ! algexpr/const
    procedure divcexpr ! const/algexpr
    procedure divcc    ! const/const
    procedure divcn    ! const/real
    procedure divnc    ! real/const
  end interface operator (/)
  
  interface operator (**)
    procedure powexprn ! algexpr**real
    procedure powcn    ! const**real
  end interface operator (**)
  
contains

  ! Functions implementing the operators
  type(algexpr) function addexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(addexpr%operand1)
    allocate(addexpr%operand2)
    
    addexpr%operand1 = e1
    addexpr%operand2 = e2
    addexpr%oprtr    = '+'
  end function addexpr
  
  type(algexpr) function addexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    character(23)             :: num
    
    addexprn = e1
    addexprn%add = addexprn%add + n2
    write(num,'(E14.6)') n2
    if (allocated(addexprn%addcode)) then
      addexprn%addcode = addexprn%addcode // " + " // trim(adjustl(num))
    else
      addexprn%addcode = trim(adjustl(num))
    end if
  end function addexprn
  
  type(algexpr) function addnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    character(23)             :: num
    
    addnexpr = e2
    addnexpr%add = addnexpr%add + n1
    write(num,'(E14.6)') n1
    if (allocated(addnexpr%addcode)) then
      addnexpr%addcode = addnexpr%addcode // " + " // trim(adjustl(num))
    else
      addnexpr%addcode = trim(adjustl(num))
    end if
  end function addnexpr
  
  type(algexpr) function addexprc(e1,c2)
    implicit none
    type(algexpr), intent(in) :: e1
    type(const),   intent(in) :: c2
    
    addexprc = e1
    addexprc%add = addexprc%add + c2%value
    if (allocated(addexprc%addcode)) then
      addexprc%addcode = addexprc%addcode // " + " // c2%token
    else
      addexprc%addcode = c2%token
    end if
  end function addexprc
  
  type(algexpr) function addcexpr(c1,e2)
    implicit none
    type(const),   intent(in) :: c1
    type(algexpr), intent(in) :: e2
    
    addcexpr = e2
    addcexpr%add = addcexpr%add + c1%value
    if (allocated(addcexpr%addcode)) then
      addcexpr%addcode = addcexpr%addcode // " + " // c1%token
    else
      addcexpr%addcode = c1%token
    end if
  end function addcexpr
  
  type(const) function addcc(c1,c2)
    implicit none
    type(const), intent(in) :: c1, c2
    
    addcc%value = c1%value + c2%value
    addcc%token = c1%token // " + " // c2%token
  end function addcc
  
  type(const) function addcn(c1,n2)
    implicit none
    type(const), intent(in) :: c1
    real*8,      intent(in) :: n2
    character(23)           :: num
    
    addcn%value = c1%value + n2
    write(num,'(E14.6)') n2
    addcn%token = c1%token // " + " // trim(adjustl(num))
  end function addcn
  
  type(const) function addnc(n1,c2)
    implicit none
    real*8,      intent(in) :: n1
    type(const), intent(in) :: c2
    character(23)           :: num
    
    addnc%value = n1 + c2%value
    write(num,'(E14.6)') n1
    addnc%token = trim(adjustl(num)) // " + " // c2%token
  end function addnc
  
  type(algexpr) function subexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(subexpr%operand1)
    allocate(subexpr%operand2)
    
    subexpr%operand1 = e1
    subexpr%operand2 = e2
    subexpr%oprtr    = '-'
  end function subexpr
  
  type(algexpr) function subexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    character(23)             :: num
    
    subexprn = e1
    subexprn%add = subexprn%add - n2
    write(num,'(E14.6)') n2
    if (allocated(subexprn%addcode)) then
      subexprn%addcode = subexprn%addcode // " - " // trim(adjustl(num))
    else
      subexprn%addcode = "(-" // trim(adjustl(num)) // ")"
    end if
  end function subexprn
  
  type(algexpr) function subnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    character(23)             :: num
    
    subnexpr = e2
    subnexpr%factor = -subnexpr%factor
    subnexpr%add = n1 - subnexpr%add
    write(num,'(E14.6)') n1
    if (allocated(subnexpr%factcode)) then
      subnexpr%factcode = "(-" // subnexpr%factcode // ")"
    else
      subnexpr%factcode = "(-1)"
    end if
    if (allocated(subnexpr%addcode)) then
      subnexpr%addcode = trim(adjustl(num)) // " - (" // subnexpr%addcode // ")"
    else
      subnexpr%addcode = trim(adjustl(num))
    end if
  end function subnexpr
  
  type(algexpr) function subexprc(e1,c2)
    implicit none
    type(algexpr), intent(in) :: e1
    type(const),   intent(in) :: c2
    
    subexprc = e1
    subexprc%add = subexprc%add - c2%value
    if (allocated(subexprc%addcode)) then
      subexprc%addcode = subexprc%addcode // " - (" // c2%token // ")"
    else
      subexprc%addcode = "(-(" // c2%token // "))"
    end if
  end function subexprc
  
  type(algexpr) function subcexpr(c1,e2)
    implicit none
    type(const),   intent(in) :: c1
    type(algexpr), intent(in) :: e2
    
    subcexpr = e2
    subcexpr%factor = -subcexpr%factor
    subcexpr%add = c1%value - subcexpr%add
    if (allocated(subcexpr%factcode)) then
      subcexpr%factcode = "(-" // subcexpr%factcode // ")"
    else
      subcexpr%factcode = "(-1)"
    end if
    if (allocated(subcexpr%addcode)) then
      subcexpr%addcode = c1%token // " - (" // subcexpr%addcode // ")"
    else
      subcexpr%addcode = c1%token
    end if
  end function subcexpr
  
  type(const) function subcc(c1,c2)
    implicit none
    type(const), intent(in) :: c1, c2
    
    subcc%value = c1%value - c2%value
    subcc%token = c1%token // " - (" // c2%token // ")"
  end function subcc
  
  type(const) function subcn(c1,n2)
    implicit none
    type(const), intent(in) :: c1
    real*8,      intent(in) :: n2
    character(23)           :: num
    
    subcn%value = c1%value - n2
    write(num,'(E14.6)') n2
    subcn%token = c1%token // " - " // trim(adjustl(num))
  end function subcn
  
  type(const) function subnc(n1,c2)
    implicit none
    real*8,      intent(in) :: n1
    type(const), intent(in) :: c2
    character(23)           :: num
    
    subnc%value = n1 - c2%value
    write(num,'(E14.6)') n1
    subnc%token = trim(adjustl(num)) // " - (" // c2%token // ")"
  end function subnc
  
  type(algexpr) function negate(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    
    negate = expr
    negate%factor = -negate%factor
    negate%add = -negate%add
    if (allocated(negate%factcode)) then
      negate%factcode = "(-" // negate%factcode // ")"
    else
      negate%factcode = "(-1)"
    end if
    if (allocated(negate%addcode)) negate%addcode = "-(" // negate%addcode // ")"
  end function negate
  
  type(const) function negatec(c)
    implicit none
    type(const), intent(in) :: c
    
    negatec%value = -c%value
    negatec%token = "(-(" // c%token // "))"
  end function negatec
  
  type(algexpr) function multexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(multexpr%operand1)
    allocate(multexpr%operand2)
    
    multexpr%operand1 = e1
    multexpr%operand2 = e2
    multexpr%oprtr    = '*'
  end function multexpr
  
  type(algexpr) function multexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    character(23)             :: num
    
    multexprn = e1
    multexprn%factor = n2*multexprn%factor
    multexprn%add = n2*multexprn%add
    write(num,'(E14.6)') n2
    if (allocated(multexprn%factcode)) then
      multexprn%factcode = trim(adjustl(num)) // "*" // multexprn%factcode
    else
      multexprn%factcode = trim(adjustl(num))
    end if
    if (allocated(multexprn%addcode)) multexprn%addcode = trim(adjustl(num)) // "*(" // multexprn%addcode // ")"
  end function multexprn
  
  type(algexpr) function multnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    character(23)             :: num
    
    multnexpr = e2
    multnexpr%factor = n1*multnexpr%factor
    multnexpr%add = n1*multnexpr%add
    write(num,'(E14.6)') n1
    if (allocated(multnexpr%factcode)) then
      multnexpr%factcode = trim(adjustl(num)) // "*" // multnexpr%factcode
    else
      multnexpr%factcode = trim(adjustl(num))
    end if
    if (allocated(multnexpr%addcode)) multnexpr%addcode = trim(adjustl(num)) // "*(" // multnexpr%addcode // ")"
  end function multnexpr
  
  type(algexpr) function multexprc(e1,c2)
    implicit none
    type(algexpr), intent(in) :: e1
    type(const),   intent(in) :: c2
    
    multexprc = e1
    multexprc%factor = c2%value*multexprc%factor
    multexprc%add = c2%value*multexprc%add
    if (allocated(multexprc%factcode)) then
      multexprc%factcode = "(" // c2%token // ")*" // multexprc%factcode
    else
      multexprc%factcode = "(" // c2%token // ")"
    end if
    if (allocated(multexprc%addcode)) multexprc%addcode = "(" // c2%token // ")*(" // multexprc%addcode // ")"
  end function multexprc
  
  type(algexpr) function multcexpr(c1,e2)
    implicit none
    type(const),   intent(in) :: c1
    type(algexpr), intent(in) :: e2
    
    multcexpr = e2
    multcexpr%factor = c1%value*multcexpr%factor
    multcexpr%add = c1%value*multcexpr%add
    if (allocated(multcexpr%factcode)) then
      multcexpr%factcode = "(" // c1%token // ")*" // multcexpr%factcode
    else
      multcexpr%factcode = "(" // c1%token // ")"
    end if
    if (allocated(multcexpr%addcode)) multcexpr%addcode = "(" // c1%token // ")*(" // multcexpr%addcode // ")"
  end function multcexpr
  
  type(const) function multcc(c1,c2)
    implicit none
    type(const), intent(in) :: c1, c2
    
    multcc%value = c1%value*c2%value
    multcc%token = "(" // c1%token // ")*(" // c2%token // ")"
  end function multcc
  
  type(const) function multcn(c1,n2)
    implicit none
    type(const), intent(in) :: c1
    real*8,      intent(in) :: n2
    character(23)           :: num
    
    multcn%value = c1%value*n2
    write(num,'(E14.6)') n2
    multcn%token = "(" // c1%token // ")*" // trim(adjustl(num))
  end function multcn
  
  type(const) function multnc(n1,c2)
    implicit none
    real*8,      intent(in) :: n1
    type(const), intent(in) :: c2
    character(23)           :: num
    
    multnc%value = n1*c2%value
    write(num,'(E14.6)') n1
    multnc%token = trim(adjustl(num)) // "*(" // c2%token // ")"
  end function multnc
  
  type(algexpr) function divexpr(e1,e2)
    implicit none
    type(algexpr), intent(in) :: e1, e2
    
    allocate(divexpr%operand1)
    allocate(divexpr%operand2)
    
    divexpr%operand1 = e1
    divexpr%operand2 = e2
    divexpr%oprtr    = '/'
  end function divexpr
  
  type(algexpr) function divexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    real*8,        intent(in) :: n2
    character(23)             :: num
    
    divexprn = e1
    divexprn%factor = divexprn%factor/n2
    divexprn%add = divexprn%add/n2
    write(num,'(E14.6)') n2
    if (allocated(divexprn%factcode)) then
      divexprn%factcode = "(" // divexprn%factcode // "/" // trim(adjustl(num)) // ")"
    else
      divexprn%factcode = "(1/" // trim(adjustl(num)) // ")"
    end if
    if (allocated(divexprn%addcode)) divexprn%addcode = "((" // divexprn%addcode // ")/" // trim(adjustl(num)) // ")"
  end function divexprn
  
  type(algexpr) function divnexpr(n1,e2)
    implicit none
    real*8,        intent(in) :: n1
    type(algexpr), intent(in) :: e2
    character(23)             :: num
    
    allocate(divnexpr%operand1)
    allocate(divnexpr%operand2)
    
    divnexpr%operand1     = one
    divnexpr%operand1%add = n1
    divnexpr%operand2     = e2
    divnexpr%oprtr        = '/'
    write(num,'(E14.6)') n1
    divnexpr%operand1%addcode = trim(adjustl(num))
  end function divnexpr
  
  type(algexpr) function divexprc(e1,c2)
    implicit none
    type(algexpr), intent(in) :: e1
    type(const),   intent(in) :: c2
    
    divexprc = e1
    divexprc%factor = divexprc%factor/c2%value
    divexprc%add = divexprc%add/c2%value
    if (allocated(divexprc%factcode)) then
      divexprc%factcode = "(" // divexprc%factcode // "/(" // c2%token // "))"
    else
      divexprc%factcode = "(1/(" // c2%token // "))"
    end if
    if (allocated(divexprc%addcode)) divexprc%addcode = "((" // divexprc%addcode // ")/(" // c2%token // "))"
  end function divexprc
  
  type(algexpr) function divcexpr(c1,e2)
    implicit none
    type(const),   intent(in) :: c1
    type(algexpr), intent(in) :: e2
    
    allocate(divcexpr%operand1)
    allocate(divcexpr%operand2)
    
    divcexpr%operand1        = one
    divcexpr%operand1%factor = c1%value
    divcexpr%operand2        = e2
    divcexpr%oprtr           = '/'
    divcexpr%operand1%addcode = c1%token
  end function divcexpr
  
  type(const) function divcc(c1,c2)
    implicit none
    type(const), intent(in) :: c1, c2
    
    divcc%value = c1%value/c2%value
    divcc%token = "((" // c1%token // ")/(" // c2%token // "))"
  end function divcc
  
  type(const) function divcn(c1,n2)
    implicit none
    type(const), intent(in) :: c1
    real*8,      intent(in) :: n2
    character(23)           :: num
    
    divcn%value = c1%value/n2
    write(num,'(E14.6)') n2
    divcn%token = "((" // c1%token // ")/" // trim(adjustl(num)) // ")"
  end function divcn
  
  type(const) function divnc(n1,c2)
    implicit none
    real*8,      intent(in) :: n1
    type(const), intent(in) :: c2
    character(23)           :: num
    
    divnc%value = n1/c2%value
    write(num,'(E14.6)') n1
    divnc%token = "(" // trim(adjustl(num)) // "/(" // c2%token // "))"
  end function divnc
  
  type(algexpr) function powexprn(e1,n2)
    implicit none
    type(algexpr), intent(in) :: e1
    integer,       intent(in) :: n2
    integer                   :: i
    
    if (n2 .eq. 0) then
      powexprn = one
    else
      powexprn = e1
      do i=2,abs(n2)
        powexprn = powexprn*e1
      end do
      if (n2 .lt. 0) powexprn = 1.0/powexprn
    end if
  end function powexprn
  
  type(const) function powcn(c1,n2)
    implicit none
    type(const), intent(in) :: c1
    real*8,      intent(in) :: n2
    character(23)           :: num
    
    powcn%value = c1%value**n2
    write(num,'(E14.6)') n2
    powcn%token = "((" // c1%token // ")**(" // trim(adjustl(num)) // "))"
  end function powcn
  
  ! Functions applying differential operators to the expression
  type(algexpr) function dx(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    dx    = expr
    dx%dx = dx%dx + 1
  end function dx
  
  type(algexpr) function dy(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    dy    = expr
    dy%dy = dy%dy + 1
  end function dy
  
  type(algexpr) function dp(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    dp    = expr
    dp%dp = dp%dp + 1
  end function dp
  
  ! This function makes a deep copy of the entire expression tree
  type(algexpr) recursive function deepcopy(expr) result(res)
    implicit none
    type(algexpr), intent(in) :: expr
    
    if (expr%basic) then
      res = expr
    else
      res%oprtr  = expr%oprtr
      res%dx     = expr%dx
      res%dy     = expr%dy
      res%dp     = expr%dp
      res%var    = expr%var
      res%factor = expr%factor
      res%add    = expr%add
      
      if (allocated(expr%factcode)) res%factcode = expr%factcode
      if (allocated(expr%addcode)) res%addcode  = expr%addcode
      
      allocate(res%operand1)
      allocate(res%operand2)
      res%operand1 = deepcopy(expr%operand1)
      res%operand2 = deepcopy(expr%operand2)
    end if
  end function deepcopy
  

  type(algexpr) recursive function Dexpand(expr) result(res)
    implicit none
    type(algexpr), target, intent(in) :: expr
    type(algexpr), pointer            :: oldop1, oldop2, rescp
    
    if (.not. expr%basic) then
      res = deepcopy(expr)
      if (expr%dx .ne. 0) then
        select case (expr%oprtr)
          case ('+')
            res = dx(res%operand1) + dx(res%operand2)
          case ('-')
            res = dx(res%operand1) - dx(res%operand2)
          case ('*')
            res = dx(res%operand1)*res%operand2 + res%operand1*dx(res%operand2)
          case ('/')
            res = (res%operand2*dx(res%operand1) - res%operand1*dx(res%operand2))/(res%operand2*res%operand2)
        end select
        res%dx = expr%dx - 1
        res%dy = expr%dy
        res%dp = expr%dp
        res%add = 0
      else if (expr%dy .ne. 0) then
        select case (expr%oprtr)
          case ('+')
            res = dy(res%operand1) + dy(res%operand2)
          case ('-')
            res = dy(res%operand1) - dy(res%operand2)
          case ('*')
            res = dy(res%operand1)*res%operand2 + res%operand1*dy(res%operand2)
          case ('/')
            res = (res%operand2*dy(res%operand1) - res%operand1*dy(res%operand2))/(res%operand2*res%operand2)
        end select
        res%dx = expr%dx
        res%dy = expr%dy - 1
        res%dp = expr%dp
        res%add = 0
      else if (expr%dp .ne. 0) then
        select case (expr%oprtr)
          case ('+')
            res = dp(res%operand1) + dp(res%operand2)
          case ('-')
            res = dp(res%operand1) - dp(res%operand2)
          case ('*')
            res = dp(res%operand1)*res%operand2 + res%operand1*dp(res%operand2)
          case ('/')
            res = (res%operand2*dp(res%operand1) - res%operand1*dp(res%operand2))/(res%operand2*res%operand2)
        end select
        res%dx = expr%dx
        res%dy = expr%dy
        res%dp = expr%dp - 1
        res%add = 0
      end if
      res%factor = expr%factor
      if (allocated(expr%factcode)) res%factcode = expr%factcode
      
      if (res%dx .ne. 0 .or. res%dy .ne. 0 .or. res%dp .ne. 0) then
        allocate(rescp)
        rescp = res
        res = Dexpand(rescp)
      else
        oldop1 => res%operand1
        oldop2 => res%operand2
        allocate(res%operand1)
        allocate(res%operand2)
        res%operand1 = Dexpand(oldop1)
        res%operand2 = Dexpand(oldop2)
      end if
    else
      res = expr
      if (res%dx .ne. 0 .or. res%dy .ne. 0 .or. res%dp .ne. 0) res%add = 0
    end if
  end function Dexpand
  
  ! Print an expression using basic text
  recursive subroutine printexpr(expr)
    implicit none
    type(algexpr), intent(in) :: expr
    
    if (expr%factor .ne. 1.d0) write(*,"(F10.3,A)",advance='no') expr%factor,"*"
    
    if (expr%dx .ne. 0) write(*,"(A,I1)",advance='no') "dx",expr%dx
    if (expr%dy .ne. 0) write(*,"(A,I1)",advance='no') "dy",expr%dy
    if (expr%dp .ne. 0) write(*,"(A,I1)",advance='no') "dp",expr%dp
    
    if (expr%basic) then
      write(*,"(A,I1,A)",advance='no') "[var",expr%var,"]"
    else
      write(*,"(A)",advance='no') "("
      call printexpr(expr%operand1)
      write(*,"(A)",advance='no') expr%oprtr
      call printexpr(expr%operand2)
      write(*,"(A)",advance='no') ")"
    end if
    
    if (expr%add .ne. 0.d0) write(*,"(A,F10.3)",advance='no') " + ",expr%add
  end subroutine printexpr
  
  ! Print LaTeX code for an expression
  ! varsymb: a string containing one-character symbols (in proper order) for each variable
  recursive subroutine printlatex(expr,varsymb)
    implicit none
    type(algexpr),    intent(in) :: expr
    character(len=*), intent(in) :: varsymb
    
    if (expr%factor .ne. 1.d0) write(*,"(F10.3,A)",advance='no') expr%factor,"*"
    
    if (expr%dx .eq. 1) then
      write(*,"(A)",advance='no') "\frac{\partial}{\partial x}"
    else if (expr%dx .gt. 1) then
      write(*,"(A,I1,A,I1,A)",advance='no') "\frac{\partial^",expr%dx,"}{\partial x^",expr%dx,"}"
    end if
    
    if (expr%dy .eq. 1) then
      write(*,"(A)",advance='no') "\frac{\partial}{\partial y}"
    else if (expr%dy .gt. 1) then
      write(*,"(A,I1,A,I1,A)",advance='no') "\frac{\partial^",expr%dy,"}{\partial y^",expr%dy,"}"
    end if
    
    if (expr%dp .eq. 1) then
      write(*,"(A)",advance='no') "\frac{\partial}{\partial p}"
    else if (expr%dp .gt. 1) then
      write(*,"(A,I1,A,I1,A)",advance='no') "\frac{\partial^",expr%dp,"}{\partial p^",expr%dp,"}"
    end if
    
    if (expr%basic) then
      write(*,"(A)",advance='no') varsymb(expr%var:expr%var)
    else
      if (expr%oprtr .eq. '/') then
        write(*,"(A)",advance='no') "\frac{"
        call printlatex(expr%operand1,varsymb)
        write(*,"(A)",advance='no') "}{"
        call printlatex(expr%operand2,varsymb)
        write(*,"(A)",advance='no') "}"
      else if (expr%oprtr .eq. '*') then
        write(*,"(A)",advance='no') "\left("
        call printlatex(expr%operand1,varsymb)
        call printlatex(expr%operand2,varsymb)
        write(*,"(A)",advance='no') "\right)"
      else
        write(*,"(A)",advance='no') "\left("
        call printlatex(expr%operand1,varsymb)
        write(*,"(A)",advance='no') expr%oprtr
        call printlatex(expr%operand2,varsymb)
        write(*,"(A)",advance='no') "\right)"
      end if
    end if
    
    if (expr%add .ne. 0.d0) write(*,"(A,F10.3)",advance='no') " + ",expr%add
  end subroutine printlatex

  integer recursive function countsubexprs(expr) result(res)
    implicit none
    type(algexpr), intent(in) :: expr
    
    if (expr%basic) then
      res = 0
    else
      res = 1 + countsubexprs(expr%operand1) + countsubexprs(expr%operand2)
    end if
  end function countsubexprs

  ! Generates Fortran code from an algexpr structure
  recursive function gencode(expr, varname) result(res)
    implicit none
    type(algexpr),             intent(in) :: expr
    character(:), allocatable, intent(in) :: varname
    character(:), allocatable             :: res
    character(12)                         :: indices
    
    if (.not. expr%basic .and. (expr%dx .ne. 0 .or. expr%dy .ne. 0 .or. expr%dp .ne. 0)) then
      write(*,*)
      write(*,*) ">>>>> Cannot generate code from unexpanded algexpr: ABORTING COMPILATION <<<<<"
      write(*,*)
      stop
    end if
    
    if (expr%basic) then
#if JOREK_MODEL == 180
      write(indices,'(A,I2,A,I1,A,I1,A,I1,A)') "(", expr%var, ",", expr%dx, ",", expr%dy, ",", expr%dp, ",1)"
#else
      write(indices,'(A,I2,A,I1,A,I1,A,I1,A)') "(", expr%var, ",", expr%dx, ",", expr%dy, ",", expr%dp, ",:)"
#endif
      res = varname // indices
    else
      res = "(" // trim(gencode(expr%operand1,varname)) // expr%oprtr // trim(gencode(expr%operand2,varname)) // ")"
    end if
    
    if (allocated(expr%factcode)) res = expr%factcode // "*" // res
    if (allocated(expr%addcode)) res = res // " + " // expr%addcode
  end function gencode
end module mod_semianalytical
