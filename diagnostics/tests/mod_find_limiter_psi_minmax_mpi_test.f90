module mod_find_limiter_psi_minmax_mpi_test
use fruit
use fruit_mpi
use data_structure
implicit none
private
public :: run_fruit_find_limiter_psi_minmax_mpi
integer,parameter           :: n_tests=1
integer,parameter           :: nx=2
integer,parameter           :: ny=2
integer,parameter           :: i_harm_test=1
integer,parameter           :: i_elm_lim_expect=1
real*8,parameter            :: R_begin=5d-1
real*8,parameter            :: R_end=1.5d0
real*8,parameter            :: Z_begin=-5d-1
real*8,parameter            :: Z_end=5d-1
real*8,parameter            :: FF_0_test=1.d0
real*8,parameter            :: tol=1d-12
real*8,dimension(n_tests)   :: Rlim_expect=(/5d0/6d0/)
real*8,dimension(n_tests)   :: Zlim_expect=(/-5d-1/)
real*8,dimension(n_tests)   :: psimin_expect=(/-4d0/27d0/)
real*8,dimension(n_tests)   :: psimax_expect=(/0d0/)
type(type_node_list)        :: test_nodes
type(type_element_list)     :: test_elements
type(type_bnd_node_list)    :: test_bnd_nodes
type(type_bnd_element_list) :: test_bnd_elements
integer :: rank_loc,n_tasks_loc,ifail_loc
real*8  :: FF_0_loc
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_find_limiter_psi_minmax_mpi(&
rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: find limiter psi minmax mpi"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: find limiter psi minmax mpi"
  call run_test_case(test_find_limiter_psiminamx,'test_find_limiter_psiminamx')
  write(*,'(/A)') "  ... tearing-down: find limiter psi minmax mpi"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_find_limiter_psi_minmax_mpi

subroutine setup(rank,n_tasks,ifail)
  use phys_module, only: FF_0
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=rank; n_tasks_loc=n_tasks_loc; ifail_loc=ifail;
  FF_0_loc=FF_0; FF_0=FF_0_test;
  call init_node_list(test_nodes, n_nodes_max, test_nodes%n_dof, n_var)
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  use phys_module, only: FF_0
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc=-1; n_tasks_loc=-1; ifail=ifail_loc;
  FF_0 = FF_0_loc
  call dealloc_node_list(test_nodes)
end subroutine teardown

!> Tests ------------------------------------------
!> Test if we can find the limiter on a square grid with a
!> linear variation of poloidal flux. The grid consists only
!> of a single element
subroutine test_find_limiter_psiminamx()
  use mod_model_settings, only: var_psi
  use mod_boundary,       only: boundary_from_grid
  use phys_module,        only: FF_0
  use mod_interp,         only: interp
  use mod_projection_helpers_test_tools, only: default_square_grid
  implicit none
  integer :: ii,i_elm_lim
  real*8  :: psi_lim,R_lim,Z_lim,min_psi,max_psi
  real*8  :: R_out,Z_out,s_lim,t_lim
  real*8  :: P,P_s,P_t,P_st,P_ss,P_tt
  !> create grid
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,test_nodes,&
  test_elements,ifail_loc,R_begin,R_end,Z_begin,Z_end,&
  bnd_node_list_out=test_bnd_nodes,bnd_element_list_out=test_bnd_elements) 
  !> Fill in some analytical values, only for the values
  !> and derivatives. The corss-derivatives are not taken into
  !> account by find_limiter and psi_minmax
  call reset_node_values(test_nodes)
  call find_limiter(rank_loc,test_nodes,test_elements,test_bnd_elements,&
  psi_lim,R_lim,Z_lim)
  !> interpolat flux at limiter
  call find_RZ(test_nodes,test_elements,R_lim,Z_lim,R_out,Z_out,&
  i_elm_lim,s_lim,t_lim,ifail_loc)
  call interp(test_nodes,test_elements,i_elm_lim,var_psi,i_harm_test,&
  s_lim,t_lim,P,P_s,P_t,P_st,P_ss,P_tt)
  !> find maximum and minimum poloidal flux
  call psi_minmax(test_nodes,test_elements,i_elm_lim_expect,min_psi,max_psi)
  !> checks
  call assert_equals(Rlim_expect(1),R_lim,tol,&
  'Error find limiter psiminmax: wrong limiter R!')
  call assert_equals(Zlim_expect(1),Z_lim,tol,&
  'Error find limiter psiminmax: wrong limiter Z!')
  call assert_equals(psimin_expect(1),psi_lim,tol,&
  'Error find limiter psiminmax: wrong limiter poloidal flux!')
  call assert_equals(P,psi_lim,tol,&
  'Error find limiter psiminamx: psi at R,Z limiter and psi limiter mismatch!')
  call assert_equals(psimin_expect(1),min_psi,tol,&
  'Error find limiter psiminamx: woring minimum poloidal flux!')
  call assert_equals(psimax_expect(1),max_psi,tol,&
  'Error find limiter psiminamx: woring maximum poloidal flux!')
end subroutine test_find_limiter_psiminamx

!> Tools ------------------------------------------
subroutine reset_node_values(node_list)
  implicit none
  type(type_node_list),intent(inout) :: node_list
  integer :: ii
  do ii=1,node_list%n_nodes
    node_list%node(ii)%values(1,:,1) = 0.d0
  enddo
  node_list%node(1)%values(1,2,1) = -1.d0
end subroutine reset_node_values
!> ------------------------------------------------
end module mod_find_limiter_psi_minmax_mpi_test
