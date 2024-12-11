module mod_penning_case_jorek_spec_mpi_test
use fruit
use fruit_mpi
use data_structure
implicit none
private
public :: run_fruit_penning_case_jorek_spec_mpi
!> Variables --------------------------------------
integer,parameter             :: message_len=100
integer,parameter             :: nx=30
integer,parameter             :: ny=30
integer,parameter             :: nrad=30
integer,parameter             :: npol=32
integer,parameter             :: master_rank=0
real*8,parameter              :: time_sol=0.d0
real*8,parameter              :: tor_angle_sol=0.d0
real*8,parameter              :: central_mass_sol=2.d0
real*8,parameter              :: central_density_sol=1.d0
logical,parameter             :: apply_dirichlet_bnd=.false.
real*8,dimension(2),parameter :: st_sol=(/5.d-1,5.d-1/)
real*8,dimension(3),parameter :: tolsE=(/1.d-10,3.d-11,0.d0/)
real*8,dimension(3),parameter :: tolsB=(/2.d-11,2.d-11,2.d0/)
integer                       :: rank_loc,n_tasks_loc,ifail_loc
real*8                        :: central_mass_restore
real*8                        :: central_density_restore
contains
!> Fruit basket -----------------------------------
subroutine run_fruit_penning_case_jorek_spec_mpi(rank,n_tasks,ifail)
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  write(*,'(/A)') "  ... setting-up: penning case jorek spec mpi"
  call setup(rank,n_tasks,ifail)
  write(*,'(/A)') "  ... running: penning case jorek spec mpi"
  call run_test_case(test_penning_case_jorek_square_10_10,'test_penning_case_jorek_square_10_10')
  call run_test_case(test_penning_case_jorek_polar_30_32,'test_penning_case_jorek_polar_30_32')
  write(*,'(/A)') "  ... tearing-down: penning case jorek spec mpi"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_penning_case_jorek_spec_mpi
!> Set-ups and teard-downs ------------------------
subroutine setup(rank,n_tasks,ifail)
  use phys_module, only: central_mass,central_density
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  rank_loc = rank; n_tasks_loc = n_tasks; ifail_loc = ifail;
  central_mass_restore    = central_mass
  central_density_restore = central_density
  central_mass            = central_mass_sol
  central_density         = central_density_sol
end subroutine setup

subroutine teardown(rank,n_tasks,ifail)
  use phys_module, only: central_mass,central_density
  implicit none
  integer,intent(inout) :: ifail
  integer,intent(in)    :: rank,n_tasks
  ifail=ifail_loc; rank_loc=-1; n_tasks_loc=0;
  central_mass    = central_mass_restore
  central_density = central_density_restore
end subroutine teardown

!> Tests ------------------------------------------
subroutine test_penning_case_jorek_square_10_10()
  use mod_fields_linear,                 only: jorek_fields_interp_linear
  use mod_projection_helpers_test_tools, only: default_square_grid
  use mod_penning_case_jorek,            only: jorek_penning_fields
  implicit none
  type(jorek_fields_interp_linear) :: fields 
  character(len=message_len),parameter :: message='penning case jorek square 10 10 test'
  allocate(fields%node_list,fields%element_list)
  call init_node_list(fields%node_list, n_nodes_max, fields%node_list%n_dof, n_var)
  call default_square_grid(rank_loc,n_tasks_loc,nx,ny,fields%node_list,&
  fields%element_list,ifail_loc)
  call jorek_penning_fields(fields%node_list,fields%element_list,&
  apply_dirichlet_in=apply_dirichlet_bnd,rank_in=rank_loc,&
  master_in=master_rank,ifail_out=ifail_loc)
  call verify_solution(fields,trim(message))
  deallocate(fields%node_list,fields%element_list);
end subroutine test_penning_case_jorek_square_10_10

subroutine test_penning_case_jorek_polar_30_32()
  use mod_fields_linear,                 only: jorek_fields_interp_linear
  use mod_projection_helpers_test_tools, only: default_polar_grid
  use mod_penning_case_jorek,            only: jorek_penning_fields
  implicit none
  type(jorek_fields_interp_linear) :: fields
  character(len=message_len),parameter :: message='penning case jorek polar 30 32 test'
  allocate(fields%node_list,fields%element_list)
  call init_node_list(fields%node_list, n_nodes_max, fields%node_list%n_dof, n_var)
  call default_polar_grid(rank_loc,n_tasks_loc,npol,nrad,&
  fields%node_list,fields%element_list,ifail_loc)
  call jorek_penning_fields(fields%node_list,fields%element_list,&
  apply_dirichlet_in=apply_dirichlet_bnd,rank_in=rank_loc,&
  master_in=master_rank,ifail_out=ifail_loc)
  call verify_solution(fields,trim(message))
  deallocate(fields%node_list,fields%element_list);
end subroutine test_penning_case_jorek_polar_30_32

!> Tools ------------------------------------------
!> Given a fields type, verify the solution at time 0 
!> against the reference case
subroutine verify_solution(fields,message_in)
  use phys_module,       only: F0
  use mod_interp,        only: interp_RZ
  use mod_fields_linear, only: jorek_fields_interp_linear
  use mod_penning_case,  only: case_penning_cylindrical
  implicit none
  type(case_penning_cylindrical),parameter    :: ref = case_penning_cylindrical()
  type(jorek_fields_interp_linear),intent(in) :: fields
  character(len=*)             :: message_in
  integer                      :: i_elm
  real*8                       :: psi,U,R,Z
  real*8,dimension(3)          :: E,B,E_ref,B_ref
  character(len=2*message_len) :: message
  !> loop on the mesh elements
  do i_elm=1,fields%element_list%n_elements
    call interp_RZ(fields%node_list,fields%element_list,&
    i_elm,st_sol(1),st_sol(2),R,Z)
    E_ref = ref%E([R,Z,tor_angle_sol],time_sol)
    B_ref = ref%B([R,Z,tor_angle_sol],time_sol)
    call fields%calc_EBpsiU(time_sol,i_elm,st_sol,tor_angle_sol,E,B,psi,U)
    write(message,'(A,A,A)') 'Error ',trim(message_in),': unexpected E_R!'
    call assert_equals(E_ref(1),E(1),tolsE(1),message)
    write(message,'(A,A,A)') 'Error ',trim(message_in),': unexpected E_Z!'
    call assert_equals(E_ref(2),E(2),tolsE(2),message)
    write(message,'(A,A,A)') 'Error ',trim(message_in),': unexpected E_phi!'
    call assert_equals(E_ref(3),E(3),tolsE(3),message)
    write(message,'(A,A,A)') 'Error ',trim(message_in),': unexpected B_R!'
    call assert_equals(B_ref(1),B(1),tolsB(1),message)
    write(message,'(A,A,A)') 'Error ',trim(message_in),': unexpected B_Z!'
    call assert_equals(B_ref(2),B(2),tolsB(2),message)
    write(message,'(A,A,A)') 'Error ',trim(message_in),': unexpected B_phi!'
    call assert_equals(B_ref(3),B(3),F0*tolsB(3),message)
  enddo
end subroutine verify_solution
!> ------------------------------------------------
end module mod_penning_case_jorek_spec_mpi_test
