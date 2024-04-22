!> the mod_light_vertices_test module contains variables
!> and methods for testing all procedures contained in
!> mod_light_vertices.f90. Due to the fact that light
!> vertices is an abstract class, the synchrotron vertex
!> type is used instead.
module mod_light_vertices_test
use fruit
use mod_particle_sim, only: particle_sim
use mod_full_synchrotron_light_dist_vertices, only: full_synchrotron_light_dist
use mod_particle_types, only: particle_fieldline_id,particle_gc_id
use mod_particle_types, only: particle_gc_vpar_id,particle_gc_Qin_id
use mod_particle_types, only: particle_kinetic_id,particle_kinetic_leapfrog_id
use mod_particle_types, only: particle_kinetic_relativistic_id
use mod_particle_types, only: particle_gc_relativistic_id
implicit none

private
public :: run_fruit_light_vertices

!> Variables ---------------------------------------------------------
character(len=12),parameter               :: input_file='light_inputs'
integer,parameter                         :: read_unit=43
integer,parameter                         :: n_x=3
integer,parameter                         :: n_properties_sol=11
integer,parameter                         :: n_particle_types=2
integer,parameter                         :: n_times_sol=3
integer,parameter                         :: n_lights_sol=1242
integer,dimension(2),parameter            :: n_input_params_sol=(/2,0/)
integer,dimension(n_times_sol),parameter  :: n_groups_per_sim=(/3,1,2/)
integer,parameter                         :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_particle_types)       :: particle_types_to_test=(/&
        particle_kinetic_relativistic_id,particle_kinetic_leapfrog_id/)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
        reshape((/55,30,44,27,0,0,37,27,0/),shape(n_particles_per_group))
integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types_sol=&
        reshape((/particle_gc_vpar_id,particle_kinetic_relativistic_id,&
        particle_kinetic_relativistic_id,particle_kinetic_relativistic_id,0,0,&
        particle_kinetic_leapfrog_id,particle_kinetic_relativistic_id,0/),&
        shape(particle_types_sol))
integer,parameter                          :: n_particles_max=maxval(n_particles_per_group)
real*8,parameter                           :: survival_threshold=0.21
real*8,parameter                           :: tol_real8=5.d-16
type(full_synchrotron_light_dist)                               :: vertex_sol
type(particle_sim),dimension(n_times_sol)                       :: sims_particles
integer,dimension(n_times_sol)                                  :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol)                     :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol)     :: active_particle_ids_sol
real*8,dimension(n_times_sol)                                   :: time_vector_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol)  :: x_cart_sol

!> Interfaces --------------------------------------------------------
contains
!> Fruit basket ------------------------------------------------------
!> fruit basket containing all set-up, test and tear-down methods
subroutine run_fruit_light_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: light vertices tests"
  call prepare_input_file
  call setup
  write(*,'(/A)') "  ... running: light vertices tests"
  call run_test_case(test_deallocate_light_vertices,'test_deallocate_light_vertices')
  call run_test_case(test_ligth_read_inputs,'test_ligth_read_inputs')
  call run_test_case(test_find_all_active_particles_ids,&
  'test_find_all_active_particles_ids')
  call run_test_case(test_find_all_active_particles_ids_types,&
  'test_find_all_active_particles_ids_types')
  call run_test_case(test_store_light_from_particle_id,&
  'test_store_light_from_particle_id')
  call run_test_case(test_extract_all_particle_types,&
  'test_extract_all_particle_types')
  call run_test_case(test_extract_all_n_particles,'test_extract_all_n_particles')
  call run_test_case(test_extract_all_n_groups,'test_extract_all_n_groups')
  call run_test_case(test_fill_time_vector,'test_fill_time_vector')
  write(*,'(/A)') "  ... tearing-down: light vertices tests"
  call teardown
end subroutine run_fruit_light_vertices

!> Set-up and teard-down ---------------------------------------------
!> create an input file for testing the reading procedure
subroutine prepare_input_file()
  implicit none
  !> variables
  integer :: ifail
  !> write the file
  open(read_unit,file=input_file,status='unknown',action='write',iostat=ifail)
  write(read_unit,'(/A)') '&light_in'
  write(read_unit,'(/A,I6)') 'n_times = ',n_times_sol
  write(read_unit,'(/A,I6)') 'n_lights = ',n_lights_sol
  write(read_unit,'(/A)') '/'
  close(read_unit)
end subroutine prepare_input_file

!> set-up features common to all unit test
subroutine setup()
  use mod_gnu_rng,                          only: gnu_rng_interval
  use mod_particle_common_test_tools,       only: sim_time_interval
  use mod_particle_common_test_tools,       only: allocate_one_particle_list_type
  use mod_particle_common_test_tools,       only: fill_particles
  use mod_particle_common_test_tools,       only: invalidate_particles
  use mod_particle_common_test_tools,       only: obtain_active_particle_ids
  use mod_light_vertices_common_test_tools, only: compute_x_cart_particles
  implicit none
  integer :: ii,ifail
  !> initialisation
  ifail = 0; n_active_particles_sol = 0;
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  vertex_sol%n_property_vertex = n_properties_sol; 
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)

  !> allocate and initialise particle lists
  active_particle_ids_sol = 0;
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_particles(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
  enddo

  !> compute and store the cartesian particle position for all particles
  call compute_x_cart_particles(n_times_sol,n_groups_max,n_particles_max,&
  sims_particles,x_cart_sol)
end subroutine setup

subroutine teardown()
  implicit none
  vertex_sol%n_mhd = 0; vertex_sol%n_particle_types = 0;
  if(allocated(vertex_sol%particle_types)) deallocate(vertex_sol%particle_types)
  call vertex_sol%deallocate_vertices
  !> remove input file
  call system("rm "//input_file)
end subroutine teardown

!> Tests -------------------------------------------------------------
!> test deallocate light vertices
subroutine test_deallocate_light_vertices()
  type(full_synchrotron_light_dist) :: vertex_test
  !> allocate test vertices
  vertex_test%n_property_vertex = n_properties_sol; 
  call vertex_test%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)
  call vertex_test%setup_light_class
  call vertex_test%deallocate_light_vertices
  !> check deallocation of main variables 
  call assert_equals(0,vertex_test%n_mhd,&
  "Error deallocate light vertices: n_mhd not 0!")
  call assert_equals(0,vertex_test%n_particle_types,&
  "Error deallocate light vertices: n_particle_types not 0!")
  call assert_false(allocated(vertex_test%particle_types),&
  "Error deallocate light vertices: particle_types not deallocated!")
  call assert_equals(0,vertex_test%n_times,&
  "Error deallocate light vertices: n_times not reset!")
  call assert_equals(0,vertex_test%n_vertices,&
  "Error deallocate light vertices: n_vertices not reset!")
  call assert_false(allocated(vertex_test%n_active_vertices),&
  "Error deallocate light vertices: n_active_vertices allocated!")
  call assert_false(allocated(vertex_test%times),&
  "Error deallocate light vertices: times allocated!")
  call assert_false(allocated(vertex_test%x),&
  "Error deallocate light vertices: x allocated!")
  call assert_false(allocated(vertex_test%properties),&
  "Error deallocate light vertices: properties allocated!")
end subroutine test_deallocate_light_vertices

!> test the procedure used for reading light input files
subroutine test_ligth_read_inputs()
  implicit none
  !> variables
  integer :: rank,ifail
  integer,dimension(2) :: n_input_params
  integer,dimension(:),allocatable :: int_param
  real*8,dimension(:),allocatable  :: real_param
  !> initialisation 
  rank = 0
  !> read input file
  n_input_params = vertex_sol%return_n_light_inputs()
  open(read_unit,file=input_file,status='old',action='read',iostat=ifail)
  call vertex_sol%read_light_inputs(rank,read_unit,int_param,real_param)
  close(read_unit)
  !> checks
  call assert_equals(n_input_params_sol,n_input_params,2,&
  "Error light read inputs: N# input parameters mismatch!")
  call assert_equals((/n_times_sol,n_lights_sol/),int_param,2,&
  "Error light read inputs: int parameters mismatch!")
  call assert_false(allocated(real_param),&
  "Error light read inputs: real parameters allocated!")
  !> cleanup
  if(allocated(int_param))  deallocate(int_param)
  if(allocated(real_param)) deallocate(real_param)
end subroutine test_ligth_read_inputs

!> procedure for testing the find active particles for all particle types
subroutine test_find_all_active_particles_ids()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_groups_max,n_times_sol) :: n_active_particles
  integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids
  !> find active particles
  call vertex_sol%find_active_particles_id_time(n_groups_max,n_particles_max,&
  n_groups_per_sim,n_particles_per_group,sims_particles,n_active_particles,&
  active_particle_ids)
  !> check solutions
  call assert_equals(n_active_particles_sol,n_active_particles,n_groups_max,n_times_sol,&
  "Error light vertices find all active particle ids: N active particles mismatch!")
  call assert_equals(n_active_vertices_sol,vertex_sol%n_active_vertices,n_times_sol,&
  "Error light vertices find all active particle ids: N active vertices mismatch!")
  do ii=1,n_times_sol
    call assert_equals(active_particle_ids_sol(:,:,ii),active_particle_ids(:,:,ii),&
    n_particles_max,n_groups_max,&
    "Error light vertices find all active particle ids: active particle ids mismatch!")
  enddo
end subroutine test_find_all_active_particles_ids

!> procedure for testing the find active particles for all particle types
subroutine test_find_all_active_particles_ids_types()
  implicit none
  integer :: ii,jj
  integer,dimension(n_times_sol)              :: n_active_vertices_loc
  integer,dimension(n_groups_max,n_times_sol) :: n_active_particles,n_active_particles_loc
  integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids
  integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_loc

  !> find active particle for a specific type
  if(allocated(vertex_sol%particle_types)) deallocate(vertex_sol%particle_types);
  vertex_sol%n_particle_types = n_particle_types
  allocate(vertex_sol%particle_types(vertex_sol%n_particle_types))
  vertex_sol%particle_types = particle_types_to_test;
  call vertex_sol%find_active_particles_id_time(n_groups_max,n_particles_max,&
  n_groups_per_sim,n_particles_per_group,sims_particles,n_active_particles,&
  active_particle_ids)
  !> compute solution
  n_active_particles_loc = n_active_particles_sol;
  active_particle_ids_loc = active_particle_ids_sol;
  do ii=1,n_times_sol
    do jj=1,n_particle_types
      if(all(particle_types_sol(jj,ii).ne.particle_types_to_test).and.&
      (particle_types_sol(jj,ii).ne.0)) then
        active_particle_ids_loc(:,ii,jj) = 0
        n_active_particles_loc(ii,jj) = 0
      endif
    enddo
  enddo
  n_active_vertices_loc = sum(n_active_particles_loc,dim=1)
  !> check solutions
  call assert_equals(n_active_particles_loc,n_active_particles,n_groups_max,n_times_sol,&
  "Error light vertices find active particle ids type: N active particles mismatch!")
  call assert_equals(n_active_vertices_loc,vertex_sol%n_active_vertices,n_times_sol,&
  "Error light vertices find active particle ids type: N active vertices mismatch!")
  do jj=1,n_times_sol
    call assert_equals(active_particle_ids_loc(:,:,jj),active_particle_ids(:,:,jj),&
    n_particles_max,n_groups_max,&
    "Error light vertices find active particle ids type: active particle ids mismatch!")
  enddo
  ! restore previous conditions
  call vertex_sol%setup_light_class
end subroutine test_find_all_active_particles_ids_types

!> test store light from particle id
subroutine test_store_light_from_particle_id()
  implicit none
  integer :: kk,jj,ii,counter
  !> loop for storing the particles
  do kk = 1,n_times_sol
    counter = 0
    do jj = 1,n_groups_per_sim(kk)
      do ii=1,n_particles_per_group(jj,kk)
        counter = counter + 1
        call vertex_sol%store_light_x_from_particle_id(counter,kk,&
        sims_particles(kk)%groups(jj)%particles(ii))
      enddo
    enddo
  enddo
  !> check solution
  do ii=1,n_times_sol
    call assert_equals(x_cart_sol(:,:,ii),vertex_sol%x(:,:,ii),n_particles_max,n_groups_max,&
    tol_real8,"Error light vertices store light x from particles: positions mismatch!")
  enddo
end subroutine test_store_light_from_particle_id

!> test extract all particle types
subroutine test_extract_all_particle_types()
  implicit none
  integer,dimension(n_groups_max,n_times_sol) :: particle_types
  call vertex_sol%extract_particle_types_all_particle_sims(&
  sims_particles,n_groups_max,particle_types)
  call assert_equals(particle_types_sol,particle_types,n_groups_max,n_times_sol,&
  "Error light vertices extract all particle types: particle types mismatch!")
end subroutine test_extract_all_particle_types

!> test extract all number of particles
subroutine test_extract_all_n_particles()
  implicit none
  integer,dimension(n_groups_max,n_times_sol) :: n_particles
  call vertex_sol%extract_n_particles_all_particle_sims(sims_particles,n_groups_max,n_particles)
  call assert_equals(n_particles_per_group,n_particles,n_groups_max,n_times_sol,&
  "Error light vertices extract all particle numbers: N particles mismatch!")
end subroutine test_extract_all_n_particles

!> test extract all number of groups
subroutine test_extract_all_n_groups()
  implicit none
  integer,dimension(n_times_sol) :: n_groups
  call vertex_sol%extract_n_groups_all_particle_sims(sims_particles,n_groups)
  call assert_equals(n_groups_per_sim,n_groups,n_times_sol,&
  "Error light vertices extract all group numbers: N groups mismatch!")
end subroutine test_extract_all_n_groups

!> test fill time vector
subroutine test_fill_time_vector()
  implicit none
  call vertex_sol%fill_time_vector_particle_sims(sims_particles)
  call assert_equals(time_vector_sol,vertex_sol%times,n_times_sol,tol_real8,&
  "Error light vertices fill time vector: time vector mismatch!")
end subroutine test_fill_time_vector

!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_light_vertices_test
