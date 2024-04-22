!> mod_particle_types_test contains variables and 
!> procedure for testing the methods coded in
!> mod_particle_types
module mod_particle_types_test
use fruit
use mod_particle_types
use mod_particle_sim, only: particle_group
use mod_particle_common_test_tools, only: n_particle_types
implicit none

private
public :: run_fruit_particle_types

!> Variables ------------------------------------
integer,parameter :: n_particles=1000
real*8,parameter  :: survival_threshold=3.25d-1
integer,dimension(n_particle_types) :: n_active_particles_sol
integer,dimension(n_particle_types) :: particle_type_list_sol
integer*1,dimension(n_particles,n_particle_types) :: particle_charge_list_sol
integer,dimension(n_particles,n_particle_types)   :: active_particle_ids_sol
type(particle_group),dimension(n_particle_types)  :: groups_sol
!> Interfaces -----------------------------------

contains

!> Fruit basket ---------------------------------
!> fruit basket containing all set-up, test and
!> tear-down procedures
subroutine run_fruit_particle_types
  implicit none
  write(*,'(/A)') "  ... setting-up: particle types tests"
  call setup
  write(*,'(/A)') "  ... running: particle types tests"
  call run_test_case(test_particle_copy,'test_particle_copy')
  call run_test_case(test_individual_particle_copy,&
  'test_individual_particle_copy')
  call run_test_case(test_particle_get_q,'test_particle_get_q')
  call run_test_case(test_codify_single_particle_type,&
  'test_codify_single_particle_type')
  call run_test_case(test_codify_particle_list,&
  'test_codify_particle_list')
#ifdef UNIT_TESTS
  call run_test_case(test_find_active_particle_id_seq,&
  'test_find_active_particle_id_seq')
  call run_test_case(test_find_active_particle_id_omp,&
  'test_find_active_particle_id_omp')
#endif
  call run_test_case(test_find_active_particle_id,&
  'test_find_active_particle_id')
  call run_test_case(test_find_active_particle_id_type,&
  'test_find_active_particle_id_type')
  write(*,'(/A)') "  ... tearing-down: particle types tests"
end subroutine run_fruit_particle_types

!> Set-up and tear-down -------------------------
!> set-up particle types test features
subroutine setup
  use mod_particle_common_test_tools, only: fill_particles
  use mod_particle_common_test_tools, only: fill_groups
  use mod_particle_common_test_tools, only: invalidate_particles
  use mod_particle_common_test_tools, only: obtain_active_particle_ids
  use mod_particle_common_test_tools, only: obtain_particle_charges
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  implicit none
  !> variables
  integer :: ifail
  
  !> allocate the particle lists
  ifail = 0
  call allocate_one_particle_list_type(n_particle_types,n_particles,groups_sol,ifail)
  call assert_true(ifail.eq.0,"Error particle_types test setup: particle list not allocated!")

  !> fill-up the group and particle base variables
  call fill_groups(n_particle_types,groups_sol)
  call fill_particles(n_particle_types,groups_sol)

  !> invalidate particles in particle lists
  call invalidate_particles(n_particle_types,n_particles,survival_threshold,&
  n_active_particles_sol,groups_sol)
  call obtain_active_particle_ids(n_particle_types,n_particles,&
  active_particle_ids_sol,groups_sol)

  !> initialise the particle type list
  particle_type_list_sol = (/particle_fieldline_id,particle_gc_id,particle_gc_vpar_id,&
  particle_kinetic_id,particle_kinetic_leapfrog_id,particle_kinetic_relativistic_id,&
  particle_gc_relativistic_id,particle_gc_Qin_id/)
  
  !> get particle charges
  call obtain_particle_charges(n_particle_types,n_particles,particle_charge_list_sol,groups_sol)
end subroutine setup

!> Tests ----------------------------------------
!> test particle copy function
subroutine test_particle_copy
  use mod_particle_sim, only: particle_group
  use mod_particle_assert_equal, only: assert_equal_particle
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  implicit none
  !> variables
  type(particle_group),dimension(n_particle_types) :: group_particles
  integer :: ii,jj,ifail
  !> allocate particle lists
  ifail = 0
  call allocate_one_particle_list_type(n_particle_types,n_particles,group_particles,ifail)
  call assert_true(ifail.eq.0,"Error particle_types test copy: particle list not allocated!")
  !> copy particles
  !$omp parallel do default(shared) private(ii,jj) collapse(2)
  do jj=1,n_particle_types
    do ii=1,n_particles
      group_particles(jj)%particles(ii) = groups_sol(jj)%particles(ii)
    enddo
  enddo
  !$omp end parallel do
  !> compare particle list
  do ii=1,n_particle_types
    call assert_equal_particle(n_particles,groups_sol(ii)%particles,group_particles(ii)%particles)
  enddo
end subroutine test_particle_copy

!> test indiviual particle copy functions
!> it assumes that the standard particle
!> copy function works
subroutine test_individual_particle_copy
  use mod_particle_assert_equal, only: assert_equal_particle
  implicit none
  type(particle_kinetic_leapfrog) :: p_leapfrog_in,p_leapfrog_out
  integer :: ii

  !> initialise particle in
  !> copy particle types assuming that the standard particle copy function works
  do ii=1,n_particle_types
    select type (p_in=>groups_sol(ii)%particles(1))
    type is (particle_kinetic_leapfrog)
      p_leapfrog_in = p_in
    end select
  enddo
  
  !> apply individual copy function
  call copy_particle_kinetic_leapfrog(p_leapfrog_in,p_leapfrog_out)

  !> test equality
  call assert_equal_particle(p_leapfrog_out,p_leapfrog_in)
end subroutine test_individual_particle_copy

!> test return charge function
subroutine test_particle_get_q
  implicit none
  !> variables
  integer :: ii,jj
  integer*1,dimension(n_particles,n_particle_types) :: q_array
  !> extract particle charge and test
  !$omp parallel do default(shared) private(ii,jj) collapse(2)
  do jj=1,n_particle_types
    do ii=1,n_particles
      q_array(ii,jj) = particle_get_q(groups_sol(jj)%particles(ii))
    enddo
  enddo
  !$omp end parallel do
  call assert_equals(int(particle_charge_list_sol),int(q_array),n_particles,&
  n_particle_types,"Errpr particle_type test get q: particle charge mismatch!")
end subroutine test_particle_get_q

!> test return particle type code for each single particle
subroutine test_codify_single_particle_type()
  implicit none
  !> variables
  integer :: ii,jj
  integer,dimension(n_particles) :: particle_code

  !> extract the particle code for every particle and test it
  do jj=1,n_particle_types 
    particle_code = -1
    !$omp parallel do default(shared) firstprivate(jj) private(ii)
    do ii=1,n_particles
      particle_code(ii) = codify_particle_type(groups_sol(jj)%particles(ii))
    enddo
    !$omp end parallel do
    call assert_true(all(particle_code.eq.particle_type_list_sol(jj)),&
    "Error in particle_types codify single particle type: particle types mismatch!")
  enddo
end subroutine test_codify_single_particle_type

!> test return particle type code for the whole particle list
subroutine test_codify_particle_list()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: list_code

  !> extract particle list code 
  do ii=1,n_particle_types
    list_code(ii) = codify_particle_type(groups_sol(ii)%particles)
  enddo
  call assert_equals(particle_type_list_sol,list_code,n_particle_types,&
  "Error in particle_types codify particle list type: list types mismatch!")
end subroutine test_codify_particle_list

#ifdef UNIT_TESTS
!> test find active particle id, sequential version
subroutine test_find_active_particle_id_seq()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: n_active_particles
  integer,dimension(n_particles,n_particle_types) :: active_particle_ids
  
  !> compute active particles and their id
  n_active_particles = -1; active_particle_ids = -1;
  do ii=1,n_particle_types
    call find_active_particle_id_seq(n_particles,groups_sol(ii)%particles,&
    n_active_particles(ii),active_particle_ids(:,ii))
  enddo
  !> check results
  call assert_equals(n_active_particles_sol,n_active_particles,n_particle_types,&
  "Error particle_types find active particle sequential: n_active_particles mismatch!")
  call assert_equals(active_particle_ids_sol,active_particle_ids,n_particles,n_particle_types,&
  "Error particle_types find active particle sequential: active_particle_ids mismatch!")
end subroutine test_find_active_particle_id_seq

!> test find active particle id, parallel version
subroutine test_find_active_particle_id_omp()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: n_active_particles
  integer,dimension(n_particles,n_particle_types) :: active_particle_ids

  !> compute active particles and their id
  n_active_particles = -1; active_particle_ids = -1;
  do ii=1,n_particle_types
    call find_active_particle_id_openmp(n_particles,groups_sol(ii)%particles,&
    n_active_particles(ii),active_particle_ids(:,ii))
  enddo 
  !> check results
  call assert_equals(n_active_particles_sol,n_active_particles,n_particle_types,&
  "Error particle_types find active particle openmp: n_active_particles mismatch!")
  call assert_equals(active_particle_ids_sol,active_particle_ids,n_particles,n_particle_types,&
  "Error particle_types find active particle openmp: active_particle_ids mismatch!")
end subroutine test_find_active_particle_id_omp
#endif

!> test find active particle id interface notype
subroutine test_find_active_particle_id()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: n_active_particles
  integer,dimension(n_particles,n_particle_types) :: active_particle_ids

  !> compute active particles and their id
  n_active_particles = -1; active_particle_ids = -1;
  do ii=1,n_particle_types
    call find_active_particle_id(n_particles,groups_sol(ii)%particles,&
    n_active_particles(ii),active_particle_ids(:,ii))
  enddo 
  !> check results
  call assert_equals(n_active_particles_sol,n_active_particles,n_particle_types,&
  "Error particle_types find active particle notype: n_active_particles mismatch!")
  call assert_equals(active_particle_ids_sol,active_particle_ids,n_particles,n_particle_types,&
  "Error particle_types find active particle notype: active_particle_ids mismatch!")
end subroutine test_find_active_particle_id

!> test find active particle id interface notype
subroutine test_find_active_particle_id_type()
  implicit none
  !> variables
  integer :: ii
  integer,dimension(n_particle_types) :: n_active_particles
  integer,dimension(n_particles,n_particle_types) :: active_particle_ids

  !> compute active particles and their id
  n_active_particles = -1; active_particle_ids = -1;
  do ii=1,n_particle_types
    call find_active_particle_id(particle_type_list_sol(ii),&
    n_particles,groups_sol(ii)%particles,&
    n_active_particles(ii),active_particle_ids(:,ii))
  enddo 
  !> check results
  call assert_equals(n_active_particles_sol,n_active_particles,n_particle_types,&
  "Error particle_types find active particle type: n_active_particles mismatch!")
  call assert_equals(active_particle_ids_sol,active_particle_ids,n_particles,n_particle_types,&
  "Error particle_types find active particle type: active_particle_ids mismatch!")
end subroutine test_find_active_particle_id_type

!>-----------------------------------------------
end module mod_particle_types_test
