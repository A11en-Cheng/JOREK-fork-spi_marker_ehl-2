!> the mod_particle_type_openmp module tests the
!> the compatibility of coupling strategies
!> between the template-like particle type and
!> the openmp directives
module mod_particle_type_openmp
use fruit
use mod_particle_types, only: particle_base
use mod_particle_types, only: particle_kinetic_relativistic
use mod_particle_types, only: particle_gc
use mod_particle_sim, only: particle_sim
implicit none

private
public :: run_fruit_particle_type_openmp

!> Parameters ----------------------------------------
!> variables *_rel_fo: relativistic full orbits
!> variables *_gc: guiding center
integer :: n_fields_rel_fo=8
integer :: n_fields_gc=7
integer*1,parameter :: q_rel_fo=1
integer*1,parameter :: q_gc=3
integer,parameter   :: i_elm_rel_fo=512
integer,parameter   :: i_elm_gc=1024
real*8,parameter    :: E_gc=1.d2
real*8,parameter    :: mu_gc=1.d-10
real*8,dimension(3),parameter :: x_rel_fo=(/3.2d0,1.d0,7.d-1/)
real*8,dimension(3),parameter :: x_gc=(/1.3d0,-5.d-2,3.d0/)
real*8,dimension(3),parameter :: p_rel_fo=(/2.d0,1.d-5,4.d1/)
real*8,dimension(2),parameter :: st_rel_fo=(/1.d-1,9.d-1/)
real*8,dimension(2),parameter :: st_gc=(/5.d-2,5.5d-1/)

!> Variables -----------------------------------------
type(particle_sim) :: sim_test,sim_sol
integer :: n_groups,n_particles
!> solutions for particle lists
integer*1,dimension(:),allocatable :: int1_rel_fo_sol,int1_gc_sol
integer,dimension(:),allocatable   :: int_rel_fo_sol,int_gc_sol
real*8,dimension(:,:),allocatable  :: real8_rel_fo_sol,real8_gc_sol

contains
!> Test baskets --------------------------------------
subroutine run_fruit_particle_type_openmp()
  implicit none

  write(*,"(/A)") "  ... set-up particle type openmp tests"
  call setup
  write(*,"(/A)") "  ... run particle type openmp tests"
  call run_test_case(test_particle_to_array,'test_particle_to_array')
  call run_test_case(test_particle_copy,'test_particle_copy')
  write(*,"(/A)") "  ... tear-down particle type openmp tests"
  call teardown
  
end subroutine run_fruit_particle_type_openmp

!> Set-up tear-down ----------------------------------
!> initialise tests
subroutine setup()
  implicit none

  !> variables
  integer :: ii

  !> variables
  n_groups = 2
  n_particles = 1000

  !> allocate groups
  allocate(sim_sol%groups(n_groups))
  allocate(sim_test%groups(n_groups))
  !> groups are initialised one by one because the particle
  !> type must be specified
  allocate(particle_kinetic_relativistic::sim_sol%groups(1)%particles(n_particles))
  allocate(particle_gc::sim_sol%groups(2)%particles(n_particles))
  allocate(particle_kinetic_relativistic::sim_test%groups(1)%particles(n_particles))
  allocate(particle_gc::sim_test%groups(2)%particles(n_particles))
  allocate(int1_rel_fo_sol(n_particles)); 
  allocate(int1_gc_sol(n_particles));
  allocate(int_rel_fo_sol(n_particles)); 
  allocate(int_gc_sol(n_particles));
  allocate(real8_rel_fo_sol(n_fields_rel_fo,n_particles)); 
  allocate(real8_gc_sol(n_fields_gc,n_particles));

  !> initialisation
  !> fill the constant type particle arrays
  int1_rel_fo_sol = q_rel_fo; int1_gc_sol = q_gc
  int_rel_fo_sol = i_elm_rel_fo; int_gc_sol = i_elm_gc;
  do ii=1,n_particles
    !> copy relativistic full orbit double variables
    real8_rel_fo_sol(1:3,ii) = x_rel_fo
    real8_rel_fo_sol(4:5,ii) = st_rel_fo
    real8_rel_fo_sol(6:8,ii) = p_rel_fo
    select type (p=>sim_sol%groups(1)%particles(ii))
    type is (particle_kinetic_relativistic)
      p%x = x_rel_fo; p%st = st_rel_fo; p%i_elm = i_elm_rel_fo;
      p%p = p_rel_fo; p%q = q_rel_fo;
    end select
  enddo
  do ii=1,n_particles
    !> copy gc double varibales
    real8_gc_sol(1:3,ii) = x_gc
    real8_gc_sol(4:5,ii) = st_gc
    real8_gc_sol(6:7,ii) = (/E_gc,mu_gc/)   
    select type (p=>sim_sol%groups(2)%particles(ii))
    type is (particle_gc)
      p%x = x_gc; p%st = st_gc; p%i_elm = i_elm_gc;
      p%E = E_gc; p%mu = mu_gc; p%q = q_gc;
    end select
  enddo

end subroutine setup

!> clean up simulation varibales
subroutine teardown()
  implicit none
  !> deallocate simulations and all their allocatables
  deallocate(sim_sol%groups); deallocate(sim_test%groups);
  deallocate(int1_rel_fo_sol); deallocate(int1_gc_sol);
  deallocate(int_rel_fo_sol);  deallocate(int_gc_sol);
  deallocate(real8_rel_fo_sol); 
  deallocate(real8_gc_sol);
end subroutine teardown

!> Tests ---------------------------------------------

!> test the copy of particle data to array using openmp
!> and type select clauses
subroutine test_particle_to_array()
  implicit none
  !> variables
  integer :: ii,jj
  integer*1,dimension(n_particles) :: int1_rel_fo,int1_gc
  integer,dimension(n_particles)   :: int_rel_fo,int_gc 
  real*8,dimension(n_fields_rel_fo,n_particles) :: real8_rel_fo
  real*8,dimension(n_fields_gc,n_particles)     :: real8_gc


  !> try to use openmp directly in select type clauses <= against standard
  !> try to use select type clauses in omp loop with private pointer
  !$omp parallel do default(private) shared(n_groups,&
  !$omp n_particles,sim_sol,int1_rel_fo,int_rel_fo,real8_rel_fo,&
  !$omp int1_gc,int_gc,real8_gc) collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
      select type (p=>sim_sol%groups(jj)%particles(ii))
      type is (particle_kinetic_relativistic)
        int1_rel_fo(ii) = p%q
        int_rel_fo(ii)  = p%i_elm
        real8_rel_fo(1:3,ii) = p%x
        real8_rel_fo(4:5,ii) = p%st
        real8_rel_fo(6:8,ii) = p%p
      type is (particle_gc)
        int1_gc(ii) = p%q
        int_gc(ii)  = p%i_elm
        real8_gc(1:3,ii) = p%x
        real8_gc(4:5,ii) = p%st
        real8_gc(6:7,ii) = (/p%E,p%mu/)
      end select
    enddo
  enddo
  !$omp end parallel do

  ! checks
  call assert_equals(int(int1_rel_fo_sol),int(int1_rel_fo),n_particles,&
  "Error particle to array: relativistic full orbit q do not match")
  call assert_equals(int(int1_gc_sol),int(int1_gc),n_particles,&
  "Error particle to array: guiding center q do not match")
  call assert_equals(int_rel_fo_sol,int_rel_fo,n_particles,&
  "Error particle to array: relativistic full orbit i_elm do not match")
  call assert_equals(int_gc_sol,int_gc,n_particles,&
  "Error particle to array: guiding center i_elm do not match")
  call assert_equals(real8_rel_fo_sol,real8_rel_fo,n_fields_rel_fo,&
  n_particles,"Error particle array: relativistic full orbit x, st, p do not match")
  call assert_equals(real8_gc_sol,real8_gc,n_fields_gc,n_particles,&
  "Error particle array: guiding center x, st, E and mu do not match")

end subroutine test_particle_to_array

!> test the copy of particle data using openmp
subroutine test_particle_copy()
  implicit none
  !> variables
  integer :: ii,jj
  integer*1,dimension(n_particles) :: int1_rel_fo,int1_gc
  integer,dimension(n_particles)   :: int_rel_fo,int_gc 
  real*8,dimension(n_fields_rel_fo,n_particles) :: real8_rel_fo
  real*8,dimension(n_fields_gc,n_particles)     :: real8_gc

  !> perform copy
  !$omp parallel do default(private) shared(sim_test,sim_sol) collapse(2)
  do jj=1,n_groups
    do ii=1,n_particles
      sim_test%groups(jj)%particles(ii) = sim_sol%groups(jj)%particles(ii)
    enddo
  enddo
  !$omp end parallel do

  !> extract arrays
  call particle_rel_fo_to_array(n_groups,n_particles,sim_test,&
  int1_rel_fo,int_rel_fo,real8_rel_fo)
  call particle_gc_to_array(n_groups,n_particles,sim_test,&
  int1_gc,int_gc,real8_gc)

  ! checks  
  call assert_equals(int(int1_rel_fo_sol),int(int1_rel_fo),n_particles,&
  "Error particle copy: relativistic full orbit q do not match")
  call assert_equals(int(int1_gc_sol),int(int1_gc),n_particles,&
  "Error particle copy: guiding center q do not match")
  call assert_equals(int_rel_fo_sol,int_rel_fo,n_particles,&
  "Error particle copy: relativistic full orbit i_elm do not match")
  call assert_equals(int_gc_sol,int_gc,n_particles,&
  "Error particle copy: guiding center i_elm do not match")
  call assert_equals(real8_rel_fo_sol,real8_rel_fo,n_fields_rel_fo,&
  n_particles,"Error particle copy: relativistic full orbit x, st, p do not match")
  call assert_equals(real8_gc_sol,real8_gc,n_fields_gc,n_particles,&
  "Error particle copy: guiding center x, st, E and mu do not match")

end subroutine test_particle_copy

!> Tools ---------------------------------------------

!> copy particle relativistic kinetic (rel fo) q, i_elm,
!> x, st and p fields into array
!> inputs:
!>   n_gr:   (integer) number of groups
!>   n_part: (integer) number of particles
!>   p_sim:  (particle_sim) particle simulation
!> outputs:
!>   int1_rel_fo:  (integer1) array of charges q
!>   int_rel_fo:   (integer) array of elements i_elm
!>   real8_rel_fo: (real8)(8,n_part) arrays of:
!>             1:3 -> global rel fo position x
!>             4:5 -> local rel fo position st
!>             6:8 -> rel fo mumentum p
subroutine particle_rel_fo_to_array(n_gr,n_part,p_sim,&
int1_rel_fo,int_rel_fo,real8_rel_fo)
  !> inputs:
  integer,intent(in) :: n_gr,n_part
  type(particle_sim),intent(in) :: p_sim
  !> outputs:
  integer*1,dimension(n_part),intent(out) :: int1_rel_fo
  integer,dimension(n_part),intent(out)   :: int_rel_fo
  real*8,dimension(8,n_part),intent(out)  :: real8_rel_fo
  !> variables
  integer :: ii,jj

  do jj=1,n_gr
    select type (p=>p_sim%groups(jj)%particles)
    type is (particle_kinetic_relativistic)
      do ii=1,n_part
        int1_rel_fo(ii) = p(ii)%q
        int_rel_fo(ii)  = p(ii)%i_elm
        real8_rel_fo(1:3,ii) = p(ii)%x
        real8_rel_fo(4:5,ii) = p(ii)%st
        real8_rel_fo(6:8,ii) = p(ii)%p
      enddo
    end select
  enddo

end subroutine particle_rel_fo_to_array

!> copy particle gc q, i_elm, x, st, E, mu fields
!> into array
!> inputs:
!>   n_gr:   (integer) number of groups
!>   n_part: (integer) number of particles
!>   p_sim:  (particle_sim) particle simulation
!> outputs:
!>   int1_gc:  (integer1) array of charges q
!>   int_gc:   (integer) array of elements i_elm
!>   real8_gc: (real8)(7,n_part) arrays of:
!>             1:3 -> global gc position x
!>             4:5 -> local gc position st
!>             6:7 -> energy and magnetic moment E,mu
subroutine particle_gc_to_array(n_gr,n_part,p_sim,&
int1_gc,int_gc,real8_gc)
  !> inputs:
  integer,intent(in) :: n_gr,n_part
  type(particle_sim),intent(in) :: p_sim
  !> outputs:
  integer*1,dimension(n_part),intent(out) :: int1_gc
  integer,dimension(n_part),intent(out)   :: int_gc
  real*8,dimension(7,n_part),intent(out)  :: real8_gc
  !> variables
  integer :: ii,jj

  do jj=1,n_gr
    select type(p=>p_sim%groups(jj)%particles)
    type is (particle_gc)
      do ii=1,n_part
        int1_gc(ii) = p(ii)%q
        int_gc(ii)  = p(ii)%i_elm
        real8_gc(1:3,ii) = p(ii)%x
        real8_gc(4:5,ii) = p(ii)%st
        real8_gc(6:7,ii) = (/p(ii)%E,p(ii)%mu/)
      enddo
    end select
  enddo

end subroutine particle_gc_to_array

end module mod_particle_type_openmp
