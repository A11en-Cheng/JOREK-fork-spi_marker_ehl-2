!> mod_light_vertices_common_test_tools contains testing procedures
!> commont to all light_vertices test
module mod_light_vertices_common_test_tools
use fruit
implicit none

private 
public :: compute_x_cart_particles

!> Variables --------------------------------------------------
integer,parameter :: n_x=3

!> Interfaces -------------------------------------------------
interface compute_x_cart_particles
  module procedure compute_x_cart_all_particles
end interface compute_x_cart_particles

contains
!> Procedures -------------------------------------------------
!> compute and store particle positions in cartesian coordinates
!> inputs:
!>   n_times:         (integer) number of time slices
!>   n_groups_max:    (integer) maximum number of groups
!>   n_particles_max: (integer) maximum number of particles
!>   sims_particles:  (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   sims_particles:  (particle_sim)(n_times) array of particle simulations
!>   x_cart:          (n_x,n_particles,n_times) position in 
!>                    cartesian coordinates
subroutine compute_x_cart_all_particles(n_times,n_groups_max,&
n_particles_max,sims_particles,x_cart)
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_particle_sim,          only: particle_sim
  implicit none
  !> inputs
  integer,intent(in) :: n_times,n_groups_max,n_particles_max
  !> inputs-outputs
  type(particle_sim),dimension(n_times),intent(inout) :: sims_particles
  !> outputs
  real*8,dimension(n_x,n_groups_max*n_particles_max,n_times) :: x_cart
  !> variables
  integer :: ii,jj,kk,counter

  !> initialisation
  x_cart = 0.d0;
  !> compute and store particle positions in cartesian coordinates
  do kk=1,n_times
    counter = 0
    do jj=1,size(sims_particles(kk)%groups)
      do ii=1,size(sims_particles(kk)%groups(jj)%particles)
        counter = counter + 1
        x_cart(:,counter,kk) = cylindrical_to_cartesian(&
        sims_particles(kk)%groups(jj)%particles(ii)%x)
      enddo
    enddo
  enddo
end subroutine compute_x_cart_all_particles
!>-------------------------------------------------------------

end module mod_light_vertices_common_test_tools
