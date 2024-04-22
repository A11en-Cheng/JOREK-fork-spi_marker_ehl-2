!> The mod_light_vertices module contains variables
!> and procedures used for defining and defining actions
!> of the light points
module mod_light_vertices
use mod_vertices, only: vertices
implicit none

private
public :: light_vertices

!> Variables --------------------------------------------
type,abstract,extends(vertices) :: light_vertices
  integer                          :: n_mhd            !< size of the JOREK MHD field array
  integer                          :: n_particle_types !< number of particles types
  integer,dimension(:),allocatable :: particle_types   !< list of particles type ids
  contains
  procedure,pass(light_vert)                                    :: return_n_light_inputs
  procedure,pass(light_vert)                                    :: read_light_inputs
  procedure,pass(light_vert)                                    :: fill_time_vector_particle_sims
  procedure,pass(light_vert)                                    :: extract_n_groups_all_particle_sims 
  procedure,pass(light_vert)                                    :: extract_n_particles_all_particle_sims
  procedure,pass(light_vert)                                    :: extract_particle_types_all_particle_sims
  procedure,pass(light_vert)                                    :: store_light_x_from_particle_id
  procedure,pass(light_vert)                                    :: find_active_particles_id_time
  procedure,pass(light_vert)                                    :: init_lights_from_particles
  procedure,pass(light_vert)                                    :: fill_lights_from_particles
  procedure,pass(light_vert)                                    :: deallocate_light_vertices
  procedure(direct_funct),pass(light_vert),deferred             :: directionality_funct
  procedure(spect_irradiance),pass(light_vert),deferred         :: spectral_irradiance
  procedure(comp_light_prop),pass(light_vert),deferred          :: compute_light_properties
  procedure(comp_mhd_fields),pass(light_vert),deferred          :: compute_mhd_fields
  procedure(setup_light),pass(light_vert),deferred              :: setup_light_class
  procedure(comp_particle_from_light),pass(light_vert),deferred :: compute_particle_from_light
  procedure(check_x_shaded),nopass,deferred                     :: check_x_shaded_in_emission_zone 
  procedure(check_angles_shaded),nopass,deferred                :: check_angles_shaded_in_emission_zone
end type light_vertices

!> Interfaces -------------------------------------------

interface

  !> computes the directionality function for a given point
  !> in space (cartesian coordinate)-time and a given light
  !> for all spectra wavelengths
  !> inputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_id:   (integer) id of the light to use
  !>   x_shaded:   (real8)(n_x) point illuminated by the light
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_dstb: (real8)(n_points,n_spectra) light intensity distribution
  !>               from the light_id light to the x_shaded point for all
  !>               spectra points and all spectra
  subroutine direct_funct(light_vert,spectra,time_id,light_id,x_shaded,light_dstb)
    use mod_spectra,  only: spectrum_base
    IMPORT :: light_vertices
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert
    class(spectrum_base),intent(inout)  :: spectra
    !> inputs
    integer,intent(in)                          :: time_id,light_id
    real*8,dimension(light_vert%n_x),intent(in) :: x_shaded
    !> outputs
    real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_dstb
  end subroutine direct_funct

  !> compute the spectral irradiance of a light source for a given point
  !> in space (cartesian coordinates)-time and a given light source
  !> for all wavelengths and spectra 
  !> inputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_id:   (integer) id of the light to use
  !>   x_shaded:   (real8)(n_x) point illuminated by the light
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  !>   spectra:    (spectrum_base) initilises spectra
  !>   light_spec_irradiance: (real8)(n_points,n_spectra) spectral irradiance
  !>                          from the light_id light at ethe time time_id to 
  !>                          the x_shaded point for all spectra points and all spectra
  subroutine spect_irradiance(light_vert,spectra,time_id,light_id,x_shaded,light_spec_irradiance)
    use mod_spectra,  only: spectrum_base
    IMPORT :: light_vertices
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert
    class(spectrum_base),intent(inout)  :: spectra
    !> inputs
    integer,intent(in)                  :: time_id,light_id
    real*8,dimension(light_vert%n_x),intent(in)    :: x_shaded
    !> outputs
    real*8,dimension(spectra%n_points,spectra%n_spectra),intent(out) :: light_spec_irradiance
  end subroutine spect_irradiance

  !> compute the particle light sources properties from particle simulations
  !> inputs:
  !>   light_vert:   (light_vertices) initialised light vertices
  !>   property_id:  (integer) id of the property to be initialised
  !>   time_id:      (integer) time id
  !>   particle:     (particle_base) particle from which computing the light properties
  !>   mass:         (real8) particle mass
  !>   mhd_fields:   (real8)(n_mhd) JOREK MHD field array
  !> outputs:
  !>   light_vert: (light_vertices) initialised light vertices
  subroutine comp_light_prop(light_vert,property_id,time_id,particle_in,mass,mhd_fields)
    use mod_particle_types, only: particle_base
    IMPORT :: light_vertices
    implicit none
    !> inputs-outputs
    class(light_vertices),intent(inout)           :: light_vert
    !> inputs
    class(particle_base),intent(in)               :: particle_in
    integer,intent(in)                            :: property_id,time_id
    real*8,intent(in)                             :: mass
    real*8,dimension(light_vert%n_mhd),intent(in) :: mhd_fields
  end subroutine comp_light_prop

  !> compute the JOREK MHD fields at a given location
  !> inputs:
  !>   light_vert:   (light_vertices) light vertices
  !>   fields:       (fields_base) JOREK MHD fields data structure
  !>   particle:     (particle_base) JOREK particle base 
  !>   time_id:      (integer) id of the particle simulation time
  !>   mass:         (real8) particle mass
  !> outputs:
  !>   mhd_fields:   (real8)(n_mhd) interpolated JOREK MHD fields
  subroutine comp_mhd_fields(light_vert,fields,particle_in,time_id,mass,mhd_fields)
    use mod_fields,         only: fields_base
    use mod_particle_types, only: particle_base
    IMPORT :: light_vertices
    implicit none
    !> inputs
    class(light_vertices),intent(in)          :: light_vert
    class(fields_base),intent(in)             :: fields
    class(particle_base),intent(in)           :: particle_in
    integer,intent(in)                        :: time_id
    real*8,intent(in)                         :: mass
    !> outputs
    real*8,dimension(light_vert%n_mhd),intent(out) :: mhd_fields
  end subroutine comp_mhd_fields

  !> setup the light class parameters and parameter arrays
  !> light vertex class
  !> inputs:
  !>   light_vert: (light_vertices) light vertices
  !> outputs:
  !>   light_vert: (light_vertices) light vertices
  subroutine setup_light(light_vert)
    IMPORT :: light_vertices
    implicit none
    !> inputs-outpus
    class(light_vertices),intent(inout) :: light_vert   
  end subroutine setup_light

  !> Reconstruct a particle given a light property
  !> inputs: 
  !>   light_vert:   (light_vertices) light vertices
  !>   fields:       (fields_base) JOREK MHD fields data structure
  !>   light_id:     (integer) index of the light to treat 
  !>   time_id:      (integer) index of the time light to treat
  !>   mass:         (real8) particle mass
  !> outputs:
  !>   particle_out: (particle_base) reconstructed particle 
  subroutine comp_particle_from_light(light_vert,&
  fields,light_id,time_id,mass,particle_out)
  use mod_fields,         only: fields_base
  use mod_particle_types, only: particle_base
  IMPORT :: light_vertices
  implicit none
  !> inputs
  class(light_vertices),intent(in) :: light_vert
  class(fields_base),intent(in)    :: fields
  integer,intent(in)               :: light_id,time_id
  real*8,intent(in)                :: mass
  !> outputs:
  class(particle_base),intent(out) :: particle_out
  end subroutine comp_particle_from_light

  !> check if a gather point is within the emission cone of a light source
  !> inputs:
  !>   n_x:          (integer) size of the spatial coordinate vector
  !>   x_shaded:     (real8)(n_x) x,y,z coordinates of the shaded (gather) vertex
  !>   x_light:      (real8)(n_x) x,y,z coordinates of the light vertex
  !>   n_int_param:  (integer) size of the integer parameter array 
  !>   n_real_param: (integer) size of the real parameter array
  !>   int_param:    (integer)(n_integer_param) integer parameter array
  !>   int_real:     (real8)(n_real_param) real parameter array
  !> outputs:
  !>   in_range:     (logical) true is the gather point is shaded by the light
  function check_x_shaded(n_x,x_shaded,x_light,n_int_param,n_real_param,&
  int_param,real_param) result(in_range)
    implicit none
    !> inputs:
    integer,intent(in)                        :: n_x,n_int_param,n_real_param
    integer,dimension(n_int_param),intent(in) :: int_param
    real*8,dimension(n_x),intent(in)          :: x_shaded,x_light
    real*8,dimension(n_real_param),intent(in) :: real_param
    !> outputs:
    logical :: in_range
  end function check_x_shaded

  !> check if a gather point is within the emission cone of a light source 
  !> given the gather-light point angles
  !> inputs:
  !>   n_x:          (integer) size of the spatial coordinate vector
  !>   x_shaded:     (real8)(n_x) x,y,z coordinates of the shaded (gather) vertex
  !>   x_light:      (real8)(n_x) x,y,z coordinates of the light vertex
  !>   n_int_param:  (integer) size of the integer parameter array 
  !>   n_real_param: (integer) size of the real parameter array
  !>   int_param:    (integer)(n_integer_param) integer parameter array
  !>   int_reak:     (real8)(n_real_param) real parameter array
  !> outputs:
  !>   in_range:     (logical) true is the gather point is shaded by the light
  function check_angles_shaded(n_angles,angles,n_int_param,n_real_param,&
  int_param,real_param) result(in_range)
    implicit none
    !> inputs:
    integer,intent(in)                        :: n_angles,n_int_param,n_real_param
    integer,dimension(n_int_param),intent(in) :: int_param
    real*8,dimension(n_angles),intent(in)     :: angles
    real*8,dimension(n_real_param),intent(in) :: real_param
    !> outputs:
    logical :: in_range
  end function check_angles_shaded
end interface

contains

!> Procedures ------------------------------------------- 
!> Deallocate the vertices and light vertices structures
subroutine deallocate_light_vertices(light_vert)
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vert
  !> clean-up light vertices data structure
  if(allocated(light_vert%particle_types)) deallocate(light_vert%particle_types)
  light_vert%n_particle_types = 0; light_vert%n_mhd = 0;
  !> clean-up general vertices data structure
  call light_vert%deallocate_vertices
end subroutine deallocate_light_vertices

!> Computes the position and properties of the particle light vertices 
!> for each particle and stores them in the position proprerties arrays
!> inputs:
!>   light_vert:     (light_vertices) empty light vertices
!>   n_times:        (integer) number of simulation times
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_lights_in:    (integer)(optional) number of requested lights
!> outputs:
!>   light_vert:     (light_vertices) initialised light vertices
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
subroutine init_lights_from_particles(light_vert,n_times,&
sims_particles,n_lights_in)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout)                 :: light_vert
  type(particle_sim),dimension(n_times),intent(inout) :: sims_particles
  !> inputs
  integer,intent(in)                                  :: n_times
  integer,intent(in),optional                         :: n_lights_in
  !> variables
  integer :: ii,jj
  integer :: n_lights,n_groups_max,n_particles_max
  integer,dimension(n_times)           :: n_groups,n_particles_required
  integer,dimension(:,:),allocatable   :: n_particles,particle_types,n_active_particles
  integer,dimension(:,:,:),allocatable :: active_particle_id

  !> setup-light vertices class
  call light_vert%setup_light_class
  !> initialise time vector
  call light_vert%allocate_time_vector(n_times)
  call light_vert%fill_time_vector_particle_sims(sims_particles)
  !> allocate and extract number of particles and particles type
  call light_vert%extract_n_groups_all_particle_sims(sims_particles,n_groups)
  n_groups_max = maxval(n_groups)
  allocate(n_particles(n_groups_max,light_vert%n_times)) 
  allocate(particle_types(n_groups_max,light_vert%n_times))
  call light_vert%extract_n_particles_all_particle_sims(sims_particles,&
  n_groups_max,n_particles)
  call light_vert%extract_particle_types_all_particle_sims(sims_particles,&
  n_groups_max,particle_types)
  !> compute the number of required particles per each time
  n_particles_required = 0
  do ii=1,light_vert%n_times
    do jj=1,light_vert%n_particle_types
      n_particles_required(ii) = n_particles_required(ii) + sum(n_particles(:,ii),&
      mask=particle_types(:,ii).eq.light_vert%particle_types(jj))
    enddo
  enddo
  n_particles_max = maxval(n_particles_required)
  n_lights = n_particles_max
  if(present(n_lights_in)) then
    if(n_lights.lt.n_lights_in) then
      n_lights =  n_lights_in   
    else
      write(*,*) "Error initialise light vertices from particles"
      write(*,*) "Requested number of lights < number of particles,use: ",n_lights
    endif
  endif
  !> allocate vertices
  call light_vert%allocate_x_properties(n_lights)
  !> allocate active particle arrays
  allocate(n_active_particles(n_groups_max,light_vert%n_times));
  allocate(active_particle_id(n_particles_max,n_groups_max,light_vert%n_times));
  !> find active particles for all groups and times
  call light_vert%find_active_particles_id_time(n_groups_max,n_particles_max,&
  n_groups,n_particles,sims_particles,n_active_particles,active_particle_id)
  !> fill the light data structure
  call light_vert%fill_lights_from_particles(sims_particles,&
  n_groups_max,n_particles_max,n_groups,particle_types,&
  n_active_particles,active_particle_id)
  !> cleanup 
  deallocate(n_particles); deallocate(particle_types);
  deallocate(n_active_particles); deallocate(active_particle_id)
end subroutine init_lights_from_particles


!> Fill the light vertex x and properties arrays from
!> particle lists (basic and simple openmp parallelisation)
!> inputs:
!>   light_vert:         (light_vertices) empty particle light vertices
!>   sims_particles:     (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:       (integer) maximum size of groups
!>   n_particles_max:    (integer) maximum number of particles
!>   n_groups:           (integer)(n_times) size of each group
!>   particle_types:     (integer)(n_groups_max,n_times) particle types for
!>                       each group and time 
!>   n_active_particles: (integer)(n_group_max,n_times) number of active particles
!>                       per group and per time
!>   active_particle_id: (integer)(n_particle_max,n_group_max,n_times) indices
!>                       of the active particles
!> outputs:
!>   light_vert:         (light_vertices) empty particle light vertices
subroutine fill_lights_from_particles(light_vert,&
sims_particles,n_groups_max,n_particles_max,n_groups,&
particle_types,n_active_particles,active_particles_id)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vert
  !> inputs:
  type(particle_sim),dimension(light_vert%n_times),intent(in)   :: sims_particles
  integer,intent(in)                                            :: n_groups_max
  integer,intent(in)                                            :: n_particles_max
  integer,dimension(light_vert%n_times),intent(in)              :: n_groups
  integer,dimension(n_groups_max,light_vert%n_times),intent(in) :: n_active_particles
  integer,dimension(n_groups_max,light_vert%n_times),intent(in) :: particle_types
  integer,dimension(n_particles_max,n_groups_max,light_vert%n_times),intent(in)::active_particles_id
  !> variables
  integer :: ii,jj,kk,pp
  real*8  :: psi,U
  real*8,dimension(light_vert%n_mhd) :: mhd_fields
  mhd_fields = 0d0
  !> compute synchrotron light properties from particle simulations
  do ii=1,light_vert%n_times
    pp = 0
    do jj=1,n_groups(ii)
        if(all(particle_types(jj,ii).ne.light_vert%particle_types)) cycle
        !$omp parallel do default(private) firstprivate(ii,jj,pp,n_active_particles,&
        !$omp mhd_fields) shared(sims_particles,active_particles_id,light_vert)
        do kk=1,n_active_particles(jj,ii)
          call light_vert%store_light_x_from_particle_id(pp+kk,ii,&
          sims_particles(ii)%groups(jj)%particles(active_particles_id(kk,jj,ii))) !< store position
          !> compute MHD fields
          call light_vert%compute_mhd_fields(sims_particles(ii)%fields,&
          sims_particles(ii)%groups(jj)%particles(active_particles_id(kk,jj,ii)),&
          ii,sims_particles(ii)%groups(jj)%mass,mhd_fields)
          !> compute light properties
          call light_vert%compute_light_properties(pp+kk,ii,&
          sims_particles(ii)%groups(jj)%particles(active_particles_id(kk,jj,ii)),&
          sims_particles(ii)%groups(jj)%mass,mhd_fields)
        enddo
        !$omp end parallel do
        pp = pp + n_active_particles(jj,ii) 
    enddo
  enddo
end subroutine fill_lights_from_particles

!> read_light_inputs read all the inputs required
!> for initialising lights from a opened file
!> inputs:
!>   light_vert: (light_vertices) light vertices
!>   my_id:      (integer) mpi rank
!>   r_unit:     (integer) reading unit id
!> outputs:
!>   light_vert: (light_vertices) light vertices
!>   int_param:  (integer)(:) integer input parameters
!>   real_param: (real8)(:) real8 input parameters
subroutine read_light_inputs(light_vert,my_id,r_unit,&
int_param,real_param)
  implicit none
  !> inputs-outpus
  class(light_vertices),intent(inout) :: light_vert
  !> inputs:
  integer,intent(in) :: my_id,r_unit
  !> outputs:
  integer,dimension(:),allocatable,intent(out) :: int_param
  real*8,dimension(:),allocatable,intent(out)  :: real_param
  !> variables
  integer :: n_times,n_lights
  integer,dimension(2) :: n_params
  !> definition and initialisation
  namelist /light_in/ n_times, n_lights
  !> if master, read the input file
  if(my_id.eq.0) then
    n_params = light_vert%return_n_light_inputs()
    if(allocated(int_param)) deallocate(int_param)
    allocate(int_param(n_params(1)))
    read(r_unit,light_in); int_param = (/n_times,n_lights/)
  endif
end subroutine read_light_inputs

!> return the size of the integer and real parameters to be read
!> inputs:
!>   light_vert: (light_vertices) light vertices
!> outputs:
function return_n_light_inputs(light_vert) result(n_inputs)
  implicit none
  !> inputs
  class(light_vertices),intent(in) :: light_vert
  !> outputs:
  real*8,dimension(2) :: n_inputs
  n_inputs = (/2,0/)
end function return_n_light_inputs

!> find active particles for all particle lists, particle
!> groups and simulation times
!> inputs:
!>   light_vert: (light_vertices) initialised light vertices
!>   n_groups_max: (integer) maximum number of groups
!>   n_particles_max: (integer) maximum number of particles
!>   n_groups:        (integer)(n_times) number of groups per time
!>   n_particles:     (integer)(n_groups,n_times) number of particles
!>                    per group per time
!>   sims_particles:  (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   light_vert:         (light_vertices) initialised light vertices
!>   sims_particles:     (particle_sim)(n_times) array of particle simulations
!>   n_active_particles: (integer)(n_groups,n_times) number of active
!>                       particles for each group and time
!>   active_particle_id: (integer)(n_particles_max,n_groups_max,n_times)
!>                       particle list index of active particles
subroutine find_active_particles_id_time(light_vert,n_groups_max,n_particles_max,&
n_groups,n_particles,sims_particles,n_active_particles,active_particle_id)
  use mod_particle_sim,only: particle_sim
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> inputs
  integer,intent(in) :: n_groups_max,n_particles_max
  integer,dimension(light_vert%n_times),intent(in) :: n_groups
  integer,dimension(n_groups_max,light_vert%n_times),intent(in) :: n_particles
  !> outputs
  integer,dimension(n_groups_max,light_vert%n_times),intent(out) :: n_active_particles
  integer,dimension(n_particles_max,n_groups_max,light_vert%n_times),intent(out)::active_particle_id
  !> variables
  integer :: ii
  n_active_particles = 0; active_particle_id = 0;
  if(allocated(light_vert%particle_types)) then
    if(size(light_vert%particle_types).eq.light_vert%n_particle_types) then
      do ii=1,light_vert%n_times
        call sims_particles(ii)%find_active_particles_groups(n_groups(ii),n_particles_max,&
        n_particles(:,ii),n_active_particles(:,ii),active_particle_id(:,:,ii),&
        light_vert%n_particle_types,light_vert%particle_types)
        light_vert%n_active_vertices(ii) = sum(n_active_particles(:,ii))
      enddo
      return
    endif
  endif
  do ii=1,light_vert%n_times
    call sims_particles(ii)%find_active_particles_groups(n_groups(ii),n_particles_max,&
    n_particles(:,ii),n_active_particles(:,ii),active_particle_id(:,:,ii))
    light_vert%n_active_vertices(ii) = sum(n_active_particles(:,ii))    
  enddo
end subroutine find_active_particles_id_time

!> fill the time vector from particle simulations
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
subroutine fill_time_vector_particle_sims(light_vert,sims_particles)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> variables
  integer :: ii
  do ii=1,light_vert%n_times
    light_vert%times(ii) = sims_particles(ii)%time
  enddo
end subroutine fill_time_vector_particle_sims

!> Extract the number of groups for all times (simulations)
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_groups:       (integer)(n_times) number of groups
subroutine extract_n_groups_all_particle_sims(light_vert,sims_particles,n_groups)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outpus:
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> outputs:
  integer,dimension(light_vert%n_times),intent(out) :: n_groups
  !> variables
  integer :: ii
  !> extract number of groups for all simulations
  n_groups = 0
  do ii=1,light_vert%n_times
    n_groups(ii) = sims_particles(ii)%compute_group_size()
  enddo
end subroutine extract_n_groups_all_particle_sims

!> extract the number of particles for all groups and times (simulations)
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:   (integer) maximum number of among all simulations
!>   n_groups:       (integer)(n_times) number of groups per simulation
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_particles:    (integer)(n_groups_max,n_times) number of particles
!>                   per group and per simulation
subroutine extract_n_particles_all_particle_sims(light_vert,sims_particles,&
n_groups_max,n_particles)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-ouputs:
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> inputs:
  integer,intent(in)                                             :: n_groups_max
  !> outouts:
  integer,dimension(n_groups_max,light_vert%n_times),intent(out) :: n_particles
  !> variables
  integer :: ii,n_groups
  !> extract number of particles
  n_particles = 0; n_groups = n_groups_max;
  do ii=1,light_vert%n_times
    call sims_particles(ii)%compute_particle_sizes(n_groups,n_particles(:,ii))
  enddo
end subroutine extract_n_particles_all_particle_sims

!> extract the particle types for all groups and simulations
!> inputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   n_groups_max:   (integer) maximum number of among all simulations
!> outputs:
!>   light_vert:     (light_vertices) light vertices type
!>   sims_particles: (particle_sim)(n_times) array of particle simulations
!>   particle_types: (integer)(n_groups_max,n_times) particle types
!>                   per group and per simulation. CODEX:
!>                   1 -> particle_fieldline
!>                   2 -> particle_gc
!>                   3 -> particle_gc_vpar
!>                   4 -> particle_gc_Qin
!>                   5 -> particle_kinetic
!>                   6 -> particle_kinetic_leapfrog
!>                   7 -> particle_kinetic_relativistic
!>                   8 -> particle_gc_relativistic
subroutine extract_particle_types_all_particle_sims(light_vert,sims_particles,&
n_groups_max,particle_types)
  use mod_particle_sim, only: particle_sim
  implicit none
  !> inputs-outputs:
  class(light_vertices),intent(inout)                            :: light_vert
  type(particle_sim),dimension(light_vert%n_times),intent(inout) :: sims_particles
  !> inputs:
  integer,intent(in)                                          :: n_groups_max
  !> outputs:
  integer,dimension(n_groups_max,light_vert%n_times),intent(out) :: particle_types
  !> variables
  integer :: ii,n_groups
  particle_types = 0; n_groups = n_groups_max;
  do ii=1,light_vert%n_times
    call sims_particles(ii)%find_particle_types(n_groups,particle_types(:,ii))
  enddo
end subroutine extract_particle_types_all_particle_sims

!> store the light position in cartesian coordinate give a particle,
!> the light and time index
!> inputs:
!>   light_vert: (light_vertices) light vertices structure
!>   light_id:   (integer) index of the light in the x table
!>   time_id:    (integer) index of the time in the x table
!>   particle:   (particle_base) particle structure
!> outputs:
!>   light_vert: (light_vertices) light vertices with new x entry
subroutine store_light_x_from_particle_id(light_vert,light_id,time_id,particle)
  use mod_particle_types,       only: particle_base
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  implicit none
  !> inputs-outputs
  class(light_vertices),intent(inout) :: light_vert
  !> inputs
  class(particle_base),intent(in) :: particle
  integer,intent(in)              :: light_id,time_id
  light_vert%x(:,light_id,time_id) = cylindrical_to_cartesian(particle%x)
end subroutine store_light_x_from_particle_id

!> Tools ------------------------------------------------
!> compute particle array size per requested type

end module mod_light_vertices
