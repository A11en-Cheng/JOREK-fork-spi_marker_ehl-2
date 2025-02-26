!> mod_gyroaverage_synchrotron_light_dist_test contains all variables and
!> procedures for testing the gyroaveraged synchrotron light vertices
module mod_gyroaverage_synchrotron_light_dist_vertices_test
use fruit
use mod_particle_types,                              only: particle_gc_relativistic_id
use mod_particle_types,                              only: particle_kinetic_relativistic_id
use mod_particle_types,                              only: particle_kinetic_id
use mod_particle_types,                              only: particle_gc_vpar_id
use mod_particle_sim,                                only: particle_sim
use mod_spectra_deterministic,                       only: spectrum_integrator_2nd
use mod_gyroaverage_synchrotron_light_dist_vertices, only: gyroaverage_synchrotron_light_dist
implicit none

private
public :: run_fruit_gyroaverage_synchrotron_light_dist_vertices

!> Variables ---------------------------------------------------------
!> general parameters
real*8,parameter :: onethird=1d0/3d0
real*8,parameter :: twothirds=2d0/3d0
real*8,parameter :: tol_real8=5d-12
real*8,parameter :: tol2_real8=5d-1
real*8,parameter :: tol3_real8=5d-8
real*8,parameter :: tol4_real8=5d1
real*8,parameter :: mass_RE=5.48579909065d-4
!> parameters for generating synchrotron lights
integer,parameter :: n_mhd_sol=16
integer,parameter :: n_properties=10
integer,parameter :: n_x=3
integer,parameter :: fill_type_base=1 !< use cylindrical initialisation
integer,parameter :: n_times_sol=2
integer,parameter :: n_particle_types_check_sol=1
integer,dimension(n_times_sol),parameter              :: n_groups_per_sim=(/4,3/)
integer,parameter                                     :: n_groups_max=maxval(n_groups_per_sim)
integer,dimension(n_groups_max,n_times_sol),parameter :: n_particles_per_group=&
           reshape((/132,324,42,10,237,143,23,34/),shape(n_particles_per_group))
integer,parameter :: n_particles_max=maxval(n_particles_per_group)
integer,dimension(n_groups_max,n_times_sol),parameter :: particle_types_sol=&
           reshape((/particle_kinetic_relativistic_id,particle_gc_relativistic_id,&
           particle_kinetic_id,particle_gc_vpar_id,particle_gc_relativistic_id,&
           particle_gc_vpar_id,particle_kinetic_id,0/),shape(particle_types_sol))
real*8,parameter                                      :: survival_threshold=0.45
!> parameters for generating spectra
integer,parameter                                     :: n_spectra=2
integer,parameter                                     :: n_lines_per_spectrum=16
real*8,dimension(n_spectra),parameter                 :: min_wlen=(/3d-6,2.5d-7/)
real*8,dimension(n_spectra),parameter                 :: max_wlen=(/3.5d-6,4.2d-7/)
!> parameters for generating shadowed points
integer,parameter                                     :: n_shaded_points=57
integer,parameter                                     :: n_shaded_points_per_particle=7
real*8,dimension(2),parameter                         :: length_shadowed=(/2d-2,7d1/) 
!> variables for generating synchrotron lights
type(gyroaverage_synchrotron_light_dist)              :: vertex_sol
type(particle_sim),dimension(n_times_sol)             :: sims_particles
integer                                               :: n_gc_RE_max
integer,dimension(n_particle_types_check_sol)         :: particle_types_check_sol
integer,dimension(n_times_sol)                        :: n_active_vertices_sol
integer,dimension(n_groups_max,n_times_sol)           :: n_active_particles_sol
integer,dimension(n_particles_max,n_groups_max,n_times_sol) :: active_particle_ids_sol
real*8,dimension(n_times_sol)                         :: time_vector_sol
real*8,dimension(n_particles_max*n_groups_max,n_times_sol)              :: rel_fact_parallel_sol
real*8,dimension(n_particles_max*n_groups_max,n_times_sol)              :: normB_sol
real*8,dimension(n_particles_max*n_groups_max,n_times_sol)              :: weight_sol
real*8,dimension(n_x,n_particles_max*n_groups_max,n_times_sol)          :: x_cart_sol
real*8,dimension(n_properties,n_particles_max*n_groups_max,n_times_sol) :: properties_sol
!> variables for generating spectra
type(spectrum_integrator_2nd)                         :: spectrum
!> variables for generating shadowed points
real*8,dimension(:,:,:,:),allocatable                 :: x_shadowed

!> Interfaces --------------------------------------------------------
contains
!> Fruit test basket -------------------------------------------------
!> fruit basket having all set-up, tests and tearing-down procedures
subroutine run_fruit_gyroaverage_synchrotron_light_dist_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: gyroaverage synchrotron light vertices tests"
  call setup
  write(*,'(/A)') "  ... running: gyroaverage synchrotron light vertices tests"
  call run_test_case(test_setup_gyroaverage_synchrotron_radiation_class,&
  'test_setup_gyroaverage_synchrotron_radiation_class')
  call run_test_case(test_compute_gyroaverage_synchrotron_mhd_fields,&
  'test_compute_gyroaverage_synchrotron_mhd_fields')
  call run_test_case(test_compute_gyroaverage_synchrotron_light_properties,&
  'test_compute_gyroaverage_synchrotron_light_properties')
  call run_test_case(test_init_gyroaverage_synchrotron_lights_from_gc,&
  'test_init_gyroaverage_synchrotron_lights_from_gc')
  call run_test_case(test_check_shaded_angles_in_gyroaverage_synchrotron_cone,&
  'test_check_shaded_angles_in_gyroaverage_synchrotron_cone')
  call run_test_case(test_gyroaverage_synchrotron_irradiance_directionality_funct,&
  'test_gyroaverage_synchrotron_irradiance_directionality_funct')
  call run_test_case(test_gyroaverage_synchrotron_irradiance_dir_funct_taskloop,&
  'test_gyroaverage_synchrotron_irradiance_dir_funct_taskloop')
  call run_test_case(test_particle_from_gyroaverage_synchrotron_lights,&
  'test_particle_from_gyroaverage_synchrotron_lights')
  write(*,'(/A)') "  ... tearing-down: gyroaverage synchrotron light vertices tests"
  call teardown
end subroutine run_fruit_gyroaverage_synchrotron_light_dist_vertices

!> Set-up and tear-down procedures------------------------------------
!> allocate and initialise the unit test features
subroutine setup()
  use mod_rng,                        only: type_rng
  use mod_pcg32_rng,                  only: pcg32_rng
  use mod_gnu_rng,                    only: gnu_rng_interval
  use mod_common_test_tools,          only: omp_initialize_rngs
  use mod_particle_common_test_tools, only: sim_time_interval
  use mod_particle_common_test_tools, only: allocate_one_particle_list_type
  use mod_particle_common_test_tools, only: fill_groups,fill_mass_RE 
  use mod_particle_common_test_tools, only: fill_particles_tokamak
  use mod_particle_common_test_tools, only: invalidate_particles
  use mod_particle_common_test_tools, only: obtain_active_particle_ids
  use mod_fields_linear,              only: jorek_fields_interp_linear
  !$ use omp_lib
  implicit none
  !> variables
  integer :: ii,jj,ifail,n_gc_RE_max_loc,n_threads
  class(type_rng),dimension(:),allocatable :: rngs
  !> initialisation
  ifail = 0; particle_types_check_sol = (/particle_gc_relativistic_id/);
  vertex_sol%n_mhd = n_mhd_sol; vertex_sol%n_property_vertex = n_properties;
  n_gc_RE_max = 0; n_active_particles_sol = 0; n_threads = 1;
  !$ n_threads = omp_get_max_threads()
  call gnu_rng_interval(n_times_sol,sim_time_interval,time_vector_sol)
  call vertex_sol%allocate_vertices(n_times_sol,n_particles_max*n_groups_max)
  !> allocate and initialise particle lists
  do ii=1,n_times_sol
    sims_particles(ii)%time = time_vector_sol(ii)
    allocate(jorek_fields_interp_linear::sims_particles(ii)%fields)
    allocate(sims_particles(ii)%groups(n_groups_per_sim(ii)))
    call allocate_one_particle_list_type(n_groups_per_sim(ii),&
    n_particles_per_group(1:n_groups_per_sim(ii),ii),&
    particle_types_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups,ifail)
    call fill_groups(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_mass_RE(n_groups_per_sim(ii),sims_particles(ii)%groups)
    call fill_particles_tokamak(n_groups_per_sim(ii),sims_particles(ii)%groups,fill_type_base)
    call invalidate_particles(n_groups_per_sim(ii),n_particles_max,survival_threshold,&
    n_active_particles_sol(1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    call obtain_active_particle_ids(n_groups_per_sim(ii),n_particles_max,&
    active_particle_ids_sol(:,1:n_groups_per_sim(ii),ii),sims_particles(ii)%groups)
    n_active_vertices_sol(ii) = sum(n_active_particles_sol(:,ii))
    where(particle_types_sol(:,ii).ne.particle_gc_relativistic_id)
      n_active_particles_sol(:,ii) = 0
    endwhere
    n_gc_RE_max_loc = 0
    do jj=1,n_groups_per_sim(ii)
      if(particle_types_sol(jj,ii).ne.particle_gc_relativistic_id) cycle
      n_gc_RE_max_loc = n_gc_RE_max_loc + n_particles_per_group(jj,ii)
    enddo
      n_gc_RE_max = max(n_gc_RE_max,n_gc_RE_max_loc)
  enddo
  !> initialise positions and properties tables
  call compute_gyroavg_synch_x_properties_ana()
  !> initialise deterministic spectra
  spectrum = spectrum_integrator_2nd(n_lines_per_spectrum,n_spectra,min_wlen,max_wlen)
  call spectrum%generate_spectrum()
  !> generate shadowed point positions
  allocate(x_shadowed(n_x,n_shaded_points_per_particle,n_gc_RE_max,n_times_sol))
  call compute_x_shadowed_gc
end subroutine setup 

!> destroy all test features
subroutine teardown()
  implicit none
  call vertex_sol%deallocate_vertices; call spectrum%deallocate_spectrum;
  if(allocated(x_shadowed)) deallocate(x_shadowed)
end subroutine teardown

!> Tests -------------------------------------------------------------
!> Test setup the gyroaverage synchrotron radiation class
subroutine test_setup_gyroaverage_synchrotron_radiation_class()
  implicit none
  !> setup the gyroaverage synchrotron light class
  call vertex_sol%setup_light_class
  call assert_equals(n_properties,vertex_sol%n_property_vertex,&
  "Error check setup gyroaverage synchrotron light class: wrong size of the vertex properties array!")
  call assert_equals(n_mhd_sol,vertex_sol%n_mhd,&
  "Error check setup gyroaverage synchrotron light class: wrong size of the mhd array!")
  call assert_equals(n_particle_types_check_sol,vertex_sol%n_particle_types,&
  "Error check setup gyroaverage synchrotron light class: wrong size of the particle types array!")
  call assert_equals(particle_types_check_sol,vertex_sol%particle_types,&
  n_particle_types_check_sol,"Error check setup gyroaverage synchrotron light class: wrong particle types list!")
end subroutine test_setup_gyroaverage_synchrotron_radiation_class

!> test the method used for computing the MHD fields required by GCs
subroutine test_compute_gyroaverage_synchrotron_mhd_fields()
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  real*8,dimension(n_mhd_sol,n_particles_max*n_groups_max) :: mhd_fields_test
  real*8,dimension(n_mhd_sol,n_particles_max*n_groups_max) :: mhd_fields_sol
  !> loop over all active particles for performing the tests
  do kk=1,n_times_sol
    counter = 0; mhd_fields_test = 0d0; mhd_fields_sol = 0d0;
    do jj=1,n_groups_per_sim(kk)
      do ii=1,n_particles_per_group(jj,kk)
        if(sims_particles(kk)%groups(jj)%particles(ii)%i_elm.le.0) cycle
        counter = counter + 1
        !> compute the solution mhd fields
        call compute_mhd_fields_gc_cyl(sims_particles(kk)%groups(jj)%particles(ii),&
        mhd_fields_sol(1:3,counter),mhd_fields_sol(4:6,counter),mhd_fields_sol(7,counter),&
        mhd_fields_sol(8:10,counter),mhd_fields_sol(11:13,counter),mhd_fields_sol(14:16,counter))
        !> compute the test mhd fields
        call vertex_sol%compute_mhd_fields(sims_particles(kk)%fields,&
        sims_particles(kk)%groups(jj)%particles(ii),kk,&
        sims_particles(kk)%groups(jj)%mass,mhd_fields_test(:,counter))
      enddo
    enddo
    !> check the results
    call assert_equals(mhd_fields_sol,mhd_fields_test,n_mhd_sol,n_particles_max*n_groups_max,&
    tol_real8,"Error gyroaverage synchrotron light compute MHD fields: too large errors!")
  enddo
end subroutine test_compute_gyroaverage_synchrotron_mhd_fields

!> test the method for computing the gyroaverage synchrotron radiation properties
subroutine test_compute_gyroaverage_synchrotron_light_properties()
  use constants, only: PI,ATOMIC_MASS_UNIT,EL_CHG,EPS_ZERO,SPEED_OF_LIGHT
  use mod_particle_types, only: particle_gc_relativistic
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  real*8  :: factor
  real*8,dimension(16) :: mhd_fields
  real*8,dimension(n_particles_max*n_groups_max)              :: error_coeff,zeros_1
  real*8,dimension(n_properties,n_particles_max*n_groups_max) :: error,zeros
  !> loop for computing the properties and do the testing
  zeros = 0d0; zeros_1 = 0d0; error_coeff = 0d0; error = 0d0;
  do kk=1,n_times_sol
    counter = 0; error = 0d0;
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is (particle_gc_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.le.0) cycle
          counter = counter + 1
          !> compute the analytical MHD fields for computing the GC velocity
          call compute_mhd_fields_gc_cyl(p_list(ii),mhd_fields(1:3),mhd_fields(4:6),&
          mhd_fields(7),mhd_fields(8:10),mhd_fields(11:13),mhd_fields(14:16))
          !> compute the light properties
          call vertex_sol%compute_light_properties(counter,kk,p_list(ii),&
          sims_particles(kk)%groups(jj)%mass,mhd_fields)
          !> check if the pre-multiplicative coefficients match
          factor = p_list(ii)%weight*((9d0*((EL_CHG*abs(real(p_list(ii)%q,kind=8)))**5)*&
          (properties_sol(4,counter,kk)**9)*(properties_sol(5,counter,kk)**2)*(mhd_fields(7)**3))/&
          (2.56d2*(PI**3)*EPS_ZERO*(SPEED_OF_LIGHT**2)*(rel_fact_parallel_sol(counter,kk)**2)*&
          ((sims_particles(kk)%groups(jj)%mass*ATOMIC_MASS_UNIT)**3)))
          error_coeff(counter) = abs(vertex_sol%properties(10,counter,kk)*&
          vertex_sol%properties(9,counter,kk) - factor)
          if(factor.ne.0d0) error_coeff(counter) = error_coeff(counter)/factor
        enddo
      end select
    enddo
    where(properties_sol(:,:,kk).ne.0d0) error = abs((vertex_sol%properties(:,:,kk) - &
    properties_sol(:,:,kk))/properties_sol(:,:,kk))
    !> check if the properties arrays are equal
    call assert_equals(zeros,error,n_properties,n_particles_max*n_groups_max,&
    tol_real8,"Error gyroaverage synchrotron light compute properties: too large errors!")
    call assert_equals(zeros_1,error_coeff,n_particles_max*n_groups_max,&
    tol_real8,"Error gyroaverage synchrotron light compute properties: normalisation mismatch!")
  enddo
end subroutine test_compute_gyroaverage_synchrotron_light_properties

!> test the initialisation of gyroaverage synchrotron lights form relativistic gc
subroutine test_init_gyroaverage_synchrotron_lights_from_gc()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer,parameter                                      :: n_synch_fail=1
  integer                                                :: ii,n_particles_time
  real*8,dimension(n_x,n_gc_RE_max,n_times_sol)          :: x_cart_loc
  real*8,dimension(n_properties,n_gc_RE_max,n_times_sol) :: properties_loc
  !> test datastructure initialisation with given size
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles,&
  n_particles_max*n_groups_max)
  call assert_equals_rel_error(n_x,n_particles_max*n_groups_max,n_times_sol,&
  x_cart_sol,vertex_sol%x,tol_real8,&
  "Error init gyroaverage synchrotron lights from particles set n lights large: positions errors too large!")
  call assert_equals_rel_error(n_properties,n_particles_max*n_groups_max,&
  n_times_sol,properties_sol,vertex_sol%properties,tol_real8,&
  "Error init gyroaverage synchrotron lights from particles set n lights large: properties errors too large!")

  !> copy values of x and properties
  do ii=1,n_times_sol
    n_particles_time                        = sum(n_active_particles_sol(:,ii))
    x_cart_loc(:,1:n_particles_time,ii)     = x_cart_sol(:,1:n_particles_time,ii)
    properties_loc(:,1:n_particles_time,ii) = properties_sol(:,1:n_particles_time,ii)
  enddo
  !> test lights initialisation without inputs
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles)
  call assert_equals(shape(x_cart_loc),shape(vertex_sol%x),3,&
  "Error init gyroaverage synchrotron lights from particles: positions shape mismatch!")
  call assert_equals(shape(properties_loc),shape(vertex_sol%properties),3,&
  "Error init gyroaverage synchrotron lights from particles: properties shape mismatch!")
  call assert_equals_rel_error(n_x,n_gc_RE_max,n_times_sol,&
  x_cart_loc,vertex_sol%x,tol_real8,&
  "Error init gyroaverage synchrotron lights from particles: positions errors too large!")
  call assert_equals_rel_error(n_properties,n_gc_RE_max,&
  n_times_sol,properties_loc,vertex_sol%properties,tol_real8,&
  "Error init gyroaverage synchrotron lights from particles: properties errors too large!")

  !> test lights initialisation with too small number of vertices
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles,n_synch_fail)
  call assert_equals(shape(x_cart_loc),shape(vertex_sol%x),3,&
  "Error init gyroaverage synchrotron lights from particles set n lights small: positions shape mismatch!")
  call assert_equals(shape(properties_loc),shape(vertex_sol%properties),3,&
  "Error init gyroaverage synchrotron lights from particles set n lights small: properties shape mismatch!")
  call assert_equals_rel_error(n_x,n_gc_RE_max,n_times_sol,&
  x_cart_loc,vertex_sol%x,tol_real8,&
  "Error init gyroaverage synchrotron lights from particles set n lights small: positions errors too large!")
  call assert_equals_rel_error(n_properties,n_gc_RE_max,&
  n_times_sol,properties_loc,vertex_sol%properties,tol_real8,&
  "Error init gyroaverage synchrotron lights from particles set n lights small: properties errors too large!")  

end subroutine test_init_gyroaverage_synchrotron_lights_from_gc

!> test the check if a point is within the radiation cone of 
!> a synchrotron light or outside (basically the test of the 
!> same routine of the full synchrotron model) using their angles.
subroutine test_check_shaded_angles_in_gyroaverage_synchrotron_cone()
  use constants,                 only: PI,TWOPI
  use mod_coordinate_transforms, only: cylindrical_to_cartesian_velocity
  use mod_coordinate_transforms, only: vectors_to_orthonormal_basis
  use mod_sampling,              only: sample_uniform_sphere_corona_rcosphi
  use mod_sampling,              only: sample_uniform_cone
  use mod_particle_types,        only: particle_gc_relativistic
  use mod_gc_relativistic,       only: compute_relativistic_factor
  implicit none
  !> variables
  integer                            :: ii,n_int_param,n_real_param
  integer                            :: n_angles,n_real_param_2
  integer,dimension(0)               :: int_param
  integer,dimension(n_x)             :: particle_id
  real*8                             :: normB,costheta,sintheta,cospsi
  real*8                             :: sinpitch_angle,cospitch_angle
  real*8                             :: pitch_angle,p_perp,cosmu
  real*8,dimension(4)                :: real_param
  real*8,dimension(n_x)              :: E_field,b_field,gradB,curlb,dbdt
  real*8,dimension(n_x)              :: rnd3,x_shaded,x_light,tang,nor,binor
  logical                            :: ifail
  logical,dimension(n_shaded_points) :: in_cone_points,out_cone_points
  !> initialisation: extract a random light
  n_angles = 2; n_int_param = size(int_param); n_real_param = size(real_param);
  n_real_param_2 = 1; in_cone_points = .false.; out_cone_points = .false.;
  x_light = 0d0; real_param = 0d0; ifail = .true.;
  do while(ifail)
    call random_number(rnd3); 
   particle_id(1) = min(1+floor(real(n_times_sol,kind=8)*rnd3(1)),n_times_sol)
   particle_id(2) = min(1+floor(real(n_groups_per_sim(particle_id(1)),kind=8)*rnd3(2)),&
                    n_groups_per_sim(particle_id(1)))
   particle_id(3) = min(1+floor(real(n_particles_per_group(particle_id(2),particle_id(1)),&
                    kind=8)*rnd3(3)),n_particles_per_group(particle_id(2),particle_id(1)))
   select type(p=>sims_particles(particle_id(1))%groups(particle_id(2))%particles(particle_id(3)))
     type is (particle_gc_relativistic)
     x_light = p%x; call compute_mhd_fields_gc_cyl(p,E_field,b_field,normB,gradB,curlb,dbdt)
     call compute_gc_velocity_cartesian(p,sims_particles(particle_id(1))%groups(particle_id(2))%mass,&
     E_field,b_field,normB,gradB,curlb,dbdt,real_param(1:3))
     real_param(1:3) = real_param(1:3)/norm2(real_param(1:3))
     real_param(4) = compute_relativistic_factor(p,&
     sims_particles(particle_id(1))%groups(particle_id(2))%mass,normB)
     p_perp = sqrt(2d0*sims_particles(particle_id(1))%groups(particle_id(2))%mass*p%p(2)*normB); 
     pitch_angle = atan2(p_perp,abs(p%p(1)));
     if(pitch_angle.lt.0d0) pitch_angle = TWOPI + pitch_angle
     ifail = .false.
   end select  
  enddo

  !> generate a velocity based reference system mu = pitch_angle + 1/gamma
  call random_number(rnd3);
  call vectors_to_orthonormal_basis(real_param(1:3),rnd3,tang,nor,binor)
  cospitch_angle = cos(pitch_angle); sinpitch_angle = sin(pitch_angle);
  costheta = cos((asin(1d0/real_param(4))+pitch_angle))
  !> identify in_cone shaded points
  do ii=1,n_shaded_points
    call random_number(rnd3)
    x_shaded = sample_uniform_cone([costheta,cospitch_angle],rnd3,real_param(1:3),x_light,length_shadowed)
    cosmu  = dot_product((x_shaded-x_light),real_param(1:3))/norm2(x_shaded-x_light)
    cospsi = cosmu*cospitch_angle + sinpitch_angle*sqrt(1d0-cosmu**2) 
    in_cone_points(ii) = vertex_sol%check_angles_shaded_in_emission_zone(&
    n_angles,[sqrt(1d0-cospsi**2),cospsi],n_int_param,n_real_param_2,int_param,[real_param(4)])
  enddo
  !> identify out_cone shaded points
  do ii=1,n_shaded_points
    call random_number(rnd3)
    x_shaded = sample_uniform_sphere_corona_rcosphi(length_shadowed**3,[-1d0,costheta],[0d0,TWOPI],rnd3)
    sintheta = sqrt(1d0-x_shaded(2)**2); x_shaded = x_light+x_shaded(1)*&
    (tang*x_shaded(2)+sintheta*(nor*cos(x_shaded(3))+binor*sin(x_shaded(3))))
    cosmu  = dot_product((x_shaded-x_light),real_param(1:3))/norm2(x_shaded-x_light)
    cospsi = cosmu*cospitch_angle + sinpitch_angle*sqrt(1d0-cosmu**2) 
    out_cone_points(ii) = vertex_sol%check_angles_shaded_in_emission_zone(&
    n_angles,[sqrt(1d0-cospsi**2),cospsi],n_int_param,n_real_param_2,int_param,[real_param(4)])
  enddo
  !> check results
  call assert_true(all(in_cone_points),&
  "Error check shaded angles in gyroaverage synchrotron cone: false negative detected!")
  call assert_true(all(.not.out_cone_points),&
  "Error check shaded angles in gyroaverage synchrotron cone: false positive detected!")
end subroutine test_check_shaded_angles_in_gyroaverage_synchrotron_cone

!> test the gyroaverage synchrotron radiation irradiance and directionality functions
subroutine test_gyroaverage_synchrotron_irradiance_directionality_funct()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer                                                :: ii,jj,kk
  integer,dimension(n_times_sol)                         :: n_particles_time
  real*8,dimension(n_x,n_gc_RE_max,n_times_sol)          :: x_cart_loc
  real*8,dimension(n_properties,n_gc_RE_max,n_times_sol) :: properties_loc
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shaded_points_per_particle) :: dir_fun,irradiance
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shaded_points_per_particle) :: dir_fun_sol
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shaded_points_per_particle) :: irradiance_sol
  !>  
  !> initialise the gyroaverage synchrotron lights
  x_cart_loc = 0d0; properties_loc = 0d0; n_particles_time = 0;
  do ii=1,n_times_sol
    n_particles_time(ii) = sum(n_active_particles_sol(:,ii))
    x_cart_loc(:,1:n_particles_time(ii),ii) = x_cart_sol(:,1:n_particles_time(ii),ii)
    properties_loc(:,1:n_particles_time(ii),ii) = properties_sol(:,1:n_particles_time(ii),ii)
  enddo
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles)
  !> compute and check the irradiance and directionality functions
  do kk=1,n_times_sol
    do jj=1,n_particles_time(kk)
      do ii=1,n_shaded_points_per_particle
        call compute_gyroavg_synch_directionality_irradiance(normB_sol(jj,kk),&
        rel_fact_parallel_sol(jj,kk),weight_sol(jj,kk),1d0,mass_RE,x_shadowed(:,ii,jj,kk),&
        x_cart_loc(:,jj,kk),properties_loc(:,jj,kk),dir_fun_sol(:,:,ii),irradiance_sol(:,:,ii))
        call vertex_sol%directionality_funct(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),dir_fun(:,:,ii))
        call vertex_sol%spectral_irradiance(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),irradiance(:,:,ii))
      enddo
      !> check the solution via relative error
      call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
      n_shaded_points_per_particle,dir_fun_sol,dir_fun,tol2_real8,&
      "Error gyroaverage synchrotron directionality function: directionality function mismatch!")
      call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
      n_shaded_points_per_particle,irradiance_sol,irradiance,tol3_real8,&
      "Error gyroaverage synchrotron directionality function: irradiance mismatch!")
    enddo
  enddo
end subroutine test_gyroaverage_synchrotron_irradiance_directionality_funct

!> test the gyroaverage synchrotron radiation irradiance and directionality functions
!> using taskloop parallelisation
subroutine test_gyroaverage_synchrotron_irradiance_dir_funct_taskloop()
  use mod_assert_equals_tools, only: assert_equals_rel_error
  implicit none
  !> variables
  integer                                                :: ii,jj,kk
  integer,dimension(n_times_sol)                         :: n_particles_time
  real*8,dimension(n_x,n_gc_RE_max,n_times_sol)          :: x_cart_loc
  real*8,dimension(n_properties,n_gc_RE_max,n_times_sol) :: properties_loc
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shaded_points_per_particle) :: dir_fun,irradiance
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shaded_points_per_particle) :: dir_fun_sol
  real*8,dimension(spectrum%n_points,spectrum%n_spectra,n_shaded_points_per_particle) :: irradiance_sol
 
  !> initialise the gyroaverage synchrotron lights
  x_cart_loc = 0d0; properties_loc = 0d0; n_particles_time = 0;
  do ii=1,n_times_sol
    n_particles_time(ii) = sum(n_active_particles_sol(:,ii))
    x_cart_loc(:,1:n_particles_time(ii),ii) = x_cart_sol(:,1:n_particles_time(ii),ii)
    properties_loc(:,1:n_particles_time(ii),ii) = properties_sol(:,1:n_particles_time(ii),ii)
  enddo
  call vertex_sol%init_lights_from_particles(n_times_sol,sims_particles)
  !> compute and check the irradiance and directionality functions
  do kk=1,n_times_sol
    do jj=1,n_particles_time(kk)
      do ii=1,n_shaded_points_per_particle
        call compute_gyroavg_synch_directionality_irradiance(normB_sol(jj,kk),&
        rel_fact_parallel_sol(jj,kk),weight_sol(jj,kk),1d0,mass_RE,x_shadowed(:,ii,jj,kk),&
        x_cart_loc(:,jj,kk),properties_loc(:,jj,kk),dir_fun_sol(:,:,ii),irradiance_sol(:,:,ii))
        !$omp parallel default(shared) firstprivate(ii,jj,kk)
        !$omp single
        call vertex_sol%directionality_funct(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),dir_fun(:,:,ii))
        call vertex_sol%spectral_irradiance(spectrum,kk,jj,x_shadowed(:,ii,jj,kk),irradiance(:,:,ii))
        !$omp end single
        !$omp end parallel
      enddo
      !> check the solution via relative error
      call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
      n_shaded_points_per_particle,dir_fun_sol,dir_fun,tol2_real8,&
      "Error gyroaverage synchrotron directionality function taskloop: directionality function mismatch!")
      call assert_equals_rel_error(spectrum%n_points,spectrum%n_spectra,&
      n_shaded_points_per_particle,irradiance_sol,irradiance,tol3_real8,&
      "Error gyroaverage synchrotron directionality function taskloop: irradiance mismatch!")
    enddo
  enddo
end subroutine test_gyroaverage_synchrotron_irradiance_dir_funct_taskloop

!> Test the generation of particles gc relativistic from gyroaverage synchrotron lights
subroutine test_particle_from_gyroaverage_synchrotron_lights()
  use mod_particle_types,        only: particle_base,particle_gc_relativistic
  use mod_particle_assert_equal, only: assert_equal_particle
  implicit none
  !> variables
  class(particle_base),dimension(:),allocatable :: particles_test
  real*8,dimension(n_mhd_sol)                   :: mhd_fields
  real*8,dimension(15)                          :: tol_real8_particle_from_light
  integer                                       :: kk,jj,ii,counter
  !> set tolerance for particle comparison
  tol_real8_particle_from_light    = tol_real8
  tol_real8_particle_from_light(3) = tol3_real8 !< tolerance weight
  tol_real8_particle_from_light(4) = tol4_real8!< tolerance for parallel momentum/magnetic moment
  do kk=1,n_times_sol
    counter = 0;
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is (particle_gc_relativistic)
        !> allocate particle lists
        allocate(particle_gc_relativistic::particles_test(n_particles_per_group(jj,kk)))
        !> loop on the particles
        do ii=1,n_particles_per_group(jj,kk)
          !> override comparison if particle is null
          if(p_list(ii)%i_elm.le.0) then
            select type (p_out=>particles_test(ii))
              type is (particle_gc_relativistic)
              p_out = p_list(ii); 
            end select
            cycle
          endif
          counter = counter + 1
          !> compute the analytical MHD fields for computing the GC velocity
          call compute_mhd_fields_gc_cyl(p_list(ii),mhd_fields(1:3),mhd_fields(4:6),&
          mhd_fields(7),mhd_fields(8:10),mhd_fields(11:13),mhd_fields(14:16))
          !> compute the light properties
          call vertex_sol%compute_light_properties(counter,kk,p_list(ii),&
          sims_particles(kk)%groups(jj)%mass,mhd_fields)
          !> reconstruct particles the property
          call vertex_sol%compute_particle_from_light(sims_particles(kk)%fields,&
          counter,kk,sims_particles(kk)%groups(jj)%mass,particles_test(ii))
#ifdef UNIT_TESTS_AFIELDS
          particles_test(ii)%i_elm = p_list(ii)%i_elm; particles_test(ii)%st = p_list(ii)%st;
#endif
          particles_test(ii)%i_life  = p_list(ii)%i_life 
          particles_test(ii)%t_birth = p_list(ii)%t_birth
          select type (p_out=>particles_test(ii))
            type is (particle_gc_relativistic)
            p_out%q = (p_list(ii)%q/abs(p_list(ii)%q))*p_out%q;
            p_out%p(1) = (p_list(ii)%p(1)/abs(p_list(ii)%p(1)))*p_out%p(1);
          end select
        enddo
        !> check solutions
        call assert_equal_particle(n_particles_per_group(jj,kk),p_list,&
        particles_test,tol_real8_in=tol_real8_particle_from_light)
        deallocate(particles_test); !< deallocate particle lists
      end select
    enddo
  enddo
end subroutine test_particle_from_gyroaverage_synchrotron_lights

!> Tools -------------------------------------------------------------
!> compute the gyroaverage synchrotron lights irradiance and directionality functions
subroutine compute_gyroavg_synch_directionality_irradiance(normB,rel_fact_parallel,&
weight,charge,mass,x_shaded,x_light,properties,dir_funct,irradiance)
  use constants,   only: PI,TWOPI,EL_CHG,EPS_ZERO,ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_besselk, only: besselk
  implicit none
  !> inputs:
  real*8,intent(in)                         :: weight,normB,rel_fact_parallel,charge,mass
  real*8,dimension(n_x),intent(in)          :: x_shaded,x_light 
  real*8,dimension(n_properties),intent(in) :: properties 
  !> outputs:
  real*8,dimension(spectrum%n_points,spectrum%n_spectra),intent(out) :: dir_funct,irradiance
  !> variables:
  integer :: ii,jj
  integer,dimension(0) :: int_param
  real*8 :: mu_angle,thetap,psi_angle,xi,fact1,fact2,fact3,besselk13,besselk23
  !> initialisations
  irradiance = 0d0; dir_funct = 0d0;
  mu_angle = acos(dot_product((x_shaded-x_light)/norm2(x_shaded-x_light),properties(1:3)))
  thetap = atan2(properties(7),properties(6))
  if(thetap.lt.0d0) thetap = TWOPI + thetap 
  psi_angle = mu_angle - thetap
  !> check if the point is shaded by the synchrotron light
  if(.not.vertex_sol%check_angles_shaded_in_emission_zone(2,[sin(psi_angle),&
  cos(psi_angle)],0,1,int_param,[properties(4)])) return
  fact1 = (((1d0-properties(5)*cos(psi_angle))/(properties(5)*cos(psi_angle)))**2)*&
          (1d0-properties(5)*properties(6)*cos(mu_angle))
  fact2 = (5d-1*properties(5)*cos(psi_angle)*(sin(psi_angle)**2))/(1d0-properties(5)*cos(psi_angle))
  fact3 = (properties(4)**3)*sqrt(((1d0-properties(5)*cos(psi_angle))**3)/&
          (5d-1*properties(5)*cos(psi_angle)))
  !> compute the irradiance
  do ii=1,spectrum%n_spectra
    do jj=1,spectrum%n_points
      xi = fact3*(properties(8)/spectrum%points(jj,ii))
      call besselk(onethird,xi,besselk13); call besselk(twothirds,xi,besselk23);
       irradiance(jj,ii) = ((weight*(9d0*((charge*EL_CHG)**5)*(properties(5)**2)*(properties(4)**9)*&
       (normB**3)))/(2.56d2*(PI**3)*EPS_ZERO*(SPEED_OF_LIGHT**2)*(rel_fact_parallel**2)*&
       ((mass*ATOMIC_MASS_UNIT)**3)))*(((properties(8)/spectrum%points(jj,ii))**4)*fact1*&
       ((besselk23**2)+fact2*(besselk13**2)))
    enddo
  enddo
  dir_funct = irradiance/properties(10)
end subroutine compute_gyroavg_synch_directionality_irradiance

!> generate shadowed points for each light. The shadowed point position
!> is taken within the emission cone of the synchrotron radiation.
!> The emission cone half angle is approximated with 
!> cos(theta) = 1/(2*rel_fact). The relativistic factor is defined as:
!> rel_fact = sqrt(1+(p/(mass*c))**2) with p the total gc momentum and
!> c the speed of light
subroutine compute_x_shadowed_gc()
  use constants,    only: TWOPI
  use mod_sampling, only: sample_uniform_cone
  implicit none
  !> variables
  integer :: ii,jj,kk,n_gc_time
  real*8  :: cos_half_angle
  real*8,dimension(n_x) :: v_gc_dir,x_gc,rng
  !> generate shadowed points
  do kk=1,n_times_sol
    n_gc_time = sum(n_active_particles_sol(:,kk))
    do jj=1,n_gc_time
      x_gc = x_cart_sol(:,jj,kk)
      v_gc_dir = properties_sol(1:3,jj,kk)
      !> angle to be sampled 1/rel_fact + theta
      cos_half_angle = atan2(properties_sol(7,jj,kk),properties_sol(6,jj,kk))
      if(cos_half_angle.lt.0d0) cos_half_angle = cos_half_angle + TWOPI
      cos_half_angle = cos((cos_half_angle + (1d0/properties_sol(4,jj,kk))))
      do ii=1,n_shaded_points_per_particle
        call random_number(rng)
        x_shadowed(:,ii,jj,kk) = sample_uniform_cone([cos_half_angle,&
        properties_sol(6,jj,kk)],rng,v_gc_dir,x_gc,length_shadowed)
      enddo
    enddo
  enddo
end subroutine compute_x_shadowed_gc

!> compute & fill the gyroaverage synchrotron light positions and properties
subroutine compute_gyroavg_synch_x_properties_ana()
  use mod_coordinate_transforms, only: cylindrical_to_cartesian
  use mod_particle_types,        only: particle_gc_relativistic
  implicit none
  !> variables
  integer :: ii,jj,kk,counter
  !> initialise positions and properties arrays
  x_cart_sol = 0d0; properties_sol = 0d0;
  !> fill property table
  do kk=1,n_times_sol
    counter = 0
    do jj=1,n_groups_per_sim(kk)
      select type(p_list=>sims_particles(kk)%groups(jj)%particles)
        type is (particle_gc_relativistic)
        do ii=1,n_particles_per_group(jj,kk)
          if(p_list(ii)%i_elm.eq.0) cycle
          counter = counter + 1
          x_cart_sol(:,counter,kk) = cylindrical_to_cartesian(p_list(ii)%x)
          call compute_gyroavg_synch_properties_ana_1p(p_list(ii),&
          sims_particles(kk)%groups(jj)%mass,properties_sol(:,counter,kk),&
          rel_fact_parallel_sol(counter,kk),normB_sol(counter,kk),&
          weight_sol(counter,kk))
        enddo
      end select
    enddo
  enddo
end subroutine compute_gyroavg_synch_x_properties_ana

!> compute gyroaverage synchrotron light properties using the analytical
!> tokamak-like MHD fields for one guiding center
subroutine compute_gyroavg_synch_properties_ana_1p(gc_in,mass,properties,&
rel_fact_parallel,normB,weight)
  use constants,                      only: PI,TWOPI,EL_CHG,EPS_ZERO
  use constants,                      only: ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
  use mod_particle_types,             only: particle_gc_relativistic
  implicit none
  !> inputs:
  type(particle_gc_relativistic),intent(in)  :: gc_in
  real*8,intent(in)                          :: mass 
  !> outputs:
  real*8,intent(out)                        :: rel_fact_parallel,normB,weight
  real*8,dimension(n_properties),intent(out) :: properties
  !> variables:
  real*8                :: thetap,charge
  real*8,dimension(n_x) :: E_field,b_field,gradB,curlb,dbdt,v_gc_cart
  
  !> initialisation
  properties = 0d0
  !> compute the analytical MHD fields for computing the GC velocity
  call compute_mhd_fields_gc_cyl(gc_in,E_field,b_field,normB,gradB,curlb,dbdt)
  !> compute the guiding center velocity
  call compute_gc_velocity_cartesian(gc_in,mass,E_field,b_field,&
  normB,gradB,curlb,dbdt,v_gc_cart)
  properties(1:3) = v_gc_cart/norm2(v_gc_cart)
  !> compute the relativistic factor
  properties(4) = sqrt(1d0 + ((gc_in%p(1)/(mass*SPEED_OF_LIGHT))**2) +&
             ((2d0*gc_in%p(2)*normB)/(mass*(SPEED_OF_LIGHT**2))))
  rel_fact_parallel = gc_in%p(1)/(mass*SPEED_OF_LIGHT*properties(4))
  rel_fact_parallel = 1d0/(sqrt(1d0-rel_fact_parallel**2))
  !> compute the beta
  properties(5) = sqrt(1d0 - (1d0/(properties(4)**2)));
  !> compute the pitch angle
  thetap = atan2(sqrt(2d0*mass*gc_in%p(2)*normB),abs(gc_in%p(1)));
  if(thetap.lt.0.d0) thetap = TWOPI+thetap;
  properties(6:7) = abs((/cos(thetap),sin(thetap)/))
  !> critical wavelength
  charge = real(abs(gc_in%q),kind=8)*EL_CHG
  properties(8) = (4d0*PI*mass*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*rel_fact_parallel)/&
                  (3d0*charge*normB*(properties(4)**2))
  !> compute the directionality function intensity
  properties(9) = (27d0*charge*normB*(properties(4)**7))/(128d0*(PI**2)*mass*&
                  ATOMIC_MASS_UNIT*SPEED_OF_LIGHT*(rel_fact_parallel**4))
  !> compute the synchrotron power normalisation
  weight = gc_in%weight
  properties(10) = (weight*(charge**4)*((normB*properties(4)*properties(5)*rel_fact_parallel)**2))/&
                   (6d0*PI*EPS_ZERO*SPEED_OF_LIGHT*((mass*ATOMIC_MASS_UNIT)**2))
end subroutine compute_gyroavg_synch_properties_ana_1p

!> compute the MHD fields for GC in cylindrical coordinates
subroutine compute_mhd_fields_gc_cyl(particle,E_field,b_field,normB,gradB,curlb,dbdt)
  use mod_particle_types,             only: particle_base
  use mod_coordinate_transforms,      only: vector_cylindrical_to_cartesian
  use mod_particle_common_test_tools, only: compute_test_E_B_normB_gradB_curlb_Dbdt_fields
  implicit none
  !> inputs:
  class(particle_base),intent(in)   :: particle
  !> outputs:
  real*8,intent(out)                :: normB
  real*8,dimension(n_x),intent(out) :: E_field,b_field,gradB,curlb,dbdt
  !> compute the MHD fields and derivatives in cartesian coordinates
  call compute_test_E_B_normB_gradB_curlb_Dbdt_fields(particle%x,E_field,&
       b_field,normB,gradB,curlb,dbdt)
end subroutine compute_mhd_fields_gc_cyl

!> compute the relativistic gc velocity in cartesian coordinates
subroutine compute_gc_velocity_cartesian(gc_in,mass,E_field,b_field,&
  normB,gradB,curlb,dbdt,v_gc_cart)
  use mod_coordinate_transforms, only: cylindrical_to_cartesian_velocity
  use mod_particle_types,        only: particle_gc_relativistic
  use mod_gc_relativistic,       only: compute_relativistic_gc_rhs
  implicit none
  !> inputs:
  type(particle_gc_relativistic),intent(in) :: gc_in
  real*8,intent(in)                         :: mass,normB
  real*8,dimension(n_x),intent(in)  :: E_field,b_field,gradB,curlb,dbdt
  !> outputs:
  real*8,dimension(n_x),intent(out) :: v_gc_cart
  !> variables:
  real*8,dimension(4) :: xdot_gc_cart
  xdot_gc_cart = compute_relativistic_gc_rhs(int(gc_in%q,kind=4),mass,gc_in%p(2),&
  gc_in%x(1),gc_in%p(1),normB,E_field,b_field,gradB,curlb,dbdt)
  v_gc_cart = cylindrical_to_cartesian_velocity(&
  gc_in%x(1),gc_in%x(3),xdot_gc_cart(1:3)) 
end subroutine compute_gc_velocity_cartesian

!>--------------------------------------------------------------------
end module mod_gyroaverage_synchrotron_light_dist_vertices_test
 
