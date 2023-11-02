!> Example for initialising the a relativistic gc population
!> (runaway electrons) and running the simulation
!> the current density initialisation is used for generating
!> the initial gc population
program re_gc_current_density_initialisation
use constants,                only: TWOPI,PI
use phys_module,              only: xcase,xpoint,central_density,central_mass
use data_structure,           only: type_bnd_node_list,type_bnd_element_list
use mod_boundary,             only: boundary_from_grid
use mod_expression,           only: exprs_all_int,init_expr,exprs,SI_UNITS
use mod_integrals3D,          only: int3d_new
use mod_boundary,             only: boundary_from_grid
use equil_info
use mod_random_seed
use mod_particle_diagnostics, only: write_particle_diagnostics
use mod_io_actions,           only: write_action,read_action
use particle_tracer
implicit none
!> Variable declarations --------------------------------------------------------------------
type(pcg32_rng)                  :: rng_pcg32
type(event)                      :: field_reader,particle_reader
type(type_bnd_node_list)         :: bnd_node_list
type(type_bnd_element_list)      :: bnd_elm_list
type(write_particle_diagnostics) :: write_diag
type(write_action)               :: write_particles,write_restart_particles 
logical                          :: append_diag,abort_at_last_mhd_restart,file_exist
integer                          :: ii,n_variables,n_particles,n_groups,restart_id,n_gc_variables
integer                          :: fraction_of_modes,n_int_pdf_param,n_real_pdf_param
integer                          :: ifail,n_int_weight_param,n_real_weight_param
integer                          :: n_int_gdf_param,n_real_gdf_param,n_int_pdf_to_part_coord_param
integer                          :: n_real_pdf_to_part_coord_param,gc_integrator_type
integer,dimension(2)             :: diag_list
integer,dimension(:),allocatable :: int_pdf_param,int_weight_param,int_gdf_param
integer,dimension(:),allocatable :: int_pdf_to_part_coord_param
real*8                           :: diag_step,write_step,mass,charge,sim_time_interval
real*8                           :: pdf_upper_bound,gdf_upper_bound,sup_pdf_safety_factor
real*8                           :: starttime,endtime
real*8,dimension(1)              :: gc_timesteps
real*8,dimension(2)              :: Rbox,Zbox,Rbound,Zbound,phibound,Ekinbound
real*8,dimension(2)              :: Pbound,Pitchbound,Chargebound
real*8,dimension(4)              :: rk45_gc_tolerances
real*8,dimension(6,2)            :: phase_space_bounds
real*8,dimension(:),allocatable  :: real_pdf_to_part_coord_param
real*8,dimension(:),allocatable  :: real_pdf_param,real_weight_param,real_gdf_param
real*8,dimension(:),allocatable  :: DUMMY_REAL_ARRAY
character(len=3)                 :: hdf5ext 
character(len=20)                :: diag_filename,particle_restart_basename,restart_basename
procedure(real_f),pointer           :: pdf_to_use         => NULL()
procedure(real_f),pointer           :: weight_to_use      => NULL()
procedure(real_f),pointer           :: gdf_to_use         => NULL()
procedure(real_arr_inout_s),pointer :: gdf_sampler_to_use => NULL()
procedure(part_inout_s),pointer     :: pdf_to_part_coord  => NULL()

!> MPI and groups initialisation ------------------------------------------------------------
n_groups = 1; call sim%initialize(num_groups=n_groups);
!> Define inputs ----------------------------------------------------------------------------
hdf5ext                    = '.h5'            !< extension of the hdf5file
diag_filename              = 'part_diag'   !< particle diagnostics filename
particle_restart_basename  = 'part_restart'   !< particle restart file basename
restart_basename           = 'jorek'          !< JOREK restart file basename
append_diag                = .true.           !< if true, append diagnostics in the same file
abort_at_last_mhd_restart  = .false.          !< abort when the last simulation restart is reached
restart_id                 = 0                !< first JOREK restart file number
gc_integrator_type         = 1                !< type of gc integration strategy to be used
                                              !< default: Runge-Kutta 45 fixed time step
                                              !< =1 Runge-Kutta 45 with feedback time step control
                                              !< =2 Runge-Kutta 45 with feedforward time step control
n_variables                = 6                !< phase space dimensionality
n_gc_variables             = 4                !< number of variables use for integrating a gc particle 
n_particles                = 32               !< number of particles per task
fraction_of_modes          = 1                !< fraction of JOREK modes to be used for interpolation
diag_list                  = [1,6]            !< total energy, canonical toroidal momentum
gc_timesteps               = [1d-10]          !< initial gc timestep in SI units
sim_time_interval          = 1d-5             !< total simulation time interval
diag_step                  = 1d-8             !< time step at which the diagnostic files are written
write_step                 = 1d-4             !< time step at which the particle restart file is written
mass                       = 5.48579909065d-4 !< electron mass in AMU
charge                     = -1d0             !< electron charge / EL_CHG
Rbound                     = [0.d0,9.99d2]
Zbound                     = [-9.99d2,9.99d2]
Phibound                   = [0d0,TWOPI]
Ekinbound                  = [2d7-1d4,2d7+1d4] !< kinetic energy in eV
Pitchbound                 = [PI-2.95d-1,PI]
Chargebound                = -1.d0
sup_pdf_safety_factor      = 1d0
rk45_gc_tolerances         = [1d-4,1d-4,1d-4,1d1] !< RK45 tolerances for time step control
                                                  !< 1: R, 2: Z, 3: phi, 4: ppar
!> Initialisation ---------------------------------------------------------------------------
write(*,*) "Simulate runaway electrons as gc: started!"
write(*,*) "... initialise simulation parameters"
!> initialise simulation events
inquire(file=(trim(diag_filename)//trim(hdf5ext)),exist=file_exist)
if(file_exist.and.(sim%my_id.eq.0)) call system("rm "//(trim(diag_filename)//trim(hdf5ext)))
write_diag = write_particle_diagnostics(filename=(trim(diag_filename)//trim(hdf5ext)),&
             append=append_diag,only=diag_list)
write_particles = write_action(basename=trim(particle_restart_basename))
write_restart_particles = write_action(filename=(trim(particle_restart_basename)//trim(hdf5ext)))
field_reader = event(read_jorek_fields_interp_linear(basename=trim(restart_basename),i=restart_id,\
               mode_divisor=fraction_of_modes,stop_at_end=abort_at_last_mhd_restart))
call with(sim,field_reader)
events = [field_reader,event(write_diag,start=sim%time,step=diag_step),\
         event(write_particles,start=sim%time,step=write_step),\
         event(stop_action(),start=sim%time+sim_time_interval)]
!> update equilibrium state
if(sim%my_id.eq.0) call boundary_from_grid(sim%fields%node_list,sim%fields%element_list,&
bnd_node_list,bnd_elm_list,.false.)
call broadcast_boundary(sim%my_id,bnd_elm_list,bnd_node_list)
call update_equil_state(sim%my_id,sim%fields%node_list,sim%fields%element_list,bnd_elm_list,xpoint,xcase)
!> check if particle restart file exists
inquire(file=trim(particle_restart_basename)//trim(hdf5ext),exist=file_exist)
if(.not.file_exist) then
  !> initialise particle structure
  sim%groups(1:n_groups)%mass = mass
  do ii=1,n_groups
    allocate(particle_gc_relativistic::sim%groups(ii)%particles(n_particles))
  enddo
  !> difine phase space domain
  call domain_bounding_box(sim%fields%node_list,sim%fields%element_list,Rbox(1),Rbox(2),Zbox(1),Zbox(2))
  if(Rbox(1).ge.Rbound(1)) Rbound(1) = Rbox(1)
  if((Rbox(2).lt.Rbound(2)).and.((Rbox(2)-Rbound(1)).gt.0.d0)) Rbound(2) = Rbox(2)
  if((Zbox(1).gt.0.d0).and.(Zbound(1).ge.Zbox(1))) Zbound(1) = Zbox(1)
  if((Zbox(1).lt.0.d0).and.(Zbound(1).lt.Zbox(1))) Zbound(1) = Zbox(1)
  if((Zbox(2).gt.0.d0).and.(Zbound(2).ge.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
  if((Zbox(2).lt.0.d0).and.(Zbound(2).lt.Zbox(2)).and.((Zbox(2)-Zbound(1)).gt.0.d0)) Zbound(2) = Zbox(2)
  Pbound = mass*SPEED_OF_LIGHT*sqrt(((EL_CHG*Ekinbound/(ATOMIC_MASS_UNIT*mass*SPEED_OF_LIGHT**2))+1.d0)**2-1.d0)
  phase_space_bounds(:,1) = [Rbound(1),Zbound(1),Phibound(1),Pbound(1),Pitchbound(1),charge]
  phase_space_bounds(:,2) = [Rbound(2),Zbound(2),Phibound(2),Pbound(2),Pitchbound(2),charge]
  !> define particle PDF
  n_real_pdf_param = 3; allocate(real_pdf_param(n_real_pdf_param));
  real_pdf_param   = [1.d0,mass,sup_pdf_safety_factor]
  n_int_pdf_param  = 1; allocate(int_pdf_param(n_int_pdf_param));
  n_real_weight_param = 3; allocate(real_weight_param(n_real_weight_param));
  int_pdf_param(1)    = sim%my_id
  n_int_pdf_to_part_coord_param  = 0
  n_real_pdf_to_part_coord_param = 1
  allocate(real_pdf_to_part_coord_param(n_real_pdf_to_part_coord_param))
  real_pdf_to_part_coord_param(1) = mass
  pdf_to_use          => pdf_current_density_uniform_phase
  weight_to_use       => particle_weight_current_density_uniform_phase
  pdf_to_part_coord   => spherical_p_cartesian_q_to_relativistic_gc
  pdf_upper_bound     = sup_pdf_current_density_uniform_phase(n_variables,&
  phase_space_bounds(:,1),phase_space_bounds(:,2),sim%fields,&
  n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
  !> compute the integral of the current density in the volume
  allocate(DUMMY_REAL_ARRAY(n_real_weight_param+1)); call init_expr;
  call int3d_new(sim%my_id,sim%fields%node_list,sim%fields%element_list,bnd_node_list,bnd_elm_list,&
  exprs('int3d_jR_tot',1,exprs_all_int%n_coord,exprs_all_int),DUMMY_REAL_ARRAY,SI_UNITS)
  real_weight_param = [DUMMY_REAL_ARRAY(2),real(n_particles,kind=8),sim%groups(1)%mass];
  deallocate(DUMMY_REAL_ARRAY);
  call MPI_Bcast(real_weight_param,n_real_weight_param,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  !> define the sampler and samper distribution
  n_real_gdf_param = 0; n_int_gdf_param = 0;
  gdf_to_use         => gdf_uniform_phase
  gdf_sampler_to_use => gdf_uniform_sampler
  gdf_upper_bound    = sup_gdf_uniform_phase(n_variables,&
  phase_space_bounds(1:n_variables,1),phase_space_bounds(1:n_variables,2), &
  n_real_pdf_param,real_pdf_param,n_int_pdf_param,int_pdf_param)
endif
write(*,*) "... initialise simulation parameters: completed!"

!> Test particle initialisation -------------------------------------------------------------
write(*,*) "... initialising guiding center in phase space ..."
if(file_exist) then
  write(*,*) "Reading particle restart: ",trim(particle_restart_basename)//trim(hdf5ext)
  particle_reader = event(read_action(filename=trim(particle_restart_basename)//trim(hdf5ext)))
  call with(sim,particle_reader)
else
  call initialise_particles_in_phase_space(n_variables,sim%groups(1)%particles,sim%fields,rng_pcg32,&
  pdf_to_use,weight_to_use,gdf_to_use,gdf_sampler_to_use,pdf_upper_bound,gdf_upper_bound,&
  pdf_to_part_coord,sim%groups(1)%mass,sim%time,phase_space_bounds,n_real_pdf_param,real_pdf_param,&
  n_int_pdf_param,int_pdf_param,n_real_weight_param,real_weight_param,n_int_weight_param,&
  int_weight_param,n_real_gdf_param,real_gdf_param,n_int_gdf_param,int_gdf_param,&
  n_real_pdf_to_part_coord_param,real_pdf_to_part_coord_param,&
  n_int_pdf_to_part_coord_param,int_pdf_to_part_coord_param)
  call with(sim,events,at=sim%time); !< execute event at the simulation initial state
endif
write(*,*) "... initialising guiding center in phase space: completed!"

!> Test particle initialisation -------------------------------------------------------------
starttime = MPI_Wtime()
call integrate_particles(n_gc_variables,gc_integrator_type,gc_timesteps,rk45_gc_tolerances,events,sim)
endtime = MPI_WTIME()
write(*,*) "my_id: ",sim%my_id," total pusher wall time [s]: ",endtime-starttime

!> Write a particle restart file ------------------------------------------------------------
call with(sim,write_restart_particles)

!> Clean-up ---------------------------------------------------------------------------------
if(allocated(real_pdf_param))               deallocate(real_pdf_param);
if(allocated(int_pdf_param))                deallocate(int_pdf_param);
if(allocated(int_weight_param))             deallocate(int_weight_param);
if(allocated(real_weight_param))            deallocate(real_weight_param);
if(allocated(int_pdf_to_part_coord_param))  deallocate(int_pdf_to_part_coord_param);
if(allocated(real_pdf_to_part_coord_param)) deallocate(real_pdf_to_part_coord_param);
pdf_to_use => NULL(); weight_to_use      => NULL();
gdf_to_use => NULL(); gdf_sampler_to_use => NULL();
call sim%finalize()
write(*,*) "Simulate runaway electrons as gc: completed!"

contains

!> particle integrator
!> inputs:
!>   n_var:           (integer) number of variables of the Runge-Kutta integrator
!>   integrator_type: (integer) type of integrator to be used:
!>                    default: Runge-Kutta 45 fixed time step
!>                    =1 Runge-Kutta 45 feedback time step control
!>                    =2 Runge-Kutta 45 feed forward time step control based on MHD gradients
!>   time_steps:      (integer)(n_groups) initial time step per group
!>   rk_tols:         (real8)(n_var) tolerances of the runge-kutta integrator
!>   events:          (n_real8)(n_events)(n_events) events to execute at each target time
!>   sim:             (particle_sim) jorek particle simulation structure
!> outputs:
!>   events:          (n_real8)(n_events)(n_events) events to execute at each target time
!>   sim:             (particle_sim) jorek particle simulation structure
subroutine integrate_particles(n_var,integrator_type,time_steps,rk_tols,events,sim)
implicit none
!> inputs-outputs:
type(particle_sim),intent(inout)                   :: sim
type(event),dimension(:),allocatable,intent(inout) :: events
!> inputs:
integer,intent(in)                            :: n_var,integrator_type
real*8,dimension(n_var),intent(in)            :: rk_tols
real*8,dimension(size(sim%groups)),intent(in) :: time_steps 
!> variables:
real*8 :: target_time,field_time
!> loop for all the simulation time
do while (.not.sim%stop_now)
  target_time = next_event_at(sim,events);
  write(*,*) "Integrating target time: ",target_time
  if(integrator_type.eq.1) then
    call push_guiding_center_loop_rk45_feedback(n_var,size(sim%groups),time_steps,target_time,rk_tols,sim)
  elseif(integrator_type.eq.2) then
    call push_guiding_center_loop_rk45_feedforward(n_var,size(sim%groups),time_steps,target_time,sim)
  else
    call push_guiding_center_loop_rk45_fix(n_var,size(sim%groups),time_steps,target_time,sim) 
  endif
  call with(sim,events,at=sim%time);
  write(*,*) "Integration of the target time: ",target_time," completed!"
enddo
end subroutine integrate_particles   


!> guiding center loop using the fixed time step runge kutta
!> inputs:
!>   n_var:       (integer) number of variables of the Runge-Kutta integrator
!>   n_groups:    (integer) number of particle groups
!>   time_steps:  (integer)(n_groups) initial time step per group
!>   target_time: (real8) time at which the particle integration is stopped
!>   sim:         (particle_sim) jorek particle simulation structure
!> outputs:
!>   sim: (particle_sim) jorek particle simulation structure
subroutine push_guiding_center_loop_rk45_fix(n_var,n_groups,time_steps,target_time,sim)
use mod_gc_relativistic, only: runge_kutta_fixed_dt_gc_push_jorek
implicit none
!> inputs-outputs:
type(particle_sim),intent(inout)      :: sim
!> inputs:
integer,intent(in)                    :: n_var,n_groups
real*8,intent(in)                     :: target_time
real*8,dimension(n_groups),intent(in) :: time_steps
!> variables:
integer :: ii,jj
real*8  :: local_time,local_dt
!> loop on particles
do jj=1,n_groups
  !$omp parallel do default(none) firstprivate(jj,target_time,time_steps) &
  !$omp private(ii,local_time,local_dt) shared(sim)
  do ii=1,size(sim%groups(jj)%particles)
    select type (gc=>sim%groups(jj)%particles(ii))
      type is (particle_gc_relativistic)
      local_time = sim%time; local_dt = time_steps(jj)
      do while(((target_time-local_time).gt.0d0).and.(gc%i_elm.gt.0))
        call runge_kutta_fixed_dt_gc_push_jorek(&
        sim%fields,local_time,local_dt,sim%groups(jj)%mass,gc)
        local_time = local_time+local_dt
        local_dt   = min(time_steps(jj),target_time-local_time)
      enddo
    end select
  enddo
  !$omp end parallel do
enddo
sim%time = target_time
end subroutine  push_guiding_center_loop_rk45_fix

!> guiding center pushing loop using the feedforward control Runge-Kutta
!> inputs:
!>   n_var:       (integer) number of variables of the Runge-Kutta integrator
!>   n_groups:    (integer) number of particle groups
!>   time_steps:  (integer)(n_groups) initial time step per group
!>   target_time: (real8) time at which the particle integration is stopped
!>   sim:         (particle_sim) jorek particle simulation structure
!> outputs:
!>   sim: (particle_sim) jorek particle simulation structure
subroutine push_guiding_center_loop_rk45_feedforward(n_var,n_groups,&
time_steps,target_time,sim)
use mod_gc_relativistic, only: runge_kutta_adapt_dt_gc_push_jorek
implicit none
!> inputs-outputs:
type(particle_sim),intent(inout)      :: sim
!> inputs:
integer,intent(in)                    :: n_var,n_groups
real*8,intent(in)                     :: target_time
real*8,dimension(n_groups),intent(in) :: time_steps
!> variables:
integer :: ii,jj
real*8  :: local_time,local_dt
!> loop on particles
do jj=1,n_groups
  !$omp parallel do default(none) firstprivate(jj,target_time,time_steps) &
  !$omp private(local_time,local_dt,ii) shared(sim)
  do ii=1,size(sim%groups(jj)%particles)
    select type (gc => sim%groups(jj)%particles(ii))
      type is (particle_gc_relativistic)
        local_time = sim%time; local_dt = time_steps(jj);
        do while(((target_time-local_time).gt.0d0).and.(gc%i_elm.gt.0))
          call runge_kutta_adapt_dt_gc_push_jorek(sim%fields,local_time,local_dt,&
          target_time,sim%groups(jj)%mass,gc)
          local_time = local_time+local_dt
         if(local_dt.le.0d0) local_dt = time_steps(jj)
        enddo
    end select
  enddo
  !$omp end parallel do
enddo
sim%time = target_time
end subroutine push_guiding_center_loop_rk45_feedforward

!> guiding center pushing loop using the feedback control Runge-Kutta
!> inputs:
!>   n_var:       (integer) number of variables of the Runge-Kutta integrator
!>   n_groups:    (integer) number of particle groups
!>   time_steps:  (integer)(n_groups) initial time step per group
!>   target_time: (real8) time at which the particle integration is stopped
!>   rk_tols:     (real8)(n_var) tolerances of the runge-kutta integrator
!>   sim:         (particle_sim) jorek particle simulation structure
!> outputs:
!>   sim: (particle_sim) jorek particle simulation structure
subroutine push_guiding_center_loop_rk45_feedback(n_var,n_groups,time_steps,\
target_time,rk_tols,sim)
use mod_gc_relativistic, only: runge_kutta_error_control_dt_gc_push_jorek
implicit none
!> inputs-outputs:
type(particle_sim),intent(inout)      :: sim
!> inputs:
integer,intent(in)                    :: n_var,n_groups
real*8,dimension(n_groups),intent(in) :: time_steps
real*8,intent(in)                     :: target_time
real*8,dimension(n_var),intent(in)    :: rk_tols
!> variables:
integer :: jj,ii
real*8  :: local_time,local_dt,dt_try
!> loop on the particles
do jj=1,n_groups
  !$omp parallel do default(none) firstprivate(jj,target_time,time_steps,rk_tols) &
  !$omp private(ii,local_time,local_dt,dt_try) shared(sim)
  do ii=1,size(sim%groups(jj)%particles)
    select type (gc => sim%groups(jj)%particles(ii))
      type is (particle_gc_relativistic)
      local_time = sim%time; local_dt = time_steps(jj);
      do while (((target_time-local_time).gt.0d0).and.(gc%i_elm.gt.0))
        call runge_kutta_error_control_dt_gc_push_jorek(sim%fields,rk_tols,&
        local_time,local_dt,target_time,sim%groups(jj)%mass,dt_try,gc)
        local_time = local_time + local_dt
        if(local_dt.le.0d0) local_dt = time_steps(jj)
      enddo
    end select
  enddo
  !$omp end parallel do
enddo
sim%time = target_time
end subroutine push_guiding_center_loop_rk45_feedback

!> method used for transforming the momentum space from spherical
!> coordinates (p,pitch) and cartesian coordinates for the
!> charge state into particle gc relativistic coordinates
!> the order of the variables in a sample are:
!> 1: R ,2: Z, 3: phi, 4: momentum, 5: pitch angle, 6: charge.
!> we note that a spherical jacobian is still used for the 
!> momentum space despite that only the total momentum and 
!> the pitch angle are provided
!> inputs: 
!>   p_inout:      (particle_base) particle to be initialised
!>   n_x:          (integer) size of the phase space sample
!>   x:            (real8)(n_x) phase space sample in cylindrical
!>                 space-momentum coordinates cartesian charge coordinates
!>   time:         (real8) time of the simulation
!>   fields:       (fields_base) JOREK MHD fields
!>   n_real_param: (integer) number of real parameters: 1
!>   real_param:   (real8)(n_real_param) real parameters: 1:mass
!>   n_int_param:  (integer) number of integer parameters: 0
!>   int_param:    (integer)(n_real_param) integer parameters: empty
!> outputs:
!>   p_inout: (particle_base) initialised particle
subroutine spherical_p_cartesian_q_to_relativistic_gc(p_inout,&
n_x,x,time,fields,n_real_param,real_param,n_int_param,int_param)
  use mod_particle_types,        only: particle_base
  use mod_particle_types,        only: particle_gc_relativistic
  use mod_fields,                only: fields_base
  implicit none
  !> Inputs-Outputs:
  class(particle_base),intent(inout) :: p_inout
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: n_x,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(n_x),intent(in)            :: x
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> variables
  real*8              :: psi,U
  real*8,dimension(3) :: B_field,E_field
  select type (p=>p_inout)
  type is (particle_gc_relativistic)
    call fields%calc_EBpsiU(time,p%i_elm,p%st,p%x(3),&
    E_field,B_field,psi,U);
    p%p = x(4)*[cos(x(5)),&
    (x(4)*((sin(x(5)))**2))/(2d0*real_param(1)*norm2(B_field))]
    p%q = int(x(6),kind=1)
  end select
end subroutine spherical_p_cartesian_q_to_relativistic_gc

!> Phase space distribution based on the plasma current density
!> and uniform phase space distribution. The momentum distribution
!> is considered uniform for relativistic particle hence, the 
!> appearance of the relativistic factor.
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   time:         (real8) physical time at wich the particle is sampled
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution coefficient (not used)
!>                 2) particle mass in AMU
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   pdf: (real8) value of the probability density 
function pdf_current_density_uniform_phase(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(pdf)
  use constants,          only: MU_ZERO,EL_CHG,SPEED_OF_LIGHT,PI
  use mod_model_settings, only: var_zj
  use mod_interp,         only: interp_PRZ
  use mod_fields,         only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: pdf
  !> Variables:
  real*8 :: DUMMY_DOUBLE_1,DUMMY_DOUBLE_2
  real*8,dimension(1) :: jphi
  !> interpolate the jorek toroidal current density at the particle position
  call interp_PRZ(fields%node_list,fields%element_list,i_elm,[var_zj],1,&
  st(1),st(2),x(3),jphi,DUMMY_DOUBLE_1,DUMMY_DOUBLE_2)
  !> compute the pdf
  DUMMY_DOUBLE_1 = sqrt((x_max(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0)
  DUMMY_DOUBLE_2 = sqrt((x_min(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0)
  pdf = ((DUMMY_DOUBLE_1**3)-3.d0*DUMMY_DOUBLE_1) - ((DUMMY_DOUBLE_2**3)-3.d0*DUMMY_DOUBLE_2);
  pdf = pdf*(cos(x_min(5))**2 - cos(x_max(5))**2)
  pdf =(-3.d0*jphi(1))/(pdf*x(6)*PI*EL_CHG*MU_ZERO*(real_param(2)**3)*(SPEED_OF_LIGHT**4)*x(1))
end function pdf_current_density_uniform_phase

!> Upper bound phase space distribution based on the plasma current density
!> and uniform phase space distribution. The momentum distribution
!> is considered uniform for relativistic particle hence, the 
!> appearance of the relativistic factor.
!> inputs:
!>   nx:           (integer) number of variables
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   fields:       (fields_base) jorek MHD fields
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution coefficient (not used)
!>                 2) particle mass in AMU
!>                 3) safety factor: must be >1
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!>                 1) mpi rank
!> outputs:
!>   sup_pdf:      (real8) value of the probability density upper bound
function sup_pdf_current_density_uniform_phase(nx,x_min,x_max,fields,&
n_real_param,real_param,n_int_param,int_param) result(sup_pdf)
  use constants,          only: SPEED_OF_LIGHT,EL_CHG,MU_ZERO,PI
  use mod_model_settings, only: n_var,var_zj
  use mod_fields,         only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: sup_pdf
  !> Variables
  integer :: ifail
  real*8  :: density_tot,density_in,density_out
  real*8  :: pressure_tot,pressure_in,pressure_out
  real*8  :: kin_par_tot,kin_par_in,kin_par_out,mom_par_tot,mom_par_in
  real*8  :: mom_par_out
  real*8,dimension(2)     :: zj_minmax
  real*8,dimension(n_var) :: varmin,varmax
  real*8  :: max_pdf,min_pdf,sqrtpovermc2plus1_max,sqrtpovermc2plus1_min
  real*8  :: cos2pitch_max,cos2pitch_min
  !> Evalutate the upper extremum of the pdf
  call Integrals_3D(int_param(1),fields%node_list,fields%element_list,&
  density_tot,density_in,density_out,pressure_tot,pressure_in,pressure_out,&
  kin_par_tot,kin_par_in,kin_par_out,mom_par_tot,mom_par_in,mom_par_out,varmin,varmax)
  zj_minmax = [varmin(var_zj),varmax(var_zj)]
  call MPI_Bcast(zj_minmax,2,MPI_REAL8,0,MPI_COMM_WORLD,ifail)
  !> compute the maximum and the minimum of the pdf
  sqrtpovermc2plus1_max = sqrt((x_max(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0);
  cos2pitch_max = cos(x_max(5))**2;
  sqrtpovermc2plus1_min = sqrt((x_min(4)**2)/((real_param(2)*SPEED_OF_LIGHT)**2)+1.d0)
  cos2pitch_min = cos(x_min(5))**2
  max_pdf = ((sqrtpovermc2plus1_max**3)-3.d0*sqrtpovermc2plus1_max) - &
            ((sqrtpovermc2plus1_min**3)-3.d0*sqrtpovermc2plus1_min);
  max_pdf = max_pdf*(cos2pitch_min - cos2pitch_max)
  max_pdf =(-3.d0*real_param(3))/(max_pdf*PI*x_min(6)*EL_CHG*MU_ZERO*(real_param(2)**3)*&
           (SPEED_OF_LIGHT**4)*x_min(1));
  min_pdf = max_pdf*zj_minmax(1); max_pdf = max_pdf*zj_minmax(2);
  !> check which between min_pdf and max_pdf has the maximum absolute value
  if(abs(max_pdf).ge.abs(min_pdf)) then
    sup_pdf = max_pdf
  else
    sup_pdf = min_pdf
  endif
end function sup_pdf_current_density_uniform_phase

!> Compute the gdf used for sampling the particle coordinates in 
!> phase space. The gdf used here is a uniform distribution in
!> cylindrical coordinates for the spatial coordinates and 
!> uniform spherical distribution for the momentum coordinates.
!> The gyro-angle coordinate is not considered here despite
!> that the spherical topology is retained
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution weight (normally mass/volume)
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   gdf: (real8) value of the sampler probability density 
function gdf_uniform_phase(nx,x,st,time,i_elm,fields,x_min,x_max,&
n_real_param,real_param,n_int_param,int_param) result(gdf)
  use constants,  only: PI
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,i_elm,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(2),intent(in)              :: st
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: gdf
  !> Evalutate pdf
  gdf = 3.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*PI)
  if(n_real_param.gt.0) gdf = real_param(1)*gdf
end function gdf_uniform_phase

!> Upper bound of the uniform phase space sampler distribution
!> The gyro-angle coordinate is not considered here despite
!> that the spherical topology is retained
!> inputs:
!>   nx:           (integer) number of variables
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution weight (normally mass/volume)
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   sup_gdf: (real8) value of the sampler probability density upper bound
function sup_gdf_uniform_phase(nx,x_min,x_max,n_real_param,real_param,&
n_int_param,int_param) result(sup_gdf)
  use constants,  only: PI
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: sup_gdf
  !> Evalutate the upper extremum of the pdf
  sup_gdf = 3.d0/((x_max(1)**2-x_min(1)**2)*(x_max(2)-x_min(2))*&
  (x_max(3)-x_min(3))*(x_max(4)**3-x_min(4)**3)*&
  (cos(x_min(5))-cos(x_max(5)))*PI)
end function sup_gdf_uniform_phase

!> GDF uniform sampler generating particle positions in phase space
!> The gyro-angle coordinate is not considered here despite
!> that the spherical topology is retained
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random numbers in [0,1)
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   fields:       (fields_base) jorek MHD fields
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   n_real_param: (integer) N# of real input parameters of the pdf
!>   real_param:   (real8)(n_real_param) real pdf parameters
!>                 1) pdf distribution weight (normally mass/volume)
!>   n_int_param:  (integer) N# of integer input parameters of the pdf
!>   int_param:    (integer)(n_int_param) integer pdf parameters
!> outputs:
!>   i_elm:        (integer) updated jorek mesh element number
!>   st:           (real8)(2) updated local mesh coordinates
!>   x:            (real8)(nx) particle position to accept
subroutine gdf_uniform_sampler(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param)
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  integer,intent(in)                          :: nx,n_real_param
  integer,intent(in)                          :: n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(nx),intent(in)             :: x_min,x_max
  class(fields_base),intent(in)               :: fields
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Inputs-Outputs:
  integer,intent(inout)                       :: i_elm
  real*8,dimension(2),intent(inout)           :: st
  real*8,dimension(nx),intent(inout)          :: x
  !> Compute new particle position in phase space
  x(1) = sqrt(x_min(1)**2 + (x_max(1)**2 - x_min(1)**2)*x(1))
  x(2:3) = x_min(2:3) + (x_max(2:3)-x_min(2:3))*x(2:3)
  x(4) = (x_min(4)**3 + (x_max(4)**3-x_min(4)**3)*x(4))**(1d0/3d0)
  x(5) = acos(cos(x_min(5))-(cos(x_max(5))-cos(x_min(5)))*x(5))
  x(6) = x_min(6) + (x_max(6)-x_min(6))*x(6)
  !> find particle RZ coordinates
  call find_RZ(fields%node_list,fields%element_list,x(1),x(2),&
  x(1),x(2),i_elm,st(1),st(2),ifail)
end subroutine gdf_uniform_sampler

!> Method used for computing the weight of each particle. Given that the pdf
!> is directly sampled via accept-reject method, the particle weight is 
!> considered uniform for all particles and equal to the total number of
!> physical particles (from the plasma current) divided the total number of
!> simulated markers. A uniform distribution in the velocity space is used. 
!> inputs:
!>   nx:           (integer) number of variables
!>   x:            (real8)(nx) random state to accept
!>   i_elm:        (integer) jorek mesh element number
!>   st:           (real8)(2) local mesh coordinates
!>   time:         (real8) physical time at wich the particle is sampled
!>   x_min:        (real8)(nx) lower bound of the phase space interval
!>   x_max:        (real8)(nx) upper bound of the phase space interval
!>   fields:       (fields_base) jorek MHD fields
!>   n_real_param: (integer) N# of real input parameters
!>   real_param:   (real8)(n_real_param) real weight parameters
!>                 1) 3D integral of the plasma density in SI units
!>                 2) total number of particles
!>                 3) mass in AMU 
!>   n_int_param:  (integer) N# of integer input parameters
!>   int_param:    (integer)(n_int_param) integer weight parameters
!> outputs:
!>   weight:       (real) particle weight
function particle_weight_current_density_uniform_phase(nx,x,st,time,i_elm,fields,&
x_min,x_max,n_real_param,real_param,n_int_param,int_param) result(weight)
  use constants,  only: EL_CHG,SPEED_OF_LIGHT
  use mod_fields, only: fields_base
  implicit none
  !> Inputs:
  class(fields_base),intent(in)               :: fields
  integer,intent(in)                          :: nx,i_elm,n_real_param,n_int_param
  integer,dimension(:),allocatable,intent(in) :: int_param
  real*8,intent(in)                           :: time
  real*8,dimension(2),intent(in)              :: st
  real*8,dimension(nx),intent(in)             :: x,x_min,x_max
  real*8,dimension(:),allocatable,intent(in)  :: real_param
  !> Outputs:
  real*8 :: weight
  !> Variables:
  real*8 :: DUMMY_DOUBLE_1,DUMMY_DOUBLE_2
  !> Compute the particle weight
  DUMMY_DOUBLE_1 = sqrt((x_max(4)**2)/((real_param(3)*SPEED_OF_LIGHT)**2)+1.d0)
  DUMMY_DOUBLE_2 = sqrt((x_min(4)**2)/((real_param(3)*SPEED_OF_LIGHT)**2)+1.d0)
  weight = ((DUMMY_DOUBLE_1**3)-3.d0*DUMMY_DOUBLE_1) - ((DUMMY_DOUBLE_2**3)-3.d0*DUMMY_DOUBLE_2)
  weight = weight*(cos(x_min(5)+cos(x_max(5))))
  weight = (2.d0*real_param(1)*(x_max(4)**3 - x_min(4)**3))/&
  (weight*real_param(2)*x(6)*EL_CHG*(real_param(3)**3)*(SPEED_OF_LIGHT**4))
end function particle_weight_current_density_uniform_phase

end program re_gc_current_density_initialisation
