!> Program to put JOREK data into IMAS IDSs
program jorek2_IDS

#ifdef USE_IMAS
  use ids_schemas 
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get, ids_put, ids_put_slice

  use mod_jorek2IMAS 
  use constants
  use nodes_elements
  use data_structure
  use mod_import_restart
  use phys_module, only : rst_format, n_tor_restart, central_density, central_mass, t_start
  use mod_impurity, only: init_imp_adas 
  use basis_at_gaussian, only: initialise_basis
  
  implicit none
  
  character(len=200):: user, database
  character(len=64) :: file_name, name_proj
  integer :: shot_number, run_number, i_begin, i_end, i_step
  integer :: ierr, idx, stat_mhd, stat_core, stat_rad, stat_eq, n_grid, stat
  logical :: first_step, file_exists, rad_only_projections_h5
  logical :: export_MHD, export_radiation, export_core_profiles, export_equilibrium
  real*8  :: rho0, fact_time, time_SI

  type(ids_mhd), target   :: mhd_ids
  type(ids_equilibrium)   :: equilibrium_ids
  type(ids_core_profiles) :: core_profiles_ids
  type(ids_radiation)     :: radiation_ids

  namelist /imas_params/ shot_number, run_number, user, database, i_begin, i_end, &
                         export_mhd, export_radiation, export_core_profiles, n_grid, &
                         export_equilibrium, rad_only_projections_h5

  ! --- Necessary initialization ------------------
  ! --- Initialize mode and mode_type arrays
  call det_modes()
  
  ! --- Initialize the basis functions
  call initialise_basis()
  
  ! --- Preset namelist input parameters
  call preset_parameters()

  ! --- Read input file
  call initialise_parameters(0, "__NO_FILENAME__")

  ! --- Initialize ADAS
#if (defined WITH_Neutrals) || (defined WITH_Impurities)
    ! --- Read ADAS data and generate coronal equilibrium if needed
    call init_imp_adas(0)
#else
    if (use_imp_adas .and. (nimp_bg(1) > 0.d0)) then
      call init_imp_adas(0)
    endif
#endif
  ! -----------------------------------------------
  
  ! --- Preset parameters for this program
  database    = 'test'                    !< Name of the database to export the results
  shot_number = 111112;   run_number=1;   
  i_begin     = 0                         !< Starting restart file index
  i_end       = 99999                     !< Ending restart file index
  export_MHD           = .true.
  export_radiation     = .false.
  export_core_profiles = .false. 
  export_equilibrium   = .false.
  rad_only_projections_h5 = .false.    !< use only *.h5 projection files for radiation IDS (single jorek_restart.h5 still needed)
  n_grid               = 100              !< Number of points used for 1D and 2D profiles  

  call getenv('USER',user)
  
  ! --- Read parameters from namelist file 'imas.nml' if it exists
  open(42, file='imas.nml', action='read', status='old', iostat=ierr)
  if ( ierr == 0 ) then
    write(*,*) 'Reading parameters from imas.nml namelist.'
    read(42, imas_params)
    close(42)
  end if

  ! --- Try to open shot and number if it exists
  call imas_open_env( 'ids', shot_number,run_number,idx,user,database,'3',stat)! 3 is the database version  
  
  if (stat /= 0) then  ! --- Create a new shot if it doesn't exist
    write(*,*) '  Shot/run number did not exist, creating new one...'
    call imas_create_env('ids',shot_number,run_number, 0,0,idx,user,database,'3') 
  endif

  if (export_radiation)  allocate( aux_node_list )

  ! --- Time normalization
  rho0       = central_density * 1.d20 * central_mass * mass_proton
  fact_time  = sqrt( mu_zero * rho0 )

  first_step = .true.

  ! --- Loop over
  do i_step = i_begin, i_end
 
    ! --- Cycle when required files don't exist 
    if (rad_only_projections_h5 .and. export_radiation) then
      write(name_proj,'(a,i9.9,a)') 'projections000.', i_step, '.h5'  ! This formatting should be improved
      inquire (file=trim(name_proj), exist=file_exists)
      if (.not. file_exists) cycle
      file_name = 'jorek_restart' 
      fact_time = 1.d0  ! Projection times are already in SI units...
    else
      if (.not. restart_file_exists(i_step)) cycle
      write(file_name,'(a,i5.5)')   'jorek', i_step
      write(name_proj,'(a,i5.5,a)') 'projections', i_step, '.h5'
    endif

    ! --- Import restart file
    write(*,*)
    write(*,'(a,i9.9,a)') '#################### STEP ', i_step, ' ####################'
    write(*,*)
    call import_restart(node_list, element_list, file_name, rst_format, ierr)
    if (ierr /=0 ) then
       write(*,*) '  Could not read the JOREK restart file'
       stop
    endif
    time_SI = t_start * fact_time

    ! --- Check whether jorek2_IDS has been compiled with a sufficient number of harmonics
    if (n_tor < n_tor_restart) then
      write(*,*) '  You must compile jorek2_IDS at least with n_tor = ', n_tor_restart
      write(*,*) '  Otherwise you are cutting information for the IDSs'
      stop
    endif

    ! --- Fill and export an MHD IDS
    if (export_mhd)  call fill_mhd_IDS(first_step, time_SI, mhd_ids)  

    ! --- Fill and export a core_profiles IDS
    if (export_core_profiles)  call fill_core_profiles_IDS(first_step, time_SI, core_profiles_ids, n_grid)  

    ! --- Fill and export an equilibrium IDS
    if (export_equilibrium)  call fill_equilibrium_IDS(first_step, time_SI, equilibrium_ids, n_grid)

    ! --- Fill and export a radiation IDS
    if (export_radiation) then
      call import_hdf5_restart_aux(aux_node_list, name_proj, rst_format, ierr)
      if (ierr /= 0) then
        write(*,*) ' Could not open projections file where radiation is stored'
        stop
      endif
      call fill_radiation_IDS(first_step, t_start*fact_time, radiation_ids)  
    endif

    stat_mhd = 1;   stat_core = 1;   stat_rad = 1;   stat_eq = 1;

    ! --- Put IDSs into database
    if (first_step) then  
      if (export_mhd)              call ids_put(idx,'mhd',mhd_ids,stat_mhd)
      if (export_core_profiles)    call ids_put(idx,'core_profiles',core_profiles_ids,stat_core)
      if (export_equilibrium)      call ids_put(idx,'equilibrium',equilibrium_ids,stat_eq)
      if (export_radiation)        call ids_put(idx,'radiation',radiation_ids,stat_rad)
    else
      if (export_mhd)              call ids_put_slice(idx,'mhd',mhd_ids,stat_mhd)
      if (export_core_profiles)    call ids_put_slice(idx,'core_profiles',core_profiles_ids,stat_core)
      if (export_equilibrium)      call ids_put_slice(idx,'equilibrium',equilibrium_ids,stat_eq)
      if (export_radiation)        call ids_put_slice(idx,'radiation',radiation_ids,stat_rad)
    endif

    if (export_mhd           .and. (stat_mhd==0 ))   write(*,*) '    MHD IDS exported'
    if (export_core_profiles .and. (stat_core==0))   write(*,*) '    Core profiles IDS exported'
    if (export_equilibrium   .and. (stat_eq==0  ))   write(*,*) '    Equlibrium IDS exported'
    if (export_radiation     .and. (stat_rad==0 ))   write(*,*) '    Radiation IDS exported'

    if (export_mhd           .and. (stat_mhd/=0 ))   write(*,*) '    Problem saving MHD IDS'
    if (export_core_profiles .and. (stat_core/=0))   write(*,*) '    Problem saving Core profiles IDS'
    if (export_equilibrium   .and. (stat_eq/=0  ))   write(*,*) '    Problem saving Equlibrium IDS'
    if (export_radiation     .and. (stat_rad/=0 ))   write(*,*) '    Problem saving Radiation IDS'

    first_step = .false.

  enddo

  call imas_close(idx) 
 
#else

  write(*,*) 'Error: jorek2_IDS must be compiled with IMAS (USE_IMAS=1)'
  write(*,*) 'You will also have to load the IMAS module, in case you have not done it'
  write(*,*) '     module load IMAS                                                   '

#endif

end program jorek2_IDS
