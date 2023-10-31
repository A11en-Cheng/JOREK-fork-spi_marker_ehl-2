!>  Module that contais functions to fill IDSs from JOREK data
module mod_jorek2IMAS

#ifdef USE_IMAS
  use ids_schemas !, only: ids_equilibrium
  use ids_routines, only: imas_open_env, &
     imas_create_env, imas_close, ids_get

  use mod_parameters 
  use mod_new_diag
  use data_structure
  use nodes_elements
  use constants
  use mod_expression, only: exprs, exprs_all_int
  use exec_commands,  only: average, expr_list, clean_up, step_imported, qprofile, &
                            zeroD_quantities, separatrix, rectangle
  use parse_commands, only: type_command
  use settings,       only: set_setting
  
  implicit none
  
 
  public
 
  ! Transform to COCOS convention 8 --> 11
  real*8 :: fact_psi =  2.d0 * PI
  real*8 :: fact_Ip  = -1.d0

  contains


  subroutine fill_mhd_IDS(first_step, mhd_ids)  

    use phys_module, only : t_start, F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass

    implicit none

    ! --- External parameters
    logical,            intent(in)           :: first_step   ! is this the first step?
    type(ids_mhd), target,     intent(inout) :: mhd_ids

   
    ! --- Local parameters 
    integer    :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_T, fact_time, fact_v, fact_zj, rho0, fact_phi, fact_rho, fact_w
    
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: num_nodes
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
  
    if (first_step) then
      ! --- Put the grid in GGD
      allocate( mhd_ids%grid_ggd(n_grid) )
      grid => mhd_ids%grid_ggd(grid_ind)
      call grid2ggd( grid, node_list, element_list )
    endif

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * mass_proton
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )
  
    fact_time =  sqrt_mu0_rho0 
    fact_v    =  1.d0 /  sqrt_mu0_rho0 
    fact_w    = -1.d0 /  sqrt_mu0_rho0      ! Transform for COCOS convention of toroidal direction (anti-clockwise) 
    fact_phi  = -1.d0 /  sqrt_mu0_rho0 * F0 ! COCOS convection: F0 depends on phi direction
    fact_zj   = -1.d0 / mu_zero * fact_Ip   ! Last sign due to COCOS transformation
    fact_rho  =  rho0 
    fact_T    =  1.d0 / ( EL_CHG * mu_zero * central_density * 1.d20 )   
  
  
    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  mhd_ids%time(n_slice) )
    allocate(  mhd_ids%ggd(n_slice ) )

    mhd_ids%ids_properties%homogeneous_time = 1
  
    mhd_ids%time(i_slice)     = t_start * fact_time 
    mhd_ids%ggd(i_slice)%time = t_start * fact_time

    ! --- Fill MHD data
    do i=1, n_var 
  
      ! --- Poloidal magnetic flux
      if (variable_names(i) == 'Psi') then      
        allocate( mhd_ids%ggd(i_slice)%psi(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%psi(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_psi, grid_ind, grid_sub_ind, fact_psi )
      endif
  
      ! --- Electrostatic potential 
      if (variable_names(i) == 'u') then      
        allocate( mhd_ids%ggd(i_slice)%phi_potential(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%phi_potential(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_u, grid_ind, grid_sub_ind, fact_phi )
      endif
  
      ! --- Toroidal current density * R
      if (variable_names(i) == 'zj') then      
        allocate( mhd_ids%ggd(i_slice)%j_tor_r(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%j_tor_r(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_zj, grid_ind, grid_sub_ind, fact_zj )
      endif
  
      ! --- Toroidal vorticity / R 
      if (variable_names(i) == 'omega') then      
        allocate( mhd_ids%ggd(i_slice)%vorticity_over_r(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%vorticity_over_r(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_w, grid_ind, grid_sub_ind, fact_w )
      endif
  
      ! --- Mass density
      if (variable_names(i) == 'rho') then      
        allocate( mhd_ids%ggd(i_slice)%mass_density(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%mass_density(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_rho, grid_ind, grid_sub_ind, fact_rho )
      endif
  
      ! --- Total temperature 
      if (variable_names(i) == 'T') then      
        ! --- Te
        allocate( mhd_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
  
        ! --- Ti
        allocate( mhd_ids%ggd(i_slice)%t_i_average(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_T, grid_ind, grid_sub_ind, fact_T*0.5d0 )
      endif
  
      ! --- Ion temperature
      if (variable_names(i) == 'T_i') then      
        allocate( mhd_ids%ggd(i_slice)%t_i_average(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%t_i_average(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_Ti, grid_ind, grid_sub_ind, fact_T )
      endif
  
      ! --- Electron temperature
      if (variable_names(i) == 'T_e') then      
        allocate( mhd_ids%ggd(i_slice)%electrons%temperature(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%electrons%temperature(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_Te, grid_ind, grid_sub_ind, fact_T )
      endif
  
      ! --- Parallel velocity
      if (variable_names(i) == 'v_par') then      
        allocate( mhd_ids%ggd(i_slice)%velocity_parallel_over_b_field(n_grid_sub))
        ggd_scalar => mhd_ids%ggd(i_slice)%velocity_parallel_over_b_field(grid_sub_ind)
        call fill_Bezier_coefficients( ggd_scalar, node_list, var_vpar, grid_ind, grid_sub_ind, fact_v )
      endif
  
    enddo

  end subroutine fill_mhd_IDS

 






  subroutine fill_radiation_IDS(first_step, radiation_ids)  

    use phys_module, only : t_start, F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    type(ids_radiation),  intent(inout) :: radiation_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_time, rho0, fact_rad
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: num_nodes
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index
  
    if (first_step) then
      ! --- Put the grid in GGD
      allocate( radiation_ids%grid_ggd(n_grid) )
      grid => radiation_ids%grid_ggd(grid_ind)
      call grid2ggd( grid, node_list, element_list )
    endif
 
    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * mass_proton
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )
    sqrt_mu0_over_rho0 = sqrt( mu_zero / rho0 )

    fact_rad = 1.d0 / ( (gamma-1.d0) * MU_ZERO * sqrt_mu0_rho0 )
    fact_time =  sqrt_mu0_rho0 

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  radiation_ids%time(n_slice) )

    radiation_ids%ids_properties%homogeneous_time = 1
    allocate( radiation_ids%process(1))   ! --- 1 type of radiation
    allocate( radiation_ids%process(1)%ggd(n_slice) )
  
    radiation_ids%time(i_slice)                = t_start * fact_time 
    radiation_ids%process(1)%ggd(i_slice)%time = t_start * fact_time
 
    ! --- Fill radiation data 
    var_rad = 2
  
    allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1))
    allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(n_grid_sub))
    allocate( radiation_ids%process(1)%ggd(i_slice)%ion(1)%label(1) )  
    allocate( radiation_ids%process(1)%identifier%name(1) )

    radiation_ids%process(1)%identifier%name  = "Line radiation"
    radiation_ids%process(1)%identifier%index = 10

    radiation_ids%process(1)%ggd(i_slice)%ion(1)%label = imp_type(index_main_imp) 
  
    ggd_scalar => radiation_ids%process(1)%ggd(i_slice)%ion(1)%emissivity(grid_sub_ind)
    call fill_Bezier_coefficients( ggd_scalar, aux_node_list, var_rad, grid_ind, grid_sub_ind, fact_rad )

  end subroutine fill_radiation_IDS






  subroutine fill_core_profiles_IDS(first_step, core_profiles_ids, n_grid)  

    use phys_module, only : t_start, F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,      intent(in) :: first_step   ! is this the first step?
    integer,      intent(in) :: n_grid       ! Number of flux surfaces to compute average
    type(ids_core_profiles),   intent(inout) :: core_profiles_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_time, rho0
    real*8, allocatable :: result(:,:), q_prof(:), rho_tor(:)
    character(10)       :: str
    type(type_command)  :: command_tmp
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1

    allocate( core_profiles_ids%profiles_1d(n_slice) )
    allocate( core_profiles_ids%time(n_slice) )

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * mass_proton
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )

    fact_time =  sqrt_mu0_rho0 
    
    core_profiles_ids%ids_properties%homogeneous_time = 1    
    core_profiles_ids%time(i_slice) = t_start * fact_time 

    ! --- Call expressions and do a flux average
    step_imported = .true.
  
  ! --- Preset namelist input parameters for jorek2_postproc
    if (first_step) then 
      call init_new_diag(.false.)
      write(str, '(I0)') n_grid
      call set_setting('units',           '1',     ierr, 'Calculate quantities in which units (0=JOREK, 1=SI)')
      call set_setting('loop_units',      '1',     ierr, 'Use which units for time-loops (0=JOREK, 1=SI)'     )
      call set_setting('linepoints',      '200',   ierr, 'Number of points along a line e.g. for pol_line'    )
      call set_setting('tor_points',      '200',   ierr, 'Number of toroidal points e.g. for tor_line'        )
      call set_setting('surfaces',         str,    ierr, 'number for flux surfaces e.g. for qprofile'         )
      call set_setting('nsmallsteps',     '3',     ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nmaxsteps',       '2500',  ierr, 'numerical parameter for field line tracing'         )
      call set_setting('deltaphi',        '0.3',   ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_min',   '0.001', ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_max',   '0.999', ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nTht',            '32',    ierr, 'numerical parameter for field line tracing'         )
    endif

    ! --- Get average and q-profile
    command_tmp%n_args = 0
    call clean_up()
    expr_list = exprs((/'Psi_N', 'T_i', 'T_e', 'ne', 'pres', 'Phi', 'eta_T', &
                        'Jpar', 'E_||', 'Er', 'vpar', 'Vtheta_i', 'Vstar_i', 'rho', 'Psi'/), 15)
    call average(command_tmp, first_step==.true., ierr, result, .true.)
    call clean_up()
    call qprofile(command_tmp, first_step==.true., ierr, q_prof)
    ! --- Correct first and last points
    q_prof(1)      = q_prof(2)        + (q_prof(2)-q_prof(3))
    q_prof(n_grid) = q_prof(n_grid-1) + (q_prof(n_grid-1)-q_prof(n_grid-2))

    ! --- Some allocations
    allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1) )

    ! --- Fill expressions in IDSs
    do i_exp=1, expr_list%n_expr

      ! --- Psi_N
      if (expr_list%expr(i_exp)%name=='Psi_N') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%grid%rho_pol_norm(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%grid%psi_magnetic_axis = ES%Psi_axis * fact_psi
        core_profiles_ids%profiles_1d(i_slice)%grid%psi_boundary      = ES%Psi_bnd  * fact_psi
        core_profiles_ids%profiles_1d(i_slice)%grid%rho_pol_norm(:)   = sqrt(result(:,i_exp))
      endif

      ! --- Psi
      if (expr_list%expr(i_exp)%name=='Psi') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%grid%psi(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%grid%psi(:)   = result(:,i_exp) * fact_psi
      endif

      ! --- Ion temperature
      if (expr_list%expr(i_exp)%name=='T_i') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%t_i_average(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%t_i_average(:) = result(:,i_exp)
      endif

      ! --- Electron temperature
      if (expr_list%expr(i_exp)%name=='T_e') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%electrons%temperature(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%electrons%temperature(:) = result(:,i_exp)
      endif

      ! --- Electron density
      if (expr_list%expr(i_exp)%name=='ne') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%electrons%density(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%electrons%density(:) = result(:,i_exp)
      endif

      ! --- Total pressure
      if (expr_list%expr(i_exp)%name=='pres') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%pressure_thermal(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%pressure_thermal(:) = result(:,i_exp)
      endif

      ! --- Electrostatic potential
      if (expr_list%expr(i_exp)%name=='Phi') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%phi_potential(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%phi_potential(:) = result(:,i_exp)
      endif

      ! --- Parallel conductivity
      if (expr_list%expr(i_exp)%name=='eta_T') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%conductivity_parallel(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%conductivity_parallel(:) = 1.d0 / result(:,i_exp)
      endif

      ! --- Parallel current density
      if (expr_list%expr(i_exp)%name=='Jpar') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%j_total(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%j_total(:) = result(:,i_exp)
      endif

      ! --- Parallel electric field
      if (expr_list%expr(i_exp)%name=='E_||') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%e_field%parallel(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%e_field%parallel(:) = result(:,i_exp)
      endif

      ! --- Radial electric field
      if (expr_list%expr(i_exp)%name=='Er') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%e_field%radial(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%e_field%radial(:) = result(:,i_exp)
      endif

      ! --- Parallel velocity
      if (expr_list%expr(i_exp)%name=='vpar') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1)%velocity%parallel(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(1)%velocity%parallel(:) = result(:,i_exp)
      endif

      ! --- Poloidal velocity
      if (expr_list%expr(i_exp)%name=='Vtheta_i') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1)%velocity%poloidal(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(1)%velocity%poloidal(:) = result(:,i_exp)
      endif

      ! --- Diamagnetic velocity
      if (expr_list%expr(i_exp)%name=='Vstar_i') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1)%velocity%diamagnetic(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(1)%velocity%diamagnetic(:) = result(:,i_exp)
      endif

      ! --- Ion density
      if (expr_list%expr(i_exp)%name=='rho') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1)%density(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(1)%density(:) = result(:,i_exp) / (central_mass*MASS_PROTON)
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1)%element(1) )
        core_profiles_ids%profiles_1d(i_slice)%ion(1)%element(1)%a = central_mass
      endif

    end do
    
    ! --- q-profile
    allocate( core_profiles_ids%profiles_1d(i_slice)%q(n_grid) )
    core_profiles_ids%profiles_1d(i_slice)%q(:) = q_prof(:)

    ! --- Get rho_norm_tor from psi_N and q_profile
    allocate( core_profiles_ids%profiles_1d(i_slice)%grid%rho_tor_norm(n_grid) )
    allocate(rho_tor(n_grid))
    rho_tor(:) = 0.d0
    do i_psi=2, n_grid
      rho_tor(i_psi) = sum(q_prof(1:i_psi)) ! Assuming equidistant psi grid!!
    end do
    core_profiles_ids%profiles_1d(i_slice)%grid%rho_tor_norm(:) = rho_tor(:)/rho_tor(n_grid)

  end subroutine fill_core_profiles_IDS






  subroutine fill_equilibrium_IDS(first_step, equilibrium_ids, n_grid)  

    use phys_module, only : t_start, F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,      intent(in) :: first_step   ! is this the first step?
    integer,      intent(in) :: n_grid       ! Number of flux surfaces to compute average
    type(ids_equilibrium),  intent(inout)  :: equilibrium_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_time, rho0, fact_rad, R_min, Z_min, R_max, Z_max, R_node, Z_node
    real*8, allocatable :: result(:,:), res0D(:), q_prof(:), rho_tor(:), R_sep(:), Z_sep(:)
    real*8, allocatable :: result2D(:,:,:), R_vec(:), Z_vec(:)
    character(30)       :: str
    type(type_command)  :: command_tmp
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1

    allocate( equilibrium_ids%time_slice(n_slice) )
    allocate( equilibrium_ids%time(n_slice) )

    ! --- Normalization factors for IMAS
    rho0               = central_density * 1.d20 * central_mass * mass_proton
    sqrt_mu0_rho0      = sqrt( mu_zero * rho0 )

    fact_time =  sqrt_mu0_rho0 
    
    equilibrium_ids%ids_properties%homogeneous_time = 1    
    equilibrium_ids%time(i_slice) = t_start * fact_time 

    ! --- Call expressions and do a flux average
    step_imported = .true.
  
  ! --- Preset namelist input parameters for jorek2_postproc
    if (first_step) then 
      call init_new_diag(.false.)
      write(str, '(I0)') n_grid
      call set_setting('units',           '1',     ierr, 'Calculate quantities in which units (0=JOREK, 1=SI)')
      call set_setting('loop_units',      '1',     ierr, 'Use which units for time-loops (0=JOREK, 1=SI)'     )
      call set_setting('linepoints',      '200',   ierr, 'Number of points along a line e.g. for pol_line'    )
      call set_setting('tor_points',      '200',   ierr, 'Number of toroidal points e.g. for tor_line'        )
      call set_setting('surfaces',         str,    ierr, 'number for flux surfaces e.g. for qprofile'         )
      call set_setting('nsmallsteps',     '3',     ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nmaxsteps',       '2500',  ierr, 'numerical parameter for field line tracing'         )
      call set_setting('deltaphi',        '0.3',   ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_min',   '0.001', ierr, 'numerical parameter for field line tracing'         )
      call set_setting('rad_range_max',   '0.999', ierr, 'numerical parameter for field line tracing'         )
      call set_setting('nTht',            '32',    ierr, 'numerical parameter for field line tracing'         )
    endif

    ! --- Get average and q-profile
    command_tmp%n_args = 0
    call clean_up()
    expr_list = exprs((/'Psi', 'pres', 'FFprime_loc', 'p_prime_loc', 'Jpar'/), 5)
    call average(command_tmp, first_step==.true., ierr, result, .true.)
    call clean_up()
    call qprofile(command_tmp, first_step==.true., ierr, q_prof)
    ! --- Correct first and last points
    q_prof(1)      = q_prof(2)        + (q_prof(2)-q_prof(3))
    q_prof(n_grid) = q_prof(n_grid-1) + (q_prof(n_grid-1)-q_prof(n_grid-2))

    ! --- Fill profiles
    do i_exp=1, expr_list%n_expr

      ! --- psi
      if (expr_list%expr(i_exp)%name=='Psi') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%psi(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%psi(:)   = result(:,i_exp) * fact_psi
      endif

      ! --- pressure
      if (expr_list%expr(i_exp)%name=='pres') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%pressure(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%pressure(:)   = result(:,i_exp) 
      endif

      ! --- p'
      if (expr_list%expr(i_exp)%name=='p_prime_loc') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%dpressure_dpsi(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%dpressure_dpsi(:)   = result(:,i_exp) / fact_psi
      endif

      ! --- FF'
      if (expr_list%expr(i_exp)%name=='FFprime_loc') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%f_df_dpsi(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%f_df_dpsi(:) = result(:,i_exp) / fact_psi
      endif

      ! --- Parallel current
      if (expr_list%expr(i_exp)%name=='Jpar') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%j_parallel(n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_1d%j_parallel(:) = result(:,i_exp) 
      endif

    end do
    
    ! --- q-profile
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%q(n_grid) )
    equilibrium_ids%time_slice(i_slice)%profiles_1d%q(:) = q_prof(:)

    ! --- Get rho_norm_tor from psi_N and q_profile
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_1d%rho_tor_norm(n_grid) )
    allocate(rho_tor(n_grid))
    rho_tor(:) = 0.d0
    do i_psi=2, n_grid
      rho_tor(i_psi) = sum(q_prof(1:i_psi)) ! Assuming equidistant psi grid!!
    end do
    equilibrium_ids%time_slice(i_slice)%profiles_1d%rho_tor_norm(:) = rho_tor(:)/rho_tor(n_grid)

    ! --- Information about the toroidal field
    equilibrium_ids%vacuum_toroidal_field%r0 = R_geo
    allocate(equilibrium_ids%vacuum_toroidal_field%b0(n_slice))
    equilibrium_ids%vacuum_toroidal_field%b0(i_slice) = F0/R_geo * fact_Ip
    
    ! --- Fill global quantities (call mod_integrals3D)
    equilibrium_ids%time_slice(i_slice)%global_quantities%psi_axis        = ES%Psi_axis * fact_psi
    equilibrium_ids%time_slice(i_slice)%global_quantities%psi_boundary    = ES%Psi_bnd  * fact_psi
    equilibrium_ids%time_slice(i_slice)%global_quantities%magnetic_axis%r = ES%R_axis
    equilibrium_ids%time_slice(i_slice)%global_quantities%magnetic_axis%z = ES%Z_axis

    command_tmp%n_args = 0
    call clean_up()
    call zeroD_quantities(command_tmp, first_step==.true., ierr, res0D)
    
    do i_exp=1, exprs_all_int%n_expr

      ! --- Beta poloidal
      if (exprs_all_int%expr(i_exp)%name=='beta_p') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%beta_pol   = res0D(i_exp)
      endif

      ! --- Beta poloidal
      if (exprs_all_int%expr(i_exp)%name=='beta_t') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%beta_tor   = res0D(i_exp)
      endif

      ! --- Normalized beta
      if (exprs_all_int%expr(i_exp)%name=='beta_n') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%beta_normal = abs(res0D(i_exp))
      endif

      ! --- Total current
      if (exprs_all_int%expr(i_exp)%name=='Ip_tot') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%ip = res0D(i_exp) * fact_Ip
      endif

      ! --- li(3)
      if (exprs_all_int%expr(i_exp)%name=='li3') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%li_3 = res0D(i_exp)
      endif

      ! --- Volume
      if (exprs_all_int%expr(i_exp)%name=='volume') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%volume = res0D(i_exp)
      endif

      ! --- Area of poloidal cross section inside LCFS
      if (exprs_all_int%expr(i_exp)%name=='area') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%area = res0D(i_exp)
      endif

      ! --- Current centre - R
      if (exprs_all_int%expr(i_exp)%name=='R_curr_cent') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%current_centre%r = res0D(i_exp)
      endif

      ! --- Current centre - Z
      if (exprs_all_int%expr(i_exp)%name=='Z_curr_cent') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%current_centre%z = res0D(i_exp)
      endif

      ! --- q_axis
      if (exprs_all_int%expr(i_exp)%name=='q02') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%q_axis = res0D(i_exp)
      endif

      ! --- q_95
      if (exprs_all_int%expr(i_exp)%name=='q95') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%q_95 = res0D(i_exp)
      endif

      ! --- Thermal energy
      if (exprs_all_int%expr(i_exp)%name=='Thermal_tot') then
        equilibrium_ids%time_slice(i_slice)%global_quantities%energy_mhd = res0D(i_exp)
      endif

    end do

    ! --- Shaping parameters, T. Luce, PPCF 55 (2013) 095009, equations (1-6)
    if (ES%limiter_plasma) then
      equilibrium_ids%time_slice(i_slice)%boundary_separatrix%type = 0
    else
      equilibrium_ids%time_slice(i_slice)%boundary_separatrix%type = 1
    endif
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%psi                    = ES%Psi_bnd * fact_psi
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%minor_radius           = ES%LCFS_a
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%elongation             = ES%LCFS_kappa
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%triangularity_upper    = ES%LCFS_deltaU
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%triangularity_lower    = ES%LCFS_deltaL
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%geometric_axis%r       = ES%LCFS_Rgeo
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%geometric_axis%z       = ES%LCFS_Zgeo
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%active_limiter_point%r = ES%R_lim
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%active_limiter_point%z = ES%Z_lim

    allocate(equilibrium_ids%time_slice(i_slice)%boundary_separatrix%x_point(2))
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%x_point(:)%r = ES%R_xpoint(:)
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%x_point(:)%z = ES%Z_xpoint(:)

    ! --- Export separatrix
    call separatrix(command_tmp, ierr, R_sep, Z_sep)
    allocate(equilibrium_ids%time_slice(i_slice)%boundary_separatrix%outline%r(size(R_sep)))
    allocate(equilibrium_ids%time_slice(i_slice)%boundary_separatrix%outline%z(size(R_sep)))
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%outline%r(:) = R_sep(:)
    equilibrium_ids%time_slice(i_slice)%boundary_separatrix%outline%z(:) = Z_sep(:)

    ! --- Export 2D quantities on a rectangular grid
    call clean_up()

    ! --- Find out range for the rectangular grid
    R_min = 1.d99;  R_max = -1.d99
    Z_min = 1.d99;  Z_max = -1.d99
    do i=1, node_list%n_nodes
      R_node = node_list%node(i)%x(1,1,1)
      Z_node = node_list%node(i)%x(1,1,2)
      if (R_node < R_min) R_min = R_node
      if (Z_node < Z_min) Z_min = Z_node
      if (R_node > R_max) R_max = R_node
      if (Z_node > Z_max) Z_max = Z_node
    enddo

    command_tmp%n_args = 7
    write(str, '(F16.12)') R_min
    command_tmp%args(1) = str  ! Rmin
    write(str, '(F16.12)') R_max
    command_tmp%args(2) = str  ! Rmax
    write(str, '(I0)') n_grid
    command_tmp%args(3) = str  ! nR
    write(str, '(F16.12)') Z_min
    command_tmp%args(4) = str  ! Zmin
    write(str, '(F16.12)') Z_max
    command_tmp%args(5) = str  ! Zmax
    write(str, '(I0)') n_grid
    command_tmp%args(6) = str  ! nZ
    command_tmp%args(7) = '0'  ! phi

    allocate(R_vec(n_grid), Z_vec(n_grid))
    R_vec = [(R_min + float((i-1)) * (R_max-R_min) / float((n_grid - 1)), i = 1, n_grid)]
    Z_vec = [(Z_min + float((i-1)) * (Z_max-Z_min) / float((n_grid - 1)), i = 1, n_grid)]
    
    expr_list = exprs((/'Psi', 'Jtor', 'BR', 'BZ', 'Btor'/), 5)
    call rectangle(command_tmp, first_step, ierr, only_n0=.true., res2D_out=result2D)
    
    allocate(equilibrium_ids%time_slice(i_slice)%profiles_2d(1))
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%type%index      = 0
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid_type%index = 1 ! --- Rectangular
    
    allocate(equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim1(n_grid))
    allocate(equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim2(n_grid))
    
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim1(:) = R_vec(:)
    equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%grid%dim2(:) = Z_vec(:)
    
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%r(n_grid, n_grid) )
    allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%z(n_grid, n_grid) )
    do i=1, n_grid
      do j=1, n_grid
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%r(i,j) = R_vec(i)  ! --- Som plotting tools use this field as well
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%z(i,j) = Z_vec(j)
      enddo
    enddo

    ! --- Fill profiles
    do i_exp=1, expr_list%n_expr

      ! --- psi
      if (expr_list%expr(i_exp)%name=='Psi') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%psi(n_grid, n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%psi(:,:) = result2D(:,:,i_exp) * fact_psi
      endif

      ! --- Jtor
      if (expr_list%expr(i_exp)%name=='Jtor') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%j_tor(n_grid, n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%j_tor(:,:) = result2D(:,:,i_exp) * fact_Ip
      endif

      ! --- B_R
      if (expr_list%expr(i_exp)%name=='BR') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_r(n_grid, n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_r(:,:) = result2D(:,:,i_exp)
      endif

      ! --- B_Z
      if (expr_list%expr(i_exp)%name=='BZ') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_z(n_grid, n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_z(:,:) = result2D(:,:,i_exp)
      endif

      ! --- B_tor
      if (expr_list%expr(i_exp)%name=='Btor') then
        allocate( equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_tor(n_grid, n_grid) )
        equilibrium_ids%time_slice(i_slice)%profiles_2d(1)%b_field_tor(:,:) = result2D(:,:,i_exp) * fact_Ip
      endif

    enddo

  end subroutine fill_equilibrium_IDS







  ! --- Fills Bezier coefficients in GGD
  subroutine fill_Bezier_coefficients( ggd_scalar, node_list, var_index, grid_ind, grid_sub_ind, res_fact )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_scalar),  intent(inout) ::  ggd_scalar
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index, grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
  
    ! --- Local parameters
    integer :: num_nodes, idof, inode, icoeff
    
    num_nodes = node_list%n_nodes
  
    if ( associated(ggd_scalar%coefficients) ) then
      deallocate( ggd_scalar%coefficients )
      allocate( ggd_scalar%coefficients(n_tor, num_nodes*(n_order+1)) ) 
    else
      allocate( ggd_scalar%coefficients(n_tor, num_nodes*(n_order+1)) ) 
    endif
  
    do inode=1, num_nodes
      do idof=1, n_order+1
  
        icoeff = inode + (idof-1)*num_nodes
  
        ggd_scalar%coefficients(:,icoeff)=node_list%node(inode)%values(:,idof,var_index)
  
      enddo
    enddo
  
    ggd_scalar%grid_index = grid_ind
    ggd_scalar%grid_subset_index = grid_sub_ind
 
   ! --- Re-scale coefficients
    ggd_scalar%coefficients = ggd_scalar%coefficients * res_fact
  
  end subroutine fill_Bezier_coefficients
  
  
  
  
  !< Fills JOREK grid into GGD
  subroutine grid2ggd( grid, node_list, element_list )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_aos3_root),     pointer   :: grid
    type (type_node_list),    intent(in)            :: node_list
    type (type_element_list), intent(in)            :: element_list
  
  
    ! --- Local parameters
    type(ids_generic_grid_dynamic_space), pointer   ::  space_RZ
    type(ids_generic_grid_dynamic_space), pointer   ::  space_fourier
    type(ids_generic_grid_dynamic_space_dimension), pointer :: ids_cells
    
    integer:: idx, shot_number, run_number, num_nodes, num_cells   
    integer :: gs_index, i, j
    integer, allocatable :: vertex_elm_array(:,:)
    real*8,  allocatable :: RZ(:,:)
    
    integer :: itor, idof, n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub
  
    ! --- create vertex - elements array
    allocate(  vertex_elm_array(n_vertex_max, element_list%n_elements)  )
    do i=1, element_list%n_elements
      vertex_elm_array(:,i) = element_list%element(i)%vertex
    enddo 
  
    ! Get values of R and Z at the nodes
    allocate( RZ(node_list%n_nodes, 2) )
    do i=1, node_list%n_nodes
      RZ(i,:) = node_list%node(i)%x(1,1,:)  
    enddo
  
    ! Write grid geometry
    allocate(  grid%space(2)                           )
    allocate(  grid%space(1)%objects_per_dimension(4)  )
    allocate(  grid%space(1)%coordinates_type(2)       )
  
    ! Set coordinates type to [R, Z]
    grid%space(1)%coordinates_type = (/ 4, 3 /)
  
    allocate(grid%identifier%description(1))
    allocate(grid%identifier%name(1))
    grid%identifier%description(1) = "Mesh JOREK output HDF5 file grid with quantities"
    grid%identifier%name = "JOREK output HDF5 file grid with quantities"
    grid%identifier%index = 1
  
    num_nodes = size(RZ,1)
    space_RZ  => grid%space(1)
  
    ! Fill simplified grid nodes (uses geometry instead of geometry_2D and misses Bezier representation) 
    allocate( space_RZ%objects_per_dimension(1)%object(num_nodes) )  ! Allocate to number of nodes
    do i=1, num_nodes 
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry(2) ) ! Allocate dimensions per each node
      space_RZ%objects_per_dimension(1)%object(i)%geometry(:) = RZ(i,:)
    enddo
  
    ! Fill dummy variables for 1D elements (edges)
    allocate( space_RZ%objects_per_dimension(2)%object(1) )          ! Allocate just one edge    
    allocate( space_RZ%objects_per_dimension(2)%object(1)%nodes(1))  
    space_RZ%objects_per_dimension(2)%object(1)%nodes(1) = 0
  
    ! Fill JOREK 2D elements (or cells)
    num_cells = element_list%n_elements
    ids_cells => space_RZ%objects_per_dimension(3)
    allocate(    ids_cells%object(num_cells) )
    do i=1, num_cells
      allocate(  ids_cells%object(i)%nodes(n_vertex_max)  )
      ids_cells%object(i)%nodes(:) = vertex_elm_array(:,i) 
    enddo
  
    ! Writing grid_subsets
    allocate(grid%grid_subset(2))  ! 2 grid subsets 
  
    ! Subset for points
    gs_index = 1
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name(1)        = "Nodes"
    grid%grid_subset(gs_index)%identifier%index          = gs_index 
    grid%grid_subset(gs_index)%identifier%description(1) = "All points/nodes/vertices/0D objects in the domain."
    grid%grid_subset(gs_index)%dimension                 = 1
   
    allocate( grid%grid_subset(gs_index)%element(num_nodes) )
   
    do i=1, num_nodes  
      allocate(  grid%grid_subset(gs_index)%element(i)%object(1)   )
      grid%grid_subset(gs_index)%element(i)%object(1)%space     = 1
      grid%grid_subset(gs_index)%element(i)%object(1)%index     = i 
      grid%grid_subset(gs_index)%element(i)%object(1)%dimension = 1
    enddo 
  
    ! Subset for cells
    gs_index = 2
  
    allocate(grid%grid_subset(gs_index)%identifier%name(1))
    allocate(grid%grid_subset(gs_index)%identifier%description(1))
  
    grid%grid_subset(gs_index)%identifier%name(1)        = "2D cells "
    grid%grid_subset(gs_index)%identifier%index          = gs_index 
    grid%grid_subset(gs_index)%identifier%description(1) = "All points/nodes/vertices/0D objects in the domain."
    grid%grid_subset(gs_index)%dimension                 = 3   ! A bit confusing, but 1 means point, 2 line and 3 surface
   
    allocate( grid%grid_subset(gs_index)%element(num_cells) )
  
    do i=1, num_cells
      allocate(  grid%grid_subset(gs_index)%element(i)%object(1)   )
      grid%grid_subset(gs_index)%element(i)%object(1)%space     = 1
      grid%grid_subset(gs_index)%element(i)%object(1)%index     = i 
      grid%grid_subset(gs_index)%element(i)%object(1)%dimension = 3
    enddo 
  
    ! Fill toroidal space (must be adapted for multiple time slices?)
    space_fourier  => grid%space(2)
    allocate(    space_fourier%coordinates_type(1)    )
    allocate(    space_fourier%identifier%description(1)  )
    space_fourier%coordinates_type(1) = 5          ! The coordinate type is 5, phi angle
    space_fourier%geometry_type%index = n_period   ! Fourier periodicity
    space_fourier%identifier%description(1) = "Toroidal space"             
  
    allocate(  space_fourier%objects_per_dimension(1)               )  ! We have only one dimension of
    allocate(  space_fourier%objects_per_dimension(1)%object(n_tor) )  ! toroidal harmonics
  
    do i=1, n_tor
      allocate(  space_fourier%objects_per_dimension(1)%object(i)%geometry(1) )  ! toroidal harmonics
      space_fourier%objects_per_dimension(1)%object(i)%geometry(1) = i
    enddo
  
    ! Fill in grid Bezier coefficients
    ! Needs generalization for JOREK 3D STELLERATOR EXTENSION!!
    space_RZ%geometry_type%index = 0  ! Standard geometry (non Fourier)
    do i=1, num_nodes
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,n_order+1) )
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(1,:) = node_list%node(i)%x(1,:,1 )   ! R dofs
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,:) = node_list%node(i)%x(1,:,2 )   ! Z dofs
    enddo
  
    ! JOREK element sizes
    do i=1, num_cells
      allocate( space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(n_order+1, n_vertex_max) )
      do j=1, n_vertex_max
        space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(:,j) = element_list%element(i)%size(j,:)   
      enddo
    enddo
  
  end subroutine grid2ggd




  ! --- Checks if a restart file exists in the current directory
  logical function restart_file_exists(i_step)

    use phys_module, only : rst_hdf5

    implicit none

    integer,  intent(in) :: i_step
    character(len=64)    :: file_name

    write(file_name,'(a,i5.5)') 'jorek', i_step
    if ( rst_hdf5 .ne. 0 ) then
      inquire (file=trim(file_name)//'.h5', exist=restart_file_exists)
    else
      inquire (file=trim(file_name)//'.rst', exist=restart_file_exists)
    end if

  end function restart_file_exists



#endif

end module mod_jorek2IMAS
