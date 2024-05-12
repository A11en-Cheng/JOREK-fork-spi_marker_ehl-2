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
                            zeroD_quantities, separatrix, rectangle, boundary_quantities
  use parse_commands, only: type_command
  use settings,       only: set_setting
  
  implicit none
  
 
  public
 
  ! Transform to COCOS convention 8 --> 11
  real*8 :: fact_psi =  2.d0 * PI
  real*8 :: fact_Ip  = -1.d0

    
  ! *******************************************************************************************************
  ! *************** Data structures needed to read geometry of the STARWALL coils *************************
  ! *******************************************************************************************************
  ! --- Constants for STARWALL coils 
  character(len=12), parameter    :: AXISYM_THICK = 'axisym_thick' !< Axisymmetric coils with dR,dZ extent
  character(len=12), parameter    :: AXISYM_FILA  = 'axisym_fila ' !< Axisymmetric coils specified by filaments
  character(len=12), parameter    :: GENERAL_THIN = 'general_thin' !< 3D coils specified by a list of points
  integer, parameter              :: N_MAX_FILA  = 2001       !< Maximum number of filaments per c
  integer, parameter              :: N_MAX_PARTS = 2001       !< Maximum number of parts per coil 
  integer, parameter              :: N_MAX_PTS   = 2001       !< Maximum number of points per coil

  ! --- Data structure for a single STARWALL coil of arbitrary type
  type :: t_coil_starwall
    ! --- For all coil types
    character(len=128)   :: name     = 'NONE'        !< Name of the coil
    character(len=32)    :: coil_type                !< Coil type according to constants specified above
    real*8               :: resist   = -999.         !< Total coil resistivity in Ohm
    real*8               :: nturns   = -1d4          !< Total number of turns in the coil
    integer              :: ntorpts  = 0             !< Number of toroidal points (axisym coils)
    ! --- For axisymmetric coils with dR and dZ extent
    real*8               :: R(N_MAX_PARTS)        = -999.      !< Position in R (center of coil)
    real*8               :: Z(N_MAX_PARTS)        = -999.      !< Position in Z (center of coil)
    real*8               :: dR(N_MAX_PARTS)       = -999.      !< Full width in R
    real*8               :: dZ(N_MAX_PARTS)       = -999.      !< Full width in Z
    integer              :: nbands_R = 0                       !< #radial bands for each coil part
    integer              :: nbands_Z = 0                       !< #toroidal bands for each coil part
    real*8               :: n_thick_turns(N_MAX_PARTS) = -1.d4 !< Number of turns in each coil part
    integer              :: nparts_coil  = 1                  !< Number of bands
    ! --- For coils specified by filaments
    integer              :: n_fila   = 0                  !< Number of filaments for the coil
    real*8               :: n_fila_turns(N_MAX_FILA) = 1. !< Number of turns in each filament
    real*8               :: R_fila(N_MAX_FILA) = -999.    !< R-position of filaments
    real*8               :: Z_fila(N_MAX_FILA) = -999.    !< Z-position of filaments
    ! --- For 3D coil specified by a list of points
    integer              :: n_pts    = 0             !< Number of points along the coil
    real*8               :: width    = -999.         !< Width of coil
    real*8               :: xpts(N_MAX_PTS) = -999.  !< x-Position of points
    real*8               :: ypts(N_MAX_PTS) = -999.  !< y-Position of points
    real*8               :: zpts(N_MAX_PTS) = -999.  !< z-Position of points
    real*8               :: dir3d(3*N_MAX_PTS) = -999. !< optional, band direction for each point
  end type t_coil_starwall

  ! --- Data structure for a coil set
  type :: t_coil_set_starwall
    character(len=32)    :: name         = 'NONE'    !< Name of coil set (should not be changed by user)
    character(len=256)   :: description  = 'NONE'    !< Human readable description
    integer              :: ncoil        = 0         !< Number of coils in the set
    integer              :: index_start  = -999      !< Index of first coil of this set in all_coil_set
    type(t_coil_starwall), allocatable :: coil(:)             !< Coils in this coil set
  end type t_coil_set_starwall
  ! *******************************************************************************************************


  contains


  subroutine fill_mhd_IDS(first_step, time_SI, mhd_ids)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass

    implicit none

    ! --- External parameters
    logical,   intent(in) :: first_step   ! is this the first step?
    real*8,    intent(in) :: time_SI
    
    type(ids_mhd), target,     intent(inout) :: mhd_ids

   
    ! --- Local parameters 
    integer    :: i, j, k, m, etype, irst, int, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: fact_T, fact_v, fact_zj, rho0, fact_phi, fact_rho, fact_w
    
    
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
  
    mhd_ids%time(i_slice)     = time_SI 
    mhd_ids%ggd(i_slice)%time = time_SI

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




  ! --- Fill a wall IDS from STARWALL data, needs to be adapted to CARIDDI!
  subroutine fill_wall_IDS(first_step, time_SI, wall_thickness, wall_ids)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    use vacuum
    use vacuum_response,  only: reconstruct_triangle_potentials

    implicit none

    ! --- External parameters
    logical,                 intent(in)    :: first_step      ! is this the first step?
    real*8,                  intent(in)    :: time_SI
    real*8,                  intent(in)    :: wall_thickness  !< effective thin wall thickness
    type(ids_wall),  target, intent(inout) :: wall_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id=0, ierr, i_tri
    integer    :: i_phi, i_pol, n_phi, n_pol, ipol1, ipol2, ipol3, ipol4, itor1, itor2, itor3, itor4
    integer    :: i1, i2, i3, i4, i_exp
    real*8     :: rho0, fact_rad, co, si, Rtri
    real*8     :: phi1, phi2, phi3, r1(3), r2(3), r3(3), r21(3), r32(3), r13(3), r21_cross_r32(3)
    real*8     :: j_lin(3), r_mid(3), Iw_net_tor
    real*8, allocatable :: tripot_w(:)
    real*8, allocatable :: result(:,:,:)
    character(10)       :: str
    type(type_command)  :: command_tmp
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    type(ids_generic_grid_scalar),      pointer :: ggd_scalar
    type(ids_generic_grid_aos3_root),   pointer :: grid
    
    integer:: n_wall_nodes, n_wall_triangles, n_walls, i_vv, i_fw
    
    integer :: n_slice, i_slice, grid_ind, grid_sub_ind, n_grid_sub, n_grid
    ! **********************************************************************************
  
    ! --- Number of grids and grid subsets
    n_grid       = 1
    n_grid_sub   = 1
    grid_ind     = 1  ! Index
    grid_sub_ind = 1  ! Index    
    n_walls      = 2  ! 1 for vessel current, 2 for first wall fluxes
    i_vv         = 1
    i_fw         = 2

    wall_ids%ids_properties%homogeneous_time = 1

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  wall_ids%time(n_slice) )

    n_wall_triangles = sr%ntri_w

    if (first_step) then
      allocate( wall_ids%description_ggd(n_walls))
    endif

    ! *******************************************************************************
    ! ************* Export STARWALL wall currents to wall_IDS ***********************
    ! *******************************************************************************
    if (first_step) then

      ! --- Put the wall grid in GGD
      allocate( wall_ids%description_ggd(i_vv)%grid_ggd(n_grid) )

      wall_ids%description_ggd(i_vv)%type%index  = 2  ! For thin wall description
      
      grid => wall_ids%description_ggd(i_vv)%grid_ggd(grid_ind)
      
      grid%identifier%index = 0   ! Unspecified
      allocate( grid%identifier%description(1))
      grid%identifier%description(1) = "Thin wall described with linear triangles"

      allocate(grid%space(1))
      grid%space(1)%identifier%index = 1  ! Primary space

      grid%space(1)%geometry_type%index = 0  ! Standard (not Fourier)
      allocate(grid%space(1)%coordinates_type(3))

      grid%space(1)%coordinates_type(:) = (/ 1, 2, 3 /)  ! Identifiers for x, y, z type coordinates

      allocate(grid%space(1)%objects_per_dimension(3))

      ! --- Save wall grid nodes
      n_wall_nodes = sr%npot_w
      allocate(grid%space(1)%objects_per_dimension(1)%object(n_wall_nodes))
      grid%space(1)%objects_per_dimension(1)%geometry_content%index = 1  ! node coordinates
      do i=1, n_wall_nodes
        allocate( grid%space(1)%objects_per_dimension(1)%object(i)%geometry(3) ) ! Allocate dimensions for each node
        grid%space(1)%objects_per_dimension(1)%object(i)%geometry(:) = sr%xyzpot_w(i,:)
      enddo

      ! --- Save thin wall triangles 
      allocate(grid%space(1)%objects_per_dimension(3)%object(n_wall_triangles))  ! Index 3 for 2D objects (faces)
      do i=1, n_wall_triangles
        allocate( grid%space(1)%objects_per_dimension(3)%object(i)%nodes(3) ) ! 3 nodes per triangle
        grid%space(1)%objects_per_dimension(3)%object(i)%nodes(:) = sr%jpot_w(i,:)  ! The node indices of this triangle
      enddo

      ! --- Create a grid subset where the wall triangles are the elements of the subset
      ! --- 1 current density value will be assigned per triangle 
      allocate(grid%grid_subset(1))
      grid%grid_subset(1)%identifier%index = 5  ! Identifier index for 2D cells in the dictionary
      grid%grid_subset(1)%dimension        = 3  ! Index 3 means 2 dimensions in the dictionary

      ! --- Save wall thickness
      allocate(wall_ids%description_ggd(i_vv)%thickness(n_grid))
      allocate(wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(n_grid_sub))
      allocate(wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%values(n_wall_triangles))

      wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%grid_index        = grid_ind
      wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%grid_subset_index = grid_sub_ind

      wall_ids%description_ggd(i_vv)%thickness(1)%grid_subset(1)%values(:) = wall_thickness

      !--- Information about the wall component
      allocate(wall_ids%description_ggd(i_vv)%component(n_grid))
      allocate(wall_ids%description_ggd(i_vv)%component(1)%type(1))
      allocate(wall_ids%description_ggd(i_vv)%component(1)%identifiers(1))
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%grid_index        = grid_ind
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%grid_subset_index = grid_sub_ind
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%index  = 0  ! --- 0 means not specified yet
      allocate(wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%name(1))
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%name = "Vacuum Vessel"
      allocate(wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%description(1))
      wall_ids%description_ggd(i_vv)%component(1)%type(1)%identifier%description = "A single layer of the &
                                                                    vacuum vessel discretized with linear thin triangles"

    endif

    ! --- Set times
    allocate( wall_ids%description_ggd(i_vv)%ggd(n_slice) )
  
    wall_ids%time(i_slice) = time_SI 
    wall_ids%description_ggd(i_vv)%ggd(i_slice)%time = time_SI
    
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(n_grid_sub) )
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%r(n_wall_triangles) ) ! --- one value per triangle
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%z(n_wall_triangles) ) ! --- one value per triangle
    allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%toroidal(n_wall_triangles) ) ! --- one value per triangle

    wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%grid_index        = 1
    wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%grid_subset_index = 1

    ! --- Get the current density per triangle
    ! --- The triangle current density
    call reconstruct_triangle_potentials(tripot_w, wall_curr, my_id, Iw_net_tor)

    do i = 1, n_wall_triangles

      ! --- Wall potential at triangle nodes
      phi1   = tripot_w(sr%jpot_w(i,1)) + Iw_net_tor*sr%phi0_w(i,1) 
      phi2   = tripot_w(sr%jpot_w(i,2)) + Iw_net_tor*sr%phi0_w(i,2) 
      phi3   = tripot_w(sr%jpot_w(i,3)) + Iw_net_tor*sr%phi0_w(i,3) 

      ! --- Position of triangle nodes
      r1(:)  = sr%xyzpot_w(sr%jpot_w(i,1),:)
      r2(:)  = sr%xyzpot_w(sr%jpot_w(i,2),:)
      r3(:)  = sr%xyzpot_w(sr%jpot_w(i,3),:)
      r21(:) = r1(:)-r2(:)
      r32(:) = r2(:)-r3(:)
      r21_cross_r32(:) = (/ r21(2)*r32(3) - r21(3)*r32(2), r21(3)*r32(1) - r21(1)*r32(3),          &
        r21(1)*r32(2) - r21(2)*r32(1) /)

      j_lin(:) = ( phi1*(r3-r2)+phi2*(r1-r3)+phi3*(r2-r1) ) / sqrt(sum(r21_cross_r32**2) ) / mu_zero 

      r_mid(:) = (r1 + r2 + r3) / 3.d0   ! Middle point of the triangle

      Rtri = sqrt( r_mid(1)**2.d0 + r_mid(2)**2.d0 )
      co   =  r_mid(1) / Rtri
      si   = -r_mid(2) / Rtri

      wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%r(i)        =  ( j_lin(1)*co - j_lin(2)*si ) / wall_thickness
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%z(i)        =    j_lin(3)                    / wall_thickness
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%j_total(1)%toroidal(i) =  (-j_lin(1)*si - j_lin(2)*co ) / wall_thickness * fact_Ip   
    enddo

    if (first_step) then
      allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(n_grid_sub) )
      allocate( wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%values(n_wall_triangles) ) 
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%values(:) = sr%eta_thin_w * wall_thickness * wall_resistivity_fact

      wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%grid_index        = 1
      wall_ids%description_ggd(i_vv)%ggd(i_slice)%resistivity(1)%grid_subset_index = 1
    endif
    ! *******************************************************************************
   

    ! *******************************************************************************
    ! ************* Export FW and divertor fluxes and currents **********************
    ! *******************************************************************************
    ! --- Call expressions to compute boundary quantities
    step_imported = .true.

    call initialise_postproc_settings(first_step, 100)
  
    ! --- Arguments for boundary quantities
    command_tmp%n_args = 3
    command_tmp%args(1) = '0'              ! phimin
    command_tmp%args(2) = '6.28318530718'  ! phimax
    n_phi = max(32, n_tor* 3)
    write(str, '(I0)') n_phi
    command_tmp%args(3) = str  
    call clean_up()
    expr_list = exprs((/'x', 'y', 'z', 'Psi', 'A_R', 'A_Z', &
                        'BR', 'BZ', 'Btor', 'JR', 'JZ', 'Jtor',  &
                        'heatF_total', 'partF_total', 'npartF_total'/), 15)
    
    call boundary_quantities(command_tmp, first_step==.true., ierr, result)
    call clean_up()

    n_pol = size(result, 2)
    n_wall_triangles = 2 * n_phi * n_pol 
    n_wall_nodes     = n_phi * n_pol 

    ! --- Export grid boundary geometry as thin triangles (FW + divertor)
    if (first_step) then

      ! --- Put the wall grid in GGD
      allocate( wall_ids%description_ggd(i_fw)%grid_ggd(n_grid) )

      wall_ids%description_ggd(i_fw)%type%index  = 2  ! For thin wall description
      
      grid => wall_ids%description_ggd(i_fw)%grid_ggd(grid_ind)
      
      grid%identifier%index = 0   ! Unspecified
      allocate( grid%identifier%description(1))
      grid%identifier%description(1) = "Thin wall described with linear triangles"

      allocate(grid%space(1))
      grid%space(1)%identifier%index = 1  ! Primary space

      grid%space(1)%geometry_type%index = 0  ! Standard (not Fourier)
      allocate(grid%space(1)%coordinates_type(3))

      grid%space(1)%coordinates_type(:) = (/ 1, 2, 3 /)  ! Identifiers for x, y, z type coordinates

      allocate(grid%space(1)%objects_per_dimension(3))

      ! --- Create and export grid nodes and triangles
      allocate(grid%space(1)%objects_per_dimension(1)%object(n_wall_nodes))
      allocate(grid%space(1)%objects_per_dimension(3)%object(n_wall_triangles))  ! Index 3 for 2D objects (faces)
      grid%space(1)%objects_per_dimension(1)%geometry_content%index = 1  ! node coordinates
      
      i_tri = 0
      do i_phi=1, n_phi
        do i_pol=1, n_pol
          
          i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node

          ! --- Get global indices of the 4 nodes forming a quadrilateral, and make two triangles out of it
          ipol1 = i_pol;                  itor1 = i_phi;
          ipol2 = mod(i_pol,n_pol) + 1;   itor2 = i_phi;
          ipol3 = mod(i_pol,n_pol) + 1;   itor3 = mod(i_phi,n_phi) + 1;
          ipol4 = i_pol;                  itor4 = mod(i_phi,n_phi) + 1;

          i1 = ipol1 + (itor1-1)*n_pol
          i2 = ipol2 + (itor2-1)*n_pol
          i3 = ipol3 + (itor3-1)*n_pol
          i4 = ipol4 + (itor4-1)*n_pol
          
          ! --- Fill in reference node coordinates
          r1(:) = (/ result(itor1,ipol1,1), result(itor1,ipol1,2), result(itor1,ipol1,3) /) ! x, y, z coordinates
          allocate( grid%space(1)%objects_per_dimension(1)%object(i)%geometry(3) ) ! Allocate dimensions for each node
          grid%space(1)%objects_per_dimension(1)%object(i)%geometry(:) =  r1

          ! --- Fill in two triangles
          ! --- Triangle 1
          i_tri = i_tri + 1
          allocate( grid%space(1)%objects_per_dimension(3)%object(i_tri)%nodes(3) ) ! 3 nodes per triangle
          grid%space(1)%objects_per_dimension(3)%object(i_tri)%nodes(:) = (/ i1, i2, i3/)  ! The node indices of this triangle

          ! --- Triangle 2
          i_tri = i_tri + 1
          allocate( grid%space(1)%objects_per_dimension(3)%object(i_tri)%nodes(3) ) ! 3 nodes per triangle
          grid%space(1)%objects_per_dimension(3)%object(i_tri)%nodes(:) = (/i1, i3, i4 /)  ! The node indices of this triangle

        enddo
      enddo

      ! --- Create a grid subset where the wall nodes elements of the subset
      ! --- values will be assigned to the nodes
      allocate(grid%grid_subset(1))
      grid%grid_subset(1)%identifier%index = 1  ! Identifier index for 0D nodes in the dictionary
      grid%grid_subset(1)%dimension        = 1  ! Index 1 means 0 dimensions in the dictionary (for 0D nodes)

      !--- Information about the wall component
      allocate(wall_ids%description_ggd(i_fw)%component(n_grid))
      allocate(wall_ids%description_ggd(i_fw)%component(1)%type(1))
      allocate(wall_ids%description_ggd(i_fw)%component(1)%identifiers(1))
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%grid_index        = grid_ind
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%grid_subset_index = grid_sub_ind
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%index  = 0  ! --- 0 means not specified yet
      allocate(wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%name(1))
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%name = "FW + divertor"
      allocate(wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%description(1))
      wall_ids%description_ggd(i_fw)%component(1)%type(1)%identifier%description = "A single layer representing the first wall &
                                                                    + divertor surfaces discretized with linear thin triangles"

    endif

    allocate( wall_ids%description_ggd(i_fw)%ggd(n_slice) )
  
    wall_ids%time(i_slice) = time_SI 
    wall_ids%description_ggd(i_fw)%ggd(i_slice)%time = time_SI

    ! --- Fill expressions in the wall nodes
    allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(n_grid_sub) )
    wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%grid_index        = 1
    wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%grid_subset_index = 1

    do i_exp=1, expr_list%n_expr
      
      ! --- Total heatflux
      if (expr_list%expr(i_exp)%name=='heatF_total') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(n_grid_sub) )
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%values(n_wall_nodes) ) ! --- one value per node

        wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%grid_index        = 1
        wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%grid_subset_index = 1

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%power_density(1)%values(i) = result(i_phi,i_pol,i_exp)
          enddo
        enddo
      endif

      ! --- Psi
      if (expr_list%expr(i_exp)%name=='Psi') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(n_grid_sub) )
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%values(n_wall_nodes) ) ! --- one value per node

        wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%grid_index        = 1
        wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%grid_subset_index = 1

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%psi(1)%values(i) = result(i_phi,i_pol,i_exp)
          enddo
        enddo
      endif

      ! --- J_phi
      if (expr_list%expr(i_exp)%name=='Jtor') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%toroidal(n_wall_nodes) ) ! --- one value per node

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%toroidal(i) = result(i_phi,i_pol,i_exp) * fact_Ip
          enddo
        enddo
      endif

      ! --- J_R
      if (expr_list%expr(i_exp)%name=='JR') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%r(n_wall_nodes) ) ! --- one value per node

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%r(i) = result(i_phi,i_pol,i_exp) 
          enddo
        enddo
      endif

      ! --- J_Z
      if (expr_list%expr(i_exp)%name=='JZ') then
        allocate( wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%z(n_wall_nodes) ) ! --- one value per node

        do i_phi=1, n_phi
          do i_pol=1, n_pol 
            i = i_pol + (i_phi-1)*n_pol  !< Global index of refence node
            wall_ids%description_ggd(i_fw)%ggd(i_slice)%j_total(1)%z(i) = result(i_phi,i_pol,i_exp) 
          enddo
        enddo
      endif

    
    enddo
    ! *******************************************************************************
  
  end subroutine fill_wall_IDS
 




  ! --- Fill a pf_passive IDS from STARWALL data, needs to be adapted to CARIDDI!
  subroutine fill_pf_passive_IDS(first_step, time_SI, pf_passive, passive_coil_geo_file)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    use vacuum
    use vacuum_response,  only: reconstruct_coil_potentials

    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(ids_pf_passive), intent(inout) :: pf_passive
    character(len=*),        intent(in) :: passive_coil_geo_file
   
    ! --- Local parameters 
    integer :: i, j, k, m, my_id=0, ierr, i_coil
    integer :: n_slice=1, i_slice=1
    logical :: found_coil
    real*8, allocatable :: pot_c(:) 
    type(t_coil_set_starwall) :: coil_set
   
    pf_passive%ids_properties%homogeneous_time = 1

    if (sr%n_diag_coils == 0) then
      write(*,*) '  No diagnostic coils in the STARWALL response file'
      write(*,*) '  needed for pf_passive IDS'
      stop
    endif

    ! --- Read STARWALL coil geometry input file for passive conductors
    call read_coil_set_starwall(passive_coil_geo_file, coil_set)

    ! --- Set times
    allocate( pf_passive%time(n_slice) )
    pf_passive%time(i_slice) = time_SI 

    if (first_step) then
      allocate(pf_passive%loop(coil_set%ncoil))
    endif

    ! --- Reconstruct currents from wall_curr
    call reconstruct_coil_potentials(pot_c, wall_curr, my_id)

    do i_coil=1, coil_set%ncoil
       ! --- Check if the coil given in the input file exists in the JOREK restart
      found_coil = .false.
      do i=1, sr%ncoil
        if (INDEX(trim(sr%coil_name(i)),trim(coil_set%coil(i_coil)%name)) > 0) then
          found_coil = .true.
          exit
        endif
      enddo
      if ( .not. found_coil) then
        write(*,*) coil_set%coil(i_coil)%name, " is not part of the starwall_response.dat file data!"
        stop
      endif

      if (coil_set%coil(i_coil)%coil_type /= AXISYM_THICK) then
        write(*,*) coil_set%coil(i_coil)%coil_type," coil type in pf_passive not supported yet in jorek2_IDS"
        stop
      endif

      ! --- Fill coil geometry for first time step
      if (first_step) then
        allocate( pf_passive%loop(i_coil)%name(1) )
        pf_passive%loop(i_coil)%name       = trim(coil_set%coil(i_coil)%name)   
        pf_passive%loop(i_coil)%resistance = coil_set%coil(i_coil)%resist * wall_resistivity_fact
        allocate( pf_passive%loop(i_coil)%element(coil_set%coil(i_coil)%nparts_coil) )
        do k=1, coil_set%coil(i_coil)%nparts_coil
          pf_passive%loop(i_coil)%element(k)%turns_with_sign           = coil_set%coil(i_coil)%n_thick_turns(k)
          pf_passive%loop(i_coil)%element(k)%geometry.geometry_type    = 2 ! Rectangle type
          pf_passive%loop(i_coil)%element(k)%geometry.rectangle.r      = coil_set%coil(i_coil)%R(k)
          pf_passive%loop(i_coil)%element(k)%geometry.rectangle.z      = coil_set%coil(i_coil)%Z(k)
          pf_passive%loop(i_coil)%element(k)%geometry.rectangle.width  = coil_set%coil(i_coil)%dR(k)
          pf_passive%loop(i_coil)%element(k)%geometry.rectangle.height = coil_set%coil(i_coil)%dZ(k)
        enddo
      endif

      allocate( pf_passive%loop(i_coil)%current(n_slice) )

      pf_passive%loop(i_coil)%current(i_slice) = pot_c(i) / MU_ZERO * fact_Ip

    enddo

  end subroutine fill_pf_passive_IDS





   ! --- Fill a pf_active IDS from STARWALL data, needs to be adapted to CARIDDI!
  subroutine fill_pf_active_IDS(first_step, time_SI, pf_active, active_coil_geo_file)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    use vacuum
    use vacuum_response,  only: reconstruct_coil_potentials

    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(ids_pf_active),  intent(inout) :: pf_active
    character(len=*),        intent(in) :: active_coil_geo_file
   
    ! --- Local parameters 
    integer :: i, j, k, m, my_id=0, ierr, i_coil
    integer :: n_slice=1, i_slice=1
    logical :: found_coil
    real*8, allocatable :: pot_c(:) 
    type(t_coil_set_starwall) :: coil_set
   
    pf_active%ids_properties%homogeneous_time = 1

    if (sr%n_pol_coils == 0) then
      write(*,*) '  No poloidal active coils in the STARWALL response file'
      write(*,*) '  needed for pf_active IDS'
      stop
    endif

    ! --- Read STARWALL coil geometry input file for active conductors
    call read_coil_set_starwall(active_coil_geo_file, coil_set)

    ! --- Set times
    allocate( pf_active%time(n_slice) )
    pf_active%time(i_slice) = time_SI 

    if (first_step) then
      allocate(pf_active%coil(coil_set%ncoil))
    endif

    ! --- Reconstruct currents from wall_curr
    call reconstruct_coil_potentials(pot_c, wall_curr, my_id)

    do i_coil=1, coil_set%ncoil
       ! --- Check if the coil given in the input file exists in the JOREK restart
      found_coil = .false.
      do i=1, sr%ncoil
        if (INDEX(trim(sr%coil_name(i)),trim(coil_set%coil(i_coil)%name)) > 0) then
          found_coil = .true.
          exit
        endif
      enddo
      if ( .not. found_coil) then
        write(*,*) coil_set%coil(i_coil)%name, " is not part of the starwall_response.dat file data!"
        stop
      endif

      if (coil_set%coil(i_coil)%coil_type /= AXISYM_THICK) then
        write(*,*) coil_set%coil(i_coil)%coil_type," coil type in pf_active not supported yet in jorek2_IDS"
        stop
      endif

      ! --- Fill coil geometry for first time step
      if (first_step) then
        allocate( pf_active%coil(i_coil)%name(1) )
        pf_active%coil(i_coil)%name       = trim(coil_set%coil(i_coil)%name)   
        pf_active%coil(i_coil)%resistance = coil_set%coil(i_coil)%resist * wall_resistivity_fact
        allocate( pf_active%coil(i_coil)%element(coil_set%coil(i_coil)%nparts_coil) )
        do k=1, coil_set%coil(i_coil)%nparts_coil
          pf_active%coil(i_coil)%element(k)%turns_with_sign           = coil_set%coil(i_coil)%n_thick_turns(k)
          pf_active%coil(i_coil)%element(k)%geometry.geometry_type    = 2 ! Rectangle type
          pf_active%coil(i_coil)%element(k)%geometry.rectangle.r      = coil_set%coil(i_coil)%R(k)
          pf_active%coil(i_coil)%element(k)%geometry.rectangle.z      = coil_set%coil(i_coil)%Z(k)
          pf_active%coil(i_coil)%element(k)%geometry.rectangle.width  = coil_set%coil(i_coil)%dR(k)
          pf_active%coil(i_coil)%element(k)%geometry.rectangle.height = coil_set%coil(i_coil)%dZ(k)
        enddo
      endif

      allocate( pf_active%coil(i_coil)%current%data(n_slice) )

      pf_active%coil(i_coil)%current%data(i_slice) = pot_c(i) / MU_ZERO * fact_Ip

    enddo

  end subroutine fill_pf_active_IDS






  subroutine fill_radiation_IDS(first_step, time_SI, radiation_ids)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,                 intent(in) :: first_step   ! is this the first step?
    real*8,                  intent(in) :: time_SI
    type(ids_radiation),  intent(inout) :: radiation_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: rho0, fact_rad
    
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

    ! --- Set times
    n_slice = 1  
    i_slice = 1
    allocate(  radiation_ids%time(n_slice) )

    radiation_ids%ids_properties%homogeneous_time = 1
    allocate( radiation_ids%process(1))   ! --- 1 type of radiation
    allocate( radiation_ids%process(1)%ggd(n_slice) )
  
    radiation_ids%time(i_slice)                = time_SI 
    radiation_ids%process(1)%ggd(i_slice)%time = time_SI
 
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






  subroutine fill_core_profiles_IDS(first_step, time_SI, core_profiles_ids, n_grid)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,      intent(in) :: first_step   ! is this the first step?
    real*8,       intent(in) :: time_SI
    integer,      intent(in) :: n_grid       ! Number of flux surfaces to compute average
    type(ids_core_profiles),   intent(inout) :: core_profiles_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    integer    :: a_imp, z_imp, i_ion_main, i_ion_imp
    real*8     :: rho0
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
    
    core_profiles_ids%ids_properties%homogeneous_time = 1    
    core_profiles_ids%time(i_slice) = time_SI 

    ! --- Call expressions and do a flux average
    step_imported = .true.

    call initialise_postproc_settings(first_step, n_grid)
  
    ! --- Get average and q-profile
    command_tmp%n_args = 0
    call clean_up()
    expr_list = exprs((/'Psi_N', 'T_i', 'T_e', 'ne', 'pres', 'Phi', 'eta_T', &
                        'Jpar', 'E_||', 'Er', 'vpar', 'Vtheta_i', 'Vstar_i', 'rho', 'Psi', &
                        'Z_eff', 'nimp', 'ni_main', 'nn_main'/), 19)
    
    ! --- If loss of LCFS, abort core_profiles
    if ( .not. ES%LCFS_is_lost ) then
      call average(command_tmp, first_step==.true., ierr, result, .true.)
      call clean_up()
      call qprofile(command_tmp, first_step==.true., ierr, q_prof)
    else 
      if(allocated(result)) deallocate(result)
      if(allocated(q_prof)) deallocate(q_prof)
      allocate(result(n_grid, expr_list%n_expr))
      allocate(q_prof(n_grid))
      result = -1.d99;   q_prof = -1.d99;
      write(*,*) '  core_profiles cannot be produced without closed flux surfaces'
    endif
    ! --- Correct first and last points
    q_prof(1)      = q_prof(2)        + (q_prof(2)-q_prof(3))
    q_prof(n_grid) = q_prof(n_grid-1) + (q_prof(n_grid-1)-q_prof(n_grid-2))

    ! --- Some allocations
    allocate( core_profiles_ids%profiles_1d(i_slice)%ion(1+n_adas) ) ! First index is for main ions
    allocate( core_profiles_ids%profiles_1d(i_slice)%neutral(1+n_adas) )
    i_ion_main = 1;  i_ion_imp = 2;

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
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%parallel(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%parallel(:) = result(:,i_exp)
      endif

      ! --- Poloidal velocity
      if (expr_list%expr(i_exp)%name=='Vtheta_i') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%poloidal(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%poloidal(:) = result(:,i_exp)
      endif

      ! --- Diamagnetic velocity
      if (expr_list%expr(i_exp)%name=='Vstar_i') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%diamagnetic(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%velocity%diamagnetic(:) = result(:,i_exp)
      endif

      ! --- Z_eff
      if (expr_list%expr(i_exp)%name=='Z_eff') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%zeff(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%zeff(:) = result(:,i_exp)
      endif

      ! --- Ion density
      if (expr_list%expr(i_exp)%name=='ni_main') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%density(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%density(:) = result(:,i_exp) 
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%element(1) )
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%element(1)%a   = central_mass
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_main)%element(1)%z_n = 1   ! Main ions have Z=1
      endif

      ! --- Neutral density (of main ions)
      if (expr_list%expr(i_exp)%name=='nn_main') then
        allocate( core_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%density(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%density(:) = result(:,i_exp) 
        allocate( core_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%element(1) )
        core_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%element(1)%a   = central_mass
        core_profiles_ids%profiles_1d(i_slice)%neutral(i_ion_main)%element(1)%z_n = 1   ! Main ions have Z=1
      endif

      ! --- Main impurity density
      if (expr_list%expr(i_exp)%name=='nimp') then   ! ion index 2 is for the main impurity species
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%density(n_grid) )
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%density(:) = result(:,i_exp)
        allocate( core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%element(1) )
        select case ( trim(imp_type(index_main_imp)) )
        case('D2')
           a_imp  = 2
           z_imp  = 1
        case('Ar')
           a_imp  = 40
           z_imp  = 18
        case('Ne')
           a_imp  = 20
           z_imp  = 10
        case('Be')
            a_imp  = 9
            z_imp  = 4
        case('W')
            a_imp  = 184
            z_imp  = 74
        case default
            write(*,*) '!! Impurity type "', trim(imp_type(index_main_imp)), '" unknown !!'
            write(*,*) '=> We assume D2.'
            a_imp  = 2
            z_imp  = 1
        end select
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%element(1)%a   = a_imp
        core_profiles_ids%profiles_1d(i_slice)%ion(i_ion_imp)%element(1)%z_n = z_imp
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






  ! --- This subroutine fills different IDSs, they are filled together here to calculate 0D
  ! --- quantities by calling mod_integrals3D only once and avoid duplications
  subroutine fill_IDSs_w_common_quantities(first_step, time_SI, n_grid, export_equil, export_summary, export_disruption, &
                                           equilibrium_ids, summary_ids, disruption_ids, simulation_description)

    implicit none

    ! --- External parameters
    logical,      intent(in) :: first_step   ! is this the first step?
    real*8,       intent(in) :: time_SI
    logical,      intent(in) :: export_equil, export_summary, export_disruption
    integer,      intent(in) :: n_grid       ! Number of flux surfaces to compute average
    type(ids_equilibrium),  intent(inout)  :: equilibrium_ids
    type(ids_summary),      intent(inout)  :: summary_ids
    type(ids_disruption),   intent(inout)  :: disruption_ids
    character(len=1000)                    :: simulation_description

    ! --- Local parameters
    real*8, allocatable :: res0D(:)
    integer             :: ierr
    type(type_command)  :: command_tmp

    step_imported = .true.
    call initialise_postproc_settings(first_step, n_grid)

    command_tmp%n_args = 0
    call clean_up()
    call zeroD_quantities(command_tmp, first_step==.true., ierr, res0D)

    if (export_equil)       call fill_equilibrium_IDS(first_step, time_SI, n_grid, res0D, equilibrium_ids)  
    if (export_summary)     call fill_summary_IDS(first_step, time_SI, res0D, summary_ids, simulation_description)  
    if (export_disruption)  call fill_disruption_IDS(first_step, time_SI, res0D, disruption_ids) 
  
  
  end subroutine fill_IDSs_w_common_quantities






  subroutine fill_equilibrium_IDS(first_step, time_SI, n_grid, res0D, equilibrium_ids)  

    use phys_module, only : F0, central_density, sqrt_mu0_rho0, &
                           sqrt_mu0_over_rho0, central_mass, imp_type, &
                           gamma, index_main_imp
    implicit none

    ! --- External parameters
    logical,             intent(in) :: first_step   ! is this the first step?
    real*8,              intent(in) :: time_SI
    integer,             intent(in) :: n_grid       ! Number of flux surfaces to compute average
    real*8, allocatable, intent(in) :: res0D(:)     ! List of 0D quantities defined in exprs_all_int (mod_expressions.f90)
    
    type(ids_equilibrium),  intent(inout)  :: equilibrium_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: rho0, fact_rad, R_min, Z_min, R_max, Z_max, R_node, Z_node
    real*8, allocatable :: result(:,:), q_prof(:), rho_tor(:), R_sep(:), Z_sep(:)
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
    
    equilibrium_ids%ids_properties%homogeneous_time = 1    
    equilibrium_ids%time(i_slice) = time_SI 

    ! --- Call expressions and do a flux average
    
    ! --- Get average and q-profile
    command_tmp%n_args = 0
    call clean_up()
    expr_list = exprs((/'Psi', 'pres', 'FFprime_loc', 'p_prime_loc', 'Jpar'/), 5)
    ! --- If loss of LCFS, abort profiles
    if ( .not. ES%LCFS_is_lost ) then
      call average(command_tmp, first_step==.true., ierr, result, .true.)
      call clean_up()
      call qprofile(command_tmp, first_step==.true., ierr, q_prof)
    else 
      if(allocated(result)) deallocate(result)
      if(allocated(q_prof)) deallocate(q_prof)
      allocate(result(n_grid, expr_list%n_expr))
      allocate(q_prof(n_grid))
      result = -1.d99;   q_prof = -1.d99;
      write(*,*) '  profiles cannot be produced without closed flux surfaces'
    endif
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






  subroutine fill_summary_IDS(first_step, time_SI, res0D, summary_ids, simulation_description)  

    implicit none

    ! --- External parameters
    logical,             intent(in) :: first_step   ! is this the first step?
    real*8,              intent(in) :: time_SI
    real*8, allocatable, intent(in) :: res0D(:)     ! List of 0D quantities defined in exprs_all_int (mod_expressions.f90)
    type(ids_summary),  intent(inout)  :: summary_ids
    character(len=1000)             :: simulation_description
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1
    allocate( summary_ids%time(n_slice) )
    summary_ids%ids_properties%homogeneous_time = 1   
    if (first_step) then 
      allocate(summary_ids%ids_properties%comment(1))
      summary_ids%ids_properties%comment = simulation_description
    endif
    summary_ids%time(i_slice) = time_SI 

    ! --- Information about the toroidal field
    summary_ids%global_quantities%r0%value = R_geo
    allocate(summary_ids%global_quantities%b0%value(n_slice))
    summary_ids%global_quantities%b0%value(i_slice) = F0/R_geo * fact_Ip
    
    do i_exp=1, exprs_all_int%n_expr

      ! --- Beta poloidal
      if (exprs_all_int%expr(i_exp)%name=='beta_p') then
        allocate(summary_ids%global_quantities%beta_pol%value(n_slice))
        summary_ids%global_quantities%beta_pol%value(i_slice) = res0D(i_exp)
      endif

      ! --- Beta toroidal
      if (exprs_all_int%expr(i_exp)%name=='beta_t') then
        allocate(summary_ids%global_quantities%beta_tor%value(n_slice))
        summary_ids%global_quantities%beta_tor%value(i_slice) = res0D(i_exp)
      endif

      ! --- Normalized beta
      if (exprs_all_int%expr(i_exp)%name=='beta_n') then
        allocate(summary_ids%global_quantities%beta_tor_norm%value(n_slice))
        summary_ids%global_quantities%beta_tor_norm%value(i_slice) = res0D(i_exp)
      endif

      ! --- Total current
      if (exprs_all_int%expr(i_exp)%name=='Ip_tot') then
        allocate(summary_ids%global_quantities%ip%value(n_slice))
        summary_ids%global_quantities%ip%value(i_slice) = res0D(i_exp) * fact_Ip
      endif

      ! --- li(3)
      if (exprs_all_int%expr(i_exp)%name=='li3') then
        allocate(summary_ids%global_quantities%li%value(n_slice))
        summary_ids%global_quantities%li%value(i_slice) = res0D(i_exp)
      endif

      ! --- Volume
      if (exprs_all_int%expr(i_exp)%name=='volume') then
        allocate(summary_ids%global_quantities%volume%value(n_slice))
        summary_ids%global_quantities%volume%value(i_slice) = res0D(i_exp)
      endif

      ! --- q_95
      if (exprs_all_int%expr(i_exp)%name=='q95') then
        allocate(summary_ids%global_quantities%q_95%value(n_slice))
        summary_ids%global_quantities%q_95%value(i_slice) = res0D(i_exp)
      endif

      ! --- Thermal energy
      if (exprs_all_int%expr(i_exp)%name=='Thermal_in') then
        allocate(summary_ids%global_quantities%energy_thermal%value(n_slice))
        summary_ids%global_quantities%energy_thermal%value(i_slice) = res0D(i_exp)
      endif

      ! --- Magnetic energy
      if (exprs_all_int%expr(i_exp)%name=='Wmag_in') then
        allocate(summary_ids%global_quantities%energy_b_field_pol%value(n_slice))
        summary_ids%global_quantities%energy_b_field_pol%value(i_slice) = res0D(i_exp)
      endif

      ! --- Ohmic power
      if (exprs_all_int%expr(i_exp)%name=='Ohmic_tot') then
        allocate(summary_ids%global_quantities%power_ohm%value(n_slice))
        summary_ids%global_quantities%power_ohm%value(i_slice) = res0D(i_exp)
      endif

      ! --- Radiated power
      if (exprs_all_int%expr(i_exp)%name=='Rad_tot') then
        allocate(summary_ids%global_quantities%power_radiated%value(n_slice))
        summary_ids%global_quantities%power_radiated%value(i_slice) = res0D(i_exp)
      endif

      ! --- Heating power
      if (exprs_all_int%expr(i_exp)%name=='Heat_src_tot') then
        allocate(summary_ids%heating_current_drive%power_additional%value(n_slice))
        summary_ids%heating_current_drive%power_additional%value(i_slice) = res0D(i_exp)
      endif

      ! --- R_axis
      if (exprs_all_int%expr(i_exp)%name=='R_axis') then
        allocate(summary_ids%local%magnetic_axis%position%r(n_slice))
        summary_ids%local%magnetic_axis%position%r(i_slice) = res0D(i_exp)
      endif

      ! --- Z_axis
      if (exprs_all_int%expr(i_exp)%name=='Z_axis') then
        allocate(summary_ids%local%magnetic_axis%position%z(n_slice))
        summary_ids%local%magnetic_axis%position%z(i_slice) = res0D(i_exp)
      endif

      ! --- psi_axis
      if (exprs_all_int%expr(i_exp)%name=='psi_axis') then
        allocate(summary_ids%local%magnetic_axis%position%psi(n_slice))
        summary_ids%local%magnetic_axis%position%psi(i_slice) = res0D(i_exp) * fact_psi
      endif

    end do

    ! --- Shaping parameters, T. Luce, PPCF 55 (2013) 095009, equations (1-6)
    allocate(summary_ids%boundary%type%value(n_slice))
    if (ES%limiter_plasma) then
      summary_ids%boundary%type%value(i_slice)  = 0
    else
      summary_ids%boundary%type%value(i_slice)  = 1
    endif
    allocate(summary_ids%boundary%minor_radius%value(n_slice))
    summary_ids%boundary%minor_radius%value(i_slice)           = ES%LCFS_a
    allocate(summary_ids%boundary%elongation%value(n_slice))
    summary_ids%boundary%elongation%value(i_slice)             = ES%LCFS_kappa
    allocate(summary_ids%boundary%triangularity_upper%value(n_slice))
    summary_ids%boundary%triangularity_upper%value(i_slice)    = ES%LCFS_deltaU
    allocate(summary_ids%boundary%triangularity_lower%value(n_slice))
    summary_ids%boundary%triangularity_lower%value(i_slice)    = ES%LCFS_deltaL
    allocate(summary_ids%boundary%geometric_axis_r%value(n_slice))
    summary_ids%boundary%geometric_axis_r%value(i_slice)       = ES%LCFS_Rgeo
    allocate(summary_ids%boundary%geometric_axis_z%value(n_slice))
    summary_ids%boundary%geometric_axis_z%value(i_slice)       = ES%LCFS_Zgeo

  end subroutine fill_summary_IDS






  subroutine fill_disruption_IDS(first_step, time_SI, res0D, disruption_ids)  

    implicit none

    ! --- External parameters
    logical,             intent(in) :: first_step   ! is this the first step?
    real*8,              intent(in) :: time_SI
    real*8, allocatable, intent(in) :: res0D(:)     ! List of 0D quantities defined in exprs_all_int (mod_expressions.f90)
    type(ids_disruption),  intent(inout)  :: disruption_ids
   
    ! --- Local parameters 
    integer    :: i, j, k, m, var_rad, i_var, i_tor, index, index_node, my_id, ierr
    real*8     :: vpar_power, kinpar_power
    
    ! **********************************************************************************
    ! ******************************* IMAS **********************************************
    ! **********************************************************************************
    integer :: n_slice, i_slice, i_exp, i_psi
    ! **********************************************************************************

    ! --- Set times
    n_slice = 1;   i_slice = 1
    allocate( disruption_ids%time(n_slice) )
    disruption_ids%ids_properties%homogeneous_time = 1    
    disruption_ids%time(i_slice) = time_SI 

    do i_exp=1, exprs_all_int%n_expr

      ! --- Poloidal halos (defined as positive)
      if (exprs_all_int%expr(i_exp)%name=='I_halo') then
        allocate(disruption_ids%global_quantities%current_halo_pol(n_slice))
        disruption_ids%global_quantities%current_halo_pol(i_slice) = res0D(i_exp) * 1.0d6  ! This diagnostic is in MA
      endif

      ! --- Toroidal halos
      if (exprs_all_int%expr(i_exp)%name=='Ip_out') then
        allocate(disruption_ids%global_quantities%current_halo_tor(n_slice))
        disruption_ids%global_quantities%current_halo_tor(i_slice) = res0D(i_exp) * fact_Ip
      endif

      ! --- Total ohmic power
      if (exprs_all_int%expr(i_exp)%name=='Ohmic_tot') then
        allocate(disruption_ids%global_quantities%power_ohm(n_slice))
        disruption_ids%global_quantities%power_ohm(i_slice) = res0D(i_exp) 
      endif

      ! --- Halo ohmic power
      if (exprs_all_int%expr(i_exp)%name=='Ohmic_out') then
        allocate(disruption_ids%global_quantities%power_ohm_halo(n_slice))
        disruption_ids%global_quantities%power_ohm_halo(i_slice) = res0D(i_exp) 
      endif

      ! --- Total thermal and kinetic power flowing into the wall 
      if (exprs_all_int%expr(i_exp)%name=='qn_par') then
        allocate(disruption_ids%global_quantities%power_parallel_halo(n_slice))
        vpar_power   = 0.d0
        kinpar_power = 0.d0
        ! --- Get other contributions to this power
        do j=1, exprs_all_int%n_expr
          if (exprs_all_int%expr(j)%name=='P_vn') then
            vpar_power = res0D(j) 
          endif
          if (exprs_all_int%expr(j)%name=='kinpar_flux') then
            kinpar_power = res0D(j) 
          endif
        enddo
        disruption_ids%global_quantities%power_parallel_halo(i_slice) = res0D(i_exp) + vpar_power + kinpar_power 
      endif

      ! --- Total radiated power
      if (exprs_all_int%expr(i_exp)%name=='Rad_tot') then
        allocate(disruption_ids%global_quantities%power_radiated_electrons_impurities(n_slice))
        disruption_ids%global_quantities%power_radiated_electrons_impurities(i_slice) = res0D(i_exp) 
      endif

    end do

  end subroutine fill_disruption_IDS






  ! --- Fills Bezier coefficients in GGD
  subroutine fill_Bezier_coefficients( ggd_scalar, node_list, var_index, grid_ind, grid_sub_ind, res_fact )
  
    implicit none
  
    ! --- External parameters
    type(ids_generic_grid_scalar),  intent(inout) ::  ggd_scalar
    type (type_node_list), intent(in)  :: node_list
    integer,               intent(in)  :: var_index, grid_ind, grid_sub_ind
    real*8,                intent(in)  :: res_fact
  
    ! --- Local parameters
    integer :: num_nodes, inode, inode_glob, itor
    
    num_nodes = node_list%n_nodes

    if ( associated(ggd_scalar%coefficients) ) then
      deallocate( ggd_scalar%coefficients )
      allocate( ggd_scalar%coefficients(num_nodes*n_tor, n_degrees) ) 
    else
      allocate( ggd_scalar%coefficients(num_nodes*n_tor, n_degrees) )
    endif

    if ( associated(ggd_scalar%values) ) then
      deallocate( ggd_scalar%values )
      allocate( ggd_scalar%values(num_nodes*n_tor) ) 
    else
      allocate( ggd_scalar%values(num_nodes*n_tor) )
    endif
  
    do inode=1, num_nodes    
      do itor=1, n_tor
        
        inode_glob = inode + (itor-1)*num_nodes 
  
        ggd_scalar%coefficients(inode_glob,:) = node_list%node(inode)%values(itor,:,var_index)
        ggd_scalar%values(inode_glob)         = node_list%node(inode)%values(itor,1,var_index)

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
    allocate(  grid%space(1)%objects_per_dimension(3)  )
    allocate(  grid%space(1)%coordinates_type(2)       )
  
    ! Set coordinates type to [R, Z]
    grid%space(1)%coordinates_type = (/ 4, 3 /)
  
    allocate(grid%identifier%description(1))
    allocate(grid%identifier%name(1))
    grid%identifier%description(1) = "Mesh coming from the JOREK code: combined 2D finite elements space in the poloidal plane &
                                      with Fourier space for the toroidal angle dependence"
    grid%identifier%name  = "JOREK mesh"
    grid%identifier%index = 0   ! Unspecified
  
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
    allocate(grid%grid_subset(1))  
  
    ! Subset for nodes in the combined space
    gs_index = 1
  
    allocate( grid%grid_subset(gs_index)%identifier%name(1)         )
    allocate( grid%grid_subset(gs_index)%identifier%description(1)  )
    grid%grid_subset(gs_index)%identifier%name(1)        = "nodes"
    grid%grid_subset(gs_index)%identifier%index          = 1 
    grid%grid_subset(gs_index)%identifier%description(1) = "The elements of the grid subset are the 0D nodes &
                                                         of the combined RZ x Fourier space (number of nodes &
                                                         is N_poloidal_nodes x N_fourier). "
    grid%grid_subset(gs_index)%dimension                 = 1    ! 1 is the convention for 0D nodes
  
    ! Fill toroidal space 
    space_fourier  => grid%space(2)
    allocate(    space_fourier%coordinates_type(1)    )
    allocate(    space_fourier%identifier%description(1)  )
    space_fourier%coordinates_type(1) = 5          ! The coordinate type is 5, phi angle
    space_fourier%geometry_type%index = n_period   ! Fourier periodicity
    space_fourier%identifier%description(1) = "Toroidal Fourier space"             
  
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
      allocate( space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,n_degrees) )
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(1,:) = node_list%node(i)%x(1,:,1 )   ! R dofs
      space_RZ%objects_per_dimension(1)%object(i)%geometry_2d(2,:) = node_list%node(i)%x(1,:,2 )   ! Z dofs
    enddo
  
    ! JOREK element sizes
    do i=1, num_cells
      allocate( space_RZ%objects_per_dimension(3)%object(i)%geometry_2d(n_degrees, n_vertex_max) )
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





  !> Read one coils set from a STARWALL coil file
  subroutine read_coil_set_starwall(filename, coil_set)
    
    ! --- Routine parameters
    character(len=*), intent(in)              :: filename
    type(t_coil_set_starwall), intent(inout)  :: coil_set
    
    ! --- Local variables
    character(len=256)        :: description
    integer                   :: i, ncoil, err
    integer, parameter        :: IOCH = 167
    type(t_coil_starwall), allocatable :: coil(:)
    
    ! --- Namelists
    namelist /coil_set_nml/ description, ncoil
    namelist /coils_nml   / coil
    
    write(*,*) '  reading coil file =  ', trim(filename)
      
    open(IOCH, file=trim(filename), status='old', action='read', iostat=err)
    if ( err /= 0 ) then
      write(*,*) '  ERROR opening file', trim(filename), '!'
      stop
    end if
      
    ! --- First namelist
    read(IOCH, coil_set_nml)
    coil_set%description = description
    coil_set%ncoil       = ncoil
      
    ! --- Second namelist
    allocate( coil_set%coil(coil_set%ncoil), coil(coil_set%ncoil) )
    read(IOCH, coils_nml)

    ! --- Calculate turns for AXISYM_FILA coils from n_fila_turns
    do i = 1, coil_set%ncoil
      if ( coil(i)%coil_type == AXISYM_FILA ) then
        coil(i)%nturns = sum( abs( coil(i)%n_fila_turns(1:coil(i)%n_fila) ) )
      end if
      if (  coil(i)%coil_type == AXISYM_THICK) then
        if (minval(coil(i)%n_thick_turns(1:coil(i)%nparts_coil)) > -1.d4) then
          coil(i)%nturns = sum( abs(coil(i)%n_thick_turns(1:coil(i)%nparts_coil)) )
          !write(*,*) 'Warning: As n_thick_turns was specified nturns will automatically be calculated from it.'
        else if (coil(i)%nparts_coil .eq. 1) then
          coil(i)%n_thick_turns(1) = coil(i)%nturns
          !write(*,*) 'Warning: The number of turns is now specified by n_thick_turns'
          !write(*,*) 'To model the old behavior n_thick_turns is set to nturns'
        else
          write(*,*) 'Error: When nparts_coil is bigger than 1'
          write(*,*) 'the value the turns has to be given in n_thick_turns'
          stop
        end if
        !write(*,*) 'Coil ',i,' has ', coil(i)%nturns, ' turns'
      end if
    end do
    coil_set%coil(:) = coil(:)
    deallocate(coil)
    
    close(IOCH)
      
    
  end subroutine read_coil_set_starwall






  subroutine initialise_postproc_settings(first_step, n_grid)
    
    implicit none

    logical, intent(in) :: first_step
    integer, intent(in) :: n_grid
    
    character(30)       :: str
    integer             :: ierr

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
      call set_setting('nsub_bnd',         '4',    ierr, 'numerical parameter for field line tracing'         )
    endif
  
  end subroutine initialise_postproc_settings




#endif

end module mod_jorek2IMAS
