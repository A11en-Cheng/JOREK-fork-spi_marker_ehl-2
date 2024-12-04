! Contains all the variables for the equations, to share with each equation routine
module equation_variables

  use mod_parameters
  
  implicit none

    
  ! --- Profiles and sources
  real*8 	:: zn,  dn_dpsi,  dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz
  real*8 	:: zT,  dT_dpsi,  dT_dz,  dT_dpsi2,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi_dz2,  dT_dpsi2_dz
  real*8 	:: current_source, particle_source, heat_source, total_rho_source
  real*8 	:: source_mgi, source_pellet, source_volume
  real*8 	:: Vt0_x, Vt0_y, Omega_tor0_x, Omega_tor0_y, dV_dpsi_source, dV_dz_source
  
  ! --- Diffusivities
  real*8 	:: r0_corr, r0_corr2, T0_corr
  real*8 	:: eta_T,    deta_dT,   d2eta_d2T
  real*8        :: eta_T_ohm, deta_dT_ohm
  real*8 	:: visco_T,  dvisco_dT, d2visco_dT2
  real*8 	:: visco_parr
  real*8 	:: D_prof
  real*8 	:: K_prof, K_par, dK_par
  
  ! --- Coefficients for model500 
  real*8 	:: Dn0x, Dn0y, Dn0p, S_ion, S_ion_T, phi, ksi_ion_norm
  
  ! --- Hyper diffusivities and Taylor Galerkin stabilisation
  real*8 	:: eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm
  real*8	:: TG_num1, TG_num2, TG_num5, TG_num6, TG_num7, TG_num8
  
  ! --- Diamagnetic parameters
  real*8   :: tau_IC
  real*8   :: W_dia, W_dia_rho, W_dia_T

  ! --- Neoclassical coefficients
  real*8   :: epsil, Btheta2
  real*8   :: amu_neo_prof, aki_neo_prof
  
  ! --- Bootstrap current terms
  real*8   :: jb
        
  ! --- R,Z-coords and jacobians
  real*8 	:: x_g, x_s, x_t, x_ss, x_st, x_tt
  real*8 	:: y_g, y_s, y_t, y_ss, y_st, y_tt
  real*8 	:: R, R_x
  real*8 	:: xjac, xjac_x, xjac_y
  
  ! --- Various
  real*8 	:: eps_cyl, theta, zeta
  
  ! --- Deltas
  real*8 	:: delta_g(n_var), delta_s(n_var), delta_t(n_var)
  
  ! --- Test functions
  real*8 	:: v,     v_x,	   v_y,     v_p,     v_s,     v_t,     v_ss,     v_st,     v_tt,     v_xx,     v_yy,     v_xy,     v_pp
  ! --- Variable 1
  real*8 	:: ps0,   ps0_x,   ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,   ps0_tt,   ps0_xx,   ps0_xy,   ps0_yy,   ps0_pp
  real*8 	:: psi,   psi_x,   psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,   psi_tt,   psi_xx,   psi_yy,   psi_xy,   psi_pp
  real*8 	:: delta_ps_x, delta_ps_y
  ! --- Variable 2
  real*8 	:: u0,    u0_x,    u0_y,    u0_p,    u0_s,    u0_t,    u0_ss,    u0_st,    u0_tt,    u0_xx,    u0_xy,    u0_yy,    u0_pp
  real*8 	:: u,     u_x,     u_y,     u_p,     u_s,     u_t,     u_ss,     u_st,     u_tt,     u_xx,     u_xy,     u_yy,     u_pp
  real*8 	:: vv2,  delta_u_x, delta_u_y
  ! --- Variable 3
  real*8 	:: zj0,   zj0_x,   zj0_y,   zj0_p,   zj0_s,   zj0_t,   zj0_ss,   zj0_st,   zj0_tt,   zj0_xx,   zj0_xy,   zj0_yy,   zj0_pp
  real*8 	:: zj,    zj_x,    zj_y,    zj_p,    zj_s,    zj_t,    zj_ss,    zj_st,    zj_tt,    zj_xx,    zj_xy,    zj_yy,    zj_pp
  ! --- Variable 4
  real*8 	:: w0,    w0_x,    w0_y,    w0_p,    w0_s,    w0_t,    w0_ss,    w0_st,    w0_tt,    w0_xx,    w0_xy,    w0_yy,    w0_pp
  real*8 	:: w,     w_x,     w_y,     w_p,     w_s,     w_t,     w_ss,     w_st,     w_tt,     w_xx,     w_xy,     w_yy,     w_pp
  ! --- Variable 5
  real*8 	:: r0,	  r0_x,	   r0_y,    r0_p,    r0_s,    r0_t,    r0_ss,    r0_st,    r0_tt,    r0_xx,    r0_xy,    r0_yy,    r0_pp
  real*8 	:: rho,   rho_x,   rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,   rho_tt,   rho_xx,   rho_xy,   rho_yy,   rho_pp
  real*8 	:: r0_hat,r0_x_hat,r0_y_hat
  ! --- Variable 6
  real*8 	:: T0,    T0_x,    T0_y,    T0_p,    T0_s,    T0_t,    T0_ss,    T0_st,    T0_tt,    T0_xx,    T0_xy,    T0_yy,    T0_pp
  real*8 	:: T,     T_x,     T_y,     T_p,     T_s,     T_t,     T_ss,     T_st,     T_tt,     T_xx,     T_xy,     T_yy,     T_pp
  ! --- Variable 7
  real*8 	:: Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt, Vpar0_xx, Vpar0_xy, Vpar0_yy, Vpar0_pp
  real*8 	:: Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt,  Vpar_xx,  Vpar_xy,  Vpar_yy,  Vpar_pp
  ! --- Variable 8
  real*8 	:: rn0,   rn0_x,   rn0_y,   rn0_p,   rn0_s,   rn0_t,   rn0_ss,   rn0_st,   rn0_tt,   rn0_xx,   rn0_xy,   rn0_yy,   rn0_pp
  real*8 	:: rhon,  rhon_x,  rhon_y,  rhon_p,  rhon_s,  rhon_t,  rhon_ss,  rhon_st,  rhon_tt,  rhon_xx,  rhon_xy,  rhon_yy,  rhon_pp
  ! --- Pressures
  real*8 	:: P0,	  P0_x,	  P0_y,	    P0_p,    P0_s,    P0_t,				     P0_xx,    P0_xy,    P0_yy,    P0_pp
  
  ! --- Parallel gradients
  real*8 	:: BB2, BB2_psi
  real*8 	:: Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star
  real*8 	:: Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star
  
  ! --- Equilibrium (n=0 or n_tor=1) variables
  real*8        :: r00
  real*8        :: T00
  
  ! --- Linearized equation terms
  real*8 	:: rhs_tmp(n_var),        rhs_k_tmp(n_var)
  real*8 	:: amat_tmp(n_var,n_var), amat_k_tmp(n_var,n_var), amat_n_tmp(n_var,n_var), amat_kn_tmp(n_var,n_var)
        
  ! --- !

  ! --- Declare variables as private for each thread (one module for each call to element_matrix)
  !$omp threadprivate(														           	&
  !$omp 	zn,  dn_dpsi,  dn_dz,  dn_dpsi2,  dn_dz2,  dn_dpsi_dz,  dn_dpsi3,  dn_dpsi_dz2,  dn_dpsi2_dz,					&
  !$omp 	zT,  dT_dpsi,  dT_dz,  dT_dpsi2,  dT_dz2,  dT_dpsi_dz,  dT_dpsi3,  dT_dpsi_dz2,  dT_dpsi2_dz,					&
  !$omp 	current_source, particle_source, heat_source, total_rho_source,									&
  !$omp 	source_mgi, source_pellet, source_volume,											&
  !$omp 	Vt0_x, Vt0_y, Omega_tor0_x, Omega_tor0_y, dV_dpsi_source, dV_dz_source,				                                &
  !$omp 	r0_corr, r0_corr2, T0_corr,													&
  !$omp 	eta_T,    deta_dT,   d2eta_d2T,													&
  !$omp         eta_T_ohm, deta_dT_ohm,                                                                                                         & 
  !$omp 	visco_T,  dvisco_dT, d2visco_dT2,												&
  !$omp 	visco_parr,															&
  !$omp 	D_prof,																&
  !$omp 	K_prof, K_par, dK_par,														&
  !$omp 	Dn0x, Dn0y, Dn0p, S_ion, S_ion_T, phi, ksi_ion_norm,											&
  !$omp 	eta_numm, visco_numm, visco_par_numm, D_perp_numm, K_perp_numm,									&
  !$omp         TG_num1, TG_num2, TG_num5, TG_num6, TG_num7, TG_num8,										&
  !$omp         tau_IC, W_dia, W_dia_rho, W_dia_T,												&
  !$omp         epsil, Btheta2,															&
  !$omp         amu_neo_prof, aki_neo_prof,													&
  !$omp 	jb,										   					   	&
  !$omp 	x_g, x_s, x_t, x_ss, x_st, x_tt,												&
  !$omp 	y_g, y_s, y_t, y_ss, y_st, y_tt,												&
  !$omp 	R,  R_x,															&
  !$omp 	xjac, xjac_x, xjac_y,														&
  !$omp 	eps_cyl, theta, zeta,														&
  !$omp 	delta_g, delta_s, delta_t,													&
  !$omp 	v,     v_x,	v_y,	 v_p,	  v_s,     v_t,     v_ss,     v_st,	v_tt,	  v_xx,     v_yy,     v_xy,	v_pp,		&
  !$omp 	ps0,   ps0_x,	ps0_y,   ps0_p,   ps0_s,   ps0_t,   ps0_ss,   ps0_st,	ps0_tt,   ps0_xx,   ps0_xy,   ps0_yy,	ps0_pp,		&
  !$omp 	psi,   psi_x,	psi_y,   psi_p,   psi_s,   psi_t,   psi_ss,   psi_st,	psi_tt,   psi_xx,   psi_yy,   psi_xy,	psi_pp,		&
  !$omp 	delta_ps_x, delta_ps_y,														&
  !$omp 	u0,    u0_x,	u0_y,	 u0_p,    u0_s,    u0_t,    u0_ss,    u0_st,	u0_tt,    u0_xx,    u0_xy,    u0_yy,	u0_pp,		&
  !$omp 	u,     u_x,	u_y,	 u_p,	  u_s,     u_t,     u_ss,     u_st,	u_tt,	  u_xx,     u_xy,     u_yy,	u_pp,		&
  !$omp 	vv2,  delta_u_x, delta_u_y,													&
  !$omp 	zj0,   zj0_x,	zj0_y,   zj0_p,   zj0_s,   zj0_t,   zj0_ss,   zj0_st,	zj0_tt,   zj0_xx,   zj0_xy,   zj0_yy,	zj0_pp,		&
  !$omp 	zj,    zj_x,	zj_y,	 zj_p,    zj_s,    zj_t,    zj_ss,    zj_st,	zj_tt,    zj_xx,    zj_xy,    zj_yy,	zj_pp,		&
  !$omp 	w0,    w0_x,	w0_y,	 w0_p,    w0_s,    w0_t,    w0_ss,    w0_st,	w0_tt,    w0_xx,    w0_xy,    w0_yy,	w0_pp,		&
  !$omp 	w,     w_x,	w_y,	 w_p,	  w_s,     w_t,     w_ss,     w_st,	w_tt,	  w_xx,     w_xy,     w_yy,	w_pp,		&
  !$omp 	r0,    r0_x,	r0_y,	 r0_p,    r0_s,    r0_t,    r0_ss,    r0_st,	r0_tt,    r0_xx,    r0_xy,    r0_yy,	r0_pp,		&
  !$omp 	rho,   rho_x,	rho_y,   rho_p,   rho_s,   rho_t,   rho_ss,   rho_st,	rho_tt,   rho_xx,   rho_xy,   rho_yy,	rho_pp,		&
  !$omp 	r0_hat,r0_x_hat,r0_y_hat,													&
  !$omp 	T0,    T0_x,	T0_y,	 T0_p,    T0_s,    T0_t,    T0_ss,    T0_st,	T0_tt,    T0_xx,    T0_xy,    T0_yy,	T0_pp,		&
  !$omp 	T,     T_x,	T_y,	 T_p,	  T_s,     T_t,     T_ss,     T_st,	T_tt,	  T_xx,     T_xy,     T_yy,	T_pp,		&
  !$omp 	Vpar0, Vpar0_x, Vpar0_y, Vpar0_p, Vpar0_s, Vpar0_t, Vpar0_ss, Vpar0_st, Vpar0_tt, Vpar0_xx, Vpar0_xy, Vpar0_yy, Vpar0_pp,	&
  !$omp 	Vpar,  Vpar_x,  Vpar_y,  Vpar_p,  Vpar_s,  Vpar_t,  Vpar_ss,  Vpar_st,  Vpar_tt,  Vpar_xx,  Vpar_xy,  Vpar_yy,  Vpar_pp,	&
  !$omp 	rn0,   rn0_x,	rn0_y,   rn0_p,   rn0_s,   rn0_t,   rn0_ss,   rn0_st,	rn0_tt,   rn0_xx,   rn0_xy,   rn0_yy,	rn0_pp,		&
  !$omp 	rhon,  rhon_x,  rhon_y,  rhon_p,  rhon_s,  rhon_t,  rhon_ss,  rhon_st,  rhon_tt,  rhon_xx,  rhon_xy,  rhon_yy,  rhon_pp,	&
  !$omp 	P0,    P0_x,   P0_y,	 P0_p,    P0_s,    P0_t,				  P0_xx,    P0_xy,    P0_yy,	P0_pp,		&
  !$omp 	BB2, BB2_psi,															&
  !$omp 	Bgrad_rho, Bgrad_rho_star, Bgrad_rho_k_star,											&
  !$omp 	Bgrad_T,   Bgrad_T_star,   Bgrad_T_k_star,											&
  !$omp 	r00,										   					   	&
  !$omp 	T00,										   					   	&
  !$omp 	rhs_tmp,        rhs_k_tmp,													&
  !$omp 	amat_tmp, amat_k_tmp, amat_n_tmp, amat_kn_tmp)


end module equation_variables
