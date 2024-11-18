!> Model-specific part of the routine init_live_data in communication/mod_live_data.f90
!!
!! Writes all input profiles to 'macroscopic_vars.dat'.
subroutine init_live_data_model(file_handle)
  
  use phys_module,   only: xpoint, xcase
  use diffusivities, only: get_dperp, get_zkperp, get_zk_iperp, get_zk_eperp
  use mod_sources
  use mod_equations
 
  implicit none
  
  integer, intent(in) :: file_handle
  
  integer :: i
  real*8  :: psin, FFp, dFFp_dpsi, dens, dn_dpsi, temp, dT_dpsi, S_rho, S_T, S_Ti, S_Te, d_perp, zk_perp, zk_i_perp, zk_e_perp, Phi
  real*8  :: d, d1, d2, d3, d4, d5, d6, d7, d8, d9  ! dummies
  
  if (with_TiTe) then
    write(file_handle,'(A,I5)') '@n_input_profiles: ', 9
  else
    write(file_handle,'(A,I5)') '@n_input_profiles: ', 7
  end if
  write(file_handle,'(A)') '@input_profiles_xlabel: Psi_{normalized}'
  write(file_handle,'(A)') '@input_profiles_ylabel: input profiles'
  write(file_handle,'(A)') '@input_profiles_logy: 0'
  if (with_TiTe) then
    write(file_handle,'(A)') '@input_profiles: %"psin"       "rho"    "drho/dpsin"   "S_rho"     "S_Ti"   "S_Te"   "D_perp"   "ZK_i_perp"     "ZK_e_perp"      "Phi"'
  else
    write(file_handle,'(A)') '@input_profiles: %"psin"       "rho"    "drho/dpsin"   "S_rho"     "S_T"   "D_perp"    "ZK_perp"     "Phi"'
  end if

  do i = 0, 200
    
    psin = real(i) / real(200) * 1.2d0

    call density    (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,dens,dn_dpsi,d,d1,d2,d3,d4,d5,d6)
    call potential_source   (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0, Phi, d, d1, d2, d3, d4, d5, d6, d7)
    if (with_TiTe) then
      call sources    (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,S_rho,S_Ti,S_Te)
      d_perp  = get_dperp (psin)
      zk_i_perp = get_zk_iperp(psin)
      zk_e_perp = get_zk_eperp(psin)
      write(file_handle,'(a,20es13.4e3)') '@input_profiles: ', psin, dens, dn_dpsi, S_rho, S_Ti, S_Te, d_perp, zk_i_perp, zk_e_perp, Phi
    else
      call sources    (xpoint,xcase,0.d0,(/-99.d0,-99.d0/),psin,0.d0,1.d0,S_rho,S_T)
      d_perp  = get_dperp (psin)
      zk_perp = get_zkperp(psin)
      write(file_handle,'(a,20es13.4e3)') '@input_profiles: ', psin, dens, dn_dpsi, S_rho, S_T, d_perp, zk_perp, Phi
    end if

    write(63,*) psin, S_rho
    
  end do
  write(file_handle,*)
  
end subroutine init_live_data_model
