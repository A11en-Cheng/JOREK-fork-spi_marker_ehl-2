subroutine potential_source(xpoint2,xcase2,Z,Z_xpoint,psi,psi_axis,psi_bnd,potential_profile, &
                       dPhi_dpsi,dPhi_dz,dPhi_dpsi2,dPhi_dz2,dPhi_dpsi_dz,dPhi_dpsi3,dPhi_dpsi_dz2, dPhi_dpsi2_dz)
use phys_module
use vacuum, only: current_FB_fact

implicit none

! --- Routine parameters
logical, intent(in)  :: xpoint2
integer, intent(in)  :: xcase2
real*8,  intent(in)  :: Z, Z_xpoint(2), psi, psi_axis, psi_bnd
real*8,  intent(out) :: potential_profile, dPhi_dpsi, dPhi_dz, dPhi_dpsi2, dPhi_dz2, &
                        dPhi_dpsi_dz, dPhi_dpsi3, dPhi_dpsi_dz2, dPhi_dpsi2_dz

! --- Internal variables.
real*8  :: prof0, prof1, dprof0_dpsi, dprof0_dpsi2, dprof0_dpsi3, psi_barrier
real*8  :: psi_n, psi_star, delta_psi, sig_phi, sigz, dprof1_dpsi, dprof1_dpsi2, dprof1_dpsi3
real*8  :: atn, datn, d2atn, d3atn
real*8  :: atn_z,   datn_z,   d2atn_z
real*8  :: atn_z_u, datn_z_u, d2atn_z_u, factor
real*8  :: cosh1, cosh2, cosh3, cosh3_u
real*8  :: tanh1, tanh2, tanh2_u
! for interpolating numerical profiles
integer :: left, right, mid
real*8  :: aux1, aux2, Z_star, Z_star_u

delta_psi = psi_bnd - psi_axis
psi_n     = (psi - psi_axis) / delta_psi

psi_n = max( min(psi_n, 2.), 0. )


if (.not. num_Phi) then
  prof0        = (Phi_0-phi_1)*(1.d0 + phi_coef(1) * psi_n + phi_coef(2) * psi_n**2 + phi_coef(3) * psi_n**3)
  dprof0_dpsi  = (phi_0-phi_1)*(phi_coef(1) + 2.d0 * phi_coef(2) * psi_n + 3.d0 * phi_coef(3) * psi_n**2) / delta_psi
  dprof0_dpsi2 = (phi_0-phi_1)*(2.d0 * phi_coef(2) + 6.d0 * phi_coef(3) * psi_n)                          / delta_psi**2
  dprof0_dpsi3 = (phi_0-phi_1)*(6.d0 * phi_coef(3))                                                       / delta_psi**3
  
  sig_phi     = phi_coef(4)
  psi_barrier = phi_coef(5)
  
  psi_star = (psi_n - psi_barrier)/sig_phi
  psi_star = min( max( psi_star, -40.d0), 40.d0) ! avoid floating-point exceptions
  
  tanh1 = tanh(psi_star)
  cosh1 = cosh(psi_star)
  cosh2 = cosh(2.d0*psi_star)
  
  atn   = (0.5d0 - 0.5d0*tanh1)
  datn  = - 1.d0/cosh1**2 / (2.d0 * sig_phi) / delta_psi
  d2atn =   1.d0/cosh1**2 / sig_phi**2 * tanh1 / delta_psi**2
  d3atn = - 1.d0/cosh1**4 / sig_phi**3 * (-2.d0 + cosh2) / delta_psi**3
  
  potential_profile = prof0        * atn
  dPhi_dpsi         = dprof0_dpsi  * atn +         prof0       * datn
  dPhi_dpsi2        = dprof0_dpsi2 * atn + 2.d0 * dprof0_dpsi  * datn + prof0              * d2atn
  dPhi_dpsi3        = dprof0_dpsi3 * atn + 3.d0 * dprof0_dpsi2 * datn + 3.d0 * dprof0_dpsi * d2atn + prof0 * d3atn
  dPhi_dz2          = 0.0
  dPhi_dpsi_dz      = 0.0
  dPhi_dpsi_dz2     = 0.0
  dPhi_dpsi2_dz     = 0.0

else
  left = 1
  right = num_Phi_len
  do
    if(right == left + 1) exit
    mid = (left+right)/2
    if (num_Phi_x(mid) >= psi_n) then
            right = mid
    else
            left = mid
    end if                    
  end do        
  !weights
  aux1 = (psi_n - num_Phi_x(left))/(num_Phi_x(right) - num_Phi_x(left)) 
  aux2 = (1. - aux1)
  potential_profile =  num_Phi_y0(left)*aux2 + num_Phi_y0(right)*aux1 !Phi at the desired point
  dPhi_dpsi         = (num_Phi_y1(left)*aux2 + num_Phi_y1(right)*aux1) / delta_psi !first derivative at desired point
  dPhi_dpsi2        = (num_Phi_y2(left)*aux2 + num_Phi_y2(right)*aux1) / delta_psi**2
  dPhi_dpsi3        = (num_Phi_y3(left)*aux2 + num_Phi_y3(right)*aux1) / delta_psi**3
  dPhi_dz2          = 0.0
  dPhi_dpsi_dz      = 0.0
  dPhi_dpsi_dz2     = 0.0
  dPhi_dpsi2_dz     = 0.0
end if

potential_profile = potential_profile + phi_1

return
end subroutine potential_source
