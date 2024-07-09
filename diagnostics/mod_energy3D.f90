!> This diagnostic module can be used to calculate magnetic and kinetic energies in individual 
!! mode families (see: C. Schwab 1993 Phys. Fluids B 5 3195-206 Section III). 
!!
!! This routine is intended for use with stellarator models. It can also be used as part of
!! jorek2_postproc. 
module mod_energy3D

  implicit none
  
  private
  
  public :: energy3d_new 
  
  contains


!> This subroutine calculates the energies of each individual mode
!! family in a stellarator by summing up the energies of all modes
!! and cross terms in the mode family.
subroutine energy3d_new(node_list,element_list,W_mag,W_kin)
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use mod_chi
implicit none

type(type_node_list),                       intent(in)  :: node_list
type(type_element_list),                    intent(in)  :: element_list
real*8, dimension(1+int(n_coord_period/2)), intent(out) :: W_mag, W_kin

type(type_element)      :: element
type(type_node)         :: nodes(n_vertex_max)

real*8, dimension(n_plane,n_gauss,n_gauss)                     :: x_g, x_s, x_t, x_p, y_g, y_s, y_t, y_p
real*8, dimension(n_tor,n_var,n_gauss,n_gauss)                 :: eq_g, eq_s, eq_t
real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1,n_plane) :: chi
real*8, dimension(1+int(n_coord_period/2))                     :: W_mag_local, W_kin_local

real*8  :: phi, BigR, xjac, wst, psi_x, psi_y, psi_p, BR1, Bz1, Bp1, BR2, Bz2, Bp2
real*8  :: Phi_x, Phi_y, Phi_p, vR1, vz1, vp1, vR2, vz2, vp2, Bv2(n_plane)
integer :: ife, iv, i, j, k, ms, mt, mp, in, im, imf

W_mag = 0.d0; W_kin = 0.d0

!$omp parallel default(none) &
!$omp shared(node_list, element_list, W_mag, W_kin, H, H_s, H_t, HZ, HZ_p, HZ_coord, HZ_coord_p, F0, mode) &
!$omp private(element, nodes, x_g, x_s, x_t, x_p, y_g, y_s, y_t, y_p, eq_g, eq_s, eq_t, chi, phi, BigR, xjac, wst, psi_x, psi_y, &
!$omp         psi_p, BR1, Bz1, Bp1, BR2, Bz2, Bp2, Phi_x, Phi_y, Phi_p, vR1, vz1, vp1, vR2, vz2, vp2, Bv2, ife, iv, i, j, k, ms, &
!$omp         mt, mp, in, im, imf, W_mag_local, W_kin_local)

W_mag_local = 0.d0; W_kin_local = 0.d0

!$omp do
do ife=1,element_list%n_elements
  element = element_list%element(ife)
  do iv=1,n_vertex_max
    nodes(iv) = node_list%node(element%vertex(iv))
  end do

  x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; x_p = 0.d0; y_g = 0.d0; y_s = 0.d0; y_t = 0.d0; y_p = 0.d0
  eq_g = 0.d0; eq_s = 0.d0; eq_t = 0.d0

  do i=1,n_vertex_max
    do j=1,n_order+1
      do ms=1,n_gauss
        do mt=1,n_gauss
          do mp=1,n_plane
            do in=1,n_coord_tor
              x_g(mp,ms,mt) = x_g(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H(i,j,ms,mt)*  HZ_coord(in,mp)
              x_s(mp,ms,mt) = x_s(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_s(i,j,ms,mt)*HZ_coord(in,mp)
              x_t(mp,ms,mt) = x_t(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H_t(i,j,ms,mt)*HZ_coord(in,mp)
              x_p(mp,ms,mt) = x_p(mp,ms,mt) + nodes(i)%x(in,j,1)*element%size(i,j)*H(i,j,ms,mt)*  HZ_coord_p(in,mp)

              y_g(mp,ms,mt) = y_g(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H(i,j,ms,mt)*  HZ_coord(in,mp)
              y_s(mp,ms,mt) = y_s(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_s(i,j,ms,mt)*HZ_coord(in,mp)
              y_t(mp,ms,mt) = y_t(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H_t(i,j,ms,mt)*HZ_coord(in,mp)
              y_p(mp,ms,mt) = y_p(mp,ms,mt) + nodes(i)%x(in,j,2)*element%size(i,j)*H(i,j,ms,mt)*  HZ_coord_p(in,mp)
            end do
          end do
          
          do k=1,n_var
            do in=1,n_tor
              eq_g(in,k,ms,mt) = eq_g(in,k,ms,mt) + nodes(i)%values(in,j,k)*element%size(i,j)*H(i,j,ms,mt)
              eq_s(in,k,ms,mt) = eq_s(in,k,ms,mt) + nodes(i)%values(in,j,k)*element%size(i,j)*H_s(i,j,ms,mt)
              eq_t(in,k,ms,mt) = eq_t(in,k,ms,mt) + nodes(i)%values(in,j,k)*element%size(i,j)*H_t(i,j,ms,mt)
            end do
          end do
        end do
      end do
    end do
  end do
  
  do ms=1,n_gauss
    do mt=1,n_gauss
      wst = wgauss(ms)*wgauss(mt)
      do mp=1,n_plane
        phi = 2.d0*pi*float(mp-1)/float(n_period*n_plane)
        chi(:,:,:,mp) = get_chi(x_g(mp,ms,mt),y_g(mp,ms,mt),phi,node_list,element_list,i,xgauss(ms),xgauss(mt),1)
        Bv2(mp) = chi(1,0,0,mp)**2 + chi(0,1,0,mp)**2 + chi(0,0,1,mp)**2/x_g(mp,ms,mt)**2
      end do
      
      do in=1,n_tor
        do im=1,n_tor
          ! Determine if modes in and im are part of the same mode family
          if (mod(mode(in)+mode(im),n_coord_period) .eq. 0 .or. mod(mode(in)-mode(im),n_coord_period) .eq. 0) then
            ! If yes, then determine the mode family imf that the modes belong to
            imf = min(mode(in),mode(im))
            do while (imf .ge. n_coord_period)
              imf = imf - n_coord_period
            end do
            if (imf .gt. int(n_coord_period/2)) imf = n_coord_period - imf
          else
            cycle
          end if
          
          do mp=1,n_plane
            BigR = x_g(mp,ms,mt)
            xjac = x_s(mp,ms,mt)*y_t(mp,ms,mt) - x_t(mp,ms,mt)*y_s(mp,ms,mt)

            psi_x = ( y_t(mp,ms,mt)*eq_s(in,var_psi,ms,mt) - y_s(mp,ms,mt)*eq_t(in,var_psi,ms,mt))*HZ(in,mp)/xjac
            psi_y = (-x_t(mp,ms,mt)*eq_s(in,var_psi,ms,mt) + x_s(mp,ms,mt)*eq_t(in,var_psi,ms,mt))*HZ(in,mp)/xjac
            psi_p = eq_g(in,var_psi,ms,mt)*HZ_p(in,mp) - psi_x*x_p(mp,ms,mt) - psi_y*y_p(mp,ms,mt)
            Phi_x = ( y_t(mp,ms,mt)*eq_s(in,var_u,ms,mt) - y_s(mp,ms,mt)*eq_t(in,var_u,ms,mt))*HZ(in,mp)/xjac
            Phi_y = (-x_t(mp,ms,mt)*eq_s(in,var_u,ms,mt) + x_s(mp,ms,mt)*eq_t(in,var_u,ms,mt))*HZ(in,mp)/xjac
            Phi_p = eq_g(in,var_u,ms,mt)*HZ_p(in,mp) - Phi_x*x_p(mp,ms,mt) - Phi_y*y_p(mp,ms,mt)
          
            BR1 =  (psi_y*chi(0,0,1,mp) - psi_p*chi(0,1,0,mp))/(F0*BigR)
            Bz1 = -(psi_x*chi(0,0,1,mp) - psi_p*chi(1,0,0,mp))/(F0*BigR)
            Bp1 =  (psi_x*chi(0,1,0,mp) - psi_y*chi(1,0,0,mp))/F0
            vR1 =  (Phi_y*chi(0,0,1,mp) - Phi_p*chi(0,1,0,mp))/(Bv2(mp)*BigR)
            vz1 = -(Phi_x*chi(0,0,1,mp) - Phi_p*chi(1,0,0,mp))/(Bv2(mp)*BigR)
            vp1 =  (Phi_x*chi(0,1,0,mp) - Phi_y*chi(1,0,0,mp))/Bv2(mp)
            
            psi_x = ( y_t(mp,ms,mt)*eq_s(im,var_psi,ms,mt) - y_s(mp,ms,mt)*eq_t(im,var_psi,ms,mt))*HZ(im,mp)/xjac
            psi_y = (-x_t(mp,ms,mt)*eq_s(im,var_psi,ms,mt) + x_s(mp,ms,mt)*eq_t(im,var_psi,ms,mt))*HZ(im,mp)/xjac
            psi_p = eq_g(im,var_psi,ms,mt)*HZ_p(im,mp) - psi_x*x_p(mp,ms,mt) - psi_y*y_p(mp,ms,mt)
            Phi_x = ( y_t(mp,ms,mt)*eq_s(im,var_u,ms,mt) - y_s(mp,ms,mt)*eq_t(im,var_u,ms,mt))*HZ(im,mp)/xjac
            Phi_y = (-x_t(mp,ms,mt)*eq_s(im,var_u,ms,mt) + x_s(mp,ms,mt)*eq_t(im,var_u,ms,mt))*HZ(im,mp)/xjac
            Phi_p = eq_g(im,var_u,ms,mt)*HZ_p(im,mp) - Phi_x*x_p(mp,ms,mt) - Phi_y*y_p(mp,ms,mt)
          
            BR2 =  (psi_y*chi(0,0,1,mp) - psi_p*chi(0,1,0,mp))/(F0*BigR)
            Bz2 = -(psi_x*chi(0,0,1,mp) - psi_p*chi(1,0,0,mp))/(F0*BigR)
            Bp2 =  (psi_x*chi(0,1,0,mp) - psi_y*chi(1,0,0,mp))/F0
            vR2 =  (Phi_y*chi(0,0,1,mp) - Phi_p*chi(0,1,0,mp))/(Bv2(mp)*BigR)
            vz2 = -(Phi_x*chi(0,0,1,mp) - Phi_p*chi(1,0,0,mp))/(Bv2(mp)*BigR)
            vp2 =  (Phi_x*chi(0,1,0,mp) - Phi_y*chi(1,0,0,mp))/Bv2(mp)
            
            W_mag_local(imf+1) = W_mag_local(imf+1) + (BR1*BR2 + Bz1*Bz2 + Bp1*Bp2)*BigR*xjac*wst*pi/n_plane
            W_kin_local(imf+1) = W_kin_local(imf+1) + eq_g(1,var_rho,ms,mt)*(vR1*vR2 + vz1*vz2 + vp1*vp2)*BigR*xjac*wst*pi/n_plane
          end do
        end do
      end do
    end do
  end do
end do
!$omp end do

!$omp critical
W_mag = W_mag + W_mag_local
W_kin = W_kin + W_kin_local
!$omp end critical
!$omp end parallel
end subroutine energy3d_new

end module mod_energy3D
