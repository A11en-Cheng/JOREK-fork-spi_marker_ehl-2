module mod_boundary_matrix_open
  implicit none
contains
subroutine boundary_matrix_open(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                psi_bnd, R_xpoint, Z_xpoint, ELM, RHS, i_tor_min, i_tor_max)
!-----------------------------------------------------------------------------------------------------
! implement consistent (i.e. nonzero) boundary conditions for zj
! the current in the grad(chi) direction on the boundary must match the corresponding quantity in GVEC
!-----------------------------------------------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use nodes_elements
use mod_chi
use mod_basisfunctions
implicit none

type(type_element)   :: element
type(type_node)      :: nodes(2)        ! the two nodes containing the boundary nodes

real*8, dimension(:,:), allocatable  :: ELM
real*8, dimension(:),   allocatable  :: RHS
integer,                intent(in)   :: i_tor_min
integer,                intent(in)   :: i_tor_max

integer    :: vertex(2), direction(2), xcase2
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
logical    :: xpoint2

integer :: direction_perp(2), i, i2, j, j2, j3, k, l, l2, ms, mp, in, im, n_tor_local, index_ij, ij3, index_kl, kl1, kl3
real*8  :: theta, zeta, zbig, element_size_ij, element_size_kl, element_size_perp, BigR, phi, grad_chi(3), Bv2
real*8  :: j_gvec(3), Bx_x_gvec, Bx_y_gvec, Bx_phi_gvec, By_x_gvec, By_y_gvec, By_phi_gvec
real*8  :: xjac, dl, v, rhs_ij_3, amat_31, amat_33, Bp_x_gvec, Bp_y_gvec, Bp_phi_gvec, zj

real*8, dimension(n_plane,n_gauss) :: x_g, x_s, x_t, x_p, y_g, y_s, y_t, y_p
real*8, dimension(n_plane,n_gauss) :: Bx_s_gvec, Bx_t_gvec, Bx_p_gvec, By_s_gvec, By_t_gvec, By_p_gvec, Bp_s_gvec, Bp_t_gvec, Bp_p_gvec

real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

type(type_node)         :: tmp_node

theta = time_evol_theta
!zeta  = time_evol_zeta
! change zeta for variable dt
zeta  = time_evol_zeta * 2.0d0 * tstep / (tstep + tstep_prev)
Zbig = 1.d12

!--------------------- reorder the nodes to have the same direction as full element (maybe not necesary)
if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 4)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 4
  vertex(2) = 3
else if ((vertex(1) .eq. 4) .and. (vertex(2) .eq. 1)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 4
else if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 2)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 2
  vertex(2) = 3
else if ((vertex(1) .eq. 2) .and. (vertex(2) .eq. 1)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 2
end if

x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; x_p = 0.d0
y_g = 0.d0; y_s = 0.d0; y_t = 0.d0; y_p = 0.d0
Bx_s_gvec = 0.d0; Bx_t_gvec = 0.d0; Bx_p_gvec = 0.d0; By_s_gvec = 0.d0; By_t_gvec = 0.d0; By_p_gvec = 0.d0
Bp_s_gvec = 0.d0; Bp_t_gvec = 0.d0; Bp_p_gvec = 0.d0

direction_perp(1) = 6/direction(2)     ! =3 if direction(2)=2, =2 if direction(2)=3
direction_perp(2) = 4

do i=1,2    ! sum over 2 verices
  do j=1,2  ! sum over two basis functions
    i2 = vertex(i)
    j2 = direction(j)
    element_size_ij = element%size(i2,j2)

    j3 = direction_perp(j)
    if (vertex(1) .eq. 1) then
      element_size_perp = 3.d0*element%size(i2,direction_perp(1))
    else
      element_size_perp = -3.d0*element%size(i2,direction_perp(1))
    end if
    
    do ms=1,n_gauss
      do mp=1,n_plane
        do in=1,n_coord_tor
          x_g(mp,ms) = x_g(mp,ms) + nodes(i)%x(in,j2,1)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)
          x_s(mp,ms) = x_s(mp,ms) + nodes(i)%x(in,j2,1)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          x_t(mp,ms) = x_t(mp,ms) + nodes(i)%x(in,j3,1)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          x_p(mp,ms) = x_p(mp,ms) + nodes(i)%x(in,j2,1)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)

          y_g(mp,ms) = y_g(mp,ms) + nodes(i)%x(in,j2,2)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)
          y_s(mp,ms) = y_s(mp,ms) + nodes(i)%x(in,j2,2)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          y_t(mp,ms) = y_t(mp,ms) + nodes(i)%x(in,j3,2)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          y_p(mp,ms) = y_p(mp,ms) + nodes(i)%x(in,j2,2)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)

          Bx_s_gvec(mp,ms) = Bx_s_gvec(mp,ms) + nodes(i)%b_field(in,j2,1)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          Bx_t_gvec(mp,ms) = Bx_t_gvec(mp,ms) + nodes(i)%b_field(in,j3,1)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          Bx_p_gvec(mp,ms) = Bx_p_gvec(mp,ms) + nodes(i)%b_field(in,j2,1)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)
          By_s_gvec(mp,ms) = By_s_gvec(mp,ms) + nodes(i)%b_field(in,j2,2)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          By_t_gvec(mp,ms) = By_t_gvec(mp,ms) + nodes(i)%b_field(in,j3,2)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          By_p_gvec(mp,ms) = By_p_gvec(mp,ms) + nodes(i)%b_field(in,j2,2)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)
          Bp_s_gvec(mp,ms) = Bp_s_gvec(mp,ms) + nodes(i)%b_field(in,j2,3)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          Bp_t_gvec(mp,ms) = Bp_t_gvec(mp,ms) + nodes(i)%b_field(in,j3,3)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          Bp_p_gvec(mp,ms) = Bp_p_gvec(mp,ms) + nodes(i)%b_field(in,j2,3)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)
        end do
      end do
    end do
  end do
end do

n_tor_local = i_tor_max - i_tor_min + 1
do ms=1,n_gauss
  do mp=1,n_plane
    BigR = x_g(mp,ms)
    phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
    chi = get_chi_domm(x_g(mp,ms),y_g(mp,ms),phi)
    grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/BigR /)
    Bv2 = dot_product(grad_chi,grad_chi)

    dl = sqrt(x_s(mp,ms)**2 + y_s(mp,ms)**2)
    xjac = x_t(mp,ms)*y_s(mp,ms) - x_s(mp,ms)*y_t(mp,ms)

    Bx_x_gvec = (Bx_t_gvec(mp,ms)*y_s(mp,ms) - Bx_s_gvec(mp,ms)*y_t(mp,ms))/xjac
    Bx_y_gvec = (Bx_s_gvec(mp,ms)*x_t(mp,ms) - Bx_t_gvec(mp,ms)*x_s(mp,ms))/xjac
    Bx_phi_gvec = Bx_p_gvec(mp,ms) - Bx_x_gvec*x_p(mp,ms) - Bx_y_gvec*y_p(mp,ms)
    By_x_gvec = (By_t_gvec(mp,ms)*y_s(mp,ms) - By_s_gvec(mp,ms)*y_t(mp,ms))/xjac
    By_y_gvec = (By_s_gvec(mp,ms)*x_t(mp,ms) - By_t_gvec(mp,ms)*x_s(mp,ms))/xjac
    By_phi_gvec = By_p_gvec(mp,ms) - By_x_gvec*x_p(mp,ms) - By_y_gvec*y_p(mp,ms)
    Bp_x_gvec = (Bp_t_gvec(mp,ms)*y_s(mp,ms) - Bp_s_gvec(mp,ms)*y_t(mp,ms))/xjac
    Bp_y_gvec = (Bp_s_gvec(mp,ms)*x_t(mp,ms) - Bp_t_gvec(mp,ms)*x_s(mp,ms))/xjac
    Bp_phi_gvec = Bp_p_gvec(mp,ms) - Bp_x_gvec*x_p(mp,ms) - Bp_y_gvec*y_p(mp,ms)
    j_gvec = (/ Bp_y_gvec - By_phi_gvec/BigR, -Bp_x_gvec + Bx_phi_gvec/BigR, By_x_gvec - Bx_y_gvec /)
    
    do i=1,2
      do j=1,2
        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)
        
        do im=i_tor_min,i_tor_max
          index_ij = n_tor_local*n_var*(n_order+1)*(vertex(i)-1) + n_tor_local * n_var * (j2-1) + im - i_tor_min + 1  ! index in the ELM matrix
          
          v = H1(i,j,ms)*element_size_ij*HZ(im,mp)
          rhs_ij_3 = -v*dot_product(grad_chi,j_gvec)*dl*zbig
          
          ij3 = index_ij + 2*n_tor_local
          
          RHS(ij3) = RHS(ij3) + rhs_ij_3*wgauss(ms)
          
          do k=1,2
            do l=1,2
              l2 = direction(l)
              element_size_kl = element%size(vertex(k),l2)
            
              do in=i_tor_min,i_tor_max
                index_kl = n_tor_local*n_var*(n_order+1)*(vertex(k)-1) + n_tor_local * n_var * (l2-1) + in - i_tor_min +1  ! index in the ELM matrix

                zj = H1(k,l,ms)*element_size_kl*HZ(in,mp)
                
                amat_33 = v*Bv2*zj*dl*zbig/F0
                
                kl3 = index_kl + 2*n_tor_local
                
                ELM(ij3,kl3) = ELM(ij3,kl3) + amat_33*wgauss(ms)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

return
end subroutine boundary_matrix_open


subroutine boundary_matrix_open_chi_correction(vertex, direction, element, nodes, xpoint2, xcase2, R_axis, Z_axis, psi_axis, &
                                               psi_bnd, R_xpoint, Z_xpoint, ELM, RHS)
!-----------------------------------------------------------------------------------------------------
! implement consistent (i.e. nonzero) boundary conditions for zj
! the current in the grad(chi) direction on the boundary must match the corresponding quantity in GVEC
!-----------------------------------------------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian
use phys_module
use nodes_elements
use mod_chi
use mod_basisfunctions
implicit none

type(type_element)   :: element
type(type_node)      :: nodes(2)        ! the two nodes containing the boundary nodes

real*8     :: ELM(n_vertex_max*(n_order+1)*n_coord_tor, n_vertex_max*(n_order+1)*n_coord_tor)
real*8     :: RHS(n_vertex_max*(n_order+1)*n_coord_tor)

integer    :: vertex(2), direction(2), xcase2
real*8     :: psi_axis, R_axis, Z_axis, psi_bnd, R_xpoint(2), Z_xpoint(2)
logical    :: xpoint2

integer :: direction_perp(2), i, i2, j, j2, j3, k, l, l2, ms, mp, in, i_tor, index_ij, index_kl, kl1, kl3
real*8  :: theta, zeta, zbig, element_size_ij, element_size_kl, element_size_perp, BigR, phi, grad_chi(3), Bv2
real*8  :: dA, v, rhs_ij, n_perp(3)

real*8, dimension(n_plane,n_gauss) :: x_g, x_s, x_t, x_p, y_g, y_s, y_t, y_p

real*8, dimension(0:n_order-1,0:n_order-1,0:n_order-1) :: chi

type(type_node)         :: tmp_node

!--------------------- reorder the nodes to have the same direction as full element (maybe not necesary)
if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 4)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 4
  vertex(2) = 3
else if ((vertex(1) .eq. 4) .and. (vertex(2) .eq. 1)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 4
else if ((vertex(1) .eq. 3) .and. (vertex(2) .eq. 2)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 2
  vertex(2) = 3
else if ((vertex(1) .eq. 2) .and. (vertex(2) .eq. 1)) then
  tmP_node = nodes(1)
  nodes(1)  = nodes(2)
  nodes(2)  = tmp_node
  vertex(1) = 1
  vertex(2) = 2
end if

x_g = 0.d0; x_s = 0.d0; x_t = 0.d0; x_p = 0.d0
y_g = 0.d0; y_s = 0.d0; y_t = 0.d0; y_p = 0.d0

direction_perp(1) = 6/direction(2)     ! =3 if direction(2)=2, =2 if direction(2)=3
direction_perp(2) = 4

do i=1,2    ! sum over 2 verices
  do j=1,2  ! sum over two basis functions
    i2 = vertex(i)
    j2 = direction(j)
    element_size_ij = element%size(i2,j2)

    j3 = direction_perp(j)
    if (vertex(1) .eq. 1) then
      element_size_perp = 3.d0*element%size(i2,direction_perp(1))
    else
      element_size_perp = -3.d0*element%size(i2,direction_perp(1))
    end if
    
    do ms=1,n_gauss
      do mp=1,n_plane
        do in=1,n_coord_tor
          x_g(mp,ms) = x_g(mp,ms) + nodes(i)%x(in,j2,1)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)
          x_s(mp,ms) = x_s(mp,ms) + nodes(i)%x(in,j2,1)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          x_t(mp,ms) = x_t(mp,ms) + nodes(i)%x(in,j3,1)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          x_p(mp,ms) = x_p(mp,ms) + nodes(i)%x(in,j2,1)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)

          y_g(mp,ms) = y_g(mp,ms) + nodes(i)%x(in,j2,2)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)
          y_s(mp,ms) = y_s(mp,ms) + nodes(i)%x(in,j2,2)*element_size_ij*H1_s(i,j,ms)*HZ_coord(in,mp)
          y_t(mp,ms) = y_t(mp,ms) + nodes(i)%x(in,j3,2)*element_size_ij*H1(i,j,ms)  *HZ_coord(in,mp)*element_size_perp
          y_p(mp,ms) = y_p(mp,ms) + nodes(i)%x(in,j2,2)*element_size_ij*H1(i,j,ms)  *HZ_coord_p(in,mp)
        end do
      end do
    end do
  end do
end do

do ms=1,n_gauss
  do mp=1,n_plane
    BigR = x_g(mp,ms)
    phi = 2.d0*pi*float(mp-1)/float(n_plane*n_period)
    chi = get_chi_domm(x_g(mp,ms),y_g(mp,ms),phi)
    grad_chi = (/ chi(1,0,0), chi(0,1,0), chi(0,0,1)/BigR /)

    n_perp = (/ -1*x_g(mp,ms)*y_s(mp,ms), x_g(mp,ms)*x_s(mp,ms), x_p(mp,ms)*y_s(mp,ms) - y_p(mp,ms)*x_s(mp,ms) /)
    dA = sqrt(sum(n_perp*n_perp))
    n_perp = n_perp / dA

    do i=1,2
      do j=1,2
        j2 = direction(j)
        element_size_ij = element%size(vertex(i),j2)
        
        do i_tor=1,n_coord_tor
          index_ij = (vertex(i)-1)*(n_order+1)*n_coord_tor + (j2-1)*n_coord_tor + i_tor
          
          v = H1(i,j,ms)*element_size_ij*HZ_coord(i_tor,mp)
          rhs_ij = -1 * v * dot_product(grad_chi, n_perp) * dA
          
          RHS(index_ij) = RHS(index_ij) + rhs_ij*wgauss(ms)
          
        end do
      end do
    end do
  end do
end do

return

end subroutine boundary_matrix_open_chi_correction

end module mod_boundary_matrix_open
