subroutine element_matrix_Poisson3D(itype,element,nodes,ivar_in,ivar_out,i_harm,ELM,RHS)
!---------------------------------------------------------------
! calculates the matrix contribution of one element to a 3D 
! Laplace equation 
!---------------------------------------------------------------
use mod_parameters
use data_structure
use gauss
use basis_at_gaussian

implicit none

type (type_element) :: element
type (type_node)    :: nodes(n_vertex_max)

real*8     :: x_g(n_plane,n_gauss,n_gauss), x_s(n_plane,n_gauss,n_gauss), x_t(n_plane,n_gauss,n_gauss), x_p(n_plane,n_gauss,n_gauss)
real*8     :: y_g(n_plane,n_gauss,n_gauss), y_s(n_plane,n_gauss,n_gauss), y_t(n_plane,n_gauss,n_gauss), y_p(n_plane,n_gauss,n_gauss)
real*8     :: ELM(n_vertex_max*(n_order+1)*n_coord_tor,n_vertex_max*(n_order+1)*n_coord_tor), RHS(n_vertex_max*(n_order+1)*n_coord_tor)

real*8     :: xjac, wst
real*8     :: v, v_x, v_y, v_p, chi, chi_x, chi_y, chi_p, rhs_ij
integer    :: ms, mt, mp, i, j, i_tor, k, l, k_tor, index_ij, index_kl, itype, ivar_in, ivar_out, i_harm

ELM=0.d0
RHS=0.d0

!---------------------------------------------------- value of (x,y) and derivatives on Gaussian points
x_g(:,:,:)   = 0.d0; x_s(:,:,:)   = 0.d0; x_t(:,:,:)   = 0.d0; x_p(:,:,:)   = 0.d0
y_g(:,:,:)   = 0.d0; y_s(:,:,:)   = 0.d0; y_t(:,:,:)   = 0.d0; y_p(:,:,:)   = 0.d0

do i=1,n_vertex_max
 do j=1,n_order+1
   do ms=1, n_gauss
     do mt=1, n_gauss
       do mp=1,n_plane
         do i_tor=1,n_coord_tor
           x_g(mp,ms,mt)  = x_g(mp,ms,mt)  + nodes(i)%x(i_tor,j,1) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord(i_tor,mp)
           x_s(mp,ms,mt)  = x_s(mp,ms,mt)  + nodes(i)%x(i_tor,j,1) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord(i_tor,mp)
           x_t(mp,ms,mt)  = x_t(mp,ms,mt)  + nodes(i)%x(i_tor,j,1) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord(i_tor,mp)
           x_p(mp,ms,mt)  = x_p(mp,ms,mt)  + nodes(i)%x(i_tor,j,1) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord_p(i_tor,mp)

           y_g(mp,ms,mt)  = y_g(mp,ms,mt)  + nodes(i)%x(i_tor,j,2) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord(i_tor,mp)
           y_s(mp,ms,mt)  = y_s(mp,ms,mt)  + nodes(i)%x(i_tor,j,2) * element%size(i,j) * H_s(i,j,ms,mt)  * HZ_coord(i_tor,mp)
           y_t(mp,ms,mt)  = y_t(mp,ms,mt)  + nodes(i)%x(i_tor,j,2) * element%size(i,j) * H_t(i,j,ms,mt)  * HZ_coord(i_tor,mp)
           y_p(mp,ms,mt)  = y_p(mp,ms,mt)  + nodes(i)%x(i_tor,j,2) * element%size(i,j) * H(i,j,ms,mt)    * HZ_coord_p(i_tor,mp)
         enddo 
       enddo
     enddo
   enddo
 enddo
enddo


!--------------------------------------------------- sum over the Gaussian integration points
do ms=1, n_gauss

 do mt=1, n_gauss
   wst = wgauss(ms)*wgauss(mt)
   
   do mp=1, n_plane  
     xjac =  x_s(mp,ms,mt)*y_t(mp,ms,mt) - x_t(mp,ms,mt)*y_s(mp,ms,mt)   
     
     do i=1,n_vertex_max
       do j=1,n_order+1
         do i_tor=1,n_coord_tor
     
           index_ij = (i-1)*(n_order+1)*n_coord_tor + (j-1)*n_coord_tor + i_tor
           
           v   = H(i,j,ms,mt)  * element%size(i,j)                                                                   * HZ_coord(i_tor,mp)
           v_x = (  y_t(mp,ms,mt) * H_s(i,j,ms,mt) - y_s(mp,ms,mt) * H_t(i,j,ms,mt) ) * element%size(i,j) / xjac     * HZ_coord(i_tor,mp)
           v_y = (- x_t(mp,ms,mt) * H_s(i,j,ms,mt) + x_s(mp,ms,mt) * H_t(i,j,ms,mt) ) * element%size(i,j) / xjac     * HZ_coord(i_tor,mp)
           v_p = H(i,j,ms,mt)   * element%size(i,j)   * HZ_coord_p(i_tor,mp) - x_p(mp,ms,mt) * v_x - y_p(mp,ms,mt) * v_y                                                               
           
           rhs_ij = 0.0 
           RHS(index_ij) = RHS(index_ij) + v * rhs_ij * x_g(mp,ms,mt) * xjac * wst
           
           do k=1,n_vertex_max
             do l=1,n_order+1
               do k_tor=1,n_coord_tor
           
                 chi   = H(k,l,ms,mt)  * element%size(k,l)                                                                 * HZ_coord(k_tor,mp)
                 chi_x = (   y_t(mp,ms,mt) * H_s(k,l,ms,mt) - y_s(mp,ms,mt) * H_t(k,l,ms,mt) ) * element%size(k,l) / xjac  * HZ_coord(k_tor,mp)
                 chi_y = ( - x_t(mp,ms,mt) * H_s(k,l,ms,mt) + x_s(mp,ms,mt) * H_t(k,l,ms,mt) ) * element%size(k,l) / xjac  * HZ_coord(k_tor,mp)
                 chi_p = H(k,l,ms,mt)  * element%size(k,l)  * HZ_coord_p(k_tor,mp) - chi_x * x_p(mp,ms,mt) - chi_y * y_p(mp,ms,mt)

                 index_kl = (k-1)*(n_order+1)*n_coord_tor + (l-1)*n_coord_tor + k_tor
                 
                 ELM(index_ij,index_kl) =  ELM(index_ij,index_kl) - (chi_x * v_x + chi_y * v_y + 1 / x_g(mp,ms,mt)**2 * chi_p * v_p) * x_g(mp,ms,mt) * xjac * wst

               enddo  ! n_coord_tor
             enddo  ! n_order + 1
           enddo  ! n_vertex_max
         enddo  ! n_coord_tor
       enddo  ! n_order + 1
     enddo  ! n_vertex_max
   enddo  ! n_plane
 enddo  ! n_gauss
enddo ! n_gauss

return
end subroutine element_matrix_poisson3D
