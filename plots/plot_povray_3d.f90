!> Plot the grid of finite elements with the correct curved edges
subroutine plot_povray_3d(node_list, element_list, boundary_list, bnd_node_list, frame, bezier)

use tr_module 
use mod_parameters
use data_structure
use basis_at_gaussian
use phys_module,          only: PI, n_tht
use mod_interp,           only: interp_RZP

implicit none

! --- Routine parameters
type (type_node_list),        intent(in)   :: node_list       !< List of grid nodes
type (type_element_list),     intent(in)   :: element_list    !< List of grid elements
type (type_bnd_element_list), intent(in)   :: boundary_list   !< List of boundary elements
type (type_bnd_node_list),    intent(in)   :: bnd_node_list   !< List of boundary nodes
logical,                      intent(in)   :: frame
logical,                      intent(in)   :: bezier

! --- Local variables
real*8             :: R1,dR1_ds,dR1_dt,dR1_dp,dR1_dst,dR1_ds2,dR1_dt2,dR1_dsp, dR1_dtp, dR1_dpp, Z1,dZ1_ds,dZ1_dt,dZ1_dp,dZ1_dst,dZ1_ds2,dZ1_dt2,dZ1_dsp, dZ1_dtp, dZ1_dpp
real*8             :: R2,dR2_ds,dR2_dt,dR2_dp,dR2_dst,dR2_ds2,dR2_dt2,dR2_dsp, dR2_dtp, dR2_dpp, Z2,dZ2_ds,dZ2_dt,dZ2_dp,dZ2_dst,dZ2_ds2,dZ2_dt2,dZ2_dsp, dZ2_dtp, dZ2_dpp
real*8             :: R3,dR3_ds,dR3_dt,dR3_dp,dR3_dst,dR3_ds2,dR3_dt2,dR3_dsp, dR3_dtp, dR3_dpp, Z3,dZ3_ds,dZ3_dt,dZ3_dp,dZ3_dst,dZ3_ds2,dZ3_dt2,dZ3_dsp, dZ3_dtp, dZ3_dpp
real*8             :: R4,dR4_ds,dR4_dt,dR4_dp,dR4_dst,dR4_ds2,dR4_dt2,dR4_dsp, dR4_dtp, dR4_dpp, Z4,dZ4_ds,dZ4_dt,dZ4_dp,dZ4_dst,dZ4_ds2,dZ4_dt2,dZ4_dsp, dZ4_dtp, dZ4_dpp
logical            :: node_1_found, node_2_found
real*8             :: phi1, phi2

real*8,allocatable :: xp(:,:), x_tot(:,:)
real*8             :: xs(2,n_dim), xx_0(n_dim), uu_0(n_dim), vv_0(n_dim), ww_0(n_dim)
real*8             :: uv_0(n_dim), uv_p(n_dim), xx_p(n_dim) ,xb(4,4,n_dim)
real*8             :: xmax, xmin, ymax, ymin, s, huv_0, huv_p, x_length
integer            :: i, i2, j, inode_0, inode_p, iplot, k, ip, np, iuv, idir_0, idir_p
integer            :: i_tor, i_plane
real*8             :: phi, phi_norm, dt_factor
real*8             :: sgn
character*3        :: label
integer            :: is_planar_elem

write(*,*) '************************************'
write(*,*) '*           plot_grid              *'
write(*,*) '************************************'
write(*,*) ' number of elements          : ',element_list%n_elements
write(*,*) ' number of boundary elements : ',boundary_list%n_bnd_elements
write(*,*) ' number of nodes             : ',node_list%n_nodes
write(*,*) ' bezier                      : ',bezier

call tr_allocate(xp,1,np,1,n_dim,"xp",CAT_GRID)
call tr_allocate(x_tot,1,np,1,n_dim,"x_tot",CAT_GRID)

iplot = 0

open(21, file='povray_input.dat') ! Open the ascii file

! Plot poloidal planes
is_planar_elem = 1
do i_plane=1, n_plane+1
  phi = 2 * PI * (i_plane-1) / float(n_plane) / float(n_period)
  do k=1, element_list%n_elements
   xb(:,:,:) = 0.0
   do i=1,4                               ! over the 4 edges
     inode_0 = element_list%element(k)%vertex(i)
     do iuv=1,4
       huv_0   = element_list%element(k)%size(i,iuv)
       
       do i_tor=1, n_coord_tor              ! toroidal harmonics
         xx_0    = node_list%node(inode_0)%x(i_tor,1,:) * HZ_coord(i_tor, mod(i_plane-1, n_plane)+1)
         uv_0    = node_list%node(inode_0)%x(i_tor,iuv,:) * HZ_coord(i_tor, mod(i_plane-1, n_plane)+1)
  
         if (iuv .eq. 1) then
           xb(i,iuv,:) = xb(i,iuv,:) + xx_0
         else
           xb(i,iuv,:) = xb(i,iuv,:) + xx_0 + uv_0 * huv_0
         endif
       enddo ! i_tor
     enddo ! iuv   
   enddo ! i vertices
  
   ! Write element control points to file
   do i=1,4                               ! over the 4 edges
     do iuv=1,4
      write(21,'(3ES16.8, 2I5)') xb(i,iuv, :), phi, is_planar_elem, i_plane ! Output the grid information to the ascii file
     enddo ! iuv   
   enddo ! i vertices
   write(21,*)
   write(21,*)

  enddo ! k elements
enddo ! i_plane

! Plot toroidal lines joining nodes across planes
dt_factor = 2 * PI / n_tht
is_planar_elem = 0
phi_norm = 2.d0*PI/float(n_plane) / float(n_period) 
do k=1, element_list%n_elements
  node_1_found = .false.
  node_2_found = .false.

  do i=1,4                               ! over the 4 edges
    inode_0 = element_list%element(k)%vertex(i)
    if (node_list%node(inode_0)%boundary .ne. 0) then
      node_1_found = .true.  

      do i2=i,4
        inode_p = element_list%element(k)%vertex(i2)
        if (node_list%node(inode_p)%boundary .ne. 0) then
          node_2_found = .true.
          
          xb(:, :, :) = 0.0

          do i_plane=1, n_plane
            phi = 2 * PI * (i_plane-1) / float(n_plane) / float(n_period)
            phi2 = 2 * PI * (i_plane) / float(n_plane) / float(n_period)
            
            ! Get derivatives from 1st vertex
            call interp_RZP(node_list,element_list,k,1.0,0.0,phi,   &
                            R1,dR1_ds,dR1_dt,dR1_dp,dR1_dst,dR1_ds2,dR1_dt2,dR1_dsp, dR1_dtp, dR1_dpp, &
                            Z1,dZ1_ds,dZ1_dt,dZ1_dp,dZ1_dst,dZ1_ds2,dZ1_dt2,dZ1_dsp, dZ1_dtp, dZ1_dpp)
            
            ! Get derivatives from 2nd vertex
            call interp_RZP(node_list,element_list,k,1.0,0.0,phi2,   &
                            R2,dR2_ds,dR2_dt,dR2_dp,dR2_dst,dR2_ds2,dR2_dt2,dR2_dsp, dR2_dtp, dR2_dpp, &
                            Z2,dZ2_ds,dZ2_dt,dZ2_dp,dZ2_dst,dZ2_ds2,dZ2_dt2,dZ2_dsp, dZ2_dtp, dZ2_dpp)
            
            ! Get derivatives from 1st vertex
            call interp_RZP(node_list,element_list,k,1.0,1.0,phi2,   &
                            R3,dR3_ds,dR3_dt,dR3_dp,dR3_dst,dR3_ds2,dR3_dt2,dR3_dsp, dR3_dtp, dR3_dpp, &
                            Z3,dZ3_ds,dZ3_dt,dZ3_dp,dZ3_dst,dZ3_ds2,dZ3_dt2,dZ3_dsp, dZ3_dtp, dZ3_dpp)
            
            ! Get derivatives from 2nd vertex
            call interp_RZP(node_list,element_list,k,1.0,1.0,phi,   &
                            R4,dR4_ds,dR4_dt,dR4_dp,dR4_dst,dR4_ds2,dR4_dt2,dR4_dsp, dR4_dtp, dR4_dpp, &
                            Z4,dZ4_ds,dZ4_dt,dZ4_dp,dZ4_dst,dZ4_ds2,dZ4_dt2,dZ4_dsp, dZ4_dtp, dZ4_dpp)

            ! Construct element control points - as in Czarny & Huysmans 2008
            xb(1, 1, 1) = R1
            xb(1, 2, 1) = R1 + 1/3.0 * dR1_dp  * phi_norm
            xb(1, 3, 1) = R1 + 1/3.0 * dR1_dt  * dt_factor
            xb(1, 4, 1) = xb(1, 3, 1) + xb(1, 2, 1) - R1 + 1/9.0 * dR1_dtp * phi_norm * dt_factor
            xb(1, 1, 2) = Z1
            xb(1, 2, 2) = Z1 + 1/3.0 * dZ1_dp  * phi_norm
            xb(1, 3, 2) = Z1 + 1/3.0 * dZ1_dt * dt_factor
            xb(1, 4, 2) = xb(1, 3, 2) + xb(1, 2, 2) - Z1 + 1/9.0 * dZ1_dtp * phi_norm * dt_factor
            
            xb(2, 1, 1) = R2
            xb(2, 2, 1) = R2 - 1/3.0 * dR2_dp  * phi_norm
            xb(2, 3, 1) = R2 + 1/3.0 * dR2_dt * dt_factor
            xb(2, 4, 1) = xb(2, 3, 1) + xb(2, 2, 1) - R2 - 1/9.0 * dR2_dtp * phi_norm * dt_factor
            xb(2, 1, 2) = Z2
            xb(2, 2, 2) = Z2 - 1/3.0 * dZ2_dp  * phi_norm
            xb(2, 3, 2) = Z2 + 1/3.0 * dZ2_dt * dt_factor
            xb(2, 4, 2) = xb(2, 3, 2) + xb(2, 2, 2) - Z2 - 1/9.0 * dZ2_dtp * phi_norm * dt_factor
            
            xb(3, 1, 1) = R3
            xb(3, 2, 1) = R3 - 1/3.0 * dR3_dp  * phi_norm
            xb(3, 3, 1) = R3 - 1/3.0 * dR3_dt * dt_factor
            xb(3, 4, 1) = xb(3, 3, 1) + xb(3, 2, 1) - R3 + 1/9.0 * dR3_dtp * phi_norm * dt_factor
            xb(3, 1, 2) = Z3
            xb(3, 2, 2) = Z3 - 1/3.0 * dZ3_dp  * phi_norm
            xb(3, 3, 2) = Z3 - 1/3.0 * dZ3_dt * dt_factor
            xb(3, 4, 2) = xb(3, 3, 2) + xb(3, 2, 2) - Z3 + 1/9.0 * dZ3_dtp * phi_norm * dt_factor
            
            xb(4, 1, 1) = R4
            xb(4, 2, 1) = R4 + 1/3.0 * dR4_dp  * phi_norm
            xb(4, 3, 1) = R4 - 1/3.0 * dR4_dt * dt_factor
            xb(4, 4, 1) = xb(4, 3, 1) + xb(4, 2, 1) - R4 - 1/9.0 * dR4_dtp * phi_norm * dt_factor
            xb(4, 1, 2) = Z4
            xb(4, 2, 2) = Z4 + 1/3.0 * dZ4_dp  * phi_norm
            xb(4, 3, 2) = Z4 - 1/3.0 * dZ4_dt * dt_factor
            xb(4, 4, 2) = xb(4, 3, 2) + xb(4, 2, 2) - Z4 - 1/9.0 * dZ4_dtp * phi_norm * dt_factor
            
            ! Write element to file - choose phi based on control point
            do j=1,4                               
              do iuv=1,4
                if ((j .eq. 1) .or. (j .eq. 4)) then
                  if ( (iuv .eq. 1) .or. (iuv .eq. 3) ) then
                    write(21,'(3ES16.8, 2I5)') xb(j,iuv, :), phi, is_planar_elem, i_plane 
                  else
                    write(21,'(3ES16.8, 2I5)') xb(j,iuv, :), phi + 1/3.0 * (phi2 - phi), is_planar_elem, i_plane
                  endif
                else
                  if ( (iuv .eq. 1) .or. (iuv .eq. 3) ) then
                    write(21,'(3ES16.8, 2I5)') xb(j,iuv, :), phi2, is_planar_elem, i_plane 
                  else
                    write(21,'(3ES16.8, 2I5)') xb(j,iuv, :), phi2 - 1/3.0 * (phi2 - phi), is_planar_elem, i_plane 
                  endif
                endif
              enddo ! iuv   
            enddo ! i vertices
            write(21,*)
            write(21,*)

          enddo ! i_plane
        endif ! node_2_found
      enddo ! i
    endif ! node_1_found

  enddo ! i

  ! Make sure routine behaves well
  if ((node_1_found) .and. (node_2_found)) then
    cycle
  else if ((.not. node_1_found) .and. (.not. node_2_found)) then
    cycle
  else
    write(*, *) "INCONSISTENT VALUES FOR NUMBER OF FOUND NODES: ", node_1_found, node_2_found
    stop
  endif

enddo

close(21) ! Close the ascii file

call tr_deallocate(xp,"xp",CAT_GRID)
call tr_deallocate(x_tot,"x_tot",CAT_GRID)

end subroutine plot_povray_3d
