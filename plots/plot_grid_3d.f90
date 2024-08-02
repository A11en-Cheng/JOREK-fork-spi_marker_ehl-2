!> Plot the grid of finite elements with the correct curved edges
subroutine plot_grid_3d(node_list, element_list, boundary_list, bnd_node_list, frame, bezier, &
  gridname)

use tr_module 
use mod_parameters
use data_structure
use basis_at_gaussian
use phys_module,          only: PI

implicit none

! --- Routine parameters
type (type_node_list),        intent(in)   :: node_list       !< List of grid nodes
type (type_element_list),     intent(in)   :: element_list    !< List of grid elements
type (type_bnd_element_list), intent(in)   :: boundary_list   !< List of boundary elements
type (type_bnd_node_list),    intent(in)   :: bnd_node_list   !< List of boundary nodes
logical,                      intent(in)   :: frame
logical,                      intent(in)   :: bezier
character(len=*),             intent(in)   :: gridname        !< Used for the name of the ascii output file

! --- Local variables
real*8,allocatable :: xp(:,:), x_tot(:,:)
real*8             :: xs(2,n_dim), xx_0(n_dim), uu_0(n_dim), vv_0(n_dim), ww_0(n_dim)
real*8             :: uv_0(n_dim), uv_p(n_dim), xx_p(n_dim) ,xb(4,n_dim)
real*8             :: xmax, xmin, ymax, ymin, s, huv_0, huv_p, x_length
integer            :: i, j, inode_0, inode_p, iplot, k, ip, np, iuv, idir_0, idir_p
integer            :: i_tor, i_plane
real*8             :: phi, phi_norm
real*8             :: sgn
character*3        :: label

write(*,*) '************************************'
write(*,*) '*           plot_grid              *'
write(*,*) '************************************'
write(*,*) ' number of elements          : ',element_list%n_elements
write(*,*) ' number of boundary elements : ',boundary_list%n_bnd_elements
write(*,*) ' number of nodes             : ',node_list%n_nodes
write(*,*) ' bezier                      : ',bezier

np = 11

call tr_allocate(xp,1,np,1,n_dim,"xp",CAT_GRID)
call tr_allocate(x_tot,1,np,1,n_dim,"x_tot",CAT_GRID)

iplot = 0

open(21, file='grid_'//trim(adjustl(gridname))//'.dat') ! Open the ascii file

! Plot poloidal planes
do i_plane=1, n_plane+1
  phi = 2 * PI * (i_plane-1) / float(n_plane) / float(n_period)
  do k=1, element_list%n_elements
   do i=1,4                               ! over the 4 edges
     do j=1,np
       xb(:,:) = 0.0
       s = (float(j-1)/float(np-1))
       
       iuv = mod(i+1,2)+1             ! the direction vector corresponding to this edge (i)
       inode_0 = element_list%element(k)%vertex(i)
       huv_0   = element_list%element(k)%size(i,iuv+1)
       
       ip      = mod(i,4)+1
       inode_p = element_list%element(k)%vertex(ip)
       huv_p   = element_list%element(k)%size(ip,iuv+1)
       do i_tor=1, n_coord_tor              ! toroidal harmonics
         xx_0    = node_list%node(inode_0)%x(i_tor,1,:) * HZ_coord(i_tor, mod(i_plane-1, n_plane)+1)
         uv_0    = node_list%node(inode_0)%x(i_tor,iuv+1,:) * HZ_coord(i_tor, mod(i_plane-1, n_plane)+1)
  
         xx_p    = node_list%node(inode_p)%x(i_tor,1,:) * HZ_coord(i_tor, mod(i_plane-1, n_plane)+1)
         uv_p    = node_list%node(inode_p)%x(i_tor,iuv+1,:) * HZ_coord(i_tor, mod(i_plane-1, n_plane)+1)
  
         xb(1,:) = xb(1,:) + xx_0
         xb(2,:) = xb(2,:) + xx_0+uv_0*huv_0
         xb(3,:) = xb(3,:) + xx_p+uv_p*huv_p
         xb(4,:) = xb(4,:) + xx_p
  
       enddo ! i_tor
       call bezier_1d(n_dim, s, xb, x_tot(j,:))
       write(21,'(3ES16.8)') x_tot(j,:), phi ! Output the grid information to the ascii file
     enddo ! j 
     write(21,*)
     write(21,*)
     
   enddo ! i edges
  enddo ! k elements
enddo ! i_plane

! Plot toroidal lines joining nodes across planes
do inode_0=1, node_list%n_nodes
  do i_plane=1, n_plane
    phi_norm = 2.d0*PI/float(n_plane) / float(n_period) 
    
    ! Interpolate edge
    do j=1,np
      s = (float(j-1)/float(np-1))
      phi = (float(i_plane - 1) + s) * phi_norm
      x_tot(j, :) = 0.0
      ! Get node control points
      do i_tor=1,n_coord_tor
        xx_0 = node_list%node(inode_0)%x(i_tor,1,:) * HZ_coord(i_tor, i_plane)
        uv_0 = node_list%node(inode_0)%x(i_tor,1,:) * HZ_coord_p(i_tor, i_plane)
        xx_p = node_list%node(inode_0)%x(i_tor,1,:) * HZ_coord(i_tor, mod(i_plane, n_plane)+1)
        uv_p = node_list%node(inode_p)%x(i_tor,1,:) * HZ_coord_p(i_tor, mod(i_plane, n_plane)+1)
      
        xb(1,:) = xx_0
        xb(2,:) = xx_0 +  uv_0 / 3.0 * phi_norm
        xb(3,:) = xx_p - uv_p / 3.0 * phi_norm
        xb(4,:) = xx_p
      
        call bezier_1d(n_dim, s, xb, xp(j,:))
        x_tot(j, :) = x_tot(j, :) + xp(j, :)
      enddo ! i_tor
      write(21,'(5ES16.8)') x_tot(j,:), phi   ! Output the grid information to the ascii file
    enddo ! end edge
    write(21,*)
    write(21,*)
  enddo ! n_plane
enddo ! n_nodes

close(21) ! Close the ascii file

call tr_deallocate(xp,"xp",CAT_GRID)
call tr_deallocate(x_tot,"x_tot",CAT_GRID)

end subroutine plot_grid_3d
