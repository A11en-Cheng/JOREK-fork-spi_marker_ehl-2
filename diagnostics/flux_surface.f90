!> This program plots s=const curves on a specified poloidal plane, where
!! s is one of the JOREK element-local coordinates. If the grid is flux-aligned,
!! these s=const curves correspond to flux surfaces in the initial equilibrium.
!!
!! Note that this routine will not work if s=const curves are plotted in the open 
!! field line region
program flux_surface
  use data_structure
  use nodes_elements
  use mod_import_restart
  use basis_at_gaussian
  use phys_module
  use mod_interp
  implicit none
  real*8  :: s, t, phi, rr, zz, R_out, Z_out
  integer :: ierr, n_lines, i_lines, curr, nr, ifail, n_sub, i_pts, j_pts, i_elm

  real*8, dimension(:), allocatable :: R_start, Z_start

  character(len=512) :: skip
  
  call det_modes
  call initialise_basis
  call initialise_parameters(0,  "__NO_FILENAME__")
  write(*,*) "****************************************"
  write(*,*) "*             flux_surface             *"
  write(*,*) "****************************************"
  write(*,*) "WARNING: This program assumes that the grid in the imported restart file is flux-aligned."
  write(*,*) "WARNING: If the grid is not flux-aligned, the output of this program will be useless."
  
  call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)
  
  open(21, file='stpts', status='old', action='read', iostat=ierr)
  if (ierr .eq. 0) then ! stpts file exists, use it.
    read(21, '(a)') skip ! read comment line (ignored)
    read(21,*) n_lines
    if (n_lines .lt. 1) then
      write(*,*) 'ERROR in stpts file: n_lines must be >= 1.'
      stop
    end if
    read(21, '(A)') skip ! read comment line (ignored)
  
    allocate(R_start(n_lines), Z_start(n_lines))
  
    curr = 0
    do
      if (curr .ge. n_lines) exit
    
      read(21, *) nr, rr, zz
    
      if (nr .eq. 1 .and. curr .eq. 0) then
        R_start(1) = rr
        Z_start(1) = zz
      else if (curr .eq. 0) then
        write(*,*) 'ERROR in stpts file: first start point must be nr=1.'
        stop
      else if (nr .lt. 1 .or. nr .gt. n_lines) then
        write(*,*) 'ERROR in stpts file: nr must be > 0 and < n_lines.'
        stop
      else if (nr .le. curr) then
        write(*,*) 'ERROR in stpts file: start points must be sorted in ascending nr-order.'
        stop
      else
        do i_lines=curr+1,nr
          R_start(i_lines) = R_start(curr) + (rr - R_start(curr))*(real(i_lines-curr)/real(nr-curr))
          Z_start(i_lines) = Z_start(curr) + (zz - Z_start(curr))*(real(i_lines-curr)/real(nr-curr))
        end do
      end if
    
      curr = nr
    end do
    close(21)
  else ! if no stpts file exists, use the following hard-coded default startpoints
    n_lines = 50
    allocate(R_start(n_lines), Z_start(n_lines))

    do i_lines=1,n_lines
      R_start(i_lines) = 1.7156 + (2.18-1.7156)*float(i_lines-1)/float(n_lines-1)
      Z_start(i_lines) = 0.12237
    end do
  end if
  
  n_sub = 4
  phi = 2.d0*pi*float(i_plane_rtree - 1)/float(n_period*n_plane)
  
  open(21,file="s_surfaces.dat",action='write',status='replace')
  do i_lines=1,n_lines
    call find_RZ(node_list,element_list,R_start(i_lines),Z_start(i_lines),R_out,Z_out,i_elm,s,t,ifail)
    
    if (ifail .ne. 0) then
      write(*,*) "Can not find RZ,", ifail
      stop
    end if
    
    if (mod(i_elm,n_tht) .eq. 0) then
      i_elm = i_elm - n_tht + 1
    else
      i_elm = int(i_elm/n_tht)*n_tht + 1
    end if
    do i_pts=1,n_tht
      do j_pts=1,n_sub
        call interp_RZP(node_list,element_list,i_elm,s,float(j_pts-1)/float(n_sub),phi,R_out,Z_out)
        write(21,'(2E18.8)') R_out, Z_out
      end do
      i_elm = i_elm + 1
    end do
    
    write(21,*)
    write(21,*)
  end do
  close(21)
end program flux_surface
