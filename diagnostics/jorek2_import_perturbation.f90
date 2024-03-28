! -----------------------------------------------------------
! --- Import an n_tor=3 perturbation into another equilibrium
! --- Author :  S.Pamela
! --- Creation date : 14/11/2013
! -----------------------------------------------------------
program jorek2_import_perturbation

  use data_structure
  use phys_module
  use basis_at_gaussian
  use nodes_elements
  use mod_import_restart
  use mod_export_restart
  implicit none

  ! Internal parameters
  integer i, j, k, l, i_tor
  integer my_id, ierr, n_mode, n_copy
  real*8, allocatable :: energies_save(:,:,:)  !< Magnetic and kinetic mode energies at timesteps.

  type(type_node_list) :: node_list2
  type(type_element_list) :: element_list2

  
  ! Name
  write(*,*)'****************************************'
  write(*,*)'****** jorek2_import_perturbation ******'
  write(*,*)'****************************************'
  
  ! Info
  write(*,*)' '
  write(*,*)'------------------------------------------------------'
  write(*,*)'This Program will import an n_tor=3 perturbation'
  write(*,*)'from the .rst|h5 file named jorek_perturbation.rst|h5 and'
  write(*,*)'export that perturbation into jorek_restart.rst|h5, whose'
  write(*,*)'final n_tor number is determined by the n_tor you used'
  write(*,*)'to compile this code (ie. in mod_parameters.f90).'
  write(*,*)'Note that jorek_restart.rst|h5 and jorek_perturbation.rst|h5'
  write(*,*)'must have exactly the same grids.'
  write(*,*)'------------------------------------------------------'
  write(*,*)' '
  
  
  ! --- Import main restart file
  call import_restart(node_list,element_list, 'jorek_restart', rst_format, ierr, .true.)

  ! --- Import rst file that contains the perturbation you want to import
  write(*,*)' '
  call import_restart(node_list2,element_list2, 'jorek_perturbation', rst_format, ierr)


  ! --- Save energies if we are importing into an equilibrium with n_tor>3
  if (n_tor .gt. 3) then
    write(*,*)' '
    write(*,*)'Saving energies of perturbation...'
    call tr_allocate(energies_save,1,3,1,2,1,index_start,"energies_save",CAT_UNKNOWN)
    do i=1,2
      do j=1,index_start
       energies_save(1,i,j) = energies(1,i,j)
       energies_save(2,i,j) = energies(2,i,j)
       energies_save(3,i,j) = energies(3,i,j)
      enddo
    enddo
  endif

  
  ! --- Import (n_tor,n_period) = (3,XX) mode only to another (n_tor,n_period) = (3,XX) equilibrium
  if (n_tor .eq. 3) then
    write(*,*)' '
    write(*,*)'Copying perturbation into node-structure...'
    do i=1,node_list%n_nodes
      do j=1,n_var
    	do k = 1, 4
    	  node_list%node(i)%values(2,k,j) = node_list2%node(i)%values(2,k,j)
    	  node_list%node(i)%deltas(2,k,j) = node_list2%node(i)%deltas(2,k,j)
    	  node_list%node(i)%values(3,k,j) = node_list2%node(i)%values(3,k,j)
    	  node_list%node(i)%deltas(3,k,j) = node_list2%node(i)%deltas(3,k,j)
    	enddo
      enddo
    enddo
  endif

  
  ! --- Import (n_tor,n_period) = (3,XX) mode only to another (n_tor,n_period) = (YY,ZZ) equilibrium
  if (n_tor .gt. 3) then
    
    ! --- Get the mode we are importing (ZZ/2)
    write(*,*)' '
    write(*,'(A,i2,A)') 'What is the n_period of the single mode that you want to import from n_tor=3 to n_tor=',n_tor,' ?'
    read (*,*) n_mode
    n_copy = n_mode/n_period
    if (n_copy .gt. (n_tor-1)/2) then
      write(*,'(A,i2)')        'Warning! You are trying to import the mode n=',n_mode,' into an equilibrium'
      write(*,'(A,i2,A,i2,A)') 'that has (n_tor,n_period) = (',n_tor,',',n_period,').'
      write(*,'(A,i2,A)')      'ie. n_max=',(n_tor-1)/2*n_period,'. Aborting...'
      stop
    endif
    
    ! --- Import variables on nodes
    write(*,*)' '
    write(*,*)'Copying perturbation into node-structure...'
    do i=1,node_list%n_nodes
      do j=1,n_var
  	do k = 1, 4
  	  node_list%node(i)%values(1         ,k,j) = node_list2%node(i)%values(1,k,j)
  	  node_list%node(i)%deltas(1         ,k,j) = node_list2%node(i)%deltas(1,k,j)
  	  node_list%node(i)%values(  n_copy+1,k,j) = node_list2%node(i)%values(2,k,j)
  	  node_list%node(i)%deltas(  n_copy+1,k,j) = node_list2%node(i)%deltas(2,k,j)
  	  node_list%node(i)%values(2*n_copy+1,k,j) = node_list2%node(i)%values(3,k,j)
  	  node_list%node(i)%deltas(2*n_copy+1,k,j) = node_list2%node(i)%deltas(3,k,j)
  	  do l = 2, n_copy
  	    node_list%node(i)%values(l,k,j) = 0.d0
  	    node_list%node(i)%deltas(l,k,j) = 0.d0
  	  enddo
  	  do l = n_copy+2, 2*n_copy
  	    node_list%node(i)%values(l,k,j) = 0.d0
  	    node_list%node(i)%deltas(l,k,j) = 0.d0
  	  enddo
  	enddo
      enddo
    enddo
    
    ! --- Import energies
    write(*,*)' '
    write(*,*)'Copying saved energies...'
    do i=1,2
      do j=1,index_start
       energies(1         ,i,j) = energies_save(1,i,j)
       energies(  n_copy+1,i,j) = energies_save(2,i,j)
       energies(2*n_copy+1,i,j) = energies_save(3,i,j)
       do l = 2, n_copy
    	 energies(l,i,j)  = 0.d0
       enddo
       do l = n_copy+2, 2*n_copy
    	 energies(l,i,j)  = 0.d0
       enddo
      enddo
    enddo
    
    ! --- Deallocate
    call tr_deallocate(energies_save,"energies_save",CAT_UNKNOWN)
  endif

  
  ! Export restart file with perturbation
  write(*,*)' '
  write(*,*)'Exporting restart...'
  call export_restart(node_list, element_list,'jorek_restart')

  write(*,*)'Finished'

end program jorek2_import_perturbation
