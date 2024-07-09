!> Fill the mode and mode_type arrays.
subroutine det_modes()
  
  use mod_parameters,  only: n_tor, n_period, n_coord_tor, n_coord_period
  use phys_module, only: mode, mode_type, mode_coord, mode_coord_type
  
  implicit none
  
  integer :: itor
  
  ! --- Fill the arrays mode (toroidal mode number n) and mode_type (cos or sin).
  do itor=1, n_tor
    mode(itor)        = int(itor / 2) * n_period
    if ( (itor==1) .or. (mod(itor,2)==0) ) then
      mode_type(itor) = 'cos'
    else
      mode_type(itor) = 'sin'
    end if
  end do
 
 do itor=1,n_coord_tor
    mode_coord(itor)     = int(itor / 2) * n_coord_period
  
    if ( (itor==1) .or. (mod(itor,2)==0) ) then
      mode_coord_type(itor) = 'cos'
    else
      mode_coord_type(itor) = 'sin'
    end if
enddo  

end subroutine det_modes
