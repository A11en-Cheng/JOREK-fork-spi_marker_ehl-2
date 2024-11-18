!> Determine the derivatives of the numerical input profiles
subroutine derive_num_profiles(my_id)
  
  use phys_module
  use profiles
  
  implicit none
  
  integer, intent(in) :: my_id
  
  if ( my_id == 0 ) then
    
    if ( num_rho ) then
      call derivProf(num_rho_x, num_rho_y0, num_rho_len, num_rho_y1)
      call derivProf(num_rho_x, num_rho_y1, num_rho_len, num_rho_y2)
      call derivProf(num_rho_x, num_rho_y2, num_rho_len, num_rho_y3)
    end if
    
    if ( num_T ) then
      call derivProf(num_T_x, num_T_y0, num_T_len, num_T_y1)
      call derivProf(num_T_x, num_T_y1, num_T_len, num_T_y2)
      call derivProf(num_T_x, num_T_y2, num_T_len, num_T_y3)
    end if

    if ( num_Ti ) then
      call derivProf(num_Ti_x, num_Ti_y0, num_Ti_len, num_Ti_y1)
      call derivProf(num_Ti_x, num_Ti_y1, num_Ti_len, num_Ti_y2)
      call derivProf(num_Ti_x, num_Ti_y2, num_Ti_len, num_Ti_y3)
    end if

    if ( num_Te ) then
      call derivProf(num_Te_x, num_Te_y0, num_Te_len, num_Te_y1)
      call derivProf(num_Te_x, num_Te_y1, num_Te_len, num_Te_y2)
      call derivProf(num_Te_x, num_Te_y2, num_Te_len, num_Te_y3)
    end if

    if ( num_Phi ) then
      call derivProf(num_Phi_x, num_Phi_y0, num_Phi_len, num_Phi_y1)
      call derivProf(num_Phi_x, num_Phi_y1, num_Phi_len, num_Phi_y2)
      call derivProf(num_Phi_x, num_Phi_y2, num_Phi_len, num_Phi_y3)
    end if
    
    if ( num_Fprofile ) then
      call derivProf(num_Fprofile_x, num_Fprofile_y0, num_Fprofile_len, num_Fprofile_y1)
      call derivProf(num_Fprofile_x, num_Fprofile_y1, num_Fprofile_len, num_Fprofile_y2)
      call derivProf(num_Fprofile_x, num_Fprofile_y2, num_Fprofile_len, num_Fprofile_y3)
    end if
    
    if ( num_ffprime ) then
      call derivProf(num_ffprime_x, num_ffprime_y0, num_ffprime_len, num_ffprime_y1)
      call derivProf(num_ffprime_x, num_ffprime_y1, num_ffprime_len, num_ffprime_y2)
    end if
    
    if ( num_rot ) then
      call derivProf(num_rot_x, num_rot_y0, num_rot_len, num_rot_y1)
      call derivProf(num_rot_x, num_rot_y1, num_rot_len, num_rot_y2)
      call derivProf(num_rot_x, num_rot_y2, num_rot_len, num_rot_y3)
    end if    
    
  end if
  
end subroutine derive_num_profiles
