!> the mod_synchrotron_light implements common
!> procedures and variables to all synchrotron
!> light distribution models
module mod_synchrotron_light_vertices
use mod_light_vertices, only: light_vertices
implicit none

private
public :: synchrotron_light

!> Variables ---------------------------------------
type,abstract,extends(light_vertices) :: synchrotron_light
  contains
  procedure,nopass :: check_x_shaded_in_emission_zone => &
                      check_shaded_x_in_synchrotron_cone
  procedure,nopass :: check_angles_shaded_in_emission_zone => &
                      check_shaded_angles_in_synchrotron_cone
end type synchrotron_light
!> Interfaces --------------------------------------

contains

!> Procedures --------------------------------------

!> check if the shaded point of the synchrotron light is inside
!> the synchrotron radiation cone of half width sin(theta) ≃ 1/gamma
!> where gamma is the relativistic factor. 
!> inputs:
!>   n_x:          (integer) size of the coordinate system
!>   x_shaded:     (real8)(n_x) position of the shaded point
!>   x_light:      (real8)(n_x) position of the point light
!>   n_int_param:  (integer) number of integer parameters: 0
!>   n_real_param: (real8) number of real parameters: 4
!>   int_param:    (integer)(n_int_param) integer parameters
!>   real_param:   (real8)(n_real_param) real_parameters:
!>                 1- x component of the emission cone (momentum) direction
!>                 2- y component of the emission cone (momentum) direction
!>                 3- z component of the emission cone (momentum) direction
!>                 4- relativistic factor
!> outouts:
!>   in_code: (logical) if true the gather point is in the synchrotron cone
function check_shaded_x_in_synchrotron_cone(n_x,x_shaded,x_light,&
n_int_param,n_real_param,int_param,real_param) result(in_range)
  !> Inputs:
  integer,intent(in)                        :: n_x,n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in) :: int_param
  real*8,dimension(n_x),intent(in)          :: x_shaded,x_light
  real*8,dimension(n_real_param),intent(in) :: real_param
  !> Outputs:
  logical :: in_range
  !> Variables:
  real*8 :: costheta
  !> check if the shaded point is in the synchrotron conede
  costheta = dot_product(x_shaded-x_light,real_param(1:3))/norm2(x_shaded-x_light)
  in_range = (costheta.ge.0).and.((sqrt(1d0-costheta**2)*real_param(4)).le.1d0)
end function check_shaded_x_in_synchrotron_cone

!> check if the shaded point of the synchrotron light is inside
!> the synchrotron radiation cone of half width sin(theta) ≃ 1/gamma
!> where gamma is the relativistic factor. The method uses the sinus
!> and cosinus of the angle between the particle and the gather point. 
!> inputs:
!>   n_angles:     (integer) size of the angles must be 2
!>   angles:       (real8)(n_angles) angles between the particle
!>                 the gather point.
!>                 1- sinus of the angle
!>                 2- cosinus of the angle
!>   n_int_param:  (integer) number of integer parameters: 0
!>   n_real_param: (real8) number of real parameters: 1
!>   int_param:    (integer)(n_int_param) integer parameters
!>   real_param:   (real8)(n_real_param) real_parameters:
!>                 1- relativistic factor
!> outouts:
!>   in_code: (logical) if true the gather point is in the synchrotron cone
function check_shaded_angles_in_synchrotron_cone(n_angles,angles,&
n_int_param,n_real_param,int_param,real_param) result(in_range)
  implicit none
  !> inputs:
  integer,intent(in)                        :: n_angles,n_int_param,n_real_param
  integer,dimension(n_int_param),intent(in) :: int_param
  real*8,dimension(n_angles),intent(in)     :: angles
  real*8,dimension(n_real_param),intent(in) :: real_param
  !> outputs:
  logical :: in_range
  !> check if in cone
  in_range = ((angles(1)*real_param(1)).le.1d0).and.(angles(2).ge.0)
end function check_shaded_angles_in_synchrotron_cone

!>-------------------------------------------------
end module mod_synchrotron_light_vertices
