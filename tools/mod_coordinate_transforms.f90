!> Module to calculate coordinate transforms between XYZ and RZPhi coordinate systems.
module mod_coordinate_transforms
  implicit none
  private
  public :: cartesian_to_cylindrical
  public :: cylindrical_to_cartesian
  public :: cylindrical_to_cartesian_velocity
  public :: cartesian_to_cylindrical_velocity
  public :: cartesian_to_spherical_latitude
  public :: spherical_latitude_to_cartesian
  public :: spherical_colatitude_to_cartesian
  public :: cartesian_to_spherical_colatitude
  public :: vector_cartesian_to_cylindrical
  public :: vector_cylindrical_to_cartesian
  public :: vector_rotation
  public :: cartesian_velocity_to_cylindrical
  public :: transform_derivatives_st_to_RZ
  public :: transform_derivatives_RZ_to_st
  public :: transform_first_derivatives_st_to_RZ
  public :: transform_second_derivatives_st_to_RZ
  public :: vectors_to_orthonormal_basis
  public :: vectors_spherical_to_cartesian
  public :: vectors_cartesian_to_spherical
  public :: rotate_vectors_cart_z
  public :: mirror_around_cart_x
  public :: mirror_around_cart_y

  !> interfaces
  interface cartesian_to_cylindrical
    module procedure cartesian_to_cylindrical_r4
    module procedure cartesian_to_cylindrical_r8
  end interface cartesian_to_cylindrical

  interface cylindrical_to_cartesian
    module procedure cylindrical_to_cartesian_real8
    module procedure cylindrical_to_cartesian_real4
  end interface cylindrical_to_cartesian

  interface cylindrical_to_cartesian_velocity
    module procedure cylindrical_to_cartesian_velocity_real4
    module procedure cylindrical_to_cartesian_velocity_real8
  end interface cylindrical_to_cartesian_velocity

  interface cartesian_to_cylindrical_velocity
    module procedure cartesian_to_cylindrical_velocity_real4
    module procedure cartesian_to_cylindrical_velocity_real8
  end interface cartesian_to_cylindrical_velocity

  !> overload speherical (latitude) to cartesian coordinate transform
  interface cartesian_to_spherical_latitude
    module procedure cartesian_to_spherical_latitude_real4
    module procedure cartesian_to_spherical_latitude_real8
  end interface cartesian_to_spherical_latitude

  !> overload cartesian to spherical coordinates (latitude) transform
  interface spherical_latitude_to_cartesian
    module procedure spherical_latitude_to_cartesian_real4
    module procedure spherical_latitude_to_cartesian_real8
  end interface spherical_latitude_to_cartesian

  !> overload speherical (colatitude) to cartesian coordinate transform
  interface cartesian_to_spherical_colatitude
    module procedure cartesian_to_spherical_colatitude_std_real4
    module procedure cartesian_to_spherical_colatitude_std_real8
  end interface cartesian_to_spherical_colatitude

  !> overload cartesian to spherical coordinates (colatitude) transform
  interface spherical_colatitude_to_cartesian
    module procedure spherical_colatitude_to_cartesian_std_real4
    module procedure spherical_colatitude_to_cartesian_std_real8
  end interface spherical_colatitude_to_cartesian

  !> overload transform derivatives from local to global coordinates
  interface transform_derivatives_st_to_RZ
     module procedure transform_first_derivatives_st_to_RZ, &
       transform_second_derivatives_st_to_RZ
  end interface transform_derivatives_st_to_RZ

  !> overload transfrom derivatives from global to local coorsinates
  interface transform_derivatives_RZ_to_st
     module procedure transform_first_derivatives_RZ_to_st, &
       transform_second_derivatives_RZ_to_st
  end interface transform_derivatives_RZ_to_st

  !> overload generation of orthonormal basis
  interface vectors_to_orthonormal_basis
    module procedure vectors_to_orthonormal_basis_3d_r4
    module procedure vectors_to_orthonormal_basis_3d_r8
  end interface vectors_to_orthonormal_basis

  !> overload the vector rotation function from a standard 
  !> spherical reference to a standard cartesian reference
  interface vectors_spherical_to_cartesian
    module procedure vectors_spherical_to_cartesian_std_r4
    module procedure vectors_spherical_to_cartesian_std_r8
  end interface vectors_spherical_to_cartesian

  !> overload the vector rotation function from a standard 
  !> cartesian reference to a standard spherical reference
  interface vectors_cartesian_to_spherical
    module procedure vectors_cartesian_to_spherical_std_r4
    module procedure vectors_cartesian_to_spherical_std_r8
  end interface vectors_cartesian_to_spherical

  !> overload method for rotating axis along the standard
  !> cartesian Z axis
  interface rotate_vectors_cart_z
    module procedure rotate_vectors_cart_z_std_r4
    module procedure rotate_vectors_cart_z_std_r8
  end interface rotate_vectors_cart_z

  !> overload method for mirroring a vector w.r.t. 
  !> to the standard x-axis
  interface mirror_around_cart_x
    module procedure mirror_around_cart_x_std_r4
    module procedure mirror_around_cart_x_std_r8
  end interface mirror_around_cart_x

  !> overload method for mirroring a vector w.r.t. 
  !> to the standard y-axis
  interface mirror_around_cart_y
    module procedure mirror_around_cart_y_std_r4
    module procedure mirror_around_cart_y_std_r8
  end interface mirror_around_cart_y

contains
  !> convert a position in xyz coordinates to RZPhi coordinates
  pure function cartesian_to_cylindrical_r4(xyz) result(cyl)
    real*4, intent(in)           :: xyz(3) !< The position in xyz coordinates
    real*4                       :: cyl(3) !< The position in RZPhi coordinates

    cyl(1) = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
    cyl(2) = xyz(3)
    cyl(3) = atan2(-xyz(2), xyz(1))
  end function cartesian_to_cylindrical_r4

  !> convert a position in xyz coordinates to RZPhi coordinates
  pure function cartesian_to_cylindrical_r8(xyz) result(cyl)
    real*8, intent(in)           :: xyz(3) !< The position in xyz coordinates
    real*8                       :: cyl(3) !< The position in RZPhi coordinates

    cyl(1) = sqrt(xyz(1)*xyz(1) + xyz(2)*xyz(2))
    cyl(2) = xyz(3)
    cyl(3) = atan2(-xyz(2), xyz(1))
  end function cartesian_to_cylindrical_r8

  !> converts a position in RZPhi coordinates to xyz coordinates
  pure function cylindrical_to_cartesian_real4(cyl) result(xyz)
    real*4, intent(in)           :: cyl(3) !< The vector components in RZPhi coordinates
    real*4                       :: xyz(3) !< The vector components in xyz coordinates

    xyz(1) = cyl(1)*cos(-cyl(3))
    xyz(2) = cyl(1)*sin(-cyl(3))
    xyz(3) = cyl(2)
  end function cylindrical_to_cartesian_real4

 !> converts a position in RZPhi coordinates to xyz coordinates
  pure function cylindrical_to_cartesian_real8(cyl) result(xyz)
    real*8, intent(in)           :: cyl(3) !< The vector components in RZPhi coordinates
    real*8                       :: xyz(3) !< The vector components in xyz coordinates

    xyz(1) = cyl(1)*cos(-cyl(3))
    xyz(2) = cyl(1)*sin(-cyl(3))
    xyz(3) = cyl(2)
  end function cylindrical_to_cartesian_real8

  !> converts the cylindrical velocity into cartesian velocity
  pure function cylindrical_to_cartesian_velocity_real4(R,phi,v_cyl) result(v_xyz)
    real*4, intent(in) :: R,phi
    real*4, intent(in) :: v_cyl(3)
    real*4             :: v_xyz(3)
    
    v_xyz(1) = v_cyl(1)*cos(-phi) + R*v_cyl(3)*sin(-phi)
    v_xyz(2) = v_cyl(1)*sin(-phi) - R*v_cyl(3)*cos(-phi)
    v_xyz(3) = v_cyl(2)
  end function cylindrical_to_cartesian_velocity_real4

  !> converts the cylindrical velocity into cartesian velocity
  pure function cylindrical_to_cartesian_velocity_real8(R,phi,v_cyl) result(v_xyz)
    real*8, intent(in) :: R,phi
    real*8, intent(in) :: v_cyl(3)
    real*8             :: v_xyz(3)
    
    v_xyz(1) = v_cyl(1)*cos(-phi) + R*v_cyl(3)*sin(-phi)
    v_xyz(2) = v_cyl(1)*sin(-phi) - R*v_cyl(3)*cos(-phi)
    v_xyz(3) = v_cyl(2)
  end function cylindrical_to_cartesian_velocity_real8

  !> converts the cartesian velocity into cylindrical velocity
  pure function cartesian_to_cylindrical_velocity_real4(x,y,v_xyz) result(v_cyl)
    real*4, intent(in) :: x,y
    real*4, intent(in) :: v_xyz(3)
    real*4             :: v_cyl(3)
     
    v_cyl(1) = (x*v_xyz(1) + y*v_xyz(2))/(sqrt(x**2 + y**2))
    v_cyl(3) = (y*v_xyz(1)-x*v_xyz(2))/(x**2 + y**2)
    v_cyl(2) = v_xyz(3)
  end function cartesian_to_cylindrical_velocity_real4

  !> converts the cartesian velocity into cylindrical velocity
  pure function cartesian_to_cylindrical_velocity_real8(x,y,v_xyz) result(v_cyl)
    real*8, intent(in) :: x,y
    real*8, intent(in) :: v_xyz(3)
    real*8             :: v_cyl(3)
     
    v_cyl(1) = (x*v_xyz(1) + y*v_xyz(2))/(sqrt(x**2 + y**2))
    v_cyl(3) = (y*v_xyz(1)-x*v_xyz(2))/(x**2 + y**2)
    v_cyl(2) = v_xyz(3)
  end function cartesian_to_cylindrical_velocity_real8

  !> return the spherical coordinates in terms of the latitude and azimutal angles
  !> (rPsiChi) given the center and the orientation of the sphere. The sphere
  !> orientation is defined by the vectors T,N,B. Single precision is used.
  !> inputs:
  !>   x:      (real4)(3) point position in cartesian coord.
  !>   origin: (real4)(3) position of the sphere origin in cartesian coord.
  !>   T:      (real4)(3) direction defining the azimuthal angle chi
  !>   N:      (real4)(3) normal direction to T
  !>   B:      (real4)(3) direction defining the latitude angle psi
  !> outputs:
  !>   rpsichi: (real4)(3) distance of x from the origin, latitude and azimuthal angles
  pure function cartesian_to_spherical_latitude_real4(x,origin,T,N,B) result(rpsichi)
    implicit none
    real*4,dimension(3),intent(in) :: x,origin,T,N,B
    real*4,dimension(3)            :: rpsichi
    real*4                         :: r_norm
    real*4,dimension(3)            :: r
    
    r = x-origin 
    r_norm = norm2(r)
    rpsichi = (/r_norm,asin(dot_product(r,B)/r_norm),&
    atan2(dot_product(r,N),dot_product(r,T))/)
  end function cartesian_to_spherical_latitude_real4

  !> return the spherical coordinates in terms of the latitude and azimutal angles
  !> (rPsiChi) given the center and the orientation of the sphere. The sphere
  !> orientation is defined by the vectors T,N,B. Double precision is used.
  !> inputs:
  !>   x:      (real8)(3) point position in cartesian coord.
  !>   origin: (real8)(3) position of the sphere origin in cartesian coord.
  !>   T:      (real8)(3) direction defining the azimuthal angle chi
  !>   N:      (real8)(3) normal direction to T
  !>   B:      (real8)(3) direction defining the latitude angle psi
  !> outputs:
  !>   rpsichi: (real8)(3) distance of x from the origin, latitude and azimuthal angles
  pure function cartesian_to_spherical_latitude_real8(x,origin,T,N,B) result(rpsichi)
    implicit none
    real*8,dimension(3),intent(in) :: x,origin,T,N,B
    real*8,dimension(3)            :: rpsichi
    real*8                         :: r_norm
    real*8,dimension(3)            :: r
    
    r = x-origin 
    r_norm = norm2(r)
    rpsichi = (/r_norm,asin(dot_product(r,B)/r_norm),&
    atan2(dot_product(r,N),dot_product(r,T))/)
  end function cartesian_to_spherical_latitude_real8


  !> return the spherical coordinates in terms of the colatitude and azimutal angles
  !> (rThetaChi) given the center and the orientation of the sphere,
  !> inputs:
  !>   x:      (real4)(3) point position in cartesian coord (X,Y,Z).
  !>   origin: (real4)(3) position of the sphere origin in cartesian coord.
  !> outputs:
  !>   rthetachi: (real4)(3) distance of x from the origin, latitude and azimuthal angles
  pure function cartesian_to_spherical_colatitude_std_real4(x,origin) result(rthetachi)
    implicit none
    real*4,dimension(3),intent(in) :: x,origin
    real*4,dimension(3)            :: rthetachi
    real*4                         :: r_norm
    real*4,dimension(3)            :: r
    
    r = x-origin 
    r_norm = norm2(r)
    rthetachi = (/r_norm,acos(r(3)/r_norm),atan2(r(2),r(1))/)
  end function cartesian_to_spherical_colatitude_std_real4

  !> return the spherical coordinates in terms of the colatitude and azimutal angles
  !> (rThetaChi) given the center and the orientation of the sphere,
  !> inputs:
  !>   x:      (real8)(3) point position in cartesian coord (X,Y,Z).
  !>   origin: (real8)(3) position of the sphere origin in cartesian coord.
  !> outputs:
  !>   rthetachi: (real8)(3) distance of x from the origin, latitude and azimuthal angles
  pure function cartesian_to_spherical_colatitude_std_real8(x,origin) result(rthetachi)
    implicit none
    real*8,dimension(3),intent(in) :: x,origin
    real*8,dimension(3)            :: rthetachi
    real*8                         :: r_norm
    real*8,dimension(3)            :: r
    
    r = x-origin 
    r_norm = norm2(r)
    rthetachi = (/r_norm,acos(r(3)/r_norm),atan2(r(2),r(1))/)
  end function cartesian_to_spherical_colatitude_std_real8

  !> Transform the position in the spherical coordinate (latitude) system of the 
  !> sphere with origin at 'origin' into global cartesian coordinates. 
  !> Single precision is used.
  !> inputs:
  !>   rpsichi: (real4)(3) distance of x from the origin, latitude and azimuthal angles
  !>   origin:  (real4)(3) position of the sphere origin in cartesian coord.
  !>   T:       (real4)(3) direction defining the azimuthal angle chi
  !>   N:       (real4)(3) normal direction to T 
  !>   B:      (real8)(3) direction defining the latitude angle psi
  !> outputs:
  !>   x: (real4)(3) point position in cartesian coordinates
  pure function spherical_latitude_to_cartesian_real4(rpsichi,origin,T,N,B) result(x)
    implicit none
    real*4,dimension(3),intent(in) :: rpsichi,origin,T,N,B
    real*4,dimension(3) :: x

    x = origin + rpsichi(1)*(cos(rpsichi(2))*&
    (T*cos(rpsichi(3))+N*sin(rpsichi(3)))+B*sin(rpsichi(2)))
  end function spherical_latitude_to_cartesian_real4

  !> Transform the position in the spherical coordinate (latitude) system of the 
  !> sphere with origin at 'origin' into global cartesian coordinates. 
  !> Double precision is used.
  !> inputs:
  !>   rpsichi: (real8)(3) distance of x from the origin, latitude and azimuthal angles
  !>   origin:  (real8)(3) position of the sphere origin in cartesian coord.
  !>   T:       (real8)(3) direction defining the azimuthal angle chi
  !>   N:       (real8)(3) normal direction to T 
  !>   B:      (real8)(3) direction defining the latitude angle psi
  !> outputs:
  !>   x: (real8)(3) point position in cartesian coordinates
  pure function spherical_latitude_to_cartesian_real8(rpsichi,origin,T,N,B) result(x)
    implicit none
    real*8,dimension(3),intent(in) :: rpsichi,origin,T,N,B
    real*8,dimension(3) :: x

    x = origin + rpsichi(1)*(cos(rpsichi(2))*&
    (T*cos(rpsichi(3))+N*sin(rpsichi(3)))+B*sin(rpsichi(2)))
  end function spherical_latitude_to_cartesian_real8

  !> Transform the position in the spherical coordinate (colatitude) system of the 
  !> sphere with origin at 'origin' into global cartesian coordinates. 
  !> Single precision is used.
  !> inputs:
  !>   rthetachi: (real4)(3) distance of x from the origin, colatitude and azimuthal angles
  !>   origin:    (real4)(3) position of the sphere origin in cartesian coord.
  !> outputs:
  !>   x: (real4)(3) point position in cartesian coordinates (X,Y,Z)
  pure function spherical_colatitude_to_cartesian_std_real4(rthetachi,origin) result(x)
    implicit none
    real*4,dimension(3),intent(in) :: rthetachi,origin
    real*4,dimension(3) :: x

    x = origin + rthetachi(1)*(/sin(rthetachi(2))*cos(rthetachi(3)),&
    sin(rthetachi(2))*sin(rthetachi(3)),cos(rthetachi(2))/)
  end function spherical_colatitude_to_cartesian_std_real4

  !> Transform the position in the spherical coordinate (colatitude) system of the 
  !> sphere with origin at 'origin' into global cartesian coordinates. 
  !> Single precision is used.
  !> inputs:
  !>   rthetachi: (real8)(3) distance of x from the origin, colatitude and azimuthal angles
  !>   origin:    (real8)(3) position of the sphere origin in cartesian coord.
  !> outputs:
  !>   x: (real4)(3) point position in cartesian coordinates (X,Y,Z)
  pure function spherical_colatitude_to_cartesian_std_real8(rthetachi,origin) result(x)
    implicit none
    real*8,dimension(3),intent(in) :: rthetachi,origin
    real*8,dimension(3) :: x

    x = origin + rthetachi(1)*(/sin(rthetachi(2))*cos(rthetachi(3)),&
    sin(rthetachi(2))*sin(rthetachi(3)),cos(rthetachi(2))/)
  end function spherical_colatitude_to_cartesian_std_real8

  pure function vector_cartesian_to_cylindrical(phi,a) result(b)
    real*8, intent(in)               :: phi !< The local toroidal angle
    real*8, dimension(3), intent(in) :: a   !< The vector components in (ex,ey,ez) basis
    real*8, dimension(3)             :: b   !< The vector components in (eR,eZ,ephi) basis
    real*8, dimension(2)             :: sincosphi

    sincosphi = (/sin(phi),cos(phi)/)
   
    b(1) = a(1)*sincosphi(2) - a(2)*sincosphi(1) 
    b(2) = a(3)
    b(3) = -1.d0*(a(1)*sincosphi(1) + a(2)*sincosphi(2))
  end function vector_cartesian_to_cylindrical  

  !> convert a vector in (eR,eZ,ephi) basis into (ex,ey,ez)  basis
  pure function vector_cylindrical_to_cartesian(phi,a) result(b)
    real*8, intent(in)               :: phi !< The local toroidal angle
    real*8, dimension(3), intent(in) :: a   !< The vector components in (eR,eZ,ephi) basis
    real*8, dimension(3)             :: b   !< The vector components in (ex,ey,ez) basis
    real*8, dimension(2)             :: sincosphi

    sincosphi = (/sin(phi),cos(phi)/)

    b(1) = a(1)*sincosphi(2) - a(3)*sincosphi(1) 
    b(2) = -1.d0*(a(1)*sincosphi(1) + a(3)*sincosphi(2))
    b(3) = a(2)
  end function vector_cylindrical_to_cartesian

  !> multiply a vector with the rotation matrix for angle phi around the z-axis in RZPhi coordinates
  pure function vector_rotation(in, phi) result(out)
    real*8, intent(in) :: in(3) !< Input vector in RZPhi coordinates
    real*8, intent(in) :: phi
    real*8             :: out(3) !< Output vector in RZPhi coordinates

    out(1) = cos(phi) * in(1) - sin(phi) * in(3)
    out(2) = in(2)
    out(3) = sin(phi) * in(1) + cos(phi) * in(3)
  end function vector_rotation

 !> Calculate the corresponding cylindrical expression for a cartesian velocity vector
  pure function cartesian_velocity_to_cylindrical(in, phi) result(out)
    real*8, intent(in) :: in(3) !< Input velocity in xyz coordinates
    real*8, intent(in) :: phi !< toroidal angle of particle (jorek coordinate 3)
    real*8             :: out(3) !< Output velocity in RZPhi coordinates

    out(1) = cos(phi) * in(1) + sin(phi) * in(2)
    out(2) = in(3)
    out(3) = -sin(phi) * in(1) + cos(phi) * in(2)
  end function cartesian_velocity_to_cylindrical

  !> Compute an orthonormal basis give two non-parallel vectors in 3D cartesian
  !> space and for real4 functions.
  !> inputs:
  !>   v1: (real4)(3) first vector
  !>   v2: (real4)(3) second vector non-parallel to v1
  !> outputs:
  !>   T: (real4)(3) first basis v1/||v1||
  !>   N: (real4)(3) second basis v2-T*(v2*T)
  !>   B: (real4)(3) third basis  T x N
  subroutine vectors_to_orthonormal_basis_3d_r4(v1,v2,T,N,B)
    use mod_math_operators, only: cross_product
    implicit none
    real*4,dimension(3),intent(in) :: v1,v2
    real*4,dimension(3),intent(out) :: T,N,B
    !> compute orthogonal basis
    T = v1/norm2(v1)
    N = v2 - T*(v2(1)*T(1)+v2(2)*T(2)+v2(3)*T(3))
    N = N/norm2(N)
    B = cross_product(T,N)
    B = B/norm2(B)
  end subroutine vectors_to_orthonormal_basis_3d_r4

  !> Compute an orthonormal basis give two non-parallel vectors in 3D cartesian
  !> space and for real8 functions.
  !> inputs:
  !>   v1: (real8)(3) first vector
  !>   v2: (real8)(3) second vector non-parallel to v1
  !> outputs:
  !>   T: (real8)(3) first basis v1/||v1||
  !>   N: (real8)(3) second basis v2-T*(v2*T)
  !>   B: (real8)(3) third basis  T x N
  subroutine vectors_to_orthonormal_basis_3d_r8(v1,v2,T,N,B)
    use mod_math_operators, only: cross_product
    implicit none
    real*8,dimension(3),intent(in) :: v1,v2
    real*8,dimension(3),intent(out) :: T,N,B
    !> compute orthogonal basis
    T = v1/norm2(v1)
    N = v2 - T*(v2(1)*T(1)+v2(2)*T(2)+v2(3)*T(3))
    N = N/norm2(N)
    B = cross_product(T,N)
    B = B/norm2(B)
  end subroutine vectors_to_orthonormal_basis_3d_r8

  !> rotate spherical vectors defined in standard spherical reference
  !> to the standard cartesian reference. Apply a first rotation along
  !> the Z axis of angle phi and a second rotation along the rotated
  !> y axis of an angle theta then multiply for the sphere radius
  !> inputs:
  !> n_v:       (integer) number of vectors to rotate
  !> rthetaphi: (real4)(3) sphere radius, colatitude (theta) and azimuth (phi)
  !> vect:      (real4)(3,n_v) vectors in the spherical system 
  !> outpus:
  !> vect:      (real4)(3,n_v) vectors in the cartesian reference
  subroutine vectors_spherical_to_cartesian_std_r4(n_v,rthetaphi,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*4,dimension(3),intent(in)        :: rthetaphi
    real*4,dimension(3,n_v),intent(inout) :: vect
    real*4,dimension(2)   :: cos_thetaphi,sin_thetaphi
    real*4,dimension(3,3) :: rot
    !> compute rotation transform
    cos_thetaphi = cos(rthetaphi(2:3));
    sin_thetaphi = sin(rthetaphi(2:3));
    rot(:,1) = (/cos_thetaphi(1)*cos_thetaphi(2),sin_thetaphi(2)*cos_thetaphi(1),-sin_thetaphi(1)/)
    rot(:,2) = (/-sin_thetaphi(2),cos_thetaphi(2),real(0d0,kind=4)/)
    rot(:,3) = (/sin_thetaphi(1)*cos_thetaphi(2),sin_thetaphi(1)*sin_thetaphi(2),cos_thetaphi(1)/)
    !> transform the vertices in the new positions
    vect = rthetaphi(1)*matmul(rot,vect)
  end subroutine vectors_spherical_to_cartesian_std_r4

  !> rotate spherical vectors defined in standard spherical reference
  !> to the standard cartesian reference. Apply a first rotation along
  !> the Z axis of angle phi and a second rotation along the rotated
  !> y axis of an angle theta then multiply for the sphere radius
  !> inputs:
  !> n_v:       (integer) number of vectors to rotate
  !> rthetaphi: (real8)(3) sphere radius, colatitude (theta) and azimuth (phi)
  !> vect:      (real8)(3,n_v) vectors in the spherical system 
  !> outpus:
  !> vect:      (real8)(3,n_v) vectors in the cartesian reference
  subroutine vectors_spherical_to_cartesian_std_r8(n_v,rthetaphi,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*8,dimension(3),intent(in)        :: rthetaphi
    real*8,dimension(3,n_v),intent(inout) :: vect
    real*8,dimension(2)   :: cos_thetaphi,sin_thetaphi
    real*8,dimension(3,3) :: rot
    !> compute rotation transform
    cos_thetaphi = cos(rthetaphi(2:3));
    sin_thetaphi = sin(rthetaphi(2:3));
    rot(:,1) = (/cos_thetaphi(1)*cos_thetaphi(2),sin_thetaphi(2)*cos_thetaphi(1),-sin_thetaphi(1)/)
    rot(:,2) = (/-sin_thetaphi(2),cos_thetaphi(2),0d0/)
    rot(:,3) = (/sin_thetaphi(1)*cos_thetaphi(2),sin_thetaphi(1)*sin_thetaphi(2),cos_thetaphi(1)/)
    !> transform the vertices in the new positions
    vect = rthetaphi(1)*matmul(rot,vect)
  end subroutine vectors_spherical_to_cartesian_std_r8

  !> rotate cartesian vectors defined in standard cartesian reference
  !> to the standard spherical reference. Apply a first rotation along
  !> the y axis of angle phi and a second rotation along the rotated
  !> z axis of an angle theta then divide by the sphere radius
  !> inputs:
  !> n_v:       (integer) number of vectors to rotate
  !> rthetaphi: (real8)(3) sphere radius, colatitude (theta) and azimuth (phi)
  !> vect:      (real8)(3,n_v) vectors in the cartesian system 
  !> outpus:
  !> vect:      (real8)(3,n_v) vectors in the spherical reference
  subroutine vectors_cartesian_to_spherical_std_r4(n_v,rthetaphi,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*4,dimension(3),intent(in)        :: rthetaphi
    real*4,dimension(3,n_v),intent(inout) :: vect
    real*4,dimension(2)   :: cos_thetaphi,sin_thetaphi
    real*4,dimension(3,3) :: rot
    !> compute rotation transform
    cos_thetaphi = cos(rthetaphi(2:3));
    sin_thetaphi = sin(rthetaphi(2:3));
    rot(:,1) = (/cos_thetaphi(1)*cos_thetaphi(2),-sin_thetaphi(2),sin_thetaphi(1)*cos_thetaphi(2)/)
    rot(:,2) = (/sin_thetaphi(2)*cos_thetaphi(1),cos_thetaphi(2),sin_thetaphi(1)*sin_thetaphi(2)/)
    rot(:,3) = (/-sin_thetaphi(1),real(0d0,kind=4),cos_thetaphi(1)/)
    !> transform the vertices in the new positions
    vect = matmul(rot,vect/rthetaphi(1))
  end subroutine vectors_cartesian_to_spherical_std_r4

  !> rotate cartesian vectors defined in standard cartesian reference
  !> to the standard spherical reference. Apply a first rotation along
  !> the y axis of angle phi and a second rotation along the rotated
  !> z axis of an angle theta then divide by the sphere radius
  !> inputs:
  !> n_v:       (integer) number of vectors to rotate
  !> rthetaphi: (real8)(3) sphere radius, colatitude (theta) and azimuth (phi)
  !> vect:      (real8)(3,n_v) vectors in the cartesian system 
  !> outpus:
  !> vect:      (real8)(3,n_v) vectors in the spherical reference
  subroutine vectors_cartesian_to_spherical_std_r8(n_v,rthetaphi,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*8,dimension(3),intent(in)        :: rthetaphi
    real*8,dimension(3,n_v),intent(inout) :: vect
    real*8,dimension(2)   :: cos_thetaphi,sin_thetaphi
    real*8,dimension(3,3) :: rot
    !> compute rotation transform
    cos_thetaphi = cos(rthetaphi(2:3));
    sin_thetaphi = sin(rthetaphi(2:3));
    rot(:,1) = (/cos_thetaphi(1)*cos_thetaphi(2),-sin_thetaphi(2),sin_thetaphi(1)*cos_thetaphi(2)/)
    rot(:,2) = (/sin_thetaphi(2)*cos_thetaphi(1),cos_thetaphi(2),sin_thetaphi(1)*sin_thetaphi(2)/)
    rot(:,3) = (/-sin_thetaphi(1),0d0,cos_thetaphi(1)/)
    !> transform the vertices in the new positions
    vect = matmul(rot,vect/rthetaphi(1))
  end subroutine vectors_cartesian_to_spherical_std_r8

  !> rotate vector around the standard cartesia z axis
  !> inputs: 
  !>   n_v:    (integer) number of vectors
  !>   angles: (real4) rotation angle
  !>   vect:   (3,n_v)(real4) vectors to be rotated
  !> outputs:
  !>   vect:   (3,n_v)(real4) rotated vectors
  subroutine rotate_vectors_cart_z_std_r4(n_v,angle,vect)
    implicit none
    integer,intent(in)                     :: n_v
    real*4,intent(in)                      :: angle
    real*4,dimension(3,n_v),intent(inout)  :: vect
    real*4,dimension(2,2)                  :: rot
    rot(1:2,1) = [cos(angle),-sin(angle)]
    rot(1:2,2) = [-rot(2,1),rot(1,1)]
    !> change the orientation of the plane points
    vect(1:2,:) = matmul(rot,vect(1:2,:))
  end subroutine rotate_vectors_cart_z_std_r4

  !> rotate vector around the standard cartesia z axis
  !> inputs: 
  !>   n_v:    (integer) number of vectors
  !>   angles: (real8) rotation angle
  !>   vect:   (3,n_v)(real8) vectors to be rotated
  !> outputs:
  !>   vect:   (3,n_v)(real8) rotated vectors
  subroutine rotate_vectors_cart_z_std_r8(n_v,angle,vect)
    implicit none
    integer,intent(in)                     :: n_v
    real*8,intent(in)                      :: angle
    real*8,dimension(3,n_v),intent(inout)  :: vect
    real*8,dimension(2,2)                  :: rot
    rot(1:2,1) = [cos(angle),-sin(angle)]
    rot(1:2,2) = [-rot(2,1),rot(1,1)]
    !> change the orientation of the plane points
    vect(1:2,:) = matmul(rot,vect(1:2,:))
  end subroutine rotate_vectors_cart_z_std_r8

  !> mirror the vectors w.r.t the standard cartesian x-axis
  !> inputs:
  !>   n_v:  (integer) number of vectors
  !>   vect: (3,n_v)(real4) vectors to mirror
  !> outputs:
  !>   vect: (3,n_v)(real4) mirrored vectors
  subroutine mirror_around_cart_x_std_r4(n_v,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*4,dimension(3,n_v),intent(inout) :: vect
    vect(2,:) = -vect(2,:)
  end subroutine mirror_around_cart_x_std_r4

  !> mirror the vectors w.r.t the standard cartesian x-axis
  !> inputs:
  !>   n_v:  (integer) number of vectors
  !>   vect: (3,n_v)(real8) vectors to mirror
  !> outputs:
  !>   vect: (3,n_v)(real8) mirrored vectors
  subroutine mirror_around_cart_x_std_r8(n_v,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*8,dimension(3,n_v),intent(inout) :: vect
    vect(2,:) = -vect(2,:)
  end subroutine mirror_around_cart_x_std_r8

  !> mirror the vectors w.r.t the standard cartesian y-axis
  !> inputs:
  !>   n_v:  (integer) number of vectors
  !>   vect: (3,n_v)(real8) vectors to mirror
  !> outputs:
  !>   vect: (3,n_v)(real8) mirrored vectors
  subroutine mirror_around_cart_y_std_r4(n_v,vect)
    implicit none
    integer,intent(in)                    :: n_v
    real*4,dimension(3,n_v),intent(inout) :: vect
    vect(1,:) = -vect(1,:)
  end subroutine mirror_around_cart_y_std_r4

  !> mirror the vectors w.r.t the standard cartesian y-axis
  !> inputs:
  !>   n_v:  (integer) number of vectors
  !>   vect: (3,n_v)(real8) vectors to mirror
  !> outputs:
  !>   vect: (3,n_v)(real8) mirrored vectors
  subroutine mirror_around_cart_y_std_r8(n_v,vect)
    implicit none
    integer,intent(in)   :: n_v
    real*8,dimension(3,n_v),intent(inout) :: vect
    vect(1,:) = -vect(1,:)
  end subroutine mirror_around_cart_y_std_r8

!--------------------------------------------------------------------------
!> This procedure expresses first order derivatives from the local (s,t)
!> to the global (R,Z) reference system.
!> inputs:
!>   n_v: (integer) number of physical quantities
!>   P_s: (real8)(n_v) physical quantities first derivative in s
!>   P_t: (real8)(n_v) physical quantities first derivatives in t
!>   R_s: (real8) major radius first derivative in s
!>   R_z: (real8) major radius first derivative in t
!>   Z_s: (real8) vertical position first derivative in s
!>   Z_t: (real8) vertical position first derivative in t
!> outputs:
!>   P_R: (real8)(n_v) physical quantities first derivative in R
!>   P_Z: (real8)(n_v) physcial quantities first derivative in Z
pure subroutine transform_first_derivatives_st_to_RZ(P_R,P_Z,n_v, &
  P_s,P_t,R_s,R_t,Z_s,Z_t)
  implicit none
  !> input variables
  integer, intent(in)                      :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_s, P_t
  real(kind=8), intent(in)                 :: R_s, R_t, Z_s, Z_t
  !> output variables
  real(kind=8), dimension(n_v), intent(out) :: P_R, P_Z
  !> internal variables
  real(kind=8) :: inverse_jacobian

  !> compute the inverse of the jacobian
  inverse_jacobian = 1.d0/(R_s*Z_t-R_t*Z_s)

  !> compute derivatives wrt. R and Z
  P_R = (P_s*Z_t-P_t*Z_s)*inverse_jacobian
  P_Z = (P_t*R_s-P_s*R_t)*inverse_jacobian
  
end subroutine transform_first_derivatives_st_to_RZ

!---------------------------------------------------------------------------

!> This procedure expresses second order derivatives from the local (s,t)
!> to the global (R,Z) reference system. Functions are assumed to be Hessian.
!> inputs:
!>   n_v:  (integer) number of physical quantities
!>   P_ss: (real8)(n_v) physical quantities second derivatives in s
!>   P_st: (real8)(n_v) physical quantities cross derivatives in s,t
!>   P_tt: (real8)(n_v) physical quantities second derivative in t
!>   P_R:  (real8)(n_v) physical quantities first derivative in R
!>   P_Z:  (real8)(n_v) physical quantities first derivative in Z
!>   R_s:  (real8) major radius first derivative in s
!>   R_t:  (real8) major radius first derivative in t
!>   R_ss: (real8) major radius second derivative in s
!>   R_st: (real8) major radius cross derivative in s,t
!>   R_tt: (real8) major radius second derivative in t
!>   Z_s:  (real8) vertical position first derivative in s
!>   Z_t:  (real8) vertical position first derivative in t
!>   Z_ss: (real8) vertical position second derivative in s
!>   Z_st: (real8) vertical position cross derivative in s,t
!>   Z_tt: (real8) vertical position second derivative in t
!> outputs:
!>   P_RR: (real8)(n_v) physical quantities second derivatives in R
!>   P_RZ: (real8)(n_v) physical quantities cross derivatives in R,Z
!>   P_ZZ: (real8)(n_v) physical quantities second derivatives in Z
pure subroutine transform_second_derivatives_st_to_RZ(P_RR,P_RZ,P_ZZ,n_v, &
  P_ss,P_st,P_tt,P_R,P_Z,R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> input variables
  integer, intent(in) :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_ss, P_st, P_tt, P_R, P_Z
  real(kind=8), intent(in) :: R_s, R_t, R_ss, R_st, R_tt
  real(kind=8), intent(in) :: Z_s, Z_t, Z_ss, Z_st, Z_tt
  !> output variables
  real(kind=8), dimension(n_v), intent(out) :: P_RR, P_RZ, P_ZZ
  !> internal variables
  integer                      :: i
  real(kind=8), dimension(10)  :: transformation_matrix !< 10:jacobian
  real(kind=8), dimension(n_v) :: RHS_RR, RHS_RZ, RHS_ZZ

  !< compute matrix elements
  transformation_matrix(1:9) = [R_s*R_s,2.d0*R_s*Z_s,Z_s*Z_s, &
     R_t*R_s,R_s*Z_t+Z_s*R_t,Z_s*Z_t,R_t*R_t,                 &
     2.d0*R_t*Z_t,Z_t*Z_t]

  !< compute the inverse of the matrix jacobian
  transformation_matrix(10) = 1.d0/(transformation_matrix(1)*( &
    transformation_matrix(5)*transformation_matrix(9)-         &
    transformation_matrix(8)*transformation_matrix(6))+        &
    transformation_matrix(2)*(transformation_matrix(7)*        &
    transformation_matrix(6)-transformation_matrix(4)*         &
    transformation_matrix(9))+transformation_matrix(3)*(       &
    transformation_matrix(4)*transformation_matrix(8)-         &
    transformation_matrix(5)*transformation_matrix(7)))
    
  !> compute RHS
  RHS_RR = P_ss - R_ss*P_R - Z_ss*P_Z
  RHS_RZ = P_st - R_st*P_R - Z_st*P_Z
  RHS_ZZ = P_tt - R_tt*P_R - Z_tt*P_Z
  
  !> compute second order derivatives
  P_RR = ((transformation_matrix(5)*transformation_matrix(9)-   &
    transformation_matrix(8)*transformation_matrix(6))*RHS_RR + &
    (transformation_matrix(8)*transformation_matrix(3)-         &
    transformation_matrix(2)*transformation_matrix(9))*RHS_RZ + &
    (transformation_matrix(2)*transformation_matrix(6)-         &
    transformation_matrix(5)*transformation_matrix(3))*         &
    RHS_ZZ)*transformation_matrix(10)
  P_RZ = ((transformation_matrix(7)*transformation_matrix(6)-   &
    transformation_matrix(4)*transformation_matrix(9))*RHS_RR + &
    (transformation_matrix(1)*transformation_matrix(9)-         &
    transformation_matrix(7)*transformation_matrix(3))*RHS_RZ + &
    (transformation_matrix(4)*transformation_matrix(3)-         &
    transformation_matrix(1)*transformation_matrix(6))*         &
    RHS_ZZ)*transformation_matrix(10)
  P_ZZ = ((transformation_matrix(4)*transformation_matrix(8)-   &
    transformation_matrix(5)*transformation_matrix(7))*RHS_RR+  &
    (transformation_matrix(2)*transformation_matrix(7)-         &
    transformation_matrix(1)*transformation_matrix(8))*RHS_RZ+  &
    (transformation_matrix(1)*transformation_matrix(5)-         &
    transformation_matrix(2)*transformation_matrix(4))*         &
    RHS_ZZ)*transformation_matrix(10) 
  
end subroutine transform_second_derivatives_st_to_RZ

!---------------------------------------------------------------------------------------

!> This procedure transforms back the first order derivatives from global RZ coordinates
!> to local s,t coordinates. This method is used mainly for tests.
!> inputs:
!>   n_v:  (integer) number of physical quantities
!>   P_R,P_Z: (real8)(n_v) physical quantities first derivatives in R and Z
!>   R_s,R_t: (real8) major radius first derivatives in s and t
!>   Z_s,Z_t: (real) vertical position first derivatives in s and t
!> outputs:
!>   P_s,P_t: (real8)(n_v) physical quantities first derivatives in s and t
pure subroutine transform_first_derivatives_RZ_to_st(P_s,P_t,n_v,P_R,P_Z, &
  R_s,R_t,Z_s,Z_t)
  implicit none
  !> input variables
  integer, intent(in)                      :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_R, P_Z
  real(kind=8), intent(in)                 :: R_s, R_t, Z_s, Z_t
  !> output varibales
  real(kind=8), dimension(n_v), intent(out) :: P_s, P_t

  P_s = P_R*R_s + P_Z*Z_s
  P_t = P_R*R_t + P_Z*Z_t
  
end subroutine transform_first_derivatives_RZ_to_st
!---------------------------------------------------------------------------------------

!> This procedure transforms back the second order derivatives from global RZ coordinates
!> to local s,t coordinates. This method is used mainly for tests.
!> inputs:
!>   n_v:            (integer) number of physical quantities
!>   P_R,P_Z:        (real8)(n_v) physical quantities first derivatives in R and Z
!>   P_RR,P_RZ,P_ZZ: (real8)(n_v) physical quantities second and cross derivatives in R and Z
!>   R_s,R_t:        (real8) major radius first derivatives in s and t
!>   Z_s,Z_t:        (real8) vertical position first derivatives in s and t
!>   R_ss,R_st,R_tt: (real8) major radius second and cross derivatives in s and t
!>   Z_ss,Z_st,Z_tt: (real8) vertical position second and cross derivatives in s and t
!> outputs:
!>   P_ss,P_st,P_tt: (real)(n_v) physical quantities second and cross derivatives in s and t
pure subroutine transform_second_derivatives_RZ_to_st(P_ss,P_st,P_tt,n_v,P_R,P_Z, &
  P_RR,P_RZ,P_ZZ,R_s,R_t,R_ss,R_st,R_tt,Z_s,Z_t,Z_ss,Z_st,Z_tt)
  implicit none
  !> input variables
  integer, intent(in)                      :: n_v
  real(kind=8), dimension(n_v), intent(in) :: P_RR, P_RZ, P_ZZ, P_R, P_Z
  real(kind=8), intent(in)                 :: R_s, R_t, R_ss, R_st, R_tt
  real(kind=8), intent(in)                 :: Z_s, Z_t, Z_ss, Z_st, Z_tt
  !> output varibales
  real(kind=8), dimension(n_v), intent(out) :: P_ss, P_st, P_tt

  P_ss = P_RR*R_s*R_s + 2.d0*P_RZ*R_s*Z_s + P_ZZ*Z_s*Z_s + P_R*R_ss + P_Z*Z_ss 
  P_st = P_RR*R_t*R_s + P_RZ*(Z_t*R_s+R_t*Z_s) + P_ZZ*Z_t*Z_s + P_R*R_st + P_Z*Z_st 
  P_tt = P_RR*R_t*R_t + 2.d0*P_RZ*R_t*Z_t + P_ZZ*Z_t*Z_t + P_R*R_tt + P_Z*Z_tt 
  
end subroutine transform_second_derivatives_RZ_to_st

!---------------------------------------------------------------------------------------

end module mod_coordinate_transforms
