!> mod_sampling_test contains all procedures used for
!> for testing the procedures contained in mod_samples
module mod_sampling_test
use fruit
use constants, only: PI,TWOPI
implicit none

private
public :: run_fruit_sampling

!> Variables -----------------------------------------
!> variables for testing the cone sampling methods
integer,parameter :: n_x=3
integer,parameter :: n_cos_half_angle=6
integer,parameter :: n_directions=4
integer,parameter :: n_origins=5
integer,parameter :: n_rays=21
integer,parameter :: n_samples=156
real*8,parameter  :: tol_real8=5.d-15
real*8,parameter  :: half_angle_lower_limit_sol=0d0
real*8,dimension(2),parameter :: colat_int=(/-PI,PI/)
real*8,dimension(2),parameter :: azimuth_int=(/0.d0,TWOPI/)
real*8,dimension(2),parameter :: half_angle_int=(/TWOPI/2.3d2,TWOPI/5.3d0/)
real*8,dimension(2),parameter :: length_int=(/1.d-1,3.5d1/)
real*8,dimension(2),parameter :: colat_sphere_int=(/PI/6.d0,2.d0*PI/3.d0/)
real*8,dimension(2),parameter :: azimuth_sphere_int=(/PI/1.1d1,TWOPI/1.5d0/)
real*8,dimension(2),parameter :: R_sphere_int=(/1.d0,3.2d2/)
real*8,dimension(3),parameter :: origin_min=(/-1.d0,2.d0,-7.d0/)
real*8,dimension(3),parameter :: origin_max=(/2.d0,5.d0,2.5d0/)
real*8,dimension(n_cos_half_angle,2) :: cos_half_angle_sol
real*8,dimension(n_x,n_directions) :: directions
real*8,dimension(n_x,n_origins)    :: origins
real*8,dimension(n_x,n_samples)    :: rand_sol

!> Interfaces ----------------------------------------

contains
!> Fruit basket --------------------------------------
!> fruit basket performing all set-ups, tests and
!> clean-up procedured
subroutine run_fruit_sampling()
  implicit none
  write(*,*) "  ... setting-up: sampling tests"
  call setup_cone
  call setup_sphere
  write(*,*) "  ... running: sampling tests"
  call run_test_case(test_sample_uniform_standard_cone,&
  'test_sample_uniform_standard_cone')
  call run_test_case(test_sample_uniform_direction_cone,&
  'test_sample_uniform_direction_cone')
  call run_test_case(test_sample_uniform_direction_length_cone,&
  'test_sample_uniform_direction_length_cone')
  call run_test_case(test_sample_uniform_sphere,&
  'test_sample_uniform_sphere')
  call run_test_case(test_sample_uniform_sphere_corona_rthetaphi,&
  'test_sample_uniform_sphere_corona_rthetaphi')
  call run_test_case(test_sample_uniform_sphere_corona_rcosphi,&
  'test_sample_uniform_sphere_corona_rcosphi')
  write(*,*) "  ... tearing-down: sampling tests"
end subroutine run_fruit_sampling

!> Set-up and tear-down ------------------------------
!> set up parameters for cone testing
subroutine setup_cone()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  integer :: ii
  real*8,dimension(n_cos_half_angle) :: half_angle_lower_limit
  real*8,dimension(n_cos_half_angle,2) :: half_angles
  real*8,dimension(n_directions)     :: colatitude,azimuth,lengths
  !> generate set of half solid angle cosines
  call gnu_rng_interval(n_cos_half_angle,half_angle_int,cos_half_angle_sol(:,1))
  half_angle_lower_limit = half_angle_lower_limit_sol
  call gnu_rng_interval(n_cos_half_angle,half_angle_lower_limit,&
  cos_half_angle_sol(:,1),cos_half_angle_sol (:,2)) 
  cos_half_angle_sol = cos(cos_half_angle_sol)
  !> generate set of origins
  do ii=1,n_origins
    call gnu_rng_interval(n_x,origin_min,origin_max,origins(:,ii))
  enddo
  !> generate set of direction vectors
  call gnu_rng_interval(n_directions,length_int,lengths)
  call gnu_rng_interval(n_directions,colat_int,colatitude)
  call gnu_rng_interval(n_directions,azimuth_int,azimuth)
  directions(1,:) = sin(colatitude)*cos(azimuth)
  directions(2,:) = sin(colatitude)*sin(azimuth)
  directions(3,:) = cos(colatitude)
  do ii=1,n_directions
    directions(:,ii) = lengths(ii)*directions(:,ii)
  enddo
end subroutine setup_cone

!> set up parameters for sphere testing
subroutine setup_sphere()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> variables
  !> generate set of random numbers
  call random_number(rand_sol)
end subroutine setup_sphere

!> Tests ---------------------------------------------
!> test the generation of random rays within a cone given
!> an origin and a direction vector and a length
subroutine test_sample_uniform_direction_length_cone()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  integer :: ii,jj,kk,pp
  real*8,dimension(3) :: u
  real*8,dimension(3) :: ray
  real*8,dimension(n_rays) :: cos_half_angle,ray_length
  !> computes rays with different origins, directions
  !> and half angles
  do pp=1,n_origins
    do kk=1,n_directions
      do jj=1,n_cos_half_angle
        do ii=1,n_rays
          call random_number(u)
          ray = sample_uniform_cone(cos_half_angle_sol(jj,:),u,&
          directions(:,kk),origins(:,pp),length_int)
          ray = ray - origins(:,pp)
          ray_length(ii) = norm2(ray)
          cos_half_angle(ii) = dot_product(ray/ray_length(ii),directions(:,kk)/norm2(directions(:,kk)))
        enddo
        call assert_true(all((cos_half_angle.ge.cos_half_angle_sol(jj,1)).and.&
        (cos_half_angle.le.cos_half_angle_sol(jj,2))),&
        "Error uniform sampling direction length cone: rays cosinus half angle out-of-bound!")
        call assert_true(all((ray_length.ge.length_int(1)).and.(ray_length.le.length_int(2))),&
        "Error uniform sampling direction length cone: rays length out-of-bound!")
      enddo
    enddo
  enddo
end subroutine test_sample_uniform_direction_length_cone

!> test the generation of random rays within a cone given
!> an origin and a direction vector
subroutine test_sample_uniform_direction_cone()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  integer :: ii,jj,kk,pp
  real*8,dimension(n_rays) :: ones=1.d0
  real*8,dimension(2) :: u
  real*8,dimension(3) :: ray
  real*8,dimension(n_rays) :: cos_half_angle,ray_length
  !> computes rays with different origins, directions
  !> and half angles
  do pp=1,n_origins
    do kk=1,n_directions
      do jj=1,n_cos_half_angle
        do ii=1,n_rays
          call random_number(u)
          ray = sample_uniform_cone(cos_half_angle_sol(jj,:),u,&
          directions(:,kk),origins(:,pp))
          ray = ray - origins(:,pp)
          ray_length(ii) = norm2(ray)
          cos_half_angle(ii) = dot_product(ray,directions(:,kk)/norm2(directions(:,kk)))
        enddo
        call assert_true(all((cos_half_angle.ge.cos_half_angle_sol(jj,1)).and.&
        (cos_half_angle.le.cos_half_angle_sol(jj,2))),&
        "Error uniform sampling direction cone: rays cosinus half angle out-of-bound!")
        call assert_equals(ones,ray_length,n_rays,tol_real8,&
        "Error uniform sampling direction cone: rays length not unitary!")
      enddo
    enddo
  enddo
end subroutine test_sample_uniform_direction_cone

!> test generation of random rays in the stantard cone
subroutine test_sample_uniform_standard_cone()
  use mod_sampling, only: sample_uniform_cone
  implicit none
  real*8,dimension(n_rays) :: ones=1.d0
  integer :: ii,jj
  real*8,dimension(3),parameter :: z_dir=(/0.d0,0.d0,1.d0/)
  real*8,dimension(2) :: u
  real*8,dimension(3) :: ray
  real*8,dimension(n_rays) :: cos_half_angle,ray_length
  !> compute rays and half angles
  do jj=1,n_cos_half_angle
    do ii=1,n_rays
      call random_number(u)
      ray = sample_uniform_cone(cos_half_angle_sol(jj,:),u)
      ray_length(ii) = norm2(ray)
      cos_half_angle(ii) = dot_product(ray,z_dir)
    enddo
    !> check correctness
    call assert_true(all((cos_half_angle.ge.cos_half_angle_sol(jj,1)).and.&
    (cos_half_angle.le.cos_half_angle_sol(jj,2))),&
    "Error uniform sampling standard cone: rays cosinus half angle out-of-bound!")
    call assert_equals(ones,ray_length,n_rays,tol_real8,&
    "Error uniform sampling standard cone: rays length not unitary!")
  enddo
end subroutine test_sample_uniform_standard_cone

!> test random sphere coordinates
subroutine test_sample_uniform_sphere()
  use mod_sampling, only: sample_uniform_sphere
  implicit none
  integer                         :: ii
  real*8,dimension(2)             :: cos_theta
  real*8,dimension(n_x)           :: RThetaPhi
  real*8,dimension(n_x,n_samples) :: rand_test
  !> initialisation
  cos_theta = cos(colat_sphere_int)
  !> sampling the sphere and recompute the uniform random number
  do ii=1,n_samples
    RThetaPhi = sample_uniform_sphere(R_sphere_int(2),&
    cos_theta,azimuth_sphere_int,rand_sol(:,ii))
    rand_test(:,ii) = (/(RThetaPhi(1)/R_sphere_int(2))**3.d0,&
    (cos(RThetaPhi(2))-cos_theta(1))/(cos_theta(2)-cos_theta(1)),&
    (RThetaPhi(3)-azimuth_sphere_int(1))/(azimuth_sphere_int(2)-azimuth_sphere_int(1))/)
  enddo
  !> checks
  call assert_equals(rand_sol,rand_test,n_x,n_samples,tol_real8,&
  "Error uniform sampling sphere: random number mismatch!")
end subroutine test_sample_uniform_sphere

!> test random sphere corona angular coordinates
subroutine test_sample_uniform_sphere_corona_rthetaphi()
  use mod_sampling, only: sample_uniform_sphere_corona_rthetaphi
  implicit none
  !> variables
  integer                         :: ii
  real*8,dimension(2)             :: R_cube,cos_theta
  real*8,dimension(n_x)           :: RThetaPhi
  real*8,dimension(n_x,n_samples) :: rand_test
  !> initialisation
  cos_theta = cos(colat_sphere_int); R_cube = R_sphere_int**3.d0;
  !> sampling the sphere corona angular coordinates and recompute the random number
  do ii=1,n_samples
    RThetaPhi = sample_uniform_sphere_corona_rthetaphi(R_cube,cos_theta,&
    azimuth_sphere_int,rand_sol(:,ii))
    rand_test(:,ii) = (/(RThetaPhi(1)**3.d0-R_cube(1))/(R_cube(2)-R_cube(1)),&
    (cos(RThetaPhi(2))-cos_theta(1))/(cos_theta(2)-cos_theta(1)),&
    (RThetaPhi(3)-azimuth_sphere_int(1))/(azimuth_sphere_int(2)-azimuth_sphere_int(1))/)
  enddo
  !> checks
  call assert_equals(rand_sol,rand_test,n_x,n_samples,tol_real8,&
  "Error uniform sampling sphere corona angles: random number mismatch!")
end subroutine test_sample_uniform_sphere_corona_rthetaphi

!> test random sphere corona cosinus coordinates
subroutine test_sample_uniform_sphere_corona_rcosphi()
  use mod_sampling, only: sample_uniform_sphere_corona_rcosphi
  implicit none
  integer                         :: ii
  real*8,dimension(2)             :: R_cube,cos_theta
  real*8,dimension(n_x)           :: RCosPhi
  real*8,dimension(n_x,n_samples) :: rand_test
  !> initialisation
  cos_theta = cos(colat_sphere_int); R_cube = R_sphere_int**3.d0;
  !> sampling the sphere corona cosine coordinates and recompute random number
  do ii=1,n_samples
    RCosPhi = sample_uniform_sphere_corona_rcosphi(R_cube,cos_theta,&
    azimuth_sphere_int,rand_sol(:,ii))
    rand_test(:,ii) = (/(RCosPhi(1)**3.d0-R_cube(1))/(R_cube(2)-R_cube(1)),&
    (RCosPhi(2)-cos_theta(1))/(cos_theta(2)-cos_theta(1)),&
    (RCosPhi(3)-azimuth_sphere_int(1))/(azimuth_sphere_int(2)-azimuth_sphere_int(1))/)
  enddo
  !> checks
  call assert_equals(rand_sol,rand_test,n_x,n_samples,tol_real8,&
  "Error uniform sampling sphere corona cosine: random number mismatch!")
end subroutine test_sample_uniform_sphere_corona_rcosphi

!> Tools ---------------------------------------------
!>----------------------------------------------------
end module mod_sampling_test
