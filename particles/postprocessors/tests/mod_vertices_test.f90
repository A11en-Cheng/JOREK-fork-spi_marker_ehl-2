!> the mod_vertices_test module contains variables
!> and methods for testing all procedures contained
!> in mod_vertices.f90. Due to the fact that vertices
!> is an abstract class, the synchrotron vertex type
!> is used instead
module mod_vertices_test
use fruit
use mod_gyroaverage_synchrotron_light_dist_vertices, only: gyroaverage_synchrotron_light_dist
implicit none

private
public :: run_fruit_vertices

!> Variables ---------------------------------------------------------
integer,parameter :: n_x=3
integer,parameter :: code_resize_fail=11
integer,parameter :: n_times_sol=3
integer,parameter :: n_times_2_sol=5
integer,parameter :: n_vertices_sol=111
integer,parameter :: n_vertices_2_sol=222
integer,parameter :: n_properties_sol=11
integer,parameter :: n_vertex_shift=20
integer,dimension(n_times_sol),parameter   :: n_active_vertices_sol=(/50,71,33/)
integer,dimension(n_times_2_sol),parameter :: n_active_vertices_2_sol=&
                                              (/128,57,183,213,34/)
real*8,parameter :: tol_real8=5.d-16
real*8,dimension(2),parameter :: rng_interval=(/-2.d3,3.d3/)
type(gyroaverage_synchrotron_light_dist) :: vertex_sol,vertex_2_sol
integer                                  :: n_new_vertex_success,n_new_vertex_fail
integer                                  :: n_distances,ifail
!> position and property solution arrays
real*8,dimension(n_x,n_vertices_sol,n_times_sol)              :: x_sol
real*8,dimension(n_x,n_vertices_2_sol,n_times_2_sol)          :: x_2_sol
real*8,dimension(n_properties_sol,n_vertices_sol,n_times_sol) :: properties_sol
real*8,dimension(:),allocatable                               :: visibility_distances
real*8,dimension(:),allocatable                               :: visibility_distances_sol
!> Interfaces --------------------------------------------------------

contains
!> Fruit basket ------------------------------------------------------
!> fruit basked containing all set-up, test and teard-down methods
subroutine run_fruit_vertices()
  implicit none
  write(*,'(/A)') "  ... setting-up: vertices tests"
  call setup
  write(*,'(/A)') "  ... running: vertices tests"
  call run_test_case(test_de_allocate_time_vector,&
  'test_de_allocate_time_vector') 
  call run_test_case(test_de_allocates_vertex_x_properties,&
  'test_de_allocates_vertex_x_properties')
  call run_test_case(test_de_allocates_vertices,&
  'test_de_allocates_vertices')
#ifdef UNIT_TESTS
  call run_test_case(test_vertices_resize_nodataloss_seq,&
  'test_vertices_resize_nodataloss_seq')
  call run_test_case(test_vertices_resize_nodataloss_omp,&
  'test_vertices_resize_nodataloss_omp')
#endif
  call run_test_case(test_vertices_resize_nodataloss,&
  'test_vertices_resize_nodataloss')
  call run_test_case(test_fit_tables_to_active_vertices,&
  'test_fit_tables_to_active_vertices')
  write(*,'(/A)') "  ... tearing-down: vertices tests"
  call teardown
end subroutine run_fruit_vertices

!> Set-up and tear-down ----------------------------------------------
!> set-up unit test features
subroutine setup()
  use mod_gnu_rng, only: gnu_rng_interval
  implicit none
  !> variables 
  integer :: ii

  !> allocate arrays
  n_distances = sum(n_active_vertices_sol)*&
  sum(n_active_vertices_2_sol)*n_times_sol*n_times_2_sol
  allocate(visibility_distances(n_distances)); allocate(visibility_distances_sol(n_distances));
  !> initialize test variables
  ifail = 0
  vertex_sol%n_property_vertex = n_properties_sol !< set "fake" number of properties
  x_sol = 0.d0; properties_sol = 0.d0;
  n_new_vertex_success = maxval(n_active_vertices_sol)
  n_new_vertex_fail  = n_new_vertex_success - n_vertex_shift
  n_new_vertex_success = n_new_vertex_success + n_vertex_shift
  !> generate random values for position and property arrays
  do ii=1,n_times_sol
    call gnu_rng_interval(n_x,n_active_vertices_sol(ii),rng_interval,&
    x_sol(1:n_x,1:n_active_vertices_sol(ii),ii))
    call gnu_rng_interval(n_properties_sol,n_active_vertices_sol(ii),rng_interval,&
    properties_sol(1:n_properties_sol,1:n_active_vertices_sol(ii),ii))
  enddo
  !> generate random values for the second position array
  do ii=1,n_times_2_sol
    call gnu_rng_interval(n_x,n_active_vertices_2_sol(ii),rng_interval,&
    x_2_sol(1:n_x,1:n_active_vertices_2_sol(ii),ii))
  enddo
end subroutine setup

!> tear-down all unit test features
subroutine teardown()
  implicit none
  !> reset test variables
  ifail = 0; vertex_sol%n_property_vertex = 0;
  deallocate(visibility_distances); deallocate(visibility_distances_sol);
end subroutine teardown

!> Tests -------------------------------------------------------------
!> test allocate and deallocate time vector
subroutine test_de_allocate_time_vector()
   use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> test allocate time vector from unallocated from allocated vector
  call vertex_sol%allocate_time_vector(n_times_sol)
  call assert_equals(vertex_sol%n_times,n_times_sol,&
  "Error allocate vertices time vector from unallocated: n_times mismatch!")
  call assert_equals_allocatable_arrays(n_times_sol,vertex_sol%n_active_vertices,&
  0,"Error allocate vertices time vector from unallocated: n_active_vertices")
  call assert_equals_allocatable_arrays(n_times_sol,0.d0,vertex_sol%times,&
  "Error allocate vertices time vector from unallocated: times")
  !> test allocated time vector from allocated vector
  call vertex_sol%allocate_time_vector(n_times_2_sol) 
  call assert_equals(vertex_sol%n_times,n_times_2_sol,&
  "Error allocate vertices time vector from allocated: n_times mismatch!")
  call assert_equals_allocatable_arrays(n_times_2_sol,vertex_sol%n_active_vertices,&
  0,"Error allocate vertices time vector from allocated: n_active_vertices")
  call assert_equals_allocatable_arrays(n_times_2_sol,0.d0,vertex_sol%times,&
  "Error allocate vertices time vector from allocated: times")
  !> test deallocate time vector
  call vertex_sol%deallocate_time_vector
  call assert_equals(0,vertex_sol%n_times,&
  "Error deallocate vertices time vector: n_times not reset!")
  call assert_false(allocated(vertex_sol%n_active_vertices),&
  "Error deallocate vertices time vector: n_active_vertices allocated!")
  call assert_false(allocated(vertex_sol%times),&
  "Error deallocate vertices time vector: times allocated!")
end subroutine test_de_allocate_time_vector

!> test the allocation / deallocation of vertex position and properties
subroutine test_de_allocates_vertex_x_properties()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> initialise the time structure first
  call vertex_sol%allocate_time_vector(n_times_sol) 
  !> allocate x and properties from unallocated vertex
  call vertex_sol%allocate_x_properties(n_vertices_sol)
  call assert_equals(n_vertices_sol,vertex_sol%n_vertices,&
  "Error allocate vertices x-property vector from unallocated: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_x,n_vertices_sol,n_times_sol,0.d0,&
  vertex_sol%x,"Error allocate vertices  x-property vector from unallocated: x")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_sol,n_times_sol,&
  0.d0,vertex_sol%properties,&
  "Error allocate vertices  x-property vector from unallocated: properties")
  !> allocate x and properties from allocated vertes
  call vertex_sol%allocate_x_properties(n_vertices_2_sol)
  call assert_equals(n_vertices_2_sol,vertex_sol%n_vertices,&
  "Error allocate vertices x-property vector from allocated: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_x,n_vertices_2_sol,n_times_sol,0.d0,&
  vertex_sol%x,"Error allocate vertices  x-property vector from allocated: x")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_2_sol,n_times_sol,&
  0.d0,vertex_sol%properties,&
  "Error allocate vertices  x-property vector from allocated: properties")
  !> deallocate x and properties
  call vertex_sol%deallocate_x_properties
  call assert_equals(0,vertex_sol%n_vertices,&
  "Error deallocate vertices x-property vector: n_vertices not reset!")
  call assert_false(allocated(vertex_sol%x),&
  "Error deallocate vertices x-property vector: x allocated!")
  call assert_false(allocated(vertex_sol%properties),&
  "Error deallocate vertices x-property vector: properties allocated!") 
end subroutine test_de_allocates_vertex_x_properties

!> check allocation and deallocation of vertices
subroutine test_de_allocates_vertices()
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> test vertex allocation from unallocated
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  call assert_equals(n_times_sol,vertex_sol%n_times,&
  "Error allocate vertices from unallocated: n_times mismatch!")
  call assert_equals(n_vertices_sol,vertex_sol%n_vertices,&
  "Error allocate vertices from unallocated: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_times_sol,&
  vertex_sol%n_active_vertices,0,"Error allocate vertices from unallocated: n_active_vertices")
  call assert_equals_allocatable_arrays(n_times_sol,vertex_sol%times,&
  "Error allocate vertices from unallocated: times")
  call assert_equals_allocatable_arrays(n_x,n_vertices_sol,n_times_sol,&
  vertex_sol%x,"Error allocate vertices from unallocated: x size")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_sol,n_times_sol,&
  0.d0,vertex_sol%properties,"Error allocate vertices from unallocated: properties size")
  !> test vertex allocation from allocated
  call vertex_sol%allocate_vertices(n_times_2_sol,n_vertices_2_sol)
  call assert_equals(n_times_2_sol,vertex_sol%n_times,&
  "Error allocate vertices from allocated: n_times mismatch!")
  call assert_equals(n_vertices_2_sol,vertex_sol%n_vertices,&
  "Error allocate vertices from allocated: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_times_2_sol,vertex_sol%n_active_vertices,&
  0,"Error allocate vertices from allocated: n_active_vertices")
  call assert_equals_allocatable_arrays(n_times_2_sol,0.d0,vertex_sol%times,&
  "Error allocate vertices from allocated: times")
  call assert_equals_allocatable_arrays(n_x,n_vertices_2_sol,n_times_2_sol,0.d0,&
  vertex_sol%x,"Error allocate vertices from allocated: x size")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_2_sol,n_times_2_sol,&
  0.d0,vertex_sol%properties,&
  "Error allocate vertices from allocated: properties size")
  !> test vertex deallocation
  call vertex_sol%deallocate_vertices
  call assert_equals(0,vertex_sol%n_times,&
  "Error deallocate vertices: n_times not reset!")
  call assert_equals(0,vertex_sol%n_vertices,&
  "Error deallocate vertices: n_vertices not reset!")
  call assert_false(allocated(vertex_sol%n_active_vertices),&
  "Error deallocate vertices: n_active_vertices allocated!")
  call assert_false(allocated(vertex_sol%times),&
  "Error deallocate vertices: times allocated!")
  call assert_false(allocated(vertex_sol%x),&
  "Error deallocate vertices: x allocated!")
  call assert_false(allocated(vertex_sol%properties),&
  "Error deallocate vertices: properties allocated!")
end subroutine test_de_allocates_vertices

#ifdef UNIT_TESTS 
!> test vertex resize without data losses: sequential
subroutine test_vertices_resize_nodataloss_seq
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_vertices,            only: resize_vertices_noloss_seq 
  implicit none
  !> variables
  integer :: ii
  !> vertices initialisation
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  ifail =0; vertex_sol%x = x_sol; vertex_sol%properties = properties_sol;
  vertex_sol%n_active_vertices = n_active_vertices_sol;
  !> test the procedure failing
  call resize_vertices_noloss_seq(vertex_sol,n_new_vertex_fail,ifail)
  call assert_equals(code_resize_fail,ifail,"Error resize vertices seq: loss of data!")
  call assert_equals(n_vertices_sol,vertex_sol%n_vertices,&
  "Error resize vertices seq: unexpected update of n_vertices!")
  !> test the table resize success
  call resize_vertices_noloss_seq(vertex_sol,n_new_vertex_success,ifail)
  call assert_equals(n_new_vertex_success,vertex_sol%n_vertices,&
  "Error resize vertices seq: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_x,n_new_vertex_success,n_times_sol,vertex_sol%x,&
  "Error resize vertices seq: unexpected positions") 
  call assert_equals_allocatable_arrays(n_properties_sol,n_new_vertex_success,n_times_sol,&
  vertex_sol%properties,"Error resize vertices seq: unexpected properties")
  do ii=1,n_times_sol
    call assert_equals(x_sol(:,1:n_new_vertex_success,ii),vertex_sol%x(:,:,ii),n_x,&
    n_new_vertex_success,tol_real8, "Error resize vertices seq: positions are not the same!")
    call assert_equals(properties_sol(:,1:n_new_vertex_success,ii),vertex_sol%properties(:,:,ii),&
    n_x,n_new_vertex_success,tol_real8,"Error resize vertices seq: properties are not the same!")
  enddo
  !> deallocate everything
  call vertex_sol%deallocate_vertices 

end subroutine test_vertices_resize_nodataloss_seq

!> test vertex resize without data losses: parallel openmp
subroutine test_vertices_resize_nodataloss_omp
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  use mod_vertices,            only: resize_vertices_noloss_omp 
  implicit none
  !> variables
  integer :: ii
  !> vertices initialisation
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  ifail =0; vertex_sol%x = x_sol; vertex_sol%properties = properties_sol;
  vertex_sol%n_active_vertices = n_active_vertices_sol;
  !> test the procedure failing
  call resize_vertices_noloss_omp(vertex_sol,n_new_vertex_fail,ifail)
  call assert_equals(code_resize_fail,ifail,"Error resize vertices omp: loss of data!")
  call assert_equals(n_vertices_sol,vertex_sol%n_vertices,&
  "Error resize vertices omp: unexpected update of n_vertices!")
  !> test the table resize success
  call resize_vertices_noloss_omp(vertex_sol,n_new_vertex_success,ifail)
  call assert_equals(n_new_vertex_success,vertex_sol%n_vertices,&
  "Error resize vertices omp: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_x,n_new_vertex_success,n_times_sol,&
  vertex_sol%x,"Error resize vertices omp: unexpected positions")
  call assert_equals_allocatable_arrays(n_properties_sol,n_new_vertex_success,n_times_sol,&
  vertex_sol%properties,"Error resize vertices omp: unexpected properties")
  do ii=1,n_times_sol
    call assert_equals(x_sol(:,1:n_new_vertex_success,ii),vertex_sol%x(:,:,ii),n_x,&
    n_new_vertex_success,tol_real8, "Error resize vertices omp: positions are not the same!")
    call assert_equals(properties_sol(:,1:n_new_vertex_success,ii),vertex_sol%properties(:,:,ii),&
    n_x,n_new_vertex_success,tol_real8,"Error resize vertices omp: properties are not the same!")
  enddo 
  !> deallocate everything
  call vertex_sol%deallocate_vertices 

end subroutine test_vertices_resize_nodataloss_omp
#endif

!> test vertex resize without data losses: interface
subroutine test_vertices_resize_nodataloss
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> variables
  integer :: ii
  !> vertices initialisation
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  ifail =0; vertex_sol%x = x_sol; vertex_sol%properties = properties_sol;
  vertex_sol%n_active_vertices = n_active_vertices_sol;
  !> test the procedure failing
  call vertex_sol%resize_vertices_noloss(n_new_vertex_fail,ifail)
  call assert_equals(code_resize_fail,ifail,"Error resize vertices: loss of data!")
  call assert_equals(n_vertices_sol,vertex_sol%n_vertices,&
  "Error resize vertices: unexpected update of n_vertices!")
  !> test the table resize success
  call vertex_sol%resize_vertices_noloss(n_new_vertex_success,ifail)
  call assert_equals(n_new_vertex_success,vertex_sol%n_vertices,&
  "Error resize vertices: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_x,n_new_vertex_success,n_times_sol,&
  vertex_sol%x,"Error resize vertices: unexpected positions")
  call assert_equals_allocatable_arrays(n_properties_sol,n_new_vertex_success,n_times_sol,&
  vertex_sol%properties,"Error resize vertices: unexpected properties")
  do ii=1,n_times_sol
    call assert_equals(x_sol(:,1:n_new_vertex_success,ii),vertex_sol%x(:,:,ii),n_x,&
    n_new_vertex_success,tol_real8, "Error resize vertices: positions are not the same!")
    call assert_equals(properties_sol(:,1:n_new_vertex_success,ii),vertex_sol%properties(:,:,ii),&
    n_x,n_new_vertex_success,tol_real8,"Error resize vertices: properties are not the same!")
  enddo
  !> deallocate everything
  call vertex_sol%deallocate_vertices 

end subroutine test_vertices_resize_nodataloss

!> test fit table to active vertices
subroutine test_fit_tables_to_active_vertices
  use mod_assert_equals_tools, only: assert_equals_allocatable_arrays
  implicit none
  !> variables
  integer :: ii,n_vertices_expected
  !> vertices initialisation
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  ifail = 0; vertex_sol%x = x_sol; vertex_sol%properties = properties_sol;
  vertex_sol%n_active_vertices = n_active_vertices_sol;
  n_vertices_expected = maxval(n_active_vertices_sol)
  !> test the table resize success
  call vertex_sol%fit_tables_to_active_vertices(ifail)
  call assert_equals(0,ifail,"Error fit table to active vertices: unexpected failure!")
  call assert_equals(n_vertices_expected,vertex_sol%n_vertices,&
  "Error fit table to active vertices: n_vertices mismatch!")
  call assert_equals_allocatable_arrays(n_x,n_vertices_expected,n_times_sol,&
  vertex_sol%x,"Error fit table to active vertices: unexpected positions")
  call assert_equals_allocatable_arrays(n_properties_sol,n_vertices_expected,n_times_sol,&
  vertex_sol%properties,"Error fit table to active vertices: unexpected properties")
  do ii=1,n_times_sol
    call assert_equals(x_sol(:,1:n_vertices_expected,ii),vertex_sol%x(:,:,ii),n_x,&
    n_vertices_expected,tol_real8, "Error fit table to active vertices: positions are not the same!")
    call assert_equals(properties_sol(:,1:n_vertices_expected,ii),vertex_sol%properties(:,:,ii),&
    n_x,n_vertices_expected,tol_real8,"Error fit table to active vertices: properties are not the same!")
  enddo
  !> deallocate everything
  call vertex_sol%deallocate_vertices 

end subroutine test_fit_tables_to_active_vertices

!> test the visibility e geometry function calculation without surfaces
subroutine test_visibility_geometry_funct_nosurfaces()
  implicit none
  integer :: ii,jj,kk,pp,counter
  !> initialisations
  counter = 0
  call vertex_sol%allocate_vertices(n_times_sol,n_vertices_sol)
  call vertex_2_sol%allocate_vertices(n_times_2_sol,n_vertices_2_sol)
  vertex_sol%n_active_vertices=n_active_vertices_sol; vertex_sol%x=x_sol;
  vertex_2_sol%n_active_vertices=n_active_vertices_2_sol; vertex_2_sol%x=x_2_sol;
  !> compute the visibilit and geometry functions
  do ii=1,vertex_sol%n_times
    do jj=1,vertex_sol%n_active_vertices(ii)
      do kk=1,vertex_2_sol%n_times
        do pp=1,vertex_2_sol%n_active_vertices(kk)
          visibility_distances_sol(counter) = 1.d0/dot_product(vertex_sol%x(:,jj,ii)-&
          vertex_2_sol%x(:,pp,kk),vertex_sol%x(:,jj,ii)-vertex_2_sol%x(:,pp,kk))
          call vertex_sol%visibility_geometry_funct(vertex_2_sol,ii,kk,jj,pp,&
          visibility_distances(counter))
        enddo
      enddo
    enddo
  enddo
  !> check solution
  call assert_equals(visibility_distances_sol,visibility_distances,n_distances,tol_real8,&
  "Error compute visibility visibility/geometry funct no surfaces: visibility/geometry mismatch!")
  !> deallocate everything
  call vertex_sol%deallocate_vertices; call vertex_2_sol%deallocate_vertices; 
end subroutine test_visibility_geometry_funct_nosurfaces

!> Tools -------------------------------------------------------------
!>--------------------------------------------------------------------
end module mod_vertices_test

