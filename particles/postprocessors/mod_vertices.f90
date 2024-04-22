!> the module mod_vertices contains all datatypes
!> and procedures common to light and gather points
module mod_vertices
implicit none

private
public :: vertices
!> public procedure only for unit testing
#ifdef UNIT_TESTS
public :: resize_vertices_noloss_seq
public :: resize_vertices_noloss_omp
#endif

!> Variable and type definitions ------------------------------------------
integer,parameter     :: n_x=3 !< number of coordinates

!> vertices: abstract class containing the basic types and procedures
!> defining light or gather points.
type,abstract :: vertices
  !> type variables
  integer               :: n_x=3             !< number of values of position vector 
  integer               :: n_times           !< number of time slices
  integer               :: n_vertices        !< total number of vertices equals per all time
  integer               :: n_property_vertex !< total number of properties per vertex
  integer,dimension(:),allocatable    :: n_active_vertices !< total number of active vertices per time
  real*8,dimension(:),allocatable     :: times             !< time of the time slices
  real*8,dimension(:,:,:),allocatable :: x                 !< position of the point (cartesian)
  real*8,dimension(:,:,:),allocatable :: properties        !< properties of the vertices
  contains
  !> type procedures
  procedure,pass(vert_inout) :: allocate_time_vector
  procedure,pass(vert_inout) :: allocate_x_properties
  procedure,pass(vert_inout) :: allocate_vertices
  procedure,pass(vert_inout) :: deallocate_time_vector
  procedure,pass(vert_inout) :: deallocate_x_properties
  procedure,pass(vert_inout) :: deallocate_vertices
  procedure,pass(vert_inout) :: resize_vertices_noloss
  procedure,pass(vert_inout) :: fit_tables_to_active_vertices
  procedure,pass(vert_1)     :: visibility_geometry_funct
end type vertices

!> Interfaces -------------------------------------------------------------
contains

!> Procedures -------------------------------------------------------------
!> allocate time vector and number of active vertices
!> inputs:
!>   vert_inout: (vertices) vertices to be allocated
!>   n_times:    (integer) number of times
!> outputs:
!>   n_vertices:  integer) number of vertices
subroutine allocate_time_vector(vert_inout,n_times)
  implicit none
  !> inputs-outputs:
  class(vertices),intent(inout) :: vert_inout
  !> inputs:
  integer,intent(in) :: n_times
  if(allocated(vert_inout%n_active_vertices)) deallocate(vert_inout%n_active_vertices)
  if(allocated(vert_inout%times)) deallocate(vert_inout%times)
  allocate(vert_inout%n_active_vertices(n_times))
  allocate(vert_inout%times(n_times))
  vert_inout%n_active_vertices = 0; vert_inout%times = 0.d0;
  vert_inout%n_times = n_times;
end subroutine allocate_time_vector

!> allocate x and properties tables and initilise them to 0.
!> Data loss is expected for preallocated tables
!> inputs:
!>   vert_inout:        (vertices) vertices to be allocated
!>   n_times:           (integer) number of times
!>   n_vertices:        (integer) number of vertices
!> outputs:
!>   vert_inout: (vertices) allocated vertices
subroutine allocate_x_properties(vert_inout,n_vertices)
  implicit none
  !> inputs
  integer,intent(in) :: n_vertices
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  !> check,allocate and initialize to 0 x and properties array
  if(allocated(vert_inout%x)) deallocate(vert_inout%x)
  if(allocated(vert_inout%properties)) deallocate(vert_inout%properties)
  allocate(vert_inout%x(vert_inout%n_x,n_vertices,vert_inout%n_times))
  allocate(vert_inout%properties(vert_inout%n_property_vertex,n_vertices,vert_inout%n_times))
  vert_inout%x = 0.d0; vert_inout%properties = 0.d0;
  vert_inout%n_vertices = n_vertices;
end subroutine allocate_x_properties

!> allocate all vertices and initialise them to zero
!> Data loss is expected for preallocated tables
!> inputs:
!>   vert_inout:        (vertices) vertices to be allocated
!>   n_times:           (integer) number of times
!>   n_vertices:        (integer) number of vertices
!> outputs:
!>   vert_inout: (vertices) allocated vertices
subroutine allocate_vertices(vert_inout,n_times,n_vertices)
  implicit none
  !> inputs
  integer,intent(in) :: n_times,n_vertices
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  call vert_inout%allocate_time_vector(n_times)
  call vert_inout%allocate_x_properties(n_vertices)
end subroutine allocate_vertices

!> deallocate time vector tables. Data loss is expected
!> inputs:
!>   vert_inout: (vertices) allocated vertices
!> outputs:
!>   vert_inout: (vertices) deallocated vertices
subroutine deallocate_time_vector(vert_inout)
  implicit none
  !> input-outputs
  class(vertices),intent(inout) :: vert_inout
  if(allocated(vert_inout%n_active_vertices)) deallocate(vert_inout%n_active_vertices)
  if(allocated(vert_inout%times)) deallocate(vert_inout%times)
  vert_inout%n_times = 0
end subroutine deallocate_time_vector

!> deallocate x and properties tables. Data loss is expected
!> inputs:
!>   vert_inout: (vertices) allocated vertices
!> outputs:
!>   vert_inout: (vertices) deallocated vertices
subroutine deallocate_x_properties(vert_inout)
  implicit none
  !> input-outputs
  class(vertices),intent(inout) :: vert_inout
  if(allocated(vert_inout%x)) deallocate(vert_inout%x)
  if(allocated(vert_inout%properties)) deallocate(vert_inout%properties)
  vert_inout%n_vertices = 0
end subroutine deallocate_x_properties

!> deallocate vertices tables. Data loss is expected
!> inputs:
!>   vert_inout: (vertices) allocated vertices
!> outputs:
!>   vert_inout: (vertices) deallocated vertices
subroutine deallocate_vertices(vert_inout)
  implicit none
  !> input-outputs
  class(vertices),intent(inout) :: vert_inout
  call vert_inout%deallocate_time_vector()
  call vert_inout%deallocate_x_properties()
end subroutine deallocate_vertices

!> resize the vertex tables without loss of data. This operation can be
!> very slow because it requires a copy of all tables data.
!> Sequetinal implementation.
!> inputs:
!>   vert_inout:   (vertices) vertex table with old sizes
!>   n_vertex_new: (integer) new size of vertex tables
!>   ifail:        (integer) failed to resize tables if =11
!> outputs:
!>  vert_inout:    (vertices) resized vertex tables
!>  ifail:         (integer) failed to resize if =11
subroutine resize_vertices_noloss(vert_inout,n_vertex_new,ifail)
  implicit none
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  integer,intent(inout)         :: ifail
  !> inputs
  integer,intent(in)            :: n_vertex_new 
  !> variables
  integer :: ii,jj,n_active_vertices_max
  real*8,dimension(:,:,:),allocatable :: x_table
  real*8,dimension(:,:,:),allocatable :: property_table

  !> call sequential or openmp vertion of the function 
#ifdef _OPENMP
  call resize_vertices_noloss_omp(vert_inout,n_vertex_new,ifail)
#else
  call resize_vertices_noloss_seq(vert_inout,n_vertex_new,ifail)
#endif

end subroutine resize_vertices_noloss

!> resize the vertex tables without loss of data. This operation can be
!> very slow because it requires a copy of all tables data.
!> Sequetinal implementation.
!> inputs:
!>   vert_inout:   (vertices) vertex table with old sizes
!>   n_vertex_new: (integer) new size of vertex tables
!>   ifail:        (integer) failed to resize tables if =11
!> outputs:
!>  vert_inout:    (vertices) resized vertex tables
!>  ifail:         (integer) failed to resize if =11
subroutine resize_vertices_noloss_seq(vert_inout,n_vertex_new,ifail)
  implicit none
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  integer,intent(inout)         :: ifail
  !> inputs
  integer,intent(in)            :: n_vertex_new 
  !> variables
  integer :: ii,jj,n_active_vertices_max
  real*8,dimension(:,:,:),allocatable :: x_table
  real*8,dimension(:,:,:),allocatable :: property_table

  !> initialisation
  n_active_vertices_max = maxval(vert_inout%n_active_vertices)
  !> check for possible data losses
  if(n_vertex_new.lt.n_active_vertices_max) then
    write(*,*) "Try to resize vertex tables but possible data loss detected!"
    ifail = 11
    return
  endif
  vert_inout%n_vertices = n_vertex_new
  allocate(x_table(vert_inout%n_x,n_active_vertices_max,vert_inout%n_times))
  allocate(property_table(vert_inout%n_property_vertex,n_active_vertices_max,vert_inout%n_times))

  !> copy data and resize tables
  do ii=1,vert_inout%n_times
    x_table(:,1:vert_inout%n_active_vertices(ii),ii) = &
    vert_inout%x(:,1:vert_inout%n_active_vertices(ii),ii)
    property_table(:,1:vert_inout%n_active_vertices(ii),ii) = &
    vert_inout%properties(:,1:vert_inout%n_active_vertices(ii),ii)
  enddo
  call vert_inout%allocate_x_properties(n_vertex_new)
  do ii=1,vert_inout%n_times
    vert_inout%x(:,1:vert_inout%n_active_vertices(ii),ii) = &
    x_table(:,1:vert_inout%n_active_vertices(ii),ii)
    vert_inout%properties(:,1:vert_inout%n_active_vertices(ii),ii) = &
    property_table(:,1:vert_inout%n_active_vertices(ii),ii)
  enddo

  !> cleanup
  deallocate(x_table); deallocate(property_table)
end subroutine resize_vertices_noloss_seq

!> resize the vertex tables without loss of data. This operation can be
!> very slow because it requires a copy of all tables data.
!> parallel implementation.
!> inputs:
!>   vert_inout:   (vertices) vertex table with old sizes
!>   n_vertex_new: (integer) new size of vertex tables
!>   ifail:        (integer) failed to resize tables if =11
!> outputs:
!>  vert_inout:    (vertices) resized vertex tables
!>  ifail:         (integer) failed to resize if =11
subroutine resize_vertices_noloss_omp(vert_inout,n_vertex_new,ifail)
  implicit none
  !> inputs-outputs
  class(vertices),intent(inout) :: vert_inout
  integer,intent(inout)         :: ifail
  !> inputs
  integer,intent(in)            :: n_vertex_new 
  !> variables
  integer :: ii,jj,n_active_vertices_max
  real*8,dimension(:,:,:),allocatable :: x_table
  real*8,dimension(:,:,:),allocatable :: property_table

  !> initialisation
  n_active_vertices_max = maxval(vert_inout%n_active_vertices)
  !> check for possible data losses
  if(n_vertex_new.lt.n_active_vertices_max) then
    write(*,*) "Try to resize vertex tables but possible data loss detected!"
    ifail = 11
    return
  endif
  vert_inout%n_vertices = n_vertex_new
  allocate(x_table(vert_inout%n_x,n_active_vertices_max,vert_inout%n_times))
  allocate(property_table(vert_inout%n_property_vertex,n_active_vertices_max,vert_inout%n_times))

  !> copy data and resize tables
  !$omp parallel default(shared) firstprivate(n_vertex_new) shared(ii,jj)
  !$omp do collapse(2)
  do ii=1,vert_inout%n_times
    do jj=1,n_active_vertices_max
      x_table(:,jj,ii) = vert_inout%x(:,jj,ii)
      property_table(:,jj,ii) = vert_inout%properties(:,jj,ii)
    enddo
  enddo
  !$omp end do
  !$omp single
  call vert_inout%allocate_x_properties(n_vertex_new)
  !$omp end single
  !$omp do collapse(2)
  do ii=1,vert_inout%n_times
    do jj=1,n_active_vertices_max
      vert_inout%x(:,jj,ii) = x_table(:,jj,ii)
      vert_inout%properties(:,jj,ii) = property_table(:,jj,ii)
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  !> cleanup
  deallocate(x_table); deallocate(property_table)
end subroutine resize_vertices_noloss_omp

!> fit the vertices table to the number of active tables.
!> This operation can be very slow because it requires 
!> a copy of all tables data
!> inputs:
!>   vert_inout: (vertices) vertex table with empty vertices
!>   ifail:      (integer) failed to resize tables if =11
!> outputs:
!>  vert_inout:  (vertices) fitted vertex tables
!>  ifail:       (integer) failed to resize if =11
subroutine fit_tables_to_active_vertices(vert_inout,ifail)
  implicit none
  class(vertices),intent(inout) :: vert_inout
  integer,intent(inout)         :: ifail
  call vert_inout%resize_vertices_noloss(maxval(vert_inout%n_active_vertices),ifail)
end subroutine fit_tables_to_active_vertices

!> compute the visibility and the geometry function in the case there
!> are no surfaces. Because there are no suraces implied, the returned
!> visibility function is always 0
!> inputs:
!>   vert_1:    (vertices) first vertex class
!>   vert_2:    (vertices) second vertex classs
!>   time_id_1: (integer) time index of the first vertex
!>   time_id_1: (integer) time index of the second vertex
!>   v_id_1:    (integer) first vertex index
!>   v_id_2:    (integer) second vertex index
!> outputs:
!>   vert_1:              (vertices) first vertex class
!>   visibility_geometry: (real8) visibility function always equal 1 divied by the
!>                        geometry function equal to the distance between two vertices
subroutine visibility_geometry_funct(vert_1,vert_2,time_id_1,time_id_2,&
v_id_1,v_id_2,visibility_geometry)
  implicit none
  class(vertices),intent(inout) :: vert_1
  class(vertices),intent(in)    :: vert_2
  integer,intent(in)            :: time_id_1,time_id_2,v_id_1,v_id_2
  real*8,intent(out)            :: visibility_geometry

  !> return the visibility function equal to 1
  visibility_geometry = 1.d0;
  visibility_geometry = visibility_geometry/&
             (((vert_1%x(1,v_id_1,time_id_1)-vert_2%x(1,v_id_2,time_id_2))**2)+ &
             ((vert_1%x(2,v_id_1,time_id_1)-vert_2%x(2,v_id_2,time_id_2))**2) + &
             ((vert_1%x(3,v_id_1,time_id_1)-vert_2%x(3,v_id_2,time_id_2))**2))
end subroutine visibility_geometry_funct

!>-------------------------------------------------------------------------
end module mod_vertices
