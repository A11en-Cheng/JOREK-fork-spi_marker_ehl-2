! The module mod_mpi_tools_test contains variables and procedures 
! used for testing the routines of the modules mod_mpi_tools
module mod_mpi_tools_mpi_test
use fruit
implicit none

private
public :: run_fruit_mpi_tools_mpi

! Variable --------------------------------------------------------
integer :: my_id,n_tasks
real*8 :: start_time,mpi_time,stop_time
contains

! Test basket ----------------------------------------------------
! run_fruit_mpi_tools executes the mpi_tools test set-up,
! tear-down and run the tests
subroutine run_fruit_mpi_tools_mpi(rank_glob,n_tasks_glob,ifail_glob)
  implicit none
  integer,intent(in)    :: rank_glob,n_tasks_glob
  integer,intent(inout) :: ifail_glob

  ! execute setup -> tests -> teardown
  write(*,"(/A)") "  ... setting-up: mpi_tools tests"
  call setup(rank_glob,n_tasks_glob,ifail_glob)
  write(*,"(/A)") "  ... running: mpi_tools tests"
  call run_test_case(test_init_mpi_threads,'test_init_mpi_threads')
  call run_test_case(test_get_mpi_wtime,'test_get_mpi_wtime')
  write(*,"(/A)") "  ... tearing-down: mpi_tools tests" 
  call teardown(rank_glob,n_tasks_glob,ifail_glob)

end subroutine run_fruit_mpi_tools_mpi

! Set-up and tear-down -------------------------------------------

! the setup procedure, initiliase the common test variables
subroutine setup(rank_glob,n_tasks_glob,ifail_glob)
  implicit none
  integer,intent(in)    :: rank_glob,n_tasks_glob
  integer,intent(inout) :: ifail_glob

  ! variables initialise to 0
  my_id = -999; n_tasks = -999; start_time = -9.d2;
  mpi_time = -9.d3; stop_time = -9.d3;
end subroutine setup

! tear-down procedure, clean-up all variables and closes
! the mpi communicator if initialised
subroutine teardown(rank_glob,n_tasks_glob,ifail_glob)
  implicit none
  integer,intent(in)    :: rank_glob,n_tasks_glob
  integer,intent(inout) :: ifail_glob
  logical               :: initialised
  integer               :: ierr
  my_id=-999; n_tasks = -999; start_time = -9.d2;
  mpi_time = -9.d3; stop_time = -9.d3;
end subroutine teardown

! Tests ----------------------------------------------------------

! test_mpi_thread_init tests the mpi thread initialisation
subroutine test_init_mpi_threads()
  use mpi
  use mod_mpi_tools, only: init_mpi_threads
  implicit none
  ! variables
  integer :: my_id_2,n_tasks_2,ierr
  logical :: is_initialised

  ! initialisation
  ierr = 0; my_id_2 = -999; n_tasks_2 = -999;

  ! initialise the mpi threads
  call init_mpi_threads(my_id,n_tasks,ierr)
  call MPI_initialized(is_initialised,ierr)

  ! checks
  call assert_true(is_initialised,"Error: init mpi threads: COMM_WORLD not initialised")
  call assert_equals(0,ierr,"Error: init mpi threads: en error occurred")
  call assert_true(my_id.ge.0,"Error: my_id must be non-negative")
  call assert_true(n_tasks.gt.0,"Error: number of tasks must be positive")

  ! extract time
  call init_mpi_threads(my_id_2,n_tasks_2,ierr,start_time)

  ! checks
  call assert_equals(0,ierr,"Error: init mpi threads: en error occurred")
  call assert_true(my_id_2.ge.0,"Error: my_id must be non-negative")
  call assert_true(n_tasks_2.gt.0,"Error: number of tasks must be positive")
  call assert_equals(my_id_2,my_id,"Error: tasks id must be equal")
  call assert_equals(n_tasks_2,n_tasks,"Error: number of tasks must be equal")
  call assert_true(start_time.ge.0.d0,"Error: the start time must be non-negative")
  
end subroutine test_init_mpi_threads

!> test getting time from mpi
subroutine test_get_mpi_wtime()
  use mod_mpi_tools, only: get_mpi_wtime
  implicit none
  mpi_time = get_mpi_wtime()
  call assert_true(mpi_time.ge.0.d0,"Error: the get time must be non negative")
end subroutine test_get_mpi_wtime

!> test finilising mpi with time
subroutine test_finalize_mpi_threads()
  use mpi
  use mod_mpi_tools, only: finalize_mpi_threads
  implicit none
  integer :: ierr
  ierr = 0; stop_time = -9.d3;
  call finalize_mpi_threads(ierr,stop_time)
  call assert_equals(0,ierr,"Error: mpi threads finalize: en error occurred")
  call assert_true(stop_time.ge.0.d0,"Error: mpi threads finalize: stop time must be non negative")
end subroutine test_finalize_mpi_threads

!-----------------------------------------------------------------

end module mod_mpi_tools_mpi_test
