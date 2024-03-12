module mod_gather_bcast_mpi_test
use fruit
use fruit_mpi
use mpi
implicit none
public :: run_fruit_gather_bcast_mpi

!> Variables and type definitions ---------------------------
integer,parameter :: master_rank=0
integer,parameter :: default_int=-1
integer           :: rank_loc,n_tasks_loc,ifail_loc
integer,dimension(:),allocatable :: mpi_task_ranks
integer,dimension(:),allocatable :: mpi_task_ranks_sol
!> Interfaces -----------------------------------------------
contains
!> Fruit basket ---------------------------------------------
!> fruit basket executing all set-ups, tests and tear-downs
subroutine run_fruit_gather_bcast_mpi(rank,n_tasks,ifail)
  use mpi
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  if(rank.eq.master_rank) write(*,*) "  ... setting-up: gather-broadcast mpi tests"
  call setup(rank,n_tasks,ifail)
  call MPI_Barrier(MPI_COMM_WORLD,ifail)
  if(rank.eq.master_rank) write(*,*) "  ... running: gather-broadcast mpi tests"
  call run_test_case(test_gather_bcast_mpi,'test_gather_bcast_mpi')
  call MPI_Barrier(MPI_COMM_WORLD,ifail)
  if(rank.eq.master_rank) write(*,*) "  ... tearing-down: gather-broadcast mpi tests"
  call teardown(rank,n_tasks,ifail)
end subroutine run_fruit_gather_bcast_mpi

!> Set-up and tear-down -------------------------------------
!> setup test features
subroutine setup(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs:
  integer,intent(in) :: rank,n_tasks
  !> variables:
  integer :: ii
  !> copy mpi data to module variables, the local copy is required because
  !> the fruit run_test_case does not accept arguments for the test method
  rank_loc = rank; n_tasks_loc = n_tasks; ifail_loc = ifail;
  !> allocate and initialise test and solution arrays
  if(allocated(mpi_task_ranks)) deallocate(mpi_task_ranks)
  allocate(mpi_task_ranks(n_tasks_loc)); mpi_task_ranks = default_int;
  if(allocated(mpi_task_ranks_sol)) deallocate(mpi_task_ranks_sol)
  allocate(mpi_task_ranks_sol(n_tasks_loc))
  mpi_task_ranks_sol = (/(ii-1,ii=1,n_tasks_loc)/)
end subroutine setup

!> teardown test features
subroutine teardown(rank,n_tasks,ifail)
  implicit none
  !> inputs-outputs:
  integer,intent(inout) :: ifail
  !> inputs: 
  integer,intent(in) :: rank,n_tasks
  !> returning the overall ifail
  ifail = ifail_loc
  !> deallocate everything
  rank_loc = 0; n_tasks_loc = 0;
  if(allocated(mpi_task_ranks))     deallocate(mpi_task_ranks)
  if(allocated(mpi_task_ranks_sol)) deallocate(mpi_task_ranks_sol)
end subroutine teardown

!> Tests ----------------------------------------------------
!> test MPI collective gather and broadcast operations
subroutine test_gather_bcast_mpi()
  use mpi
  implicit none
  !> gather all data on master and broadcast the resulting array to slaves
  call MPI_Gather(rank_loc,1,MPI_INTEGER,mpi_task_ranks,1,MPI_INTEGER,&
  master_rank,MPI_COMM_WORLD,ifail_loc)
  call MPI_Bcast(mpi_task_ranks,n_tasks_loc,MPI_INTEGER,&
  master_rank,MPI_COMM_WORLD,ifail_loc)
  !> check results
  call assert_equals(mpi_task_ranks,mpi_task_ranks_sol,n_tasks_loc,&
  "Error MPI gather and broadcast operations: task rank array mismatch!")
end subroutine test_gather_bcast_mpi
!> ----------------------------------------------------------
end module mod_gather_bcast_mpi_test
