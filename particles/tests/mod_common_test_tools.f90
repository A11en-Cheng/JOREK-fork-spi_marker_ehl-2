!> mod_common_test_tools contains and procedure common
!> to multiple unit tests
module mod_common_test_tools
implicit none

private
public :: omp_initialize_rngs,remove_file

contains
!> Procedures -------------------------------------------------------------
!> initialise the random number generators in omp loops
!> inputs:
!>   n_points_loc: (integer) number of points for rngs init
subroutine omp_initialize_rngs(n_points_loc,n_rngs,rngs,use_xor_time_pid_in)
  use mod_rng, only: type_rng
  use mod_random_seed,only: random_seed
  !$ use omp_lib
  implicit none
  !> inputs
  integer,intent(in)          :: n_points_loc,n_rngs  
  logical,intent(in),optional :: use_xor_time_pid_in 
  !> inputs-outputs
  class(type_rng),dimension(n_rngs),intent(inout) :: rngs !< random number generators
  !> variables
  integer :: n_threads,n_points_per_thread,thread_id,ifail
  logical :: use_xor_time_pid
  !> initialise the rngs using the pcg32
  n_threads = 1
  thread_id = 1
  use_xor_time_pid = .true. 
  if(present(use_xor_time_pid_in)) use_xor_time_pid = use_xor_time_pid_in
  !$omp parallel default(private) shared(rngs,n_threads,ifail) &
  !$omp firstprivate(use_xor_time_pid)
  !$ thread_id = omp_get_thread_num()+1
  !$ n_threads = omp_get_num_threads()
  n_points_per_thread = n_points_loc/n_threads
  call rngs(thread_id)%initialize(n_dims=n_points_per_thread,&
  seed=random_seed(use_xor_time_pid_in=use_xor_time_pid),&
  n_streams=n_threads,i_stream=thread_id,ierr=ifail)
  !$omp end parallel 
end subroutine omp_initialize_rngs

!> remove a file to the give path
!> inputs:
!>   filename: (character) name to the file to be removed
!>   rank:     (integer)(optional) rank of the mpi task
!>             default: 0
subroutine remove_file(filename,rank_in)
  implicit none
  character(len=*),intent(in) :: filename
  integer,intent(in),optional :: rank_in
  integer :: rank
  rank = 0; if(present(rank_in)) rank = rank_in
  if(rank.eq.0) call system("rm "//filename)
end subroutine remove_file

!>-------------------------------------------------------------------------
end module mod_common_test_tools
