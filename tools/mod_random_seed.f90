!> Module containing routines to generate a random seed based on urandom or
!> some pid and time xor magic.
module mod_random_seed
  use phys_module, only: use_manual_random_seed, manual_seed
  implicit none
  private
  public :: random_seed
contains
  !> Try some methods to get a nice random seed
  function random_seed(use_xor_time_pid_in) result(seed)
    implicit none
    logical,intent(in),optional :: use_xor_time_pid_in
    integer                     :: seed,ierr
    logical                     :: use_xor_time_pid
    ierr = 1; use_xor_time_pid = .true.;

    if (use_manual_random_seed) then
      seed = manual_seed
    else
      if(present(use_xor_time_pid_in)) use_xor_time_pid = use_xor_time_pid_in;
      if(.not.use_xor_time_pid) call read_urandom_int(seed, ierr)
      if (ierr .ne. 0) seed = xor_time_pid()
    endif
end function random_seed

  !> Read an int from /dev/urandom
  subroutine read_urandom_int(seed, ierr)
    implicit none
    integer, intent(out) :: seed
    integer, intent(out) :: ierr
    integer :: un

    open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=ierr)
    if (ierr == 0) then
      read(un) seed
      close(un)
    end if
  end subroutine read_urandom_int

  !> Fallback method to generate a random seed
  !> xor the time and PID together, and run that through a lcg
  !> Does not work for openmp, but for that streams should be used
  !> with the PCG generator
  function xor_time_pid()
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_int32_t
    implicit none
    integer(c_int64_t), save :: t = 0
    integer :: dt(8)
    integer :: xor_time_pid
    interface
      function fgetpid() bind(C)
        import c_int32_t
        integer(c_int32_t) :: fgetpid
      end function fgetpid
    end interface

    ! Reuse the value of t so that subsequent calls at the same time
    ! produce different seeds. This is a bad RNG but it's only a seed
    if (t == 0) then
      call system_clock(t)
      if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
      end if
      t = ieor(t, int(fgetpid(), kind(t)))
    end if
    xor_time_pid = lcg(t)
  end function xor_time_pid

  !> Linear congruential generator for seeding
  function lcg(s)
    use, intrinsic :: iso_fortran_env, only: int64
    integer :: lcg
    integer(int64), intent(inout) :: s
    if (s == 0) then
      s = 104729
    else
      s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
end module mod_random_seed
