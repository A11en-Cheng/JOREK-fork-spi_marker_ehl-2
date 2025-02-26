!> Particle input-output event module
!> The documentation of the module can be found at: https://www.jorek.eu/wiki/doku.php?id=jorek-particles_i_o
module mod_io_actions
use hdf5, only: H5F_ACC_TRUNC_F
use mod_particle_io
use mod_event
use mod_particle_sim
implicit none
private
public io_action, read_action, write_action, &
    get_filename

type, extends(action), abstract :: io_action
  character(len=120), public :: filename = ''         !< Filename to use
  character(len=80), public  :: basename = 'part'     !< If no filename, use basename + digits
  integer, public            :: decimal_digits = 3    !< Number of decimals before the point in timestamp
  integer, public            :: fractional_digits = 8 !< Number of decimals after the point
  character(len=5), public   :: extension = '.h5'     !< I/O file type extension
  integer, public            :: mpi_comm_io           !< mpi communicator
  integer, public            :: mpi_info_io           !< mpi information 
  logical, public            :: use_hdf5_access_properties = .false. !< if true a new hdf5 property list is created (required for MPIO), H5P_DEFAULT_F otherwise
  contains
    procedure :: get_filename
end type io_action

type, extends(io_action) :: read_action
  real*8  :: time             !< used with the formats from io_action if filename is unset
  logical :: test   = .false. !< if true avoid to read adas for unit testing
contains
  procedure :: do => do_read_action
end type read_action
interface read_action
  module procedure new_read_action
end interface read_action

type, extends(io_action) :: write_action
  integer, public :: file_access = -1 !<H5F_ACC_TRUNC_F !< define the file access behaviour of HDF5, default: truncate file if exists, create otherwise
  logical, public :: mpio_collective = .true. !< if true MPIO calls are collective (individual otherwise)
  logical, public :: use_native_hdf5_mpio = .false. !< if true native HDF5-MPIO is used, MPI_Gatherv is used otherwise
contains
  procedure :: do => do_write_action
end type write_action
interface write_action
  module procedure new_write_action
end interface write_action

contains
!> Write a filename consisting of this%basename and time as floating-point
!> number, without spaces (tricky), with this%decimal_digits before the `.` and
!> this%fractional_digits behind the `.`. Set these two to 0 to just set the name.
function get_filename(this, time) result(filename)
  class(io_action), intent(in) :: this
  real*8, intent(in)           :: time
  character(len=120)           :: filename
  character(len=20)            :: format
  if (this%decimal_digits .eq. 0 .and. this%fractional_digits .eq. 0) then
    write(filename,'(A,A)') trim(this%basename), this%extension
  else if (this%fractional_digits .eq. 0) then
    write(format,'(A,I0,A)') '(A,I', this%decimal_digits, ',A)'
    write(filename,trim(format)) trim(this%basename), int(time), this%extension
  else
    write(format,'(A,I0,A,I0,A,I0,A)') '(A,I', this%decimal_digits, '.', this%decimal_digits, &
        ',f0.', this%fractional_digits, ',A)'
    write(filename,trim(format)) trim(this%basename), int(time), time-real(int(time)), this%extension
  end if
end function get_filename

!> Constructor for read_action.
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_read_action(filename, basename, decimal_digits, fractional_digits, &
extension,use_hdf5_access_properties_in, mpi_comm_in, mpi_info_in,test_in)
  use mpi
  implicit none
  type(read_action) :: new_read_action
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer,          intent(in), optional :: decimal_digits
  integer,          intent(in), optional :: fractional_digits
  character(len=*), intent(in), optional :: extension
  integer,          intent(in), optional :: mpi_comm_in,mpi_info_in
  logical,          intent(in), optional :: test_in,use_hdf5_access_properties_in
  new_read_action%mpi_comm_io = MPI_COMM_WORLD
  new_read_action%mpi_info_io = MPI_INFO_NULL
  if (present(filename)) new_read_action%filename = filename
  if (present(basename)) new_read_action%basename = basename
  if (present(decimal_digits)) new_read_action%decimal_digits = decimal_digits
  if (present(fractional_digits)) new_read_action%fractional_digits = fractional_digits
  if (present(extension)) new_read_action%extension = extension
  if (present(use_hdf5_access_properties_in)) new_read_action%use_hdf5_access_properties = use_hdf5_access_properties_in
  if (present(mpi_comm_in)) new_read_action%mpi_comm_io = mpi_comm_in
  if (present(mpi_info_in)) new_read_action%mpi_info_io = mpi_info_in
  if (present(test_in)) new_read_action%test = test_in
  new_read_action%name = "ReadAction"
  new_read_action%log = .true.
end function new_read_action

!> Action for reading the simulation
subroutine do_read_action(this, sim, ev)
  class(read_action), intent(inout) :: this
  type(particle_sim), intent(inout) :: sim
  type(event), intent(inout), optional :: ev
  if (len_trim(this%filename) .eq. 0) then
    call read_simulation_hdf5(sim, trim(this%get_filename(sim%time)), &
    use_hdf5_access_properties=this%use_hdf5_access_properties, &
    mpi_comm_in=this%mpi_comm_io, mpi_info_in=this%mpi_info_io, test_in=this%test)
  else
    call read_simulation_hdf5(sim, trim(this%filename), &
    use_hdf5_access_properties=this%use_hdf5_access_properties, &
    mpi_comm_in=this%mpi_comm_io, mpi_info_in=this%mpi_info_io, test_in=this%test)
  endif
end subroutine do_read_action

!> Constructor for write_action
!> Be sure to use keyword arguments when initializing, to avoid confusion
function new_write_action(filename, basename, decimal_digits, fractional_digits, extension, &
file_access_in, use_native_hdf5_mpio_in, use_hdf5_access_properties_in, mpi_comm_in, mpi_info_in,&
mpio_collective_in)
  use mpi
  implicit none
  type(write_action) :: new_write_action
  character(len=*), intent(in), optional :: filename
  character(len=*), intent(in), optional :: basename
  integer,          intent(in), optional :: decimal_digits
  integer,          intent(in), optional :: fractional_digits
  character(len=*), intent(in), optional :: extension
  integer,          intent(in), optional :: file_access_in,mpi_comm_in,mpi_info_in
  logical,          intent(in), optional :: use_native_hdf5_mpio_in,mpio_collective_in
  logical,          intent(in), optional :: use_hdf5_access_properties_in
  new_write_action%name = "WriteAction"
  new_write_action%log = .true.
  new_write_action%file_access = H5F_ACC_TRUNC_F
  new_write_action%mpi_comm_io = MPI_COMM_WORLD
  new_write_action%mpi_info_io = MPI_INFO_NULL
  if (present(filename)) new_write_action%filename = filename
  if (present(basename)) new_write_action%basename = basename
  if (present(decimal_digits)) new_write_action%decimal_digits = decimal_digits
  if (present(fractional_digits)) new_write_action%fractional_digits = fractional_digits
  if (present(extension)) new_write_action%extension = extension
  if (present(file_access_in)) new_write_action%file_access = file_access_in
  if (present(use_native_hdf5_mpio_in)) new_write_action%use_native_hdf5_mpio = use_native_hdf5_mpio_in 
  if (present(use_hdf5_access_properties_in)) new_write_action%use_hdf5_access_properties = &
  use_hdf5_access_properties_in
  if (present(mpi_comm_in)) new_write_action%mpi_comm_io = mpi_comm_in
  if (present(mpi_info_in)) new_write_action%mpi_info_io = mpi_info_in
  if (present(mpio_collective_in)) new_write_action%mpio_collective = mpio_collective_in
end function new_write_action

!> Action for writing the simulation
subroutine do_write_action(this, sim, ev)
  class(write_action), intent(inout) :: this
  type(particle_sim), intent(inout)  :: sim
  type(event), intent(inout), optional :: ev
  if (len_trim(this%filename) .eq. 0) then
    call write_simulation_hdf5(sim, trim(this%get_filename(sim%time)),&
    file_access_in=this%file_access,&
    use_native_hdf5_mpio_in=this%use_native_hdf5_mpio,&
    use_hdf5_access_properties=this%use_hdf5_access_properties,&
    collective_mpio_in=this%mpio_collective,&
    mpi_comm_in=this%mpi_comm_io, mpi_info_in=this%mpi_info_io)
  else
    call write_simulation_hdf5(sim,trim(this%filename),file_access_in=this%file_access,&
    use_native_hdf5_mpio_in=this%use_native_hdf5_mpio,&
    use_hdf5_access_properties=this%use_hdf5_access_properties,&
    collective_mpio_in=this%mpio_collective,&
    mpi_comm_in=this%mpi_comm_io, mpi_info_in=this%mpi_info_io)
  endif
end subroutine do_write_action
end module mod_io_actions
