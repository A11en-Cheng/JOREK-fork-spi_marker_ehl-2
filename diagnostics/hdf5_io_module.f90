!---------------------------------------------
! file : HDF5_io.f90
! date : 25/01/2006
!  array saving and reading in HDF5 format
!
! Adapted from GYSELA routines
!---------------------------------------------
module hdf5_io_module
#ifdef USE_HDF5  
  use HDF5
  implicit none

  !******************************
  contains
  !******************************

  !---------------------------------------- 
  ! create HDF5 file 
  !----------------------------------------
  subroutine HDF5_create(filename,file_id,ierr)
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr

    integer :: ierr_HDF5

    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)

    !*** Create a new file using default properties ***
    call H5Fcreate_f(trim(filename)//char(0), &
      H5F_ACC_TRUNC_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_create


  !---------------------------------------- 
  ! open HDF5 file 
  !----------------------------------------
  subroutine HDF5_open(filename,file_id,ierr)
    character(LEN=*) , intent(in)  :: filename  ! file name
    integer(HID_T)   , intent(out) :: file_id   ! file identifier
    integer, optional, intent(out) :: ierr

    integer :: ierr_HDF5

    !*** Initialize fortran interface ***
    call H5open_f(ierr_HDF5)

    !*** open the HDF5 file ***
    call H5Fopen_f(trim(filename)//char(0), &
      H5F_ACC_RDONLY_F,file_id,ierr_HDF5)
    if (present(ierr)) ierr = ierr_HDF5
  end subroutine HDF5_open


  !---------------------------------------- 
  ! close HDF5 file 
  !----------------------------------------
  subroutine HDF5_close(file_id)
    integer(HID_T), intent(in) :: file_id   ! file identifier

    integer :: error   ! error flag

    call H5Fclose_f(file_id,error)
  end subroutine HDF5_close

  !----------------------------------------
  ! check if dataset is in HDF5 file
  !----------------------------------------
  subroutine HDF5_is_dataset_in_file(file_id,dsetname,in_file)
    use H5LT, only: h5ltfind_dataset_f
    !> inputs:
    integer(HID_T),intent(in)   :: file_id
    character(len=*),intent(in) :: dsetname
    !> outputs:
    logical,intent(out) :: in_file
    in_file = h5ltfind_dataset_f(file_id,trim(dsetname)).ge.1
  end subroutine HDF5_is_dataset_in_file

  !---------------------------------------- 
  ! extract dataset rank and shape
  !----------------------------------------
  subroutine HDF5_extract_dataset_rank_shape(file_id,rank,dims,dsetname)
    !> inputs:
    character(len=*),intent(in) ::dsetname
    integer(HID_T), intent(in)  :: file_id
    !> outputs:
    integer,intent(out) :: rank !< rank (number of indexes) of the dataset
    integer(HSIZE_T),dimension(:),allocatable,intent(out) :: dims !< shape of the dataset
    !> variables:
    integer(HID_T) :: dataset_id,dataspace_id !< datatset,dataspace identifier
    integer        :: error !< error flag
    integer(HSIZE_T),dimension(:),allocatable :: maxdims !< maximum dataset shape
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    !> allocate shapes
    allocate(dims(rank)); allocate(maxdims(rank));
    !*** get space size ***
    call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
  end subroutine HDF5_extract_dataset_rank_shape

  !----------------------------------------
  ! Open or create HDF5 file
  ! Takes a property list as optional argument, which can be used to set parameters
  ! for the file opening or access, such as block size, offsets, buffering,
  ! parallel IO and data alignment.
  ! See https://support.hdfgroup.org/HDF5/doc1.6/Files.html for more information.
  !----------------------------------------
  subroutine HDF5_open_or_create(filename,plist,file_id,ierr,file_access)
    use mpi
    character(LEN=*) , intent(in)  :: filename  !< file name
    integer(HID_T)   , intent(in), optional :: plist     !< Which features to use when opening
    integer(HID_T)   , intent(out) :: file_id   !< file identifier
    integer          , intent(in), optional :: file_access
    integer, optional, intent(out) :: ierr

    integer        :: ierr_HDF5
    integer        :: access_f
    logical        :: file_exists, is_hdf5

    !*** Initialize fortran interface ***
    access_f = H5F_ACC_EXCL_F; if(present(file_access)) access_f = file_access
    call H5open_f(ierr_HDF5)
    !*** Test if the file exists ***
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists.and.(access_f.ne.H5F_ACC_TRUNC_F)) then
      !*** Test if it is an HDF5 file ***!
      call H5Fis_hdf5_f(trim(filename)//char(0), is_hdf5, ierr_HDF5)
      if (is_hdf5) then
        !*** Try to open the HDF5 file ***!
        call H5Fopen_f(trim(filename)//char(0), &
          H5F_ACC_RDWR_F,file_id,ierr_HDF5, access_prp=plist)
        if (present(ierr)) ierr = ierr_HDF5
      else
        !*** Present an error ***!
        write(*,*) "ERROR: Invalid HDF5 file, exiting"
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
      end if
    else
      !*** Try to create an HDF5 file ***
      call H5Fcreate_f(trim(filename)//char(0), &
        access_f, file_id, ierr_HDF5, access_prp=plist)
      if (present(ierr)) ierr = ierr_HDF5
    end if
  end subroutine HDF5_open_or_create

  !*************************************************
  !  HDF5 WRITING
  !*************************************************
  !---------------------------------------- 
  ! HDF5 saving for a character string
  !----------------------------------------
  subroutine HDF5_char_saving(file_id,charvar,dsetname)
    use H5LT
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    character(LEN=*), intent(in) :: charvar
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: ierr_HDF5  ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: type_id  ! dataspace identifier

   !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    !*** Create character dataset ***
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,ierr_HDF5)
    CALL h5tset_size_f(type_id,int(len(charvar),SIZE_T),ierr_HDF5)
    CALL h5screate_f(H5S_SCALAR_F,dataspace,ierr_HDF5)
    CALL h5dcreate_f(file_id,dsetname,type_id,dataspace,dataset,ierr_HDF5)
    CALL h5dwrite_f(dataset,type_id,charvar,dim,ierr_HDF5)
    !*** Closing ***
    call H5Dclose_f(dataset,ierr_HDF5)
    call H5Sclose_f(dataspace,ierr_HDF5)
  end subroutine HDF5_char_saving


  !---------------------------------------- 
  ! HDF5 saving for a character 1D array
  !----------------------------------------
  subroutine HDF5_array1D_saving_char(file_id,array1D,dim1,dsetname)
    use H5LT
    integer(HID_T)                , intent(in) :: file_id   ! file identifier
    character(LEN=*), dimension(:), intent(in) :: array1D
    integer                       , intent(in) :: dim1
    character(LEN=*)              , intent(in) :: dsetname  ! dataset name

    integer             :: ierr_HDF5  ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,ierr_HDF5)

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_CHARACTER,&
      dataspace,dataset,ierr_HDF5)

    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    call H5Dwrite_f(dataset,H5T_NATIVE_CHARACTER,array1D,dim,ierr_HDF5)

    !*** Closing ***
    call H5Sclose_f(dataspace,ierr_HDF5)
    call H5Dclose_f(dataset,ierr_HDF5)
  end subroutine HDF5_array1D_saving_char


  !---------------------------------------- 
  ! HDF5 saving for an integer 
  !----------------------------------------
  subroutine HDF5_integer_saving(file_id,int,dsetname)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(in) :: int
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,trim(dsetname), &
      H5T_NATIVE_INTEGER,dataspace,dataset,error)

    !*** Write the integer data to the dataset ***
    !***  using default transfer properties    ***
    call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)

    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_saving


  !---------------------------------------- 
  ! HDF5 saving for a real double
  !----------------------------------------
  subroutine HDF5_real_saving(file_id,rd,dsetname)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    real*8          , intent(in) :: rd
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = 1
    rank   = 1
    call H5Screate_simple_f(rank,dim,dataspace,error)

    !*** Create integer dataset ***
    call H5Dcreate_f(file_id,dsetname, &
      H5T_NATIVE_DOUBLE,dataspace,dataset,error)

    !*** Write the integer data to the dataset ***
    !***  using default transfer properties    ***
    call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)

    !*** Closing ***
    call H5Sclose_f(dataspace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_saving


  !---------------------------------------- 
  ! HDF5 saving for a 1D array of integer
  !----------------------------------------
  subroutine HDF5_array1D_saving_int(file_id,array1D,dim1,dsetname,start,compress_level)
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    integer, dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Begin position of data
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      filespace,dataset,error,property)

    !*** Write the integer array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,shape(array1D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_int

  !---------------------------------------- 
  ! HDF5 saving for a 2D array of integer
  !----------------------------------------
  subroutine HDF5_array2D_saving_int(file_id,array2D,dim1,dim2,dsetname,start,compress_level)
    integer(HID_T)           , intent(in) :: file_id   ! file identifier
    integer, dimension(:,:)  , intent(in) :: array2D
    integer                  , intent(in) :: dim1, dim2
    character(LEN=*)         , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T)         , intent(in), optional :: start(2) !< Begin position of data
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER,&
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(2,shape(array2D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array2D,shape(array2D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_saving_int

  !---------------------------------------- 
  ! gzip HDF5 saving for a 3D array integer
  !----------------------------------------
  subroutine HDF5_array3D_saving_int(file_id,array3D, &
    dim1,dim2,dim3,dsetname,start,compress_level)
    integer(HID_T)          , intent(in) :: file_id   ! file identifier
    integer, dimension(:,:,:), intent(in) :: array3D
    integer                 , intent(in) :: dim1, dim2, dim3
    character(LEN=*)        , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(3), intent(in), optional :: start !< Offset of array to write
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    rank   = 3
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_INTEGER, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(3,shape(array3D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_INTEGER,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_saving_int


  !---------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*4
  !----------------------------------------
  subroutine HDF5_array1D_saving_r4(file_id,array1D,dim1,dsetname,start,compress_level)
    integer(HID_T)       , intent(in) :: file_id   ! file identifier
    real(4), dimension(:), intent(in) :: array1D
    integer              , intent(in) :: dim1
    character(LEN=*)     , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T)     , intent(in), optional :: start(1) !< Begin position of data
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_REAL,&
      filespace,dataset,error,property)

    !*** Write the real*4 array data to the dataset using ***
    !***  default transfer properties                     ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_REAL,array1D,shape(array1D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_REAL,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving_r4


  !-------------------------------------------- 
  ! gzip HDF5 saving for a 1D array of real*8
  !--------------------------------------------
  subroutine HDF5_array1D_saving(file_id,array1D,dim1,dsetname,start,compress_level)
    integer(HID_T)      , intent(in) :: file_id   ! file identifier
    real*8, dimension(:), intent(in) :: array1D
    integer             , intent(in) :: dim1
    character(LEN=*)    , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T)     , intent(in), optional :: start(1) !< Begin position of data
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 
    
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    rank   = 1
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE,&
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(1,shape(array1D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,shape(array1D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Sclose_f(filespace,error)
    call H5Pclose_f(property,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 2D array
  !----------------------------------------
  subroutine HDF5_array2D_saving(file_id,array2D,dim1,dim2,dsetname,start,compress_level)
    integer(HID_T)        , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:), intent(in) :: array2D
    integer               , intent(in) :: dim1, dim2
    character(LEN=*)      , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(2), intent(in), optional :: start !< Offset of array to write
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer              :: error      ! error flag
    integer              :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)       :: dim        ! dataset dimensions
    integer(HID_T)       :: dataset    ! dataset identifier
    integer(HID_T)       :: dataspace  ! dataspace identifier
    integer(HID_T)       :: filespace  ! dataspace identifier
    integer(HID_T)       :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    rank   = 2
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(2,shape(array2D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array2D,shape(array2D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 3D array
  !----------------------------------------
  subroutine HDF5_array3D_saving(file_id,array3D, &
    dim1,dim2,dim3,dsetname,start,compress_level)
    integer(HID_T)          , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:,:), intent(in) :: array3D
    integer                 , intent(in) :: dim1, dim2, dim3
    character(LEN=*)        , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(3), intent(in), optional :: start !< Offset of array to write
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    rank   = 3
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***   using default transfer properties        ***
    if (present(start)) then
      call H5Screate_simple_f(3,shape(array3D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 4D array
  !----------------------------------------
  subroutine HDF5_array4D_saving(file_id,array4d, &
    dim1,dim2,dim3,dim4,dsetname,start,compress_level)
    integer(HID_T)            , intent(in) :: file_id   ! file identifier
    real*8, dimension(:,:,:,:), intent(in) :: array4d
    integer                   , intent(in) :: dim1, dim2, dim3, dim4
    character(LEN=*)          , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(4), intent(in), optional :: start !< Offset of array to write
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier
 
    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    rank   = 4
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties ***
    if (present(start)) then
      call H5Screate_simple_f(4,shape(array4D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array4D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array4D,shape(array4D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array4D_saving


  !---------------------------------------- 
  ! gzip HDF5 saving for a 5D array
  !----------------------------------------
  subroutine HDF5_array5D_saving(file_id,array5d, &
    dim1,dim2,dim3,dim4,dim5,dsetname,start,compress_level)
    integer(HID_T)              , intent(in) :: file_id  ! file identifier
    real*8, dimension(:,:,:,:,:), intent(in) :: array5d
    integer                     , intent(in) :: dim1, dim2
    integer                     , intent(in) :: dim3, dim4, dim5
    character(LEN=*)            , intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(5), intent(in), optional :: start !< Offset of array to write
    integer, intent(in), optional :: compress_level !< if set and start is not provided compress with this level

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(5)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier
    integer(HID_T)      :: property   ! Property list identifier 

    !*** Create and initialize dataspaces for datasets ***
    dim(1) = dim1
    dim(2) = dim2
    dim(3) = dim3
    dim(4) = dim4
    dim(5) = dim5
    rank   = 5
    call H5Screate_simple_f(rank,dim,filespace,error)

    !*** Get compression property list ***
    if (present(compress_level)) then
      property = get_HDF5_plist(rank, dim, .not. present(start), compress_level)
    else
      property = get_HDF5_plist(rank, dim, .false.)
    end if

    !*** Create real dataset ***
    call H5Dcreate_f(file_id,trim(dsetname),H5T_NATIVE_DOUBLE, &
      filespace,dataset,error,property)

    !*** Write the real*8 array data to the dataset ***
    !***  using default transfer properties         ***
    if (present(start)) then
      call H5Screate_simple_f(5,shape(array5D,kind=HSIZE_T),dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array5D,kind=HSIZE_T), hdferr=error)
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array5D,shape(array5D,kind=HSIZE_T),error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dwrite_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)
    end if

    !*** Closing ***
    call H5Pclose_f(property,error)
    call H5Sclose_f(filespace,error)
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_saving


  !************************************************
  !  HDF5 READING
  !************************************************
  !---------------------------------------- 
  ! HDF5 reading for a character
  !----------------------------------------
  subroutine HDF5_char_reading(file_id,charvar,dsetname)
    use H5LT
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    character(LEN=*), intent(out) :: charvar
    character(LEN=*), intent(in)  :: dsetname  ! dataset name
    integer             :: error      ! error flag
    call H5LTread_dataset_string_f(file_id,trim(dsetname),charvar,error)
  end subroutine HDF5_char_reading


  !---------------------------------------- 
  ! HDF5 reading for a character array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading_char(file_id,array1D,dsetname)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    character       , intent(out), dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    call H5Dread_f(dataset,H5T_NATIVE_CHARACTER,array1D,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_char


  !---------------------------------------- 
  ! HDF5 reading for an integer 
  !----------------------------------------
  subroutine HDF5_integer_reading(file_id,int,dsetname)
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    integer         , intent(out) :: int
    character(LEN=*), intent(in)  :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_INTEGER,int,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_integer_reading


  !---------------------------------------- 
  ! HDF5 reading for an integer array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading_int(file_id,array1D,dsetname,start)
    integer(HID_T)  , intent(in) :: file_id   ! file identifier
    integer         , intent(out), dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_int


  !---------------------------------------- 
  ! HDF5 reading for a real double
  !----------------------------------------
  subroutine HDF5_real_reading(file_id,rd,dsetname)
    integer(HID_T)  , intent(in)  :: file_id   ! file identifier
    real*8          , intent(out) :: rd
    character(LEN=*), intent(in)  :: dsetname  ! dataset name

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: data_type

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)   

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,rd,dim,error)

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_real_reading


  !---------------------------------------- 
  ! HDF5 reading for an array 1D
  !----------------------------------------
  subroutine HDF5_array1D_reading(file_id,array1D,dsetname,start)
    integer(HID_T), intent(in)   :: file_id   ! file identifier
    real*8        , dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 1D (real*4)
  !----------------------------------------
  subroutine HDF5_array1D_reading_r4(file_id,array1D,dsetname,start)
    integer(HID_T), intent(in)   :: file_id   ! file identifier
    real*4        , dimension(:) :: array1D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(1), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(1)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***   using default transfer properties  ***
    dim(1) = size(array1D,1)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(1,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_REAL,array1D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_REAL,array1D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array1D_reading_r4

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 1D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array1D_reading(file_id,array1D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(1),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(1),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    real*8,dimension(:),allocatable,intent(inout) :: array1D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(1) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array1D)) deallocate(array1D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.1) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array1D(reqdims(1))); array1D = 0d0;
        call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array1D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array1D(dims(1))); array1D = 0d0;
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array1D,dims,error)
      endif
    else
      write(*,*) 'Read 1D allocatable array from HDF5, rank is not 1: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array1D_reading

  !---------------------------------------- 
  ! HDF5 reading for an integer allocatable array 1D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array1D_reading_int(file_id,array1D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(1),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(1),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    integer,dimension(:),allocatable,intent(inout) :: array1D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(1) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array1D)) deallocate(array1D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.1) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array1D(reqdims(1))); array1D = 0;
        call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array1D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_INTEGER,array1D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array1D(dims(1))); array1D = 0;
        call H5Dread_f(dataset_id,H5T_NATIVE_INTEGER,array1D,dims,error)
      endif
    else
      write(*,*) 'Read 1D allocatable integer array from HDF5, rank is not 1: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array1D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for an integer array 2D 
  !----------------------------------------
  subroutine HDF5_array2D_reading_int(file_id,array2D,dsetname,start)
    integer(HID_T)  , intent(in)     :: file_id   ! file identifier
    integer         , dimension(:,:) :: array2D
    character(LEN=*), intent(in)     :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(2), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array2D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(2,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_INTEGER,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for an allocatable integer array 2D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array2D_reading_int(file_id,array2D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(2),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(2),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    integer,dimension(:,:),allocatable,intent(inout) :: array2D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(2) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array2D)) deallocate(array2D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.2) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array2D(reqdims(1),reqdims(2))); array2D = 0;
        call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_INTEGER,array2D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array2D(dims(1),dims(2))); array2D = 0;
        call H5Dread_f(dataset_id,H5T_NATIVE_INTEGER,array2D,dims,error)
      endif
    else
      write(*,*) 'Read 2D integer allocatable array from HDF5, rank is not 2: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array2D_reading_int

  !---------------------------------------- 
  ! HDF5 reading for an array 2D
  !----------------------------------------
  subroutine HDF5_array2D_reading(file_id,array2D,dsetname,start)
    integer(HID_T)  , intent(in)     :: file_id   ! file identifier
    real*8          , dimension(:,:) :: array2D
    character(LEN=*), intent(in)     :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(2), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(2)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array2D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(2,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array2D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array2D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 2D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array2D_reading(file_id,array2D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(2),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(2),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    real*8,dimension(:,:),allocatable,intent(inout) :: array2D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(2) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array2D)) deallocate(array2D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.2) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array2D(reqdims(1),reqdims(2))); array2D = 0d0;
        call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array2D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array2D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array2D(dims(1),dims(2))); array2D = 0d0;
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array2D,dims,error)
      endif
    else
      write(*,*) 'Read 2D allocatable array from HDF5, rank is not 2: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array2D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 3D
  !----------------------------------------
  subroutine HDF5_array3D_reading(file_id,array3D,dsetname,&
       in1, in2, in3,start)
    integer(HID_T)     , intent(in)       :: file_id   ! file identifier
    real*8             , dimension(:,:,:) :: array3D
    character(LEN=*)   , intent(in)       :: dsetname  ! dataset name
    integer, intent(in), optional  :: in1, in2, in3
    integer(HSIZE_T), dimension(3), intent(in), optional :: start !< Offset of array to read
    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
      dimension(3)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    if (present(in1) .and. present(in2) .and. present(in3)) then
      dim = [in1, in2, in3]
    else
      dim = shape(array3D)
    end if
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(3,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array3D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array3D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 3D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array3D_reading(file_id,array3D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(3),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(3),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    real*8,dimension(:,:,:),allocatable,intent(inout) :: array3D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(3) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array3D)) deallocate(array3D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.3) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array3D(reqdims(1),reqdims(2),reqdims(3))); array3D = 0d0;
        call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array3D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array3D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array3D(dims(1),dims(2),dims(3))); array3D = 0d0;
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array3D,dims,error)
      endif
    else
      write(*,*) 'Read 3D allocatable array from HDF5, rank is not 3: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array3D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 4D
  !----------------------------------------
  subroutine HDF5_array4D_reading(file_id,array4D,dsetname,ierr,start)
    integer(HID_T), intent(in)     :: file_id   ! file identifier
    real*8, dimension(:,:,:,:)     :: array4D
    character(LEN=*), intent(in)   :: dsetname  ! dataset name
    integer, optional, intent(out) :: ierr
    integer(HSIZE_T), dimension(4), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank           ! dataset rank
    integer(HSIZE_T), &
      dimension(4)      :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)
   
    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array4D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(4,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array4D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array4D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
    if (present(ierr)) ierr = error
  end subroutine HDF5_array4D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 4D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array4D_reading(file_id,array4D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(4),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(4),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    real*8,dimension(:,:,:,:),allocatable,intent(inout) :: array4D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(4) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array4D)) deallocate(array4D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.4) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array4D(reqdims(1),reqdims(2),reqdims(3),reqdims(4))); array4D = 0d0;
        call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array4D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array4D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array4D(dims(1),dims(2),dims(3),dims(4))); array4D = 0d0;
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array4D,dims,error)
      endif
    else
      write(*,*) 'Read 4D allocatable array from HDF5, rank is not 4: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array4D_reading

  !---------------------------------------- 
  ! HDF5 reading for an array 5D
  !----------------------------------------
  subroutine HDF5_array5D_reading(file_id,array5D,dsetname,start)
    integer(HID_T), intent(in)   :: file_id   
    real*8, dimension(:,:,:,:,:) :: array5D
    character(LEN=*), intent(in) :: dsetname  ! dataset name
    integer(HSIZE_T), dimension(4), intent(in), optional :: start !< Offset of array to read

    integer             :: error      ! error flag
    integer             :: rank       ! dataset rank
    integer(HSIZE_T), &
           dimension(5) :: dim        ! dataset dimensions
    integer(HID_T)      :: dataset    ! dataset identifier
    integer(HID_T)      :: dataspace  ! dataspace identifier
    integer(HID_T)      :: filespace  ! dataspace identifier

    !*** file opening ***
    call H5Dopen_f(file_id,trim(dsetname),dataset,error)

    !*** read the integer data to the dataset ***
    !***  using default transfer properties   ***
    dim = shape(array5D)
    if (present(start)) then
      call H5Dget_space_f(dataset,filespace,error)
      call H5Screate_simple_f(5,dim,dataspace,error)
      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
      start=start, count=shape(array5D,kind=HSIZE_T), hdferr=error)
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error, &
          file_space_id=filespace, mem_space_id=dataspace)
      call H5Sclose_f(filespace,error)
      call H5Sclose_f(dataspace,error)
    else
      call H5Dread_f(dataset,H5T_NATIVE_DOUBLE,array5D,dim,error)
    end if

    !*** Closing ***
    call H5Dclose_f(dataset,error)
  end subroutine HDF5_array5D_reading

  !---------------------------------------- 
  ! HDF5 reading for an allocatable array 5D
  !---------------------------------------- 
  subroutine HDF5_allocatable_array5D_reading(file_id,array5D,dsetname,reqdims,start)
    !> inputs:
    integer(HID_T),intent(in)                         :: file_id  !< file identifier
    character(len=*),intent(in)                       :: dsetname !< dataset name
    integer(HSIZE_T),dimension(5),intent(in),optional :: start    !< offset array to read
    integer(HSIZE_T),dimension(5),intent(in),optional :: reqdims  !< requested slab dimensions
    !> inputs-outputs:
    real*8,dimension(:,:,:,:,:),allocatable,intent(inout) :: array5D
    !> variables:
    integer(HID_T) :: dataspace_id,dataspace_req_id !< dataspace identifier
    integer(HID_T) :: dataset_id   !< dataset identifier
    integer        :: rank,error   !< dataset rank and hdf5 error
    integer(HSIZE_T),dimension(5) :: dims,maxdims !< dataset dimensions

    !*** initialisation ***
    dims = 0; if(allocated(array5D)) deallocate(array5D)
    !*** get the dataset and dataspace ids ***
    call H5Dopen_f(file_id,trim(dsetname),dataset_id,error)
    call H5dget_space_f(dataset_id,dataspace_id,error)
    !*** get the space rank ***
    call H5Sget_simple_extent_ndims_f(dataspace_id,rank,error)
    if(rank.eq.5) then
      !*** get space size ***
      call H5Sget_simple_extent_dims_f(dataspace_id,dims,maxdims,error)
      if(present(reqdims).and.present(start)) then
        !*** load a data slab of size reqdims and start start ***
        !*** using default transfer properties ***
        allocate(array5D(reqdims(1),reqdims(2),reqdims(3),reqdims(4),reqdims(5))); 
        array5D = 0d0; call H5Screate_simple_f(1,reqdims,dataspace_req_id,error)
        call H5Sselect_hyperslab_f(dataset_id, H5S_SELECT_SET_F, &
        start=start, count=shape(array5D,kind=HSIZE_T), hdferr=error)
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array5D,reqdims,error, &
            file_space_id=dataspace_id, mem_space_id=dataspace_req_id)
        call H5Sclose_f(dataspace_req_id,error)       
      else
        !*** allocate and read read the dataset ***
        !*** using default transfer properties ***
        allocate(array5D(dims(1),dims(2),dims(3),dims(4),dims(5))); array5D = 0d0;
        call H5Dread_f(dataset_id,H5T_NATIVE_DOUBLE,array5D,dims,error)
      endif
    else
      write(*,*) 'Read 5D allocatable array from HDF5, rank is not 5: skip!'
    endif
    !*** Closing ***
    call H5Sclose_f(dataspace_id,error)
    call H5Dclose_f(dataset_id,error)
  end subroutine HDF5_allocatable_array5D_reading

  function get_HDF5_plist(rank, chunksize, gzip, level) result(property)
    integer, intent(in) :: rank
    integer(HSIZE_T), dimension(rank), intent(in) :: chunksize
    logical, intent(in) :: gzip
    integer, intent(in), optional :: level
    integer(HID_T) :: property

    integer        :: error, i
    integer(HSIZE_T), dimension(rank) :: chk
    integer, parameter :: cmpr_default = 6
    integer :: cmpr_level

    ! Check chunksize, if too large make it smaller.
    ! Aim for 10000 elements per chunk, by reducing the last dimension first
    chk = chunksize
    i   = rank
    do while (product(chk) .gt. 100000)
      if (chk(i) .gt. 1) then
        chk(i) = chk(i) / 2
      else
        i = i-1
      end if
    end do

    !*** Creates a new property dataset ***
    call H5Pcreate_f(H5P_DATASET_CREATE_F,property,error)

    ! If we are not doing parallel IO we can enable compression
    if (gzip) then
      if (present(level)) then
        cmpr_level = level
      else
        cmpr_level = cmpr_default
      end if
      if (cmpr_level .gt. 0) then
        call H5Pset_chunk_f(property,rank,chk,error)
        call H5Pset_deflate_f(property,cmpr_level,error)
      end if
    end if
  end function get_HDF5_plist

#endif
end module hdf5_io_module
