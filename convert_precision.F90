program convert_precision
  use mpi_f08
  implicit none
  real, parameter :: small = 1.e-5
  integer, parameter :: sp = selected_real_kind(6 , 37)
  integer, parameter :: dp = selected_real_kind(15,307)
#ifdef _SINGLE_TO_DOUBLE
  integer, parameter :: r_in  = sp
  integer, parameter :: r_out = dp
#elif  _DOUBLE_TO_SINGLE
  integer, parameter :: r_in  = dp
  integer, parameter :: r_out = sp
#else /* default double -> single */
  integer, parameter :: r_in  = dp
  integer, parameter :: r_out = sp
#endif
  character(len=*), parameter :: out_ext = '.converted'
  character(len=*), parameter :: input_file = 'files.in'
  !
  real(kind=r_in ), allocatable, dimension(:) :: data_in
  real(kind=r_out), allocatable, dimension(:) :: data_out
  !
  integer :: myid,nproc,ierr
  integer :: iunit,istatus
  character(len=1024) :: fname
  integer(kind=MPI_OFFSET_KIND) :: filesize,disp,nreals,nreals_myid
  type(MPI_DATATYPE) :: MPI_REAL_R_IN, MPI_REAL_R_OUT
  type(MPI_DATATYPE) :: MPI_INTEGER_DISP
  integer :: isize
  type(MPI_FILE) :: fh
  integer :: i,error
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  if(r_in == dp.and.r_out == sp) then
    MPI_REAL_R_IN  = MPI_DOUBLE_PRECISION
    MPI_REAL_R_OUT = MPI_REAL
    if(myid == 0) print*, 'Double to Single precision conversion'
  elseif(r_in == sp.and.r_out == dp) then
    MPI_REAL_R_IN  = MPI_REAL
    MPI_REAL_R_OUT = MPI_DOUBLE_PRECISION
    if(myid == 0) print*, 'Single to Double precision conversion'
  else
    if(myid == 0) print*, 'Error, invalid precision of input/output reals.'
    if(myid == 0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
  end if
  isize = storage_size(disp)/8
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, isize, MPI_INTEGER_DISP, ierr)
  !
  istatus = 0
  open(newunit=iunit, file = 'files.in',iostat=istatus)
  if(istatus /= 0) then
    if(myid == 0) print*, 'Error while reading file ',trim(input_file),'.'
    if(myid == 0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
  end if
  !
  do while(istatus == 0)
    read(iunit,'(A)',iostat=istatus) fname
    if(istatus /= 0) exit
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
    if(mod(filesize,int(r_in,kind(filesize))) /= 0) then
      if(myid == 0) print*, 'Error while reading file ',fname,' with size ',filesize,' not divisable by ',r_in,'.'
      if(myid == 0) print*, 'Aborting...'
      istatus = 1000
    end if
    if(istatus /= 0) exit
    ! divide reals evenly
    nreals = filesize/r_in
    if(myid+1 <= mod(nreals,int(nproc,kind(nreals)))) then
      nreals_myid = nreals/nproc+1
    else
      nreals_myid = nreals/nproc
    end if
    allocate(data_in(nreals_myid))
    !
    ! calculate displacements
    !
    if(myid == 0) then
      disp = 0
    else
      if(myid+1-1 <= mod(nreals,int(nproc,kind(nreals)))) then
        disp = nreals/nproc+1
      else
        disp = nreals/nproc
      end if
    end if
    call MPI_SCAN(MPI_IN_PLACE,disp,1,MPI_INTEGER_DISP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_FILE_SET_VIEW(fh,disp*r_in,MPI_REAL_R_IN,MPI_REAL_R_IN,'native',MPI_INFO_NULL,ierr)
    call MPI_FILE_READ(fh,data_in,int(nreals_myid),MPI_REAL_R_IN,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    allocate(data_out(nreals_myid))
    do i=1,nreals_myid
      data_out(i) = data_in(i)
    end do
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname)//trim(out_ext),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*r_out, MPI_REAL_R_OUT,MPI_REAL_R_OUT, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,data_out,int(nreals_myid),MPI_REAL_R_OUT,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    ! final step: check converted data
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname)//trim(out_ext),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
    if(filesize /= nreals*r_out) then
      if(myid == 0) print*, 'Error, unexpected size of the output file.'
      if(myid == 0) print*, 'Expected size: ', nreals*r_out, 'Actual size: ', filesize
      if(myid == 0) print*, 'Aborting...'
      exit
    end if
    call MPI_FILE_SET_VIEW(fh,disp*r_out,MPI_REAL_R_OUT,MPI_REAL_R_OUT,'native',MPI_INFO_NULL,ierr)
    call MPI_FILE_READ(fh,data_out,int(nreals_myid),MPI_REAL_R_OUT,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    error = 0
    do i=1,nreals_myid
      if(abs(real(data_out(i),kind=r_out)-real(data_in(i),kind=r_out)) > small) error = error + 1
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,error,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate(data_in,data_out)
    if(error > 0) then
      if(myid == 0) print*, 'Error, the data read form the converted file does not match the input data.'
      if(myid == 0) print*, 'Aborting...'
      exit
    end if
    if(myid == 0) print*, 'File ',trim(fname),' successfully converted from ',r_in, ' to ', r_out, ' bytes per element.'
    if(myid == 0) print*, 'New file name: ', trim(fname)//trim(out_ext),' .'
  end do
  call MPI_FINALIZE(ierr)
end program convert_precision
