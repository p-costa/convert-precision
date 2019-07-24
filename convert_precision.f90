program convert_precision
  use mpi
  implicit none
  real, parameter :: small = 1.e-5
  integer, parameter :: r_in  = 8
  integer, parameter :: r_out = 4
  character(len=*), parameter :: out_ext = '.converted'
  character(len=*), parameter :: input_file = 'file.in'
  !
  real(kind=r_in ), allocatable, dimension(:) :: data_in
  real(kind=r_out), allocatable, dimension(:) :: data_out
  !
  integer :: myid,nproc,ierr
  integer :: iunit,istatus
  character(len=100) :: fname
  integer(kind=MPI_OFFSET_KIND) :: filesize,disp
  integer :: MPI_REAL_R_IN, MPI_REAL_R_OUT
  integer :: fh
  integer :: nreals,nreals_myid,disp_myid
  integer :: i,error
  !
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  if(r_in.eq.8.and.r_out.eq.4) then
    MPI_REAL_R_IN  = MPI_DOUBLE_PRECISION
    MPI_REAL_R_OUT = MPI_REAL
  elseif(r_in.eq.4.and.r_out.eq.8) then
    MPI_REAL_R_IN  = MPI_REAL
    MPI_REAL_R_OUT = MPI_DOUBLE_PRECISION
  else
    if(myid.eq.0) print*, 'Error, invalid precision of input/output reals.'
    if(myid.eq.0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
  endif
  !
  istatus = 0
  open(newunit=iunit, file = 'files.in',iostat=istatus)
  if(istatus.ne.0) then
    if(myid.eq.0) print*, 'Error while reading file ',trim(input_file),'.'
    if(myid.eq.0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
  endif
  !
  do while(istatus.eq.0)
    read(iunit,'(A)',iostat=istatus) fname
    if(istatus.ne.0) exit
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MPI_FILE_GET_SIZE(fh,filesize,ierr)
    if(mod(filesize,r_in).ne.0) then
      if(myid.eq.0) print*, 'Error while reading file ',fname,' with size ',filesize,' not divisable by ',r_in,'.'
      if(myid.eq.0) print*, 'Aborting...'
      istatus = 1000
    endif
    if(istatus.ne.0) exit
    ! divide reals evenly
    nreals = filesize/r_in
    if(myid+1.le.mod(nreals,nproc)) then
      nreals_myid = nreals/nproc+1
    else
      nreals_myid = nreals/nproc
    endif
    allocate(data_in(nreals_myid))
    !
    ! calculate displacements
    !
    if(myid.eq.0) then
      disp_myid = 0
    else
      if(myid+1-1.le.mod(nreals,nproc)) then
        disp_myid = nreals/nproc+1
      else
        disp_myid = nreals/nproc
      endif
    endif
    call MPI_SCAN(disp_myid,disp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_FILE_SET_VIEW(fh,disp*r_in,MPI_REAL_R_IN,MPI_REAL_R_IN,'native',MPI_INFO_NULL,ierr)
    call MPI_FILE_READ(fh,data_in,nreals_myid,MPI_REAL_R_IN,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    allocate(data_out(nreals_myid))
    do i=1,nreals_myid
      data_out(i) = data_in(i)
    enddo
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname)//trim(out_ext),MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)  ! guarantee overwriting
    call MPI_FILE_SET_VIEW(fh,disp*r_out, MPI_REAL_R_OUT,MPI_REAL_R_OUT, 'native', &
         MPI_INFO_NULL, ierr)
    call MPI_FILE_WRITE(fh,data_out,nreals_myid,MPI_REAL_R_OUT,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    !
    ! final step: check converted data
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD,trim(fname)//trim(out_ext),MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    call MPI_FILE_SET_VIEW(fh,disp*r_out,MPI_REAL_R_OUT,MPI_REAL_R_OUT,'native',MPI_INFO_NULL,ierr)
    call MPI_FILE_READ(fh,data_out,nreals_myid,MPI_REAL_R_OUT,MPI_STATUS_IGNORE,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    error = 0
    do i=1,nreals_myid
      if(abs(real(data_out(i),kind=r_out)-real(data_in(i),kind=r_out)).gt.small) error = error + 1
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE,error,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate(data_in,data_out)
    if(error.gt.0) then
      if(myid.eq.0) print*, 'Error, the data read form the converted file does not match the input data.'
      if(myid.eq.0) print*, 'Aborting...'
      exit
    endif
    if(myid.eq.0) print*, 'File ',trim(fname),' successfully converted from ',r_in, ' to ', r_out, ' bytes.'
    if(myid.eq.0) print*, 'New file name: ', trim(fname)//trim(out_ext),' .'
  enddo
  call MPI_FINALIZE(ierr)
end program convert_precision
