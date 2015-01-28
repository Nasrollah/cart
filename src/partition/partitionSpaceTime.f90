program partitionSpaceTime
  !
  implicit none
  !
  include 'mpif.h'
  !
  integer :: numprocs
  integer :: numprocs_active
  integer :: numprocs_temporal
  integer :: numprocs_spatial
  integer :: myid
  integer :: myid_temporal
  integer :: myid_spatial
  integer :: spaceComm
  integer :: timeComm
  integer :: ierr
  integer :: timeGroup,spaceGroup
  integer :: globalGroup
  integer :: spaceElem,timeElem
  integer :: debug
  !
  integer, allocatable :: spaceRanks(:),timeRanks(:)
  !
  debug = 0
  !
  call mpi_init(ierr)
  call mpi_comm_group(mpi_comm_world,globalGroup,ierr)
  call mpi_comm_rank (mpi_comm_world,myid,       ierr)
  call mpi_comm_size (mpi_comm_world,numprocs,   ierr)
  !
  numprocs_temporal = 3
  numprocs_spatial  = floor(real(numprocs)/real(numprocs_temporal))
  numprocs_active   = numprocs_spatial*numprocs_temporal
  !
  if (myid.eq.0.and.debug.gt.0) then
     write(*,*) "numprocs: ", numprocs
     write(*,*) "numprocs_spatial: ", numprocs_spatial
     write(*,*) "numprocs_temporal: ", numprocs_temporal
     write(*,*) "numprocs_active: ", numprocs_active
  endif
  !
  allocate(timeRanks(numprocs_temporal))
  allocate(spaceRanks(numprocs_spatial))
  !
  !> Generate Groups for Space & Time
  !
  call getSpaceTimeRanks(myid,numprocs_spatial,numprocs_temporal,numprocs_active,spaceRanks,timeRanks,spaceElem,timeElem)
  ! Space
  call mpi_group_incl(globalGroup,spaceElem,spaceRanks(1:spaceElem),spaceGroup,ierr)
  ! Time
  call mpi_group_incl(globalGroup,timeElem,  timeRanks(1:timeElem), timeGroup, ierr)  
  !

  !
  !> Generate Communicators for Space & Time
  !
  ! Space
  call mpi_comm_create(mpi_comm_world,spaceGroup,spaceComm,ierr)
  call mpi_comm_rank(spaceComm,myid_spatial,ierr)
  call mpi_comm_size(spaceComm,numprocs_spatial,ierr) 

  ! Time
  call mpi_comm_create(mpi_comm_world,timeGroup,timeComm,ierr)
  call mpi_comm_rank(timeComm,myid_temporal,ierr)
  call mpi_comm_size(timeComm,numprocs_temporal,ierr)
  !
  if (myid.eq.6.and.debug.gt.0) then
     write(*,*) "myid: ",myid
     write(*,*) "myid_temporal: ", myid_temporal
     write(*,*) "myid_spatial: ", myid_spatial
     write(*,*) "numprocs: ", numprocs
     write(*,*) "numprocs_spatial: ", numprocs_spatial
     write(*,*) "numprocs_temporal: ", numprocs_temporal
  endif
  if(debug.gt.0) write(*,*) myid,myid_spatial,myid_temporal     
  !
  !> Deallocate Arrays
  !
  deallocate(spaceRanks,timeRanks)
  !
  !> Finalize MPI
  !
  call mpi_finalize()
  !
end program partitionSpaceTime

subroutine getSpaceTimeRanks(myid,numprocs_spatial,numprocs_temporal,numprocs_active,spaceRanks,timeRanks,spaceElem,timeElem)
  !
  implicit none
  !
  integer,                               intent(in)  :: myid
  integer,                               intent(in)  :: numprocs_spatial
  integer,                               intent(in)  :: numprocs_temporal
  integer,                               intent(in)  :: numprocs_active
  integer, dimension(numprocs_spatial),  intent(out) :: spaceRanks
  integer, dimension(numprocs_temporal), intent(out) :: timeRanks
  integer,                               intent(out) :: spaceElem,timeElem
  !
  integer :: i,j,spaceFac,timeFac
  !
  if (myid.lt.numprocs_active) then
     ! 
     !> Compute Node
     !
     spaceElem = numprocs_spatial
     spaceFac = floor(real(myid)/real(numprocs_spatial))
     do i = 0,numprocs_spatial-1
        spaceRanks(i+1) = i + spaceFac*numprocs_spatial
     enddo
     
     timeElem = numprocs_temporal
     timeFac = mod(myid,numprocs_spatial)
     do i = 0,numprocs_temporal-1
        timeRanks(i+1) = i*numprocs_spatial + timeFac
     enddo
  else
     !
     !> Non-compute Node
     !
     spaceElem = 1
     spaceRanks(1) = myid

     timeElem = 1
     timeRanks(1)  = myid
  endif
  !
  return
end subroutine getSpaceTimeRanks
