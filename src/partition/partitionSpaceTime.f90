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
  debug = 1
  !
  call mpi_init(ierr)
  call mpi_comm_group(mpi_comm_world,globalGroup,ierr)
  call mpi_comm_rank (mpi_comm_world,myid,       ierr)
  call mpi_comm_size (mpi_comm_world,numprocs,   ierr)
  !
  numprocs_temporal = 4
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
  endif
  write(*,*) "myid: ",myid," mpi_comm_world: ", mpi_comm_world
  write(*,*) "myid: ",myid," myid_t: ",myid_temporal," timeGroup: ", timeGroup
  write(*,*) "myid: ",myid," myid_s: ",myid_temporal," spaceGroup: ", spaceGroup
  if(debug.gt.0) write(*,*) myid,myid_spatial,myid_temporal     
  !
  !> Test Communicators
  !
  if (myid.lt.numprocs_active.and.debug.gt.0) then
     call testSpaceTimeComm(myid,myid_spatial,myid_temporal,spaceElem,timeElem,spaceComm,timeComm,spaceRanks,timeRanks)
  endif
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
!
subroutine testSpaceTimeComm(myid,myid_spatial,myid_temporal,spaceElem,timeElem,spaceComm,timeComm,spaceRanks,timeRanks)
  !
  implicit none
  !
  include 'mpif.h'
  !
  integer,                       intent(in) :: myid,myid_spatial,myid_temporal
  integer,                       intent(in) :: spaceElem,spaceComm
  integer,                       intent(in) :: timeElem,timeComm
  integer, dimension(spaceElem), intent(in) :: spaceRanks
  integer, dimension(timeElem),  intent(in) :: timeRanks
  !
  integer :: s1,s2,t1,t2,sumSpace,sumTime
  integer :: ierr,i
  !
  !> Test Spatial Communicator
  !
  s1 = myid
  !
  sumSpace = sum(spaceRanks)
  call mpi_reduce(s1,s2,1,MPI_INTEGER,MPI_SUM,0,spaceComm,ierr)
  if (myid_spatial.eq.0) then
     if (sumSpace.eq.s2) then
        write(*,*) "Sucessful MPI_REDUCE(s) for myid ",myid,sumSpace,s2
     else
        write(*,*) "Unsucessful MPI_REDUCE(s) for myid ",myid,sumSpace,s2
     endif
  endif
  !
  !> Test Temporal Communicator
  !
  t1 = myid
  !                                                                                                                                                                       
  sumTime = sum(timeRanks)
  call mpi_reduce(t1,t2,1,MPI_INTEGER,MPI_SUM,0,timeComm,ierr)
  if (myid_temporal.eq.0) then
     if (sumTime.eq.t2) then
        write(*,*) "Sucessful MPI_REDUCE(t) for myid ",myid,sumTime,t2
     else
        write(*,*) "Unsucessful MPI_REDUCE(t) for myid ",myid,sumTime,t2
     endif
  endif
  !
  return
end subroutine testSpaceTimeComm
!
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
