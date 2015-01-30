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
  integer :: comm1,comm2
  integer :: n1,n2
  integer :: topologySpaceTime
  integer :: myid_1,myid_2
  integer :: elem1,elem2
  integer :: group1,group2
  !
  integer, allocatable :: spaceRanks(:),timeRanks(:)
  integer, allocatable :: ranks1(:),ranks2(:)
  !
  topologySpaceTime = 1
  debug = 1
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
  n1 = numprocs_spatial
  n2 = numprocs_temporal
  if (topologySpaceTime.eq.1) then
     n1 = numprocs_temporal
     n2 = numprocs_spatial
  endif

  allocate(spaceRanks(numprocs_spatial))
  allocate(timeRanks(numprocs_temporal))
  allocate(ranks1(n1))
  allocate(ranks2(n2))
  !
  !> Generate Groups for Space & Time
  !
  call getSpaceTimeRanks(myid,n1,n2,numprocs_active,ranks1,ranks2,elem1,elem2)
  call mpi_group_incl(globalGroup,elem1,ranks1(1:elem1),group1,ierr)
  call mpi_group_incl(globalGroup,elem2,ranks2(1:elem2),group2, ierr)  
  !
  !> Generate Communicators for Space & Time
  !
  ! Direction 1
  call mpi_comm_create(mpi_comm_world,group1,comm1,ierr)
  call mpi_comm_rank(comm1,myid_1,ierr)
  call mpi_comm_size(comm1,n1,ierr) 

  ! Direction 2
  call mpi_comm_create(mpi_comm_world,group2,comm2,ierr)
  call mpi_comm_rank(comm2,myid_2,ierr)
  call mpi_comm_size(comm2,n2,ierr)
  !
  !> Assign spatial and temporal communicators 
  !
  if (topologySpaceTime.eq.0) then
     myid_spatial  = myid_1
     spaceComm     = comm1
     spaceElem     = elem1
     spaceRanks    = ranks1

     myid_temporal = myid_2
     timeComm      = comm2
     timeElem      = elem2
     timeRanks     = ranks2

  else
     myid_spatial = myid_2
     spaceComm    = comm2
     spaceElem    = elem2
     spaceRanks   = ranks2

     myid_temporal = myid_1
     timeComm      = comm1
     timeElem      = elem1
     timeRanks     = ranks1
  endif
  !
 ! if (myid.eq.6.and.debug.gt.0) then
 !    write(*,*) "myid: ",myid
 !    write(*,*) "myid_temporal: ", myid_temporal
 !    write(*,*) "myid_spatial: ", myid_spatial
 ! endif
!  write(*,*) "myid: ",myid," mpi_comm_world: ", mpi_comm_world
!  write(*,*) "myid: ",myid," myid_t: ",myid_temporal," timeGroup: ", timeGroup
!  write(*,*) "myid: ",myid," myid_s: ",myid_temporal," spaceGroup: ", spaceGroup
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
  deallocate(ranks1,ranks2,spaceRanks,timeRanks)
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
subroutine getSpaceTimeRanks(myid,n1,n2,numprocs_active,ranks1,ranks2,elem1,elem2,topologyST)
  !
  implicit none
  !
  integer,                               intent(in)  :: myid
  integer,                               intent(in)  :: n1
  integer,                               intent(in)  :: n2
  integer,                               intent(in)  :: numprocs_active
  integer, dimension(n1),                intent(out) :: ranks1
  integer, dimension(n2),                intent(out) :: ranks2
  integer,                               intent(out) :: elem1,elem2
  integer,                               intent(in)  :: topologyST
  !
  integer :: i,j,fac1,fac2
  !
  if (myid.lt.numprocs_active) then
     ! 
     !> Compute Node
     !
     elem1 = n1
     fac1 = floor(real(myid)/real(n1))
     do i = 0,n1-1
        ranks1(i+1) = i + fac1*n1
     enddo
     
     elem2 = n2
     fac2 = mod(myid,n1)
     do i = 0,n2-1
        ranks2(i+1) = i*n1 + fac2
     enddo
  else
     !
     !> Non-compute Node
     !
     elem1 = 1
     ranks1(1) = myid

     elem2 = 1
     ranks2(1)  = myid
  endif
  !
  return
end subroutine getSpaceTimeRanks
