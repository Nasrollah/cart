module spatialCommunication
  !
  use mpi
  implicit none
  !
  integer :: cartComm          !< spatial communicator
  integer :: timeComm          !< temporal communicator
  integer :: myid              !< myid in a global (space and time) sense
  integer :: myid_spatial      !< myid in local spatial communicator
  integer :: myid_temporal     !< myid in local temporal communicator
  integer :: numprocs          !< total number of processors
  integer :: numprocs_spatial  !< number of spatial procs for each time instance
  integer :: numprocs_temporal !< number of procs in spatial dimension
  integer :: numprocs_active   !< number of active processors (may be less than numprocs)
  integer :: irnum             !< request number      
  integer :: iprocs(3)         !< total number of processors in the processor map
  integer :: id(3)             !< my location in the processor map
  integer :: gdim(3)           !< global dimensions of the grid
  integer :: globalGroup       !< MPI group containing all ranks
  integer :: cartGroup         !< MPI group containing ranks in local spatial communicator
  integer :: timeGroup         !< MPI group containing ranks in local temporal communicator
  integer :: cartElem          !< number of elements in spatial communicator
  integer :: timeElem          !< number of ranks in temporal communicator
  !
  integer, allocatable :: displ1(:)  !< displacements after partition direction 1
  integer, allocatable :: displ2(:)  !< displacements after partition direction 2   
  integer, allocatable :: displ3(:)  !< displacements after partition direction 3
  !
  integer, allocatable :: pid(:,:,:)      !< processor map
  integer, allocatable :: ireq(:)         !< request handles
  integer, allocatable :: mpistatus(:,:)  !< status handles
  integer :: ierr
  !
  real*8, allocatable :: buffer(:,:)     !< buffer for communication
  !
  integer :: i,ip,ipp,idir,plane,sendid,j,k,l,m,n,iq,iloc,recvid
  integer :: idim(2,3)
  integer, dimension(3) :: stride,dim0,ldim,kdim,ipsend,iprecv
  integer :: qskip,qmult
  integer :: js,je,ks,ke,ls,le,jkmax
  integer :: jstride,dim2,dim1,buffersize,bufcount
  logical :: iskip
  !
  integer :: n1,n2
  integer :: myid_1,myid_2
  integer :: elem1,elem2
  integer :: group1,group2
  integer :: comm1,comm2
  integer :: topologySpaceTime
!  integer,allocatable :: cartRanks(:),timeRanks(:)
  integer, allocatable :: ranks1(:),ranks2(:)
  !
contains
  !
  !> Initialize the spatial communicator
  !
  subroutine initSpatialComm(comm)
    !
    implicit none
    integer, intent(in) :: comm
    !
    cartComm=comm
    !
    call mpi_comm_rank(mpi_comm_world,myid,ierr)
    call mpi_comm_rank(cartComm,myid_spatial,ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)
    call mpi_comm_size(cartComm,numprocs_spatial,ierr)
    !
    if (myid.eq.0) write(*,*) "inside initSpatialComm"
  end subroutine initSpatialComm
  !
  !> Get Proc Id
  !
  subroutine getProcID(id1)
    !
    implicit none
    !
    integer, dimension(3), intent(out) :: id1
    !
    id1 = id
    !
    return
  end subroutine getProcID
  !> release all memory and terminate
  !> the spatial communication
  !
  subroutine terminateSpatialComm
    implicit none
    if (allocated(displ1)) deallocate(displ1)
    if (allocated(displ2)) deallocate(displ2)
    if (allocated(displ3)) deallocate(displ3)
    if (allocated(pid)) deallocate(pid)
    if (allocated(ireq)) deallocate(ireq)
    if (allocated(mpistatus)) deallocate(mpistatus)
    if (allocated(buffer)) deallocate(buffer)
  end subroutine terminateSpatialComm
  !>
  !> partition the given Cartesian grid
  !> in each dimension 
  !>
  subroutine partitionGrid(jmax,kmax,lmax,dx,dy,dz,xx1,nfringe)
    !
    integer, intent(inout) :: jmax,kmax,lmax
    real*8, intent(in) :: dx,dy,dz
    real*8, intent(inout) :: xx1(3)
    integer, intent(in) :: nfringe
    !
    call factorize(numprocs_spatial,iprocs)
    !
    allocate(displ1(0:iprocs(1)))
    allocate(displ2(0:iprocs(2)))
    allocate(displ3(0:iprocs(3)))
    !
    displ1(0)=0
    displ2(0)=0
    displ3(0)=0
    !
    call divide1d(iprocs(1),jmax,displ1(1))
    call divide1d(iprocs(2),kmax,displ2(1))
    call divide1d(iprocs(3),lmax,displ3(1))
    !
    m=0
    allocate(pid(iprocs(1),iprocs(2),iprocs(3)))
    !
    do l=1,iprocs(3)
       do k=1,iprocs(2)
          do j=1,iprocs(1)
             pid(j,k,l)=m
             if (m==myid_spatial) then
                id(1)=j
                id(2)=k
                id(3)=l
             endif
	     m=m+1
          enddo
       enddo
    enddo
    !
    jmax=displ1(id(1))-displ1(id(1)-1)+2*nfringe
    kmax=displ2(id(2))-displ2(id(2)-1)+2*nfringe
    lmax=displ3(id(3))-displ3(id(3)-1)+2*nfringe
    !
    xx1(1)=(displ1(id(1)-1)-(nfringe))*dx
    xx1(2)=(displ2(id(2)-1)-(nfringe))*dy
    xx1(3)=(displ3(id(3)-1)-(nfringe))*dz
    !
  end subroutine partitionGrid
  !>
  !> initialize the data buffers for
  !> communication
  !>
  subroutine initdatabuffers(istor,jmax,kmax,lmax,nq,nfringe)
    !
    implicit none
    !
    character*(*), intent(in) :: istor         !< storage type
    integer, intent(in) :: jmax !< dimension 1
    integer, intent(in) :: kmax !< dimension 2
    integer, intent(in) :: lmax !< dimension 3
    integer, intent(in) :: nq   !< number of q-variables
    integer, intent(in) :: nfringe !< number of fringes
    !
    integer :: mdim(3)
    integer :: temp
    !
    mdim=(/jmax,kmax,lmax/)
    !
    ! inline bubble sort
    !
    if (mdim(2) > mdim(1)) then
       temp=mdim(1)
       mdim(1)=mdim(2)
       mdim(2)=temp
    endif
    if (mdim(3) > mdim(2)) then
       temp=mdim(2)
       mdim(2)=mdim(3)
       mdim(3)=temp
    endif
    !
    allocate(buffer(mdim(1)*mdim(2)*nfringe*nq,12))
    allocate(ireq(12))
    allocate(mpistatus(MPI_STATUS_SIZE,12))
    irnum=0
    !
        if (istor=='row') then
       qskip=1
       qmult=nq
    else
       qskip=jmax*kmax*lmax
       qmult=1
    endif
    !
    idim(1,1)=1
    idim(2,1)=jmax
    idim(1,2)=1
    idim(2,2)=kmax
    idim(1,3)=1
    idim(2,3)=lmax
    dim0=(/jmax,kmax,lmax/)
    jkmax=jmax*kmax
    stride=(/1,jmax,jkmax/)
    ldim=(/jkmax,jkmax,jmax/)
    kdim=(/jmax,1,1/)
    !
  end subroutine initdatabuffers
  !>
  !> Post the non-blocking sends for data to be send
  !> across each plane of the Cartesian box
  !>
  subroutine initiateDataSend(q,iperiodic,iplanes,nplanes,&
                                nfringe,jmax,kmax,lmax,nq)
    implicit none
    !
    integer, intent(in) :: nfringe !< number of fringes
    integer, intent(in) :: jmax !< dimension 1
    integer, intent(in) :: kmax !< dimension 2
    integer, intent(in) :: lmax !< dimension 3
    integer, intent(in) :: nq   !< number of q-variables
    integer, intent(in) :: nplanes             !< number of planes of data
    integer, intent(in) :: iplanes(nplanes)    !< $\in$ [1,2,3,4,5,6] 
                                               !< [1,2]=[xmin,xmax],[3,4]=[ymin,ymax]
                                               !< [5,6]=[zmin,zmax]
    integer, intent(in) :: iperiodic(3)        !< periodicity flag 1=periodic
    real*8, intent(in) :: q(jmax*kmax*lmax*nq) !< field data
    !
    irnum=0
    !
    do i=1,nplanes
       plane=iplanes(i)
       idir=(plane-1)/2+1
       !
       ! find the processor to send the data to
       !
       ipsend=id
       iskip=.false.
       !
       if (mod(plane,2)==0) then
          ipsend(idir)=ipsend(idir)+1
          if (ipsend(idir) > iprocs(idir)) then
             if (iperiodic(idir)==1) then
                ipsend(idir)=1
             else
                iskip=.true.
             endif
          endif
       else
          ipsend(idir)=ipsend(idir)-1
          if (ipsend(idir) < 1) then
             if (iperiodic(idir)==1) then
                ipsend(idir)=iprocs(idir)
             else
                iskip=.true.
             endif
          endif
       endif
       !
       if (iskip) cycle
       !
       sendid=pid(ipsend(1),ipsend(2),ipsend(3))
       !
       if (idir==2) then
          ipp=mod(idir,3)+1
          ip=mod(idir+1,3)+1
       else
          ip=mod(idir,3)+1
          ipp=mod(idir+1,3)+1
       endif
       !
       if (mod(plane,2)==1) then
          js=nfringe+1
          je=2*nfringe
       else
          js=idim(2,idir)-2*nfringe+1
          je=idim(2,idir)-nfringe
       endif
       !
       ks=idim(1,ip)
       ke=idim(2,ip)
       ls=idim(1,ipp)
       le=idim(2,ipp)
       jstride=stride(idir)   
       dim2=ldim(idir)
       dim1=kdim(idir)
       !
       buffersize=0
       !
       do l=ls,le
          do k=ks,ke
             iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
             do j=js,je
                iloc=iq*qmult+1
                do n=1,nq
                   buffersize=buffersize+1
                   buffer(buffersize,plane)=q(iloc)
                   iloc=iloc+qskip
                enddo
	        iq=iq+jstride
             enddo
          enddo
       enddo
       !       
       irnum=irnum+1       
       call mpi_isend(buffer(1,plane),buffersize,MPI_REAL8,sendId,plane,&
            cartComm,ireq(irnum),ierr)
       !
    enddo
  end subroutine initiateDataSend
  !>
  !> Post the non-blocking receives for each
  !> plane that has to receive data
  !>
  subroutine initiateDataRecv(q,iperiodic,iplanes,nplanes,&
       nfringe,jmax,kmax,lmax,nq)
    implicit none
    !
    integer, intent(in) :: jmax !< dimension 1
    integer, intent(in) :: kmax !< dimension 2
    integer, intent(in) :: lmax !< dimension 3
    integer, intent(in) :: nq   !< number of q-variables
    integer, intent(in) :: nfringe !< number of fringes
    integer, intent(in) :: iperiodic(3)        !< periodicity flag 1=periodic
    integer, intent(in) :: nplanes             !< number of planes of data
    integer, intent(in) :: iplanes(nplanes)    !< $\in$ [1,2,3,4,5,6] 
                                               !< [1,2]=[xmin,xmax],[3,4]=[ymin,ymax]
                                               !< [5,6]=[zmin,zmax]
    real*8, intent(in) :: q(jmax*kmax*lmax*nq) !< q-variables
    integer :: itag
    !
    do i=1,nplanes
       plane=iplanes(i)
       idir=(plane-1)/2+1
       !
       ! find the processor to receive the data 
       !
       iprecv=id
       !
       iskip=.false.
       !
       if (mod(plane,2)==1) then
          iprecv(idir)=iprecv(idir)-1
          if (iprecv(idir) < 1) then
             if (iperiodic(idir)==1) then
                iprecv(idir)=iprocs(idir)
             else
                iskip=.true.
             endif
          endif
          itag=plane+1
       else
          iprecv(idir)=iprecv(idir)+1
          if (iprecv(idir) > iprocs(idir)) then
             if (iperiodic(idir)==1) then
                iprecv(idir)=1
             else
                iskip=.true.
             endif
          endif
          itag=plane-1
       endif
       !
       if (iskip) cycle
       !
       recvid=pid(iprecv(1),iprecv(2),iprecv(3))
       !
       if (idir==2) then
          ipp=mod(idir,3)+1
          ip=mod(idir+1,3)+1
       else
          ip=mod(idir,3)+1
          ipp=mod(idir+1,3)+1
       endif
       !
       ks=idim(1,ip)
       ke=idim(2,ip)
       ls=idim(1,ipp)
       le=idim(2,ipp)
       !
       buffersize=(le-ls+1)*(ke-ks+1)*nfringe*nq
       irnum=irnum+1
       call mpi_irecv(buffer(1,plane+6),buffersize,MPI_REAL8,recvId,itag,&
            cartComm,ireq(irnum),ierr)
    enddo
    !
  end subroutine initiateDataRecv
  !>
  !> Finalize the communication with waitalls
  !> and data update
  !>
  subroutine finalizeCommunication(q,iperiodic,iplanes,nplanes,&
       nfringe,jmax,kmax,lmax,nq)
    implicit none
    !
    integer, intent(in) :: jmax !< dimension 1
    integer, intent(in) :: kmax !< dimension 2
    integer, intent(in) :: lmax !< dimension 3
    integer, intent(in) :: nq   !< number of q-variables
    integer, intent(in) :: nfringe !< number of fringes
    integer, intent(in) :: iperiodic(3)        !< periodicity flag 1=periodic
    integer, intent(in) :: nplanes             !< number of planes of data
    integer, intent(in) :: iplanes(nplanes)    !< $\in$ [1,2,3,4,5,6] 
                                               !< [1,2]=[xmin,xmax],[3,4]=[ymin,ymax]
                                               !< [5,6]=[zmin,zmax]
    real*8, intent(inout) :: q(jmax*kmax*lmax*nq) !< q-variables
    !
    call mpi_waitall(irnum,ireq,mpistatus,ierr)
    !
    do i=1,nplanes
       plane=iplanes(i)
       idir=(plane-1)/2+1
       !
       ! make sure to skip the boundary partitions
       ! if the BC is non-periodic
       !
       if (id(idir)==1 .and. mod(plane,2)==1 .and. iperiodic(idir)==0) cycle
       if (id(idir)==iprocs(idir) .and. mod(plane,2)==0 .and.iperiodic(idir)==0) cycle
       !
       if (idir==2) then
          ipp=mod(idir,3)+1
          ip=mod(idir+1,3)+1
       else
          ip=mod(idir,3)+1
          ipp=mod(idir+1,3)+1
       endif
       !
       if (mod(plane,2)==1) then
          js=1
          je=nfringe
       else
          js=idim(2,idir)-nfringe+1
          je=idim(2,idir)
       endif
       !
       ks=idim(1,ip)
       ke=idim(2,ip)
       ls=idim(1,ipp)
       le=idim(2,ipp)
       jstride=stride(idir)   
       dim2=ldim(idir)
       dim1=kdim(idir)
       !
       bufcount=0
       do l=ls,le
          do k=ks,ke
             iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
             do j=js,je
                iloc=iq*qmult+1
                do n=1,nq
                   bufcount=bufcount+1
                   q(iloc)=buffer(bufcount,plane+6)
                   iloc=iloc+qskip
                enddo
	        iq=iq+jstride
             enddo
          enddo
       enddo
    enddo
  end subroutine finalizeCommunication
!
end module spatialCommunication
