!
subroutine updateAllFringes(q,iperiodic,nf,jmax,kmax,lmax,nq)
  !
  use spatialCommunication, only : initSpatialComm,partitionGrid,&
	                           initDataBuffers,initiateDataSend,&
                                   initiateDataRecv,finalizeCommunication
  implicit none
  !
  integer, intent(in) :: nf   !< number of fringes
  integer, intent(in) :: jmax !< dimension 1
  integer, intent(in) :: kmax !< dimension 2
  integer, intent(in) :: lmax !< dimension 3
  integer, intent(in) :: nq   !< number of q-variables
  !< [1,2]=[xmin,xmax],[3,4]=[ymin,ymax]
  !< [5,6]=[zmin,zmax]
  integer, intent(in)   :: iperiodic(3)        !< periodicity flag 1=periodic
  real*8, intent(inout) :: q(jmax*kmax*lmax*nq) !< field data
  !
  integer :: nplanes
  integer, dimension(6) :: iplanes_send,iplanes_recv
  !
  nplanes=2
  iplanes_send=(/1,2,0,0,0,0/)
  call initiateDataSend(q,iperiodic,iplanes_send,nplanes,nf,jmax,kmax,lmax,nq)
  iplanes_recv=(/2,1,0,0,0,0/)
  call initiateDataRecv(q,iperiodic,iplanes_recv,nplanes,nf,jmax,kmax,lmax,nq)
  call finalizeCommunication(q,iperiodic,iplanes_recv,nplanes,nf,jmax,kmax,lmax,nq)    
  !
  iplanes_send=(/3,4,0,0,0,0/)
  call initiateDataSend(q,iperiodic,iplanes_send,nplanes,nf,jmax,kmax,lmax,nq)
  iplanes_recv=(/4,3,0,0,0,0/)
  call initiateDataRecv(q,iperiodic,iplanes_recv,nplanes,nf,jmax,kmax,lmax,nq)
  call finalizeCommunication(q,iperiodic,iplanes_recv,nplanes,nf,jmax,kmax,lmax,nq)    
  !
  iplanes_send=(/5,6,0,0,0,0/)
  call initiateDataSend(q,iperiodic,iplanes_send,nplanes,nf,jmax,kmax,lmax,nq)
  iplanes_recv=(/6,5,0,0,0,0/)
  call initiateDataRecv(q,iperiodic,iplanes_recv,nplanes,nf,jmax,kmax,lmax,nq)
  call finalizeCommunication(q,iperiodic,iplanes_recv,nplanes,nf,jmax,kmax,lmax,nq)    
  !
end subroutine updateAllFringes
