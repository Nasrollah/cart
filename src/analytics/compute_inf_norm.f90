subroutine compute_inf_norm(norm,rhs,jmax,kmax,lmax,nq,istor)
!
implicit none
include 'mpif.h'
!
integer, intent(in) :: jmax,kmax,lmax
real*8, intent(in) :: rhs(nq*jmax*kmax*lmax)
real*8, intent(inout) :: norm
character*(*) :: istor
real*8 :: tnorm
integer :: ierr
integer :: j,k,l
integer :: mj(3)
integer :: nq
real*8  :: maxnorm
integer :: qmult,qskip,iq,iloc,jkmax
integer :: npts
!
tnorm=0.
maxnorm=0.
jkmax=jmax*kmax
npts=jkmax*lmax
call getstride(qskip,qmult,istor,npts,nq)
!
do l=2,lmax-1
   do k=2,kmax-1
      iq=(l-1)*jkmax+(k-1)*jmax+1
      do j=2,jmax-1
         iloc=iq*qmult+1
         tnorm=tnorm+rhs(iloc)**2
         if (maxnorm < abs(rhs(iloc))) then 
           maxnorm=abs(rhs(iloc))
           mj=(/j,k,l/)
         endif
         iq=iq+1
      enddo
   enddo
enddo
!
write(6,*) 'maxnorm=',maxnorm,mj
!
!call mpi_reduce(tnorm,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
norm=tnorm
norm=sqrt(norm/(jmax*kmax*lmax))
!
return 
end subroutine compute_inf_norm
