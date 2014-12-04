subroutine compute_norm_two(norm,rhs,jmax,kmax,lmax,ninstances)
!
implicit none
include 'mpif.h'
!
integer, intent(in) :: jmax,kmax,lmax,ninstances
real*8, intent(in) :: rhs(5,jmax,kmax,lmax)
real*8, intent(inout) :: norm
real*8 :: tnorm
integer :: ierr
integer :: j,k,l
!
tnorm=0.
do l=3,lmax-2
   do k=3,kmax-2
      do j=3,jmax-2
         tnorm=tnorm+rhs(1,j,k,l)**2
      enddo
   enddo
enddo
!
call mpi_reduce(tnorm,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
norm=sqrt( norm/(jmax*kmax*lmax*ninstances) )
!
return 
end subroutine compute_norm_two

subroutine compute_norm_inf(norm,rhs,jmax,kmax,lmax,ninstances)
!
implicit none
include 'mpif.h'
!
integer, intent(in) :: jmax,kmax,lmax,ninstances
real*8, intent(in) :: rhs(5,jmax,kmax,lmax)
real*8, intent(inout) :: norm
real*8 :: tnorm
integer :: ierr
integer :: j,k,l
!
tnorm=0.
do l=3,lmax-2
   do k=3,kmax-2
      do j=3,jmax-2
         tnorm=max(abs(rhs(1,j,k,l)),tnorm)
      enddo
   enddo
enddo
!
call mpi_reduce(tnorm,norm,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpi_comm_world,ierr)
!
return 
end subroutine compute_norm_inf
