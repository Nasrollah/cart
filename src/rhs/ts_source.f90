!>
!> compute first order time-spectral derivative
!>
module temp_arrays
 integer, save :: not_allocated=1
 real*8, allocatable :: qtmp1(:,:,:,:),qtmp2(:,:,:,:)
 real*8, allocatable :: dij(:)
end module temp_arrays
! 
subroutine ts_source_term(q,rhs,vol,rank,ninstances,jmax,kmax,lmax)
!
use temp_arrays
!
implicit none
include 'mpif.h'

integer, intent(in) :: jmax,kmax,lmax
real*8, intent(in) :: q(5,jmax,kmax,lmax)
real*8, intent(in) :: vol
real*8, intent(inout) :: rhs(5,jmax,kmax,lmax)
integer, intent(in) :: rank
integer, intent(in) :: ninstances
!
integer :: i,j,k,l,n
real*8 :: pi,xx,cotan
integer :: npts,ierr
!
pi=acos(-1.)
npts=5*jmax*kmax*lmax
!
if (not_allocated.eq.1) then
 allocate(qtmp1(5,jmax,kmax,lmax))
 allocate(qtmp2(5,jmax,kmax,lmax))
 allocate(dij(ninstances))
 j=rank+1
 
 if (mod(ninstances,2).eq.0) then
    do i=1,ninstances
       xx=pi*(i-j)/ninstances
       if (i .ne. j) then
          cotan = cos(xx)/sin(xx)
          dij(i)=0.5*((-1.)**(i-j))*cotan*vol
       else
          dij(i)=0.
       endif
    enddo
 else
    do i=1,ninstances
       xx=pi*(i-j)/ninstances
       if (i .ne. j) then
          dij(i)=0.5*((-1.)**(i-j))*1.0/sin(xx)*vol
       else
          dij(i)=0.
       endif
    enddo
 endif
 !write(6,*) 'sum(dij)=',sum(dij)
 not_allocated=0
endif
!
do n=1,ninstances   
   qtmp1=dij(n)*q
   call mpi_reduce(qtmp1,qtmp2,npts,MPI_DOUBLE_PRECISION,MPI_SUM,n-1,mpi_comm_world,ierr)
enddo
!
do l=3,lmax-2
   do k=3,kmax-2
      do j=3,jmax-2
         rhs(:,j,k,l)=rhs(:,j,k,l)-qtmp2(:,j,k,l)
      enddo
   enddo
enddo
!rhs=rhs-qtmp2
!
return
end subroutine ts_source_term

