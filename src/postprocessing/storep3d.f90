!>
!> Store plot3d data
!> 
subroutine storep3d(x,q,istep,fsmach,alpha,rey,totime,jmax,kmax,lmax,nvar,nf,istor)
!
implicit none
integer, intent(in) :: istep
integer, intent(in) :: nf
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in) :: fsmach,alpha,rey,totime
real*8, intent(in) :: x(3*jmax*kmax*lmax)
real*8, intent(in) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
integer :: j,k,l,n
integer :: jkmax,npts
character*16 :: integer_string,filename
integer :: xskip,xmult,qskip,qmult
!
! write output plot3d files
!
jkmax=jmax*kmax
npts=jkmax*lmax
call getstride(xskip,xmult,istor,npts,3)
call getstride(qskip,qmult,istor,npts,nvar)
!
write(integer_string,"(I7)") 1000000+istep
filename='x'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax-2*nf,kmax-2*nf,lmax-2*nf
write(1) ((((x(((l-1)*jkmax+(k-1)*jmax+j-1)*xmult+1+(n-1)*xskip),&
     j=nf+1,jmax-nf),k=nf+1,kmax-nf),l=nf+1,lmax-nf),n=1,3)
close(1)
!
filename='q'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax-2*nf,kmax-2*nf,lmax-2*nf
!
write(1) fsmach,alpha,rey,totime
write(1) ((((q(((l-1)*jkmax+(k-1)*jmax+j-1)*qmult+1+(n-1)*qskip),&
     j=nf+1,jmax-nf),k=nf+1,kmax-nf),l=nf+1,lmax-nf),n=1,5)
close(1)

return
end subroutine storep3d

!>
!> Store plot3d data
!> 
subroutine storep3di(x,q,istep,fsmach,alpha,rey,totime,jmax,kmax,lmax,nvar,nf,istor)
!
implicit none
integer, intent(in) :: istep
integer, intent(in) :: nf
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in) :: fsmach,alpha,rey,totime
real*8, intent(in) :: x(3*jmax*kmax*lmax)
real*8, intent(in) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
integer :: j,k,l,n
integer :: jkmax,npts
character*16 :: integer_string,filename
integer :: xskip,xmult,qskip,qmult
integer, allocatable :: iblank(:,:,:)
!
! write output plot3d files
!
allocate(iblank(jmax,kmax,lmax))
iblank=1
!
do l=1,lmax
   do k=1,kmax
      do j=1,nf
         iblank(j,k,l)=-1
      end do
      do j=jmax-nf+1,jmax
         iblank(j,k,l)=-1
      enddo
   enddo
enddo
!
do l=1,lmax
   do j=1,jmax
      do k=1,nf
         iblank(j,k,l)=-1
      end do
      do k=kmax-nf+1,kmax
         iblank(j,k,l)=-1
      enddo
   enddo
enddo
!
do j=1,jmax
   do k=1,kmax
      do l=1,nf
         iblank(j,k,l)=-1
      end do
      do l=lmax-nf+1,lmax
         iblank(j,k,l)=-1
      enddo
   enddo
enddo
!
jkmax=jmax*kmax
npts=jkmax*lmax
call getstride(xskip,xmult,istor,npts,3)
call getstride(qskip,qmult,istor,npts,nvar)
!
write(integer_string,"(I7)") 1000000+istep
filename='x'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax,kmax,lmax
write(1) ((((x(((l-1)*jkmax+(k-1)*jmax+j-1)*xmult+1+(n-1)*xskip),&
     j=1,jmax),k=1,kmax),l=1,lmax),n=1,3),&
     (((iblank(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
close(1)
!
filename='q'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax,kmax,lmax
!
write(1) fsmach,alpha,rey,totime
write(1) ((((q(((l-1)*jkmax+(k-1)*jmax+j-1)*qmult+1+(n-1)*qskip),&
     j=1,jmax),k=1,kmax),l=1,lmax),n=1,5)
close(1)
!
deallocate(iblank)
return
end subroutine storep3di
