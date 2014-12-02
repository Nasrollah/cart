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
