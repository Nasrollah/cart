subroutine periodic_bc(q,nq,jmax,kmax,lmax,idir,nf,istor)
!
implicit none
!
integer, intent(in) :: jmax,kmax,lmax,nq
real*8, intent(inout) :: q(nq*jmax*kmax*lmax)
integer, intent(in) :: idir,nf
character*(*), intent(in) :: istor
!
integer :: j,k,l,n
integer :: qskip,qmult,ipp,ip,dim1,dim2,jstride,iq,iq1,iloc,iloc1,jj,jkmax
integer :: dim0(3),stride(3),ldim(3),kdim(3)
!
if (istor=='row') then
   qskip=1
   qmult=nq
else
   qskip=jmax*kmax*lmax
   qmult=1
endif
!
jkmax=jmax*kmax
dim0=(/jmax,kmax,lmax/)
stride=(/1,jmax,jkmax/)
ldim=(/jkmax,jkmax,jmax/)
kdim=(/jmax,1,1/)
!
if (idir == 2) then
   ipp=mod(idir,3)+1
   ip=mod(idir+1,3)+1
else
   ip=mod(idir,3)+1
   ipp=mod(idir+1,3)+1
endif
!
dim2=ldim(idir)
dim1=kdim(idir)
jstride=stride(idir)
!
do l=1,dim0(ipp)
   do k=1,dim0(ip)
      !
      jj=dim0(idir)-2*nf+1
      iq=(l-1)*dim2+(k-1)*dim1
      iq1=(l-1)*dim2+(k-1)*dim1+(jj-1)*jstride
      !
      do j=1,nf
         iloc=iq*qmult+1
         iloc1=iq1*qmult+1
         do n=1,nq
            q(iloc)=q(iloc1)
            iloc=iloc+qskip
            iloc1=iloc1+qskip
         enddo
         iq=iq+jstride
         iq1=iq1+jstride
      enddo
      !
      jj=dim0(idir)-nf+1
      iq=(l-1)*dim2+(k-1)*dim1+(jj-1)*jstride
      iq1=(l-1)*dim2+(k-1)*dim1+nf*jstride
      !
      do j=jj,dim0(idir)
         iloc=iq*qmult+1
         iloc1=iq1*qmult+1
         do n=1,nq
            q(iloc)=q(iloc1)
            iloc=iloc+qskip
            iloc1=iloc1+qskip
         enddo
         iq=iq+jstride
         iq1=iq1+jstride
      enddo
      !
   enddo
enddo
!
return
end subroutine periodic_bc

