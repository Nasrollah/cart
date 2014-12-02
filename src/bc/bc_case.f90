subroutine bc_case(q,nq,jmax,kmax,lmax,nf,icase,istor)
!
implicit none
integer, intent(in) :: nq,nf,jmax,kmax,lmax
real*8, intent(inout) :: q(nq*jmax*kmax*lmax)
character*(*) :: icase
character*(*) :: istor
!
if (icase=='vortex') then
   call periodic_bc(q,nq,jmax,kmax,lmax,1,nf,istor)
   call periodic_bc(q,nq,jmax,kmax,lmax,3,nf,istor)
else if (icase=='taylor-green') then
   call periodic_bc(q,nq,jmax,kmax,lmax,1,nf,istor)
   call periodic_bc(q,nq,jmax,kmax,lmax,2,nf,istor)
   call periodic_bc(q,nq,jmax,kmax,lmax,3,nf,istor)
endif
end subroutine bc_case
