!===================================================================!
!> \brief
!! This subroutine  performs a low-pass filter
!! of the given field
!! Versions:\par
!!    - Sitaraman 02/05/2015
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include filter
!!
!====================================================================!
subroutine filter(nq,nvar,jmax,kmax,lmax,q,flux_order,istor)
!
implicit none
!
integer, intent(in) :: nq                        !< number of field variables stored
integer, intent(in) :: nvar                      !< number of field variables to compute residuals for
real*8,  intent(inout) :: q(nq*jmax*kmax*lmax)      !< solution variables
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
integer, intent(in) :: flux_order                !< order of physical flux
character*(*), intent(in) :: istor               !< storage type
!
! local variables
!
integer :: j,k,l,n,mdim
integer :: js,je,ks,ke,ls,le,nf,nd
real*8, allocatable :: pressure(:)
real*8, allocatable :: f(:,:)
real*8, allocatable :: ft(:,:),qt(:,:)
real*8, allocatable :: sigma(:)
real*8 :: qcons(nvar),ds(3)
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qskip,iminus,iplus,iqp
real*8  :: faceArea,faceSpeed
real*8 :: flux(5),sflux(5),dflux(5),fluxplus(5),fluxminus(5)
real*8 ::fcoef(3),dcoef(3)
integer :: idim(2,3),stride(3),ldim(3),kdim(3)
integer :: dim0(3),dim2,dim1
!
real*8 :: gm1,ruu,vel,rvel,fdiv,eps,ediv,alphaf
!
! begin
!
mdim=max(lmax,max(jmax,kmax))
allocate(f(nq,mdim),ft(nq,mdim))
!
! use storage string to set the variable
! location order
!
if (istor=='row') then
   qskip=1
   qmult=nq
else
   qskip=jmax*kmax*lmax
   qmult=1
endif
!
nf=flux_order/2
alphaf=(-1d0)**(nf+1)*(2d0**(-2*nf))
!
js=nf+1
je=jmax-nf
ks=nf+1
ke=kmax-nf
ls=nf+1
le=lmax-nf
!
idim(1,1)=js
idim(2,1)=je
idim(1,2)=ks
idim(2,2)=ke
idim(1,3)=ls
idim(2,3)=le
jkmax=jmax*kmax
dim0=(/jmax,kmax,lmax/)
stride=(/1,jmax,jkmax/)
ldim=(/jkmax,jkmax,jmax/)
kdim=(/jmax,1,1/)
!
! coordinate 1 direction fluxes
!
do idir=1,3
   if (idir == 2) then
      ipp=mod(idir,3)+1
      ip=mod(idir+1,3)+1
   else
      ip=mod(idir,3)+1
      ipp=mod(idir+1,3)+1
   endif
   !
   js=idim(1,idir)
   je=idim(2,idir)
   ks=idim(1,ip)
   ke=idim(2,ip)
   ls=idim(1,ipp)
   le=idim(2,ipp)
   !
   do l=ls,le
      do k=ks,ke
         !
         iq=(l-1)*dim2+(k-1)*dim1
         !
         do j=1,dim0(idir)
            iloc=iq*qmult+1
            iqp=iq+1
            do n=1,nvar
               qcons(n)=q(iloc)
               iloc=iloc+qskip
            enddo
            !
            ft(:,j)=qcons
            !
            iq=iq+jstride
            !
         enddo
         !
         ! filter
         !
         do j=js,je
           if (flux_order==2) then
              f(:,j)=ft(:,j)+alphaf*(-(ft(:,j-1)+ft(:,j+1))+2d0*ft(:,j))
           else if (flux_order==4) then
              f(:,j)=ft(:,j)+alphaf*(-(ft(:,j-2)+ft(:,j+2))+&
                                            4d0*(ft(:,j-1)+ft(:,j+1))-&
                                            6d0*ft(:,j))
           else if (flux_order==6) then
               f(:,j)=ft(:,j)+alphaf*(-(ft(:,j-3)+ft(:,j+3))+&
                                             6d0*(ft(:,j-2)+ft(:,j+2))&
                                             -15d0*(ft(:,j-1)+ft(:,j+1))+&
                                             20d0*ft(:,j))
           endif
         enddo
         !
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
         !
         do j=js,je
            iloc=iq*qmult+1
            do n=1,nvar
               q(iloc)=f(n,j)
               iloc=iloc+qskip
            enddo
            iq=iq+jstride
         enddo
         !
      enddo
   enddo
enddo
!
deallocate(f,ft)
!
return
end subroutine filter
