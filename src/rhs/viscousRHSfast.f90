!===================================================================!
!> \brief
!! This subroutine computes 2nd order central difference viscous fluxes 
!! for a Cartesian grid topology. This is a low memory fast implementation
!! compared to the viscousRHS that does general 2nd, 4th and 6th order
!! viscous fluxes
!! Versions:\par
!!    - Sitaraman 09/11/2014
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include viscousRHSfast.f90
!!
!====================================================================!
subroutine viscousRHSfast(rey,pr,prtr,nq,nvar,gamma,q,s,dx,dy,dz,jmax,kmax,lmax)
!
implicit none
!
real*8, intent(in) ::  rey                       !< reynolds number (based on speed of sound)
real*8, intent(in) :: pr                         !< prandtl number (laminar diffusion)
real*8, intent(in) :: prtr                       !< turbulent prandtl number (turbulent diffusion)
integer, intent(in) :: nq                        !< number of field variables stored
integer, intent(in) :: nvar                      !< number of field variables to compute residuals for
real*8, intent(in) :: gamma                      !< ratio of specific heats
real*8,  intent(in) :: q(nq,jmax,kmax,lmax)      !< solution variables
real*8,  intent(inout) :: s(nq,jmax,kmax,lmax)   !< residual
real*8,  intent(in)    :: dx                     !< coordinate 1 spacing
real*8,  intent(in)    :: dy                     !< coordinate 2 spacing
real*8,  intent(in)    :: dz                     !< coordinate 3 spacing
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
!
! local variables
!
integer :: j,k,l,mdim,p
integer :: js,je,ks,ke,ls,le
real*8 :: mu,turmu,visfac,kcond
real*8  :: gradu(3,4)
real*8 :: qprim(nvar)
real*8  :: faceArea
!
real*8 :: gm1,ggm1
real*8 :: prfac,prtrfac,reyinv,TWOTHIRD
real*8 :: c2b=0.3678d0
integer :: idir
!
real*8, allocatable :: ft(:,:)
real*8 :: dxi,dyi,dzi
!
! begin
!
gm1=gamma-1
ggm1=gamma*gm1
prfac=1d0/(gm1*pr)
prtrfac=1d0/(gm1*prtr)
reyinv=1.d0/rey
TWOTHIRD=2.d0/3.d0
!
mdim=max(lmax,max(jmax,kmax))
allocate(ft(nq,mdim))
!
! change variables to primitive and
! compute velocity gradients
!
js=2
je=jmax-1
ls=2
le=lmax-1
ks=2
ke=kmax-1
!
call convertVariables(gamma,nq,jmax,kmax,lmax,q,0,'row')
!
! fix me need to implement this
! for row and column and unify idir=1,3
!
faceArea=dy*dz*reyinv
dxi=1./dx
dyi=0.25d0/dy
dzi=0.25d0/dz
idir=1
!
do l=ls,le
   do k=ks,ke
      do j=1,jmax-1
         !
         gradu(1,1:3)=(q(2:4,j+1,k,l)-q(2:4,j,k,l))*dxi
         gradu(2,1:3)=(q(2:4,j,k+1,l)-q(2:4,j,k-1,l)+&
                      q(2:4,j+1,k+1,l)-q(2:4,j+1,k-1,l))*dyi
         gradu(3,1:3)=(q(2:4,j,k,l+1)-q(2:4,j,k,l-1)+&
                      q(2:4,j+1,k,l+1)-q(2:4,j+1,k,l-1))*dzi
         gradu(1,4)=(q(5,j+1,k,l)-q(5,j,k,l))*dxi
         !
         qprim(1:5)=(q(2:6,j,k,l)+q(2:6,j+1,k,l))*0.5d0
         !if (j==5 .and. k==5 .and. l==5) then
         !   write(6,*) 'gradu(1,:)=',gradu(1,:)
         !   write(6,*) 'gradu(2,:)=',gradu(2,:)
         !   write(6,*) 'gradu(3,:)=',gradu(3,:)
         !   write(6,*) visfac,kcond,turmu
         !   stop
         !endif
         !
         mu=(c2b+1d0)*qprim(4)*sqrt(qprim(4))/(c2b+qprim(4)) !< sutherland's law for viscosity
         turmu=qprim(5)
         visfac=(mu+turmu)
         kcond=(mu*prfac+turmu*prtrfac)
         !
         ! arrange the interface fluxes
         !
         do p=1,3
            ft(p,j)=visfac*(gradu(p,idir)+gradu(idir,p))
         enddo
         !
         ft(1,j)=ft(1,j)-TWOTHIRD*visfac*(gradu(1,1)+gradu(2,2)+gradu(3,3))
         ft(4,j)= (ft(1,j)*qprim(1)+ft(2,j)*qprim(2)+ft(3,j)*qprim(3)) + kcond*gradu(1,4)
      enddo
      !
      do j=js,je
         s(2:5,j,k,l)=s(2:5,j,k,l)+faceArea*(ft(:,j)-ft(:,j-1))
      enddo
   enddo
enddo
!
faceArea=dx*dz*reyinv
dxi=0.25d0/dx
dyi=1d0/dy
dzi=0.25d0/dz
idir=2
!
do l=ls,le
   do j=js,je
      do k=1,kmax-1
         !
         gradu(2,1:3)=(q(2:4,j,k+1,l)-q(2:4,j,k,l))*dyi
         gradu(1,1:3)=(q(2:4,j+1,k+1,l)-q(2:4,j-1,k+1,l)+&
                      q(2:4,j+1,k,l)-q(2:4,j-1,k,l))*dxi
         gradu(3,1:3)=(q(2:4,j,k+1,l+1)-q(2:4,j,k+1,l-1)+&
                      q(2:4,j,k,l+1)-q(2:4,j,k,l-1))*dzi
         gradu(2,4)=(q(5,j,k+1,l)-q(5,j,k,l))*dyi
         !
         qprim(1:5)=(q(2:6,j,k,l)+q(2:6,j,k+1,l))*0.5d0
         !
         mu=(c2b+1d0)*qprim(4)*sqrt(qprim(4))/(c2b+qprim(4)) !< sutherland's law for viscosity
         turmu=qprim(5)
         visfac=(mu+turmu)
         kcond=(mu*prfac+turmu*prtrfac)
         !
         ! arrange the interface fluxes
         !
         do p=1,3
            ft(p,k)=visfac*(gradu(p,idir)+gradu(idir,p))
         enddo
         !
         ft(2,k)=ft(2,k)-TWOTHIRD*visfac*(gradu(1,1)+gradu(2,2)+gradu(3,3))
         ft(4,k)= (ft(1,k)*qprim(1)+ft(2,k)*qprim(2)+ft(3,k)*qprim(3)) + kcond*gradu(2,4)
      enddo
      !
      do k=ks,ke
         s(2:5,j,k,l)=s(2:5,j,k,l)+faceArea*(ft(:,k)-ft(:,k-1))
      enddo
   enddo
enddo
!
faceArea=dx*dy*reyinv
dxi=0.25d0/dx
dyi=0.25d0/dy
dzi=1d0/dz
idir=3
!
do k=ks,ke
   do j=js,je
      do l=1,lmax-1
         !
         gradu(3,1:3)=(q(2:4,j,k,l+1)-q(2:4,j,k,l))*dzi
         gradu(1,1:3)=(q(2:4,j+1,k,l+1)-q(2:4,j-1,k,l+1)+&
                      q(2:4,j+1,k,l)-q(2:4,j-1,k,l))*dxi
         gradu(2,1:3)=(q(2:4,j,k+1,l+1)-q(2:4,j,k-1,l+1)+&
                      q(2:4,j,k+1,l)-q(2:4,j,k-1,l))*dyi
         gradu(3,4)=(q(5,j,k,l+1)-q(5,j,k,l))*dzi
         !
         qprim(1:5)=(q(2:6,j,k,l)+q(2:6,j,k,l+1))*0.5d0
         !
         mu=(c2b+1d0)*qprim(4)*sqrt(qprim(4))/(c2b+qprim(4)) !< sutherland's law for viscosity
         turmu=qprim(5)
         visfac=(mu+turmu)
         kcond=(mu*prfac+turmu*prtrfac)
         !
         ! arrange the interface fluxes
         !
         do p=1,3
            ft(p,l)=visfac*(gradu(p,idir)+gradu(idir,p))
         enddo
         !
         ft(3,l)=ft(3,l)-TWOTHIRD*visfac*(gradu(1,1)+gradu(2,2)+gradu(3,3))
         ft(4,l)= (ft(1,l)*qprim(1)+ft(2,l)*qprim(2)+ft(3,l)*qprim(3)) + kcond*gradu(3,4)
      enddo
      !
      do l=ls,le
         s(2:5,j,k,l)=s(2:5,j,k,l)+faceArea*(ft(:,l)-ft(:,l-1))
      enddo
   enddo
enddo
!
! convert variables back to conservative
!
call convertVariables(gamma,nq,jmax,kmax,lmax,q,1,'row')
!
deallocate(ft)
!
return
end subroutine viscousRHSfast
