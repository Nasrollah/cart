!===================================================================!
!> \brief
!! This subroutine computes central difference fluxes for conservative variables 
!! on a Cartesian grid system. 2nd, 4th and 6 order physical fluxes and  4th and
!! 6th order dissipations are implemented
!!
!! Versions:\par
!!    - Sitaraman 07/28/2014
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include inviscidRHSupwind.f90
!!
!====================================================================!
subroutine inviscidRHSunified2(nq,nvar,gamma,q,s,spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,flux_order,&
                             diss_order,efac,istor)
!
implicit none
!
integer, intent(in) :: nq                        !< number of field variables stored
integer, intent(in) :: nvar                      !< number of field variables to compute residuals for
real*8, intent(in) :: gamma                      !< ratio of specific heats
real*8,  intent(in) :: q(nq*jmax*kmax*lmax)      !< solution variables
real*8,  intent(inout) :: s(nq*jmax*kmax*lmax)   !< residual
real*8,  intent(inout) :: spec(jmax*kmax*lmax)   !< speed of sound
real*8,  intent(in)    :: timeMetric(3)          !< grid speeds in three directions
real*8,  intent(in)    :: dx                     !< coordinate 1 spacing
real*8,  intent(in)    :: dy                     !< coordinate 2 spacing
real*8,  intent(in)    :: dz                     !< coordinate 3 spacing
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
integer, intent(in) :: flux_order                !< order of physical flux
integer, intent(in) :: diss_order                !< order of dissipation
real*8,  intent(in) :: efac                      !< scaling for dissipation
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
real*8  :: qcons(nvar),ds(3)
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qskip,iminus,iplus,iqp
real*8  :: faceArea,faceSpeed
real*8  :: flux(5),sflux(5),dflux(5),fluxplus(5),fluxminus(5)
real*8  ::fcoef(3),dcoef(3)
integer :: idim(2,3),stride(3),ldim(3),kdim(3)
integer :: dim0(3),dim2,dim1
!
real*8 :: gm1,ruu,vel,rvel,fdiv,eps,ediv
!
! begin
!
gm1=gamma-1
!
mdim=max(lmax,max(jmax,kmax))
allocate(f(nq,mdim),ft(nq,mdim),qt(nq,mdim),sigma(mdim))
allocate(pressure(jmax*kmax*lmax))
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
if (flux_order==2) then
   fdiv=0.5d0
   nf=1
elseif (flux_order==4) then
   fdiv=1.0d0/12.0d0
   nf=2
else
   fdiv=1.0d0/60.0d0
   nf=3
endif
!
if (diss_order==2) then
   ediv=0.5d0
   nf=1
elseif (diss_order==4) then
   ediv=1.0d0/12.0d0
   nf=max(2,nf)
else
   ediv=1.0d0/60.0d0
   nf=max(3,nf)
endif
!
eps=ediv*efac
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
!
! compute the pressure and speed of sound
! from conservative variables
!
iq=0
do l=1,lmax
   do k=1,kmax
      do j=1,jmax
         !
         ! get conservative variables here
         !
         iloc=iq*qmult+1
         do n=1,nvar
            qcons(n)=q(iloc)
            s(iloc)=0.0d0
            iloc=iloc+qskip
         end do
         !
         ruu=(qcons(4)**2+qcons(3)**2+qcons(2)**2)
         !
         iq=iq+1
         pressure(iq)=gm1*(qcons(5)-0.5*ruu/qcons(1))
         spec(iq)=sqrt(gamma*pressure(iq)/qcons(1))
         !
      enddo
   enddo
enddo
!
dim0=(/jmax,kmax,lmax/)
ds=(/dx,dy,dz/)
stride=(/1,jmax,jkmax/)
ldim=(/jkmax,jkmax,jmax/)
kdim=(/jmax,1,1/)
!
! coordinate 1 direction fluxes
!
do idir=1,1
   if (idir == 2) then
      ipp=mod(idir,3)+1
      ip=mod(idir+1,3)+1
   else
      ip=mod(idir,3)+1
      ipp=mod(idir+1,3)+1
   endif
   !
   faceArea=ds(ip)*ds(ipp)
   faceSpeed=timeMetric(idir)
   !
   js=idim(1,idir)
   je=idim(2,idir)
   ks=idim(1,ip)
   ke=idim(2,ip)
   ls=idim(1,ipp)
   le=idim(2,ipp)
   !
   jstride=stride(idir)   
   dim2=ldim(idir)
   dim1=kdim(idir)
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
            vel=qcons(idir+1)/qcons(1)
            rvel=vel-faceSpeed
            flux(1)=rvel*qcons(1)
            flux(2)=rvel*qcons(2)
            flux(3)=rvel*qcons(3)
            flux(4)=rvel*qcons(4)
            flux(5)=rvel*qcons(5)+pressure(iqp)*vel
            flux(idir+1)=flux(idir+1)+pressure(iqp)
            !
            ft(:,j)=flux
            qt(:,j)=qcons
            !
            sigma(j)=qcons(1)!spec(iqp)+abs(vel)
            !if (j.eq.20.and.k.eq.20.and.l.eq.20) then
               !write(*,*) "iru2 sig(20): ", sigma(j)
            !   sigma(j)=qcons(1)!spec(iqp)+abs(vel)  
            !endif

            iq=iq+jstride
            !
         enddo
         !
         do j=js-1,je
            if (flux_order==2) then
               sflux=(ft(:,j)+ft(:,j+1))*fdiv
            else if (flux_order==4) then
               sflux=fdiv*(-ft(:,j-1)+7d0*(ft(:,j+1)+ft(:,j))-ft(:,j+2))
            else
               sflux=fdiv*(ft(:,j-2)-8d0*(ft(:,j-1)+ft(:,j+2))+&
                    37d0*(ft(:,j)+ft(:,j+1))+ft(:,j+3))
            endif
            
            if (diss_order==2) then
               !dflux=eps*(sigma(j)+sigma(j+1))*(qt(:,j+1)-qt(:,j))
               dflux = (sigma(j)+sigma(j+1))*(qt(:,j+1)-qt(:,j))
            elseif (diss_order==4) then
               dflux=eps*(sigma(j)+sigma(j+1))*(qt(:,j-1)&
                     -3d0*(qt(:,j)-qt(:,j+1))-qt(:,j+2))
            else
               dflux=eps*(sigma(j)+sigma(j+1))*(-qt(:,j-2)+qt(:,j+3)+&
                    5d0*(qt(:,j-1)-qt(:,j+2))-10d0*(qt(:,j)-qt(:,j+1)))
            endif
            
            !f(:,j)=sflux-dflux
            f(:,j)= -dflux   
         enddo
         !
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
         !
         do j=js,je
            iloc=iq*qmult+1
            do n=1,nvar
               s(iloc)=s(iloc)-faceArea*(f(n,j)-f(n,j-1))
               iloc=iloc+qskip
            enddo
            iq=iq+jstride
         enddo
         !
      enddo
   enddo
enddo
!
deallocate(f,ft,qt,sigma,pressure)
!
return
end subroutine inviscidRHSunified2

         
! subroutine finiteDifferenceFlux(mdim,nq,nf,fdiv,eps,sigma,ft,qt,j,numericalflux)
! !
! integer, intent(in) :: mdim
! integer, intent(in) :: nq
! integer, intent(in) :: nf
! real*8, intent(in) :: fdiv
! real*8, intent(in) :: eps
! real*8, intent(in) :: sigma(mdim)         
! real*8, intent(in) :: ft(nq,mdim)
! real*8, intent(in) :: qt(nq,mdim)
! integer, intent(in) :: j
! real*8, intent(out) :: numericalflux(nq)
! !
! numericalFlux=fdiv*(ft(:,j-2)-8d0*(ft(:,j-1)+ft(:,j+2))+&
!      37d0*(ft(:,j)+ft(:,j+1))+ft(:,j+3))-&
!      eps*(sigma(j)+sigma(j+1))*&
!      (-qt(:,j-2)+qt(:,j+3)+&
!      5d0*(qt(:,j-1)-qt(:,j+2))-&
!      10d0*(qt(:,j)-qt(:,j+1)))
! !
! return
! end subroutine finiteDifferenceFlux
