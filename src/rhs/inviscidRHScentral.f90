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
subroutine inviscidRHScentral(nq,nvar,gamma,q,s,spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,flux_order,&
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
real*8, allocatable :: df(:,:)
real*8, allocatable :: sigma(:)
real*8 :: qcons(nvar),ds(3)
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qskip,iminus,iplus,iqp
real*8  :: faceArea,faceSpeed
real*8 :: flux(5),sflux(5),dflux(5),fluxplus(5),fluxminus(5)
real*8 ::fcoef(3),dcoef(3)
integer :: idim(2,3),stride(3),ldim(3),kdim(3)
integer :: dim0(3),dim2,dim1
!
real*8 :: gm1,ruu,vel,rvel
!
! begin
!
gm1=gamma-1
!
mdim=max(lmax,max(jmax,kmax))
allocate(f(5,mdim),df(5,mdim),sigma(mdim))
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
   fcoef=(/1d0,0d0,0d0/)*0.5d0
   nf=1
elseif (flux_order==4) then
   fcoef=(/8d0,-1d0,0d0/)*1.0d0/12.0d0
   nf=2
else
   fcoef=(/45d0,-9d0,1d0/)*1.0d0/60.0d0
   nf=3
endif
!
if (diss_order==4) then
   dcoef=efac*(/3d0,-1d0,0d0/)*1.0d0/12.0d0
   nd=2
else
   dcoef=efac*(/10d0,-5d0,1d0/)*1d0/60.0d0
   nd=3
endif
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
         iloc=iq*qmult
         do n=1,nvar
            iloc=iloc+qskip
            qcons(n)=q(iloc)
            s(iloc)=0d0
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
do idir=1,3
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
            f(:,j)=0d0       !< storage for physical flux
            df(:,j)=0d0      !< storage for dissipation flux
         enddo
         !
         do j=1,nf
            iloc=iq*qmult
            iqp=iq+1
            do n=1,nvar
               iloc=iloc+qskip
               qcons(n)=q(iloc)
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
            sigma(j)=spec(iqp)+abs(vel)
            !
            do n=1,nf
               sflux=fcoef(n)*flux
               !
               iminus=j-n
               iplus=j+n
               !
               if (iminus .ge. js) then
                  f(:,iminus)=f(:,iminus)+sflux
               endif
               !
               f(:,iplus)=f(:,iplus)-sflux
               !
            enddo
            !
            do n=1,nd
               dflux=dcoef(n)*qcons
               !
               iminus=j-n
               iplus=j+n-1
               !
               if (iminus .ge. js-1) then
                  df(:,iminus)=df(:,iminus)+dflux
               endif
               !
               df(:,iplus)=df(:,iplus)-dflux
               !
            end do
            !
            iq=iq+jstride
            !
         enddo
         !
         do j=nf+1,dim0(idir)-nf
            iloc=iq*qmult
            iqp=iq+1
            do n=1,nvar
               iloc=iloc+qskip
               qcons(n)=q(iloc)
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
            sigma(j)=spec(iqp)+abs(vel)
            !
            do n=1,nf
               sflux=fcoef(n)*flux
               !
               iminus=j-n
               iplus=j+n
               !
               f(:,iminus)=f(:,iminus)+sflux
               f(:,iplus)=f(:,iplus)-sflux
               !
            enddo
            !
            do n=1,nd
               dflux=dcoef(n)*qcons
               !
               iminus=j-n
               iplus=j+n-1
               !
               df(:,iminus)=df(:,iminus)+dflux
               df(:,iplus)=df(:,iplus)-dflux
               !
            end do
            !
            iq=iq+jstride
            !
         enddo
         !
         do j=dim0(idir)-nf+1,dim0(idir)
            iloc=iq*qmult
            iqp=iq+1
            do n=1,nvar
               iloc=iloc+qskip
               qcons(n)=q(iloc)
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
            sigma(j)=spec(iqp)+abs(vel)
            !
            do n=1,nf
               sflux=fcoef(n)*flux
               !
               iminus=j-n
               iplus=j+n
               !
               f(:,iminus)=f(:,iminus)+sflux
               if (iplus.le.je) then
                  f(:,iplus)=f(:,iplus)-sflux
               endif
               !
            enddo
            !
            do n=1,nd
               dflux=dcoef(n)*qcons
               !
               iminus=j-n
               iplus=j+n-1
               !
               df(:,iminus)=df(:,iminus)+dflux
               if (iplus.le.je) then
                  df(:,iplus)=df(:,iplus)-dflux
               endif
               !
            end do
            !
            iq=iq+jstride
            !
         enddo
         !
         do j=js-1,je
            df(:,j)=(sigma(j)+sigma(j+1))*df(:,j)
         enddo
         !
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
         !
         do j=js,je
            iloc=iq*qmult
            do n=1,nvar
               iloc=iloc+qskip
               s(iloc)=s(iloc)-faceArea*(f(n,j)-df(n,j)+df(n,j-1))
            enddo
            iq=iq+jstride
         enddo
         !
      enddo
   enddo
   !call write_s(s,jdeb,kdeb,ldeb,nvar,jmax,kmax,lmax)
enddo
!
deallocate(f,df,sigma,pressure)
!
return
end subroutine inviscidRHScentral

         
         
         


