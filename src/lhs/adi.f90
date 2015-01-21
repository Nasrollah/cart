!===================================================================!
!> \brief
!! This subroutine computes the implicit spatial update using
!! an (A)lternate-(D)irection (I)mplcit (ADI) algorithm.
!!
!! Versions:\par
!!    - Leffell 09/25/2014
!!
!! Uses:\par
!!    blockThomas.f90
!!
!! Source code:\par
!!   \include adi.f90
!!
!====================================================================!
subroutine adi(nq,nvar,gamma,q,s,spec,tscale,timeMetric,dx,dy,dz,jmax,kmax,lmax,flux_order,&
               diss_order,efac,istor)
!
implicit none
!
integer, intent(in)       :: nq                        !< number of field variables stored
integer, intent(in)       :: nvar                      !< number of field variables to compute residuals for
real*8,  intent(in)       :: gamma                     !< ratio of specific heats
real*8,  intent(in)       :: q(nq*jmax*kmax*lmax)      !< solution variables
real*8,  intent(inout)    :: s(nq*jmax*kmax*lmax)      !< residual
real*8,  intent(in)       :: spec(jmax*kmax*lmax)      !< speed of sound
real*8,  intent(in)       :: tscale(jmax*kmax*lmax)    !< time step
real*8,  intent(in)       :: timeMetric(3)             !< grid speeds in three directions
real*8,  intent(in)       :: dx                        !< coordinate 1 spacing
real*8,  intent(in)       :: dy                        !< coordinate 2 spacing
real*8,  intent(in)       :: dz                        !< coordinate 3 spacing
integer, intent(in)       :: jmax                      !< coordinate 1 dimension
integer, intent(in)       :: kmax                      !< coordinate 2 dimension
integer, intent(in)       :: lmax                      !< coordinate 3 dimension
integer, intent(in)       :: flux_order                !< order of physical flux
integer, intent(in)       :: diss_order                !< order of dissipation
real*8,  intent(in)       :: efac                      !< scaling for dissipation
character*(*), intent(in) :: istor                     !< storage type
!
! local variables
!
integer :: j,k,l,n,mdim,ii
integer :: js,je,lenj,ks,ke,ls,le,nf,nd
real*8, allocatable :: pressure(:),dtphys(:)
real*8, allocatable :: sigma(:)
real*8, allocatable :: a(:,:,:),b(:,:,:),c(:,:,:),dfdq(:,:,:)
real*8, allocatable :: dq(:,:)
real*8  :: sflux(nq,nq),dflux(nq,nq),dfdq_loc(nq,nq)
real*8  :: qcons(nvar),ds(3)
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qskip,iminus,iplus,iqp
real*8  :: faceArea,faceSpeed
real*8  :: fcoef(3),dcoef(3)
integer :: idim(2,3),stride(3),ldim(3),kdim(3)
integer :: dim0(3),dim2,dim1
real*8  :: sig,sigm,sigp
real*8  :: vol
real*8  :: h0,h0eps,i2ds,dt,dt_ivol
!
real*8 :: gm1,ruu,vel,rvel,fdiv,eps,ediv
real*8 :: acoeff,bcoeff,ccoeff
!
! begin
!
gm1=gamma-1
!
vol = dx*dy*dz
sflux = 0.0d0
dflux = 0.0d0
!
mdim=max(lmax,max(jmax,kmax))
allocate(sigma(mdim))
allocate(pressure(jmax*kmax*lmax),dtphys(mdim))
allocate(a(nq,nq,mdim),b(nq,nq,mdim),c(nq,nq,mdim),dfdq(nq,nq,mdim))
allocate(dq(nq,mdim))
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
   nf=max(1,nf)
else if (diss_order==4) then
   ediv=1.0d0/12.0d0
   nf=max(2,nf)
else
   ediv = 1.0d0/60.0d0
   nf=max(3,nf)
endif
!
eps=efac*0.5d0
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
! compute the pressure 
!
dt = tscale(1)
dt_ivol = dt/vol
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
            s(iloc)=dt_ivol*s(iloc)
            iloc=iloc+qskip
         end do
         !
         ruu=(qcons(4)**2+qcons(3)**2+qcons(2)**2)
         !
         iq=iq+1
         pressure(iq)=gm1*(qcons(5)-0.5*ruu/qcons(1))
         !spec(iq)=sqrt(gamma*pressure(iq)/qcons(1))
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
! evaluate flux jacobians to populate the tridiagonal
! system in each coordinate direction
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
   faceArea  = ds(ip)*ds(ipp)
   faceSpeed = timeMetric(idir)
   i2ds      = 0.5d0/ds(idir)
   h0        = dt*i2ds
   h0eps     = h0*eps
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
   ! Solve ADI system of equations
   !
   dfdq = 0.0d0
   dq   = 0.0d0
   a    = 0.0d0
   b    = 0.0d0
   c    = 0.0d0

   do l=ls,le
      do k=ks,ke
         
         iq=(l-1)*dim2+(k-1)*dim1

         ! Evaluate and store flux Jacobians 
         do j=1,dim0(idir)
            iloc=iq*qmult+1
            iqp=iq+1
            do n=1,nvar
               qcons(n)=q(iloc)
               iloc=iloc+qskip
            enddo
            !
            vel = qcons(idir+1)/qcons(1)
            rvel= vel-faceSpeed

            call fluxJacobian(nq,gm1,qcons,vel,rvel,pressure(iqp),&
                              flux_order,diss_order,efac,idir,dfdq(:,:,j))
            !
            sigma(j)  = spec(iqp)+abs(vel)
            !dtphys(j) = tscale(iqp)
            iq        = iq+jstride
            !
         enddo
                 
         ! Populate A,B & C matrices
         do j=js,je
            sigp           = sigma(j+1)
            sig            = sigma(j)
            sigm           = sigma(j-1)

            acoeff =  h0eps*(sigm +      sig       )
            bcoeff = -h0eps*(sigm +2.0d0*sig + sigp)
            ccoeff =  h0eps*(            sig + sigp)

            ! A
            sflux = -h0*dfdq(:,:,j-1)
            dflux = 0.0d0
            do ii = 1,nq
               dflux(ii,ii) = acoeff
            enddo
            a(:,:,j) =  sflux-dflux
            
            ! B
            sflux = 0.0d0
            dflux = 0.0d0
            do ii = 1,nq
               sflux(ii,ii) = 1.0d0
               dflux(ii,ii) = bcoeff
            enddo
            b(1:nq,1:nq,j) = sflux-dflux
            
            ! C
            sflux = h0*dfdq(1:nq,1:nq,j+1)
            dflux = 0.0d0
            do ii = 1,nq
               dflux(ii,ii) = ccoeff
            enddo
            c(1:nq,1:nq,j) =  sflux-dflux

         enddo
        
         ! Solve Implicit System in Coordinate IDIR
         
         ! Populate RHS
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride

         do j=js,je
            iloc=iq*qmult+1
            do n=1,nvar
               dq(n,j)=s(iloc)
               iloc=iloc+qskip
            enddo
            iq=iq+jstride
         enddo
       
         ! Solve block Tridiagonal System
         lenj = je - js + 1
         call blockThomas(a(:,:,js:je),b(:,:,js:je),c(:,:,js:je),dq(:,js:je),nq,lenj)

         ! Extract updated RHS
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
        
         do j=js,je
            iloc=iq*qmult+1
            do n=1,nvar
               s(iloc)=dq(n,j)
               iloc=iloc+qskip
            enddo
            iq=iq+jstride
         enddo
         !
      enddo
   enddo
enddo

! Scale by vol bc it's divided in lhs. Eventually get rid of this
s = s*vol
!
deallocate(sigma,pressure,dtphys)
deallocate(a,b,c,dfdq,dq)
!
return
end subroutine adi
