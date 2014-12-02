!===================================================================!
!> \brief
!! This subroutine computes 2nd/4th and 6th order central difference viscous fluxes 
!! for a Cartesian grid topology. Note that 6th order requires larger
!! fringe width (nf=5), which is not supported yet in the off-body mesh
!! generator
!!
!! Versions:\par
!!    - Sitaraman 08/22/2014
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include viscousRHS.f90
!!
!====================================================================!
subroutine viscousRHS(rey,pr,prtr,nq,nvar,gamma,q,qwork,s,dx,dy,dz,jmax,kmax,lmax,&
                             flux_order,istor)
!
implicit none
!
real*8, intent(in) ::  rey                       !< reynolds number (based on speed of sound)
real*8, intent(in) :: pr                         !< prandtl number (laminar diffusion)
real*8, intent(in) :: prtr                       !< turbulent prandtl number (turbulent diffusion)
integer, intent(in) :: nq                        !< number of field variables stored
integer, intent(in) :: nvar                      !< number of field variables to compute residuals for
real*8, intent(in) :: gamma                      !< ratio of specific heats
real*8,  intent(in) :: q(nq*jmax*kmax*lmax)      !< solution variables
real*8,  intent(inout) :: qwork(9*jmax*kmax*lmax) !< work array for storing velocity gradient
real*8,  intent(inout) :: s(nq*jmax*kmax*lmax)   !< residual
real*8,  intent(in)    :: dx                     !< coordinate 1 spacing
real*8,  intent(in)    :: dy                     !< coordinate 2 spacing
real*8,  intent(in)    :: dz                     !< coordinate 3 spacing
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
integer, intent(in) :: flux_order                !< order of physical viscous flux (2, 4 or 6)
character*(*), intent(in) :: istor               !< storage type 'row' or 'column'
                                                 !< row=> q(:,j,k,l), column=>q(j,k,l,:) 
!
! local variables
!
integer :: j,k,l,n,m,mdim,p,id1,id2,jj
integer :: js,je,ks,ke,ls,le,nf,nf2,iqoffset,iqstore
real*8 :: mu,turmu,visfac,kcond
real*8  :: gradu(4,3),visflux(4)
real*8 :: qprim(nvar),ds(3)
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qstride,ilocw
real*8  :: faceArea
integer :: idim0(2,3),stride(3),ldim(3),kdim(3)
integer :: dim2,dim1,dim0(3)
!
real*8 :: gm1,ggm1,dsinv
real*8 :: prfac,prtrfac,reyinv,TWOTHIRD
real*8 :: c2b=0.3678d0
real*8 :: t1,t2
!
real*8, allocatable :: ft(:,:)
real*8, allocatable :: dcof(:),icof(:)
real*8, allocatable :: cof(:)
real*8, allocatable :: qv(:,:)
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
allocate(ft(4,mdim))
allocate(qv(11,mdim))
!
! use storage string to set the variable
! location order
!
if (istor=='row') then
   qstride=1
   qmult=nq
else
   qstride=jmax*kmax*lmax
   qmult=1
endif
!
if (flux_order==2) then
   nf=1
   allocate(dcof(2*nf),icof(2*nf),cof(2*nf))
   dcof=(/-1d0,1d0/)
   icof=(/0.5d0,0.5d0/) 
elseif (flux_order==4) then
   nf=2
   allocate(dcof(2*nf),icof(2*nf),cof(2*nf))
   dcof=(/1d0,-27d0,27d0,-1d0/)/24d0
   icof=(/-1d0,9d0,9d0,-1d0/)/16d0
else
   nf=3
   allocate(dcof(2*nf),icof(2*nf),cof(2*nf))
   dcof=(/-9d0,125d0,-2250d0,2250d0,-125d0,9d0/)/1920d0
   icof=(/3d0,-25d0,75d0,75d0,-25d0,3d0/)/256d0
endif
!
nf2=2*nf
jkmax=jmax*kmax
!
idim0(1,1)=nf2
idim0(2,1)=jmax-nf2+1
idim0(1,2)=nf2
idim0(2,2)=kmax-nf2+1
idim0(1,3)=nf2
idim0(2,3)=lmax-nf2+1
!
! change variables to primitive and
! compute velocity gradients
!
call convertVariables(gamma,nq,jmax,kmax,lmax,q,0,istor)
call cpu_time(t1)
call velocityGradient(dx,dy,dz,jmax,kmax,lmax,nq,q,qwork,flux_order,istor)
call cpu_time(t2)
!write(6,*) (t2-t1)/(jmax*kmax*lmax*nvar)
!
!
ds=(/dx,dy,dz/)
stride=(/1,jmax,jkmax/)
ldim=(/jkmax,jkmax,jmax/)
kdim=(/jmax,1,1/)
dim0=(/jmax,kmax,lmax/)
!
do idir=1,3
   !
   call cpu_time(t1)
   id1=mod(idir,3)+1
   id2=mod(idir+1,3)+1
   !
   if (idir == 2) then
      ipp=id1
      ip=id2
   else
      ip=id1
      ipp=id2
   endif
   !
   ! face area here is multiplied by inverse
   ! of Reynolds number for the appropriate final
   ! scaling
   !
   faceArea=ds(ip)*ds(ipp)*reyinv
   !
   js=idim0(1,idir)
   je=idim0(2,idir)
   ks=idim0(1,ip)
   ke=idim0(2,ip)
   ls=idim0(1,ipp)
   le=idim0(2,ipp)
   !
   jstride=stride(idir)   
   dim2=ldim(idir)
   dim1=kdim(idir)
   dsinv=1.d0/ds(idir)
   iqoffset=3*mod(idir,3)+1
   cof=dcof*dsinv
   !
   do l=ls,le
      do k=ks,ke
         !
         iq=(l-1)*dim2+(k-1)*dim1
         !
         ! collect primitive variables and 
         ! cross direction gradients in "qv"
         !
         do j=1,dim0(idir)
            iloc=iq*qmult+1
            iloc=iloc+qstride !< skip over density
            do n=1,5
               qv(n,j)=q(iloc)
               iloc=iloc+qstride
            enddo
            ilocw=iq*9
            m=iqoffset
            do n=6,11
               qv(n,j)=qwork(ilocw+m)
               m=mod(m,9)+1
            enddo
            iq=iq+jstride
         enddo
         !
         iq=(l-1)*dim2+(k-1)*dim1+(js-nf-1)*jstride
         !
         do j=js-nf,je+nf-1
            !
            qprim=0d0
            gradu=0d0
            jj=j-nf
            do m=1,nf2
               gradu(:,idir)=gradu(:,idir)+cof(m)*qv(1:4,jj+m) !< derivative operator
               qprim(:)=qprim(:)+icof(m)*qv(1:5,jj+m)          !< interpolate primitive variables
               gradu(1:3,id1)=gradu(1:3,id1)+icof(m)*qv(6:8,jj+m)  !< interpolate cross-direction gradients
               gradu(1:3,id2)=gradu(1:3,id2)+icof(m)*qv(9:11,jj+m) !< interpolate cross-direction gradients           
            enddo
            !
            mu=(c2b+1d0)*qprim(4)*sqrt(qprim(4))/(c2b+qprim(4)) !< sutherland's law for viscosity
            visfac=(mu+qprim(5))                                !< add eddy viscosity to mu
            kcond=(mu*prfac+turmu*prtrfac)                      !< non-dimensional heat conduction
            !
            ! arrange the interface fluxes
            !
            do p=1,3
               ft(p,j)=visfac*(gradu(p,idir)+gradu(idir,p))
            enddo
            !
            ft(idir,j)=ft(idir,j)-TWOTHIRD*visfac*(gradu(1,1)+gradu(2,2)+gradu(3,3))
            ft(4,j)= (ft(1,j)*qprim(1)+ft(2,j)*qprim(2)+ft(3,j)*qprim(3)) + kcond*gradu(4,idir)
            !
            iq=iq+jstride
            !
         enddo
         !
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
         !
         do j=js,je
            iloc=iq*qmult+1
            visflux=0d0
            do m=1,nf2
              visflux=visflux+dcof(m)*ft(:,j-nf+m-1) !< derivative operator
            enddo
            iloc=iloc+qstride  ! skip over continuity equation residual
            do n=1,4
               s(iloc)=s(iloc)+faceArea*visflux(n)
               iloc=iloc+qstride
            enddo
            iq=iq+jstride
         enddo
         !
      enddo
   enddo
   call cpu_time(t2)
   !write(6,*) 'idir,t2-t1=',idir,(t2-t1)/(jmax*kmax*lmax*nvar)
enddo
!
! convert variables back to conservative
!
call convertVariables(gamma,nq,jmax,kmax,lmax,q,1,istor)
!
deallocate(ft,dcof,icof,cof,qv)
!
return
end subroutine viscousRHS
