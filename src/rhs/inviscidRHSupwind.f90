!===================================================================!
!> \brief
!! This subroutine computes upwind fluxes for conservative variables 
!! on a Cartesian grid system. Options for muscl and 5th order weno
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
subroutine inviscidRHSupwind(nq,nvar,gamma,q,s,spec,timeMetric,dx,dy,dz,jmax,kmax,lmax,&
                             reconstruction,istor)
!
implicit none
!
integer, intent(in) :: nq                        !< number of field variables stored
integer, intent(in) :: nvar                      !< number of field variables to compute residuals for
real*8, intent(in) :: gamma                      !< ratio of specific heats
real*8,  intent(in) :: q(nq*jmax*kmax*lmax)    !< solution variables
real*8,  intent(inout) :: s(nq*jmax*kmax*lmax) !< residual
real*8,  intent(inout) :: spec(jmax,kmax,lmax)   !< spectral radius (largest eigen value)
real*8,  intent(in)    :: timeMetric(3)          !< grid speeds in three directions
real*8,  intent(in)    :: dx                     !< coordinate 1 spacing
real*8,  intent(in)    :: dy                     !< coordinate 2 spacing
real*8,  intent(in)    :: dz                     !< coordinate 3 spacing
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
character*(*), intent(in) :: reconstruction      !< upwind reconstruction type ('muscl' or 'weno5')
character*(*), intent(in) :: istor               !< storage type
!
! local variables
!
integer :: j,k,l,m,n,mdim
integer :: js,je,ks,ke,ls,le
real*8, allocatable :: ql(:,:),qr(:,:),pressure(:,:,:)
real*8, allocatable :: faceNormal(:,:),faceSpeed(:)
real*8, allocatable :: f(:,:)
real*8, allocatable :: spc(:)
real*8 :: qcons(nvar)
real*8 :: rhoi
real*8 :: th
real*8 :: qt
real*8 :: epsj,epsk,epsl
real*8 :: fmin,fmax
integer :: ibmin,ibmax
integer :: idir
integer :: jkmax
integer :: iq,iloc,qmult,qskip
real*8  :: faceArea
real*8  :: faceSpeedScal
!
real*8 :: gm1
integer :: jdeb,kdeb,ldeb
!
jdeb=3
kdeb=21
ldeb=21
!
! begin
!
gm1=gamma-1
!
! parameters for 3rd order MUSCL
!
th=1./3
qt=0.25
!
mdim=max(jmax,kmax)
mdim=max(mdim,lmax)
!
! factors used to determine limiting
!
epsj=1e-6
epsk=1e-6
epsl=1e-6
!
js=3
je=jmax-2
ks=3
ke=kmax-2
ls=3
le=lmax-2
jkmax=jmax*kmax
!
allocate(pressure(jmax,kmax,lmax))
allocate(ql(5,mdim),qr(5,mdim))
allocate(faceNormal(3,mdim),faceSpeed(mdim))
allocate(f(5,mdim))
allocate(spc(mdim))
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
! compute the pressure field
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
            s(iloc)=0d0
            iloc=iloc+qskip
         end do
         !
         pressure(j,k,l)=gm1*(qcons(5)-0.5*(qcons(4)**2+qcons(3)**2+qcons(2)**2)/qcons(1))
         spec(j,k,l)=0.
         !
         iq=iq+1
         !
      enddo
   enddo
enddo
!
! coordinate 1 direction fluxes
!
idir=1
faceArea=dy*dz
faceSpeedScal=timeMetric(1)
!
do l=ls,le
   do k=ks,ke
      !
      iq=(l-1)*jkmax+(k-1)*jmax
      do j=1,jmax 
         iloc=iq*qmult+1
         do n=1,nvar
            qcons(n)=q(iloc)
            iloc=iloc+qskip
         enddo
         rhoi=1.0/qcons(1)
         f(1,j)=qcons(1)
         f(2,j)=qcons(2)*rhoi
         f(3,j)=qcons(3)*rhoi
         f(4,j)=qcons(4)*rhoi
         f(5,j)=pressure(j,k,l)
         iq=iq+1
      enddo
      !
      do j=1,jmax-1 !js-1,je
         faceNormal(1,j)=faceArea
         faceNormal(2,j)=0.
         faceNormal(3,j)=0.
         faceSpeed(j)=0.
      enddo
      !
      ibmin=-2
      ibmax=-2
      !
      if (reconstruction=='muscl') then
         call muscld(f,ql,qr,1,jmax,th,qt,epsj,ibmin,ibmax,mdim,nvar)
      elseif (reconstruction=='weno5') then
         call fweno(f,ql,qr,1,jmax,th,qt,epsj,ibmin,ibmax,mdim,nvar)
      endif
      !      
      call roeflxCart(gamma,f,ql,qr,idir,faceArea,faceSpeed,spc,js-1,je,mdim)
      !
      iq=(l-1)*jkmax+(k-1)*jmax+js-1
      do j=js,je
         iloc=iq*qmult+1
         do n=1,nvar
            s(iloc)=s(iloc)-f(n,j)+f(n,j-1)
            iloc=iloc+qskip
         enddo
         spec(j,k,l)=spec(j,k,l)+spc(j)+spc(j-1)
         iq=iq+1
      enddo
      !
   enddo
enddo
!call write_s(s,jdeb,kdeb,ldeb,nvar,jmax,kmax,lmax)

!
! coordinate 2 direction fluxes
!
idir=2
faceArea=dz*dx
faceSpeedScal=timeMetric(2)
!
do l=ls,le
   do j=js,je
      !
      iq=(l-1)*jkmax+j-1
      do k=1,kmax 
         iloc=iq*qmult+1
         do n=1,nvar
            qcons(n)=q(iloc)
            iloc=iloc+qskip
         enddo
         rhoi=1.0/qcons(1)
         f(1,k)=qcons(1)
         f(2,k)=qcons(2)*rhoi
         f(3,k)=qcons(3)*rhoi
         f(4,k)=qcons(4)*rhoi         
         f(5,k)=pressure(j,k,l)
         iq=iq+jmax
      enddo
      !
      do k=1,kmax-1 !js-1,je
         faceNormal(1,k)=0.
         faceNormal(2,k)=faceArea
         faceNormal(3,k)=0.
         faceSpeed(j)=0.
      enddo

      !
      ibmin=-2
      ibmax=-2
      !
      if (reconstruction=='muscl') then
         call muscld(f,ql,qr,1,kmax,th,qt,epsk,ibmin,ibmax,mdim,nvar)
      elseif (reconstruction=='weno5') then
         call fweno(f,ql,qr,1,kmax,th,qt,epsk,ibmin,ibmax,mdim,nvar)
      endif
      !
      call roeflxCart(gamma,f,ql,qr,idir,faceArea,faceSpeedScal,spc,ks-1,ke,mdim)      
      !call roeflx(f,ql,qr,faceNormal,faceSpeed,spc,js-1,je,mdim)
      !
      iq=(l-1)*jkmax+(ks-1)*jmax+j-1
      do k=ks,ke
         iloc=iq*qmult+1
         do n=1,nvar
            s(iloc)=s(iloc)-f(n,k)+f(n,k-1)
            iloc=iloc+qskip
         enddo
         spec(j,k,l)=spec(j,k,l)+spc(j)+spc(j-1)
         iq=iq+jmax
      enddo
      !
   enddo
enddo
!call write_s(s,jdeb,kdeb,ldeb,nvar,jmax,kmax,lmax)
!
! coordinate 3 direction fluxes
!
idir=3
faceArea=dx*dy
faceSpeedScal=timeMetric(3)
!
do k=ks,ke
   do j=js,je
      !
      iq=(k-1)*jmax+j-1
      do l=1,lmax 
         iloc=iq*qmult+1
         do n=1,nvar
            qcons(n)=q(iloc)
            iloc=iloc+qskip
         enddo
         rhoi=1.0/qcons(1)
         f(1,l)=qcons(1)
         f(2,l)=qcons(2)*rhoi
         f(3,l)=qcons(3)*rhoi
         f(4,l)=qcons(4)*rhoi
         f(5,l)=pressure(j,k,l)
         iq=iq+jkmax
      enddo
      !
      ibmin=-2
      ibmax=-2
      !
      if (reconstruction=='muscl') then
         call muscld(f,ql,qr,1,lmax,th,qt,epsl,ibmin,ibmax,mdim,nvar)
      elseif (reconstruction=='weno5') then
         call fweno(f,ql,qr,1,lmax,th,qt,epsl,ibmin,ibmax,mdim,nvar)
      endif
      !
      call roeflxCart(gamma,f,ql,qr,idir,faceArea,faceSpeedScal,spc,ls-1,le,mdim)
      !
      iq=(ls-1)*jkmax+(k-1)*jmax+j-1
      do l=ls,le
         iloc=iq*qmult+1
         do n=1,nvar
            s(iloc)=s(iloc)-f(n,l)+f(n,l-1)
            iloc=iloc+qskip
         enddo
         spec(j,k,l)=spec(j,k,l)+spc(j)+spc(j-1)
         iq=iq+jkmax
      enddo
      !
   enddo
enddo
!call write_s(s,jdeb,kdeb,ldeb,nvar,jmax,kmax,lmax)
!
! fix the spectral radii values
! at the boundaries
!
do j=1,js-1
  do k=ks,ke
   do l=ls,le
    spec(j,k,l)=spec(js,k,l)
   enddo
  enddo
enddo
do j=je+1,jmax
   do k=ks,ke
    do l=ls,le
     spec(j,k,l)=spec(je,k,l)
    enddo
   enddo
enddo
!
do j=1,jmax
  do k=1,ks-1
   do l=ls,le
    spec(j,k,l)=spec(j,ks,l)
   enddo
  enddo
enddo
do j=1,jmax
   do k=ke+1,kmax
    do l=ls,le
     spec(j,k,l)=spec(j,ke,l)
    enddo
   enddo
enddo
!
do j=1,jmax
  do k=1,kmax
   do l=1,ls-1
    spec(j,k,l)=spec(j,k,ls)
   enddo
  enddo
enddo
do j=1,jmax
   do k=1,kmax
    do l=le+1,lmax
     spec(j,k,l)=spec(j,k,le)
    enddo
   enddo
enddo
!
deallocate(pressure,ql,qr,f,spc)
!
return
end subroutine inviscidRHSupwind
      


         
         
         


