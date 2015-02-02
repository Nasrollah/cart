!===================================================================!
!> \brief
!! This subroutine computes the implicit spatial update using
!! Gauss-Seidel Line Relaxation
!!
!! This is a naive implemntation with no storage. A more sophisticated 
!! approach will follow.
!!
!! Versions:\par
!!    - Leffell 12/16/2014
!!
!! Uses:\par
!!    blockThomas.f90
!!
!! Source code:\par
!!   \include gslr.f90
!!
!====================================================================!
subroutine gslr(nq,nvar,gamma,q,s,spec,tscale,timeMetric,dx,dy,dz,jmax,kmax,lmax,flux_order,&
                      diss_order,efac,istor,nsweep,myid)
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
integer, intent(in)       :: nsweep                    !< number of gs sweeps in each direction
integer, intent(in) :: myid
!
! local variables
!
integer :: j,k,l,n,mdim,ii
integer :: js,je,lenj,ks,ke,ki,ls,le,li,nf,nd
real*8, allocatable :: pressure(:),dpressure(:),dtphys(:)
real*8, allocatable :: sigma(:),dsigma(:)
real*8, allocatable :: a(:,:,:),b(:,:,:),c(:,:,:),dfdq(:,:,:)
real*8, allocatable :: dq_line(:,:),dq(:)
real*8  :: sflux(nq,nq),dflux(nq,nq),dfdq_loc(nq,nq)
real*8  :: qcons(nvar),dqcons(nvar),ds(3)
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qskip,iminus,iplus,iqp
real*8  :: faceArea,faceSpeed
real*8  :: fcoef(3),dcoef(3)
integer :: idim(2,3),stride(3),ldim(3),kdim(3)
integer :: dim0(3),dim2,dim1
real*8  :: sig,sigm,sigp
real*8  :: vol
real*8  :: h0,h0eps,i2ds,dt,dt_ivol
real*8  :: dsig
real*8  :: irho,irho2,u2
integer :: isweep
!
real*8 :: gm1,ruu,vel,rvel,fdiv,eps,ediv
!real*8 :: ts,te,tc,tfj,tabc,tbt
real*8 :: acoeff,bcoeff,ccoeff
!
! begin
!
gm1=gamma-1
!
vol   = dx*dy*dz
sflux = 0.0d0
dflux = 0.0d0
!
mdim=max(lmax,max(jmax,kmax))
allocate(sigma(mdim),dsigma(jmax*kmax*lmax),dq(nq*jmax*kmax*lmax))
allocate(pressure(jmax*kmax*lmax),dpressure(jmax*kmax*lmax),dtphys(mdim))
allocate(a(nq,nq,mdim),b(nq,nq,mdim),c(nq,nq,mdim),dfdq(nq,nq,mdim))
allocate(dq_line(nq,mdim))
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
dim0=(/jmax,kmax,lmax/)
ds=(/dx,dy,dz/)
stride=(/1,jmax,jkmax/)
ldim=(/jkmax,jkmax,jmax/)
kdim=(/jmax,1,1/)
!
! evaluate the diagonal component of the entire system (dsigma)
!
dsigma = 1.0d0
do idir=1,3
   if (idir == 2) then
      ipp=mod(idir,3)+1
      ip=mod(idir+1,3)+1
   else
      ip=mod(idir,3)+1
      ipp=mod(idir+1,3)+1
   endif

   i2ds      = 0.5d0/ds(idir)
   h0        = dt*i2ds
   h0eps     = h0*eps
   
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
   ! Evaluate sigma in j-lines for each k,l pair
   !            
   do l=ls,le
      do k=ks,ke
         !                                                                    
         iq=(l-1)*dim2+(k-1)*dim1
         
         do j=1,dim0(idir)
            iloc=iq*qmult+1
            iqp=iq+1
            do n=1,nvar
               qcons(n)=q(iloc)
               iloc=iloc+qskip
            enddo
            !                                                                             
            vel = qcons(idir+1)/qcons(1)
            sigma(j)  = spec(iqp)+abs(vel)
            iq        = iq+jstride
            !                                                                             
         enddo
         
         ! Evaluate dsigma
         iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
         
         do j=js,je
            iqp = iq+1
            dsigma(iqp) = dsigma(iqp) + h0eps*(sigma(j-1) + 2.0d0*sigma(j) + sigma(j+1))
            iq=iq+jstride
         enddo
         !                                                                                
      enddo
   enddo
enddo
!
! evaluate flux jacobians to populate the tridiagonal
! system in each coordinate direction
!
dq=0.0d0
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
   dfdq    = 0.0d0
   dq_line = 0.0d0
   a       = 0.0d0
   b       = 0.0d0
   c       = 0.0d0

   do isweep = 1,nsweep

      ! Forward Sweep
      ki = 1
      li = 1
      call gslr_sweep(nq,nvar,gm1,pressure,flux_order,diss_order,efac,idir,&
           ip,ipp,dfdq,sigma,a,b,c,js,je,ks,ke,ki,ls,le,li,dim0,dim1,dim2,q,&
           dq,s,dq_line,jmax,kmax,lmax,mdim,&
           timeMetric,qskip,qmult,jstride,spec,faceSpeed,&
           h0,h0eps,eps,sflux,dflux,dsigma,myid)
      
      ! Backward Sweep
      ki = -1
      li = -1
      call gslr_sweep(nq,nvar,gm1,pressure,flux_order,diss_order,efac,idir,&
           ip,ipp,dfdq,sigma,a,b,c,js,je,ke,ks,ki,le,ls,li,dim0,dim1,dim2,q,&
           dq,s,dq_line,jmax,kmax,lmax,mdim,&
           timeMetric,qskip,qmult,jstride,spec,faceSpeed,&
           h0,h0eps,eps,sflux,dflux,dsigma,myid)
   enddo
enddo
!
! Scale by vol bc it's divided in lhs. Eventually get rid of this
!
s = vol*dq
!
deallocate(sigma,dsigma,pressure,dpressure,dtphys)
deallocate(a,b,c,dfdq,dq,dq_line)
!
return
end subroutine gslr

!
subroutine gslr_sweep(nq,nvar,gm1,pressure,flux_order,diss_order,efac,idir,&
                      ip,ipp,dfdq,sigma,a,b,c,js,je,ks,ke,ki,ls,le,li,dim0,dim1,dim2,q,&
                      dq,s,dq_line,jmax,kmax,lmax,mdim,&
                      timeMetric,qskip,qmult,jstride,spec,faceSpeed,&
                      h0,h0eps,eps,sflux,dflux,dsigma,myid)
!
implicit none
!
integer,                              intent(in)    :: nq,nvar,mdim,qskip,qmult,jstride
real*8,                               intent(in)    :: gm1
real*8,  dimension(jmax*kmax*lmax),   intent(in)    :: spec,pressure
integer,                              intent(in)    :: flux_order,diss_order
real*8,                               intent(in)    :: efac,faceSpeed
integer,                              intent(in)    :: idir,ip,ipp
real*8, dimension(mdim),              intent(inout) :: sigma
real*8, dimension(nq,nq,mdim),        intent(inout) :: dfdq,a,b,c
integer,                              intent(in)    :: js,je,ks,ke,ki,ls,le,li
integer,                              intent(in)    :: jmax,kmax,lmax
integer,                              intent(in)    :: dim1,dim2
real*8, dimension(nq*jmax*kmax*lmax), intent(inout) :: dq
real*8, dimension(nq*jmax*kmax*lmax), intent(in)    :: q,s
real*8, dimension(nq,mdim),           intent(inout) :: dq_line
real*8, dimension(3),                 intent(in)    :: timeMetric
integer, dimension(3),                intent(in)    :: dim0
real*8,                               intent(in)    :: h0,h0eps,eps
real*8, dimension(nq,nq),             intent(inout) :: sflux,dflux
real*8, dimension(jmax*kmax*lmax),    intent(in)    :: dsigma
integer,                              intent(in)    :: myid
!
integer :: j,k,l,ii,n,ik,il
integer :: iq,iloc,iqp,lenj
real*8  :: qcons(nvar),dqcons(nvar),vel,rvel
real*8  :: acoeff,bcoeff,ccoeff
real*8  :: sigm,sig,sigp
real*8  :: bdq_p(nvar),bdq_m(nvar),cdq_p(nvar),cdq_m(nvar),rhs(nvar)
real*8 :: ts,te,tfj,tabc,trhs,tblock,tupdate
!
tfj = 0.0d0
tabc = 0.0d0
trhs = 0.0d0
tblock = 0.0d0
tupdate = 0.0d0
do l=ls,le,li
   do k=ks,ke,ki
      
      !>
      ! Evaluate and store flux Jacobians in j-line 
      !>
      call cpu_time(ts)
      call store_fj(q,mdim,nq,nvar,k,l,dim0,dim1,dim2,qmult,qskip,idir,faceSpeed,gm1,efac,dfdq,&
           pressure,spec,sigma,flux_order,diss_order,jstride,jmax,kmax,lmax)
      call cpu_time(te)
      tfj = tfj + te - ts
      !>
      ! Populate A,B & C matrices in j-line
      !>
      call cpu_time(ts)
      call store_abc(nq,jmax,kmax,lmax,js,je,k,l,dfdq,mdim,a,b,c,&
	   h0,h0eps,sigma,dsigma,sflux,dflux,dim1,dim2,jstride)
      call cpu_time(te)
      tabc =tabc + te - ts
      !>                                                                                              
      !> Evaluate the matrix-vector products Adq in k- and l-coordinate directions
      !>
      call cpu_time(ts)
      do j=js,je
         
         ! k-direction
         ik = 1
         il = 0
         call rhs_gs(bdq_p,bdq_m,nq,nvar,j,k,l,ip, h0,eps,q,dq,jmax,kmax,lmax,qmult,&
	       qskip,dim1,dim2,jstride,gm1,timeMetric,spec,ik,il)
         
         ! l-direction
         ik = 0
         il = 1
         call rhs_gs(cdq_p,cdq_m,nq,nvar,j,k,l,ipp,h0,eps,q,dq,jmax,kmax,lmax,qmult,&
	    qskip,dim1,dim2,jstride,gm1,timeMetric,spec,ik,il)
    
         ! Update RHS
         rhs = bdq_m + bdq_p + cdq_m + cdq_p
         !
         call get_iloc(j,k,l,qmult,dim1,dim2,jstride,iloc,iqp)
         !
         do n=1,nvar
            dq_line(n,j)=s(iloc) - rhs(n)
            iloc=iloc+qskip
         enddo
      enddo
      call cpu_time(te)
      trhs=trhs+te-ts
      
      ! Solve block Tridiagonal System
      lenj = je - js + 1
      call cpu_time(ts)
      call blockThomas(a(:,:,js:je),b(:,:,js:je),c(:,:,js:je),dq_line(:,js:je),nq,lenj)
      call cpu_time(te)
      tblock = tblock + te-ts
      
      ! Extract updated RHS
      iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
      
      do j=js,je
         iloc=iq*qmult+1
         do n=1,nvar
            dq(iloc)=dq_line(n,j)
            iloc=iloc+qskip
         enddo
         iq=iq+jstride
      enddo
   enddo
enddo
!
!write(*,*) "tfj = ",tfj
!write(*,*) "tabc = ",tabc
!write(*,*) "trhs = ", trhs
!write(*,*) "tblock = ", tblock
!
return
end subroutine gslr_sweep

subroutine rhs_gs(adq_p,adq_m,nq,nvar,j,k,l,idir,h0,eps,q,dq,jmax,kmax,lmax,qmult,&
	          qskip,dim1,dim2,jstride,gm1,timeMetric,spec,ik,il)
!
implicit none
!
real*8, dimension(nq),                intent(inout) :: adq_p,adq_m
integer,                              intent(in)    :: nq,nvar,j,k,l,idir,jmax,kmax,lmax,qmult,qskip,dim1,dim2,jstride
real*8,                               intent(in)    :: h0,eps,gm1
real*8, dimension(3),                 intent(in)    :: timeMetric
real*8, dimension(nq*jmax*kmax*lmax), intent(in)    :: q,dq
real*8, dimension(jmax*kmax*lmax),    intent(in)    :: spec
integer,                              intent(in)    :: ik,il
!
integer :: iqp
real*8  :: vel,sigm,sig,sigp,mcoeff,pcoeff
real*8  :: qcons(nq),dqcons(nq)
!>
! (j,k,l)
!>
call get_qcons(q,qcons,nq,nvar,jmax,kmax,lmax,j,k,l,qmult,qskip,dim1,dim2,jstride,iqp)
vel    = qcons(idir+1)/qcons(1)
sig    = spec(iqp) + abs(vel)
!>
! (j,k-ik,l-il) 
!>
call get_qcons(q,  qcons,nq,nvar,jmax,kmax,lmax,j,k-ik,l-il,qmult,qskip,dim1,dim2,jstride,iqp)
call get_qcons(dq,dqcons,nq,nvar,jmax,kmax,lmax,j,k-ik,l-il,qmult,qskip,dim1,dim2,jstride,iqp)
vel    = qcons(idir+1)/qcons(1)
sigm   = spec(iqp)+abs(vel)
mcoeff = eps*(sigm+sig)
!                                                                                                                                                                                      
call fjdq(nq,idir,gm1,qcons,dqcons,timeMetric,mcoeff,adq_m)
!                                                                                                                                                                                                    
adq_m = -h0*adq_m

!>
! (j,k+ik,l+il)
!>
call get_qcons(q,  qcons,nq,nvar,jmax,kmax,lmax,j,k+ik,l+il,qmult,qskip,dim1,dim2,jstride,iqp)
call get_qcons(dq,dqcons,nq,nvar,jmax,kmax,lmax,j,k+ik,l+il,qmult,qskip,dim1,dim2,jstride,iqp)
vel    = qcons(idir+1)/qcons(1)
sigp   = spec(iqp) + abs(vel)
pcoeff = -eps*(sig+sigp)
!                                                                                                                                                                                                    
call fjdq(nq,idir,gm1,qcons,dqcons,timeMetric,pcoeff,adq_p)
!                                                                                                                                                                                                    
adq_p  = h0*adq_p
!
return
end subroutine rhs_gs

subroutine store_abc(nq,jmax,kmax,lmax,js,je,k,l,dfdq,mdim,a,b,c,h0,h0eps,sigma,dsigma,sflux,dflux,dim1,dim2,jstride)
!
implicit none
!
integer, intent(in) :: nq,jmax,kmax,lmax,mdim,js,je,k,l,dim1,dim2,jstride
real*8, dimension(nq,nq,mdim), intent(inout) :: dfdq,a,b,c
real*8, dimension(mdim), intent(in) :: sigma
real*8, dimension(jmax*kmax*lmax), intent(in) :: dsigma
real*8, dimension(nq,nq), intent(inout) :: dflux,sflux
real*8, intent(in) :: h0,h0eps
!
integer :: i,j,iq,iqp
real*8 :: acoeff,bcoeff,ccoeff
real*8 :: sigm,sig,sigp
!
iq=(l-1)*dim2+(k-1)*dim1+(js-1)*jstride
!
do j=js,je
   iqp = iq + 1
   
   sigp           = sigma(j+1)
   sig            = sigma(j)
   sigm           = sigma(j-1)
   
   acoeff =  h0eps*(sigm + sig)
   bcoeff =  dsigma(iqp)
   ccoeff =  h0eps*(sig + sigp)
   
   ! A                                                                                                                                                                                          
   sflux = -h0*dfdq(1:nq,1:nq,j-1)
   do i = 1,nq
      dflux(i,i) = acoeff
   enddo
   a(1:nq,1:nq,j) =  sflux-dflux
   
   ! B                                                                                                                                                                                          
   sflux = 0.0d0
   do i = 1,nq
      sflux(i,i) = bcoeff
   enddo
   b(1:nq,1:nq,j) = sflux
  
   ! C                                                                                                                                                                                          
   sflux = h0*dfdq(1:nq,1:nq,j+1)
   do i = 1,nq
      dflux(i,i) = ccoeff
   enddo
   c(1:nq,1:nq,j) =  sflux-dflux
   !                                                                                                                                                                                            
   iq = iq + jstride
enddo
!
return
end subroutine store_abc
