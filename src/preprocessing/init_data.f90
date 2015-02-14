!================================================================================
subroutine init_data(q,fsmach,alpha,beta,gamma,jmax,kmax,lmax,nq,istor)
!=================================================================================
! If the vortex is not used, then this subroutine
! initializes the data
!
implicit none
integer, intent(in) :: jmax,kmax,lmax,nq
real*8, intent(out) :: q(nq*jmax*kmax*lmax)
real*8, intent(in) :: fsmach
real*8, intent(in) :: alpha,beta,gamma
character*(*) :: istor
!
real*8 :: rho,p,u,v,w,e
!
integer :: j,k,l
integer :: qskip,qmult,iq,iloc,npts
!
! begin
!
! these are non-dimensional
! 
! i.e rho=RHO/rhoinf
!     [u,v,w]=[U,V,W]/a
!     e = E/(RHO*a^2)
!     p = P/pinf          
!
!
npts=jmax*kmax*lmax
call getstride(qskip,qmult,istor,npts,nq)
!
rho=1.0
u=fsmach*cos(alpha)*cos(beta)
v=fsmach*cos(alpha)*sin(beta)
w=fsmach*sin(alpha)
p=1./gamma
e=p/(gamma-1)+0.5*rho*fsmach**2
!
iq=0
!
do j=1,jmax
   do k=1,kmax
      do l=1,lmax
         iloc=iq*qmult+1
         q(iloc)=rho
         q(iloc+qskip)=rho*u
         q(iloc+2*qskip)=rho*v
         q(iloc+3*qskip)=rho*w
         q(iloc+4*qskip)=e
         iq=iq+1
      enddo
   enddo
enddo
!
return
end subroutine init_data
!
subroutine vortex_exact_cart(fsmach,gamma,q,x,js,je,ks,ke,ls,le,jmax,kmax,lmax,nvar,istor)
!========================================================================================================  
! This is the modified subroutine to initialize the lamb vortex (scully,vatistas formula)
!
implicit none
!
real*8, intent(in) :: fsmach,gamma
integer, intent(in) :: js,je,ks,ke,ls,le
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in)  :: x(3*jmax*kmax*lmax)
real*8, intent(out) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
! local variables
!
real*8 :: gm1,tinf,sinf,pi,cinf,rc,alpha,mu
real*8 :: rinf,uinf,vinf,winf,pinf,x0,y0,strnth0
real*8 :: z00,strnth,afac,bfac,scale,sigma
real*8 :: xbar,ybar,rsq,ee,eex,eey,u,v,w,v2,tpr,t,p,rho
integer :: j,k,l,n
integer :: qskip,qmult,iq,iloc,ilocx,xskip,xmult,npts
!
! begin
!
rinf=1.0
uinf=fsmach
vinf=0.0
winf=0.0
pinf=1./gamma
x0=0.0
y0=0.0
strnth0=1.
sigma=1.
!
gm1=gamma-1.
pi=4.0*atan(1.0)
!
cinf=sqrt(gamma*pinf/rinf)
tinf=pinf/rinf
sinf=pinf/rinf**gamma
!
z00=0.
strnth=strnth0*fsmach
afac=strnth/(2.*pi)
bfac=-0.5*gm1/gamma*afac**2
mu=1./10
alpha=1.25643
rc=sqrt(4*mu*0.1*alpha)
n=2
npts=jmax*kmax*lmax
!
call getstride(qskip,qmult,istor,npts,nvar)
call getstride(xskip,xmult,istor,npts,3)
!
!--- perturbation due to vortex ---
!
scale=1.0
do l=ls,le
   do k=ks,ke
      do j=js,je
         iq=(l-1)*jmax*kmax+(k-1)*jmax+(j-1)
         iloc=iq*qmult+1
         ilocx=iq*xmult+1
         xbar=x(ilocx)-x0
         ybar=x(ilocx+xskip)-y0
         rsq=xbar**2+ybar**2
         !ee=1/((rc**(2*n)+rsq**(n))**(1./n))    
         ee=exp(0.5*(1.-sigma*rsq))*scale
         !ee=1-exp(-rsq/rc**2)*scale
         u=uinf-(afac*ee)*ybar
         w=winf
         v=vinf+(afac*ee)*xbar
         v2=0.5*(u**2+v**2+w**2)
         !
         tpr=bfac*ee**2
         t=tinf+tpr
         p=(t**gamma/sinf)**(1./gm1)
         rho=p/t
         !
         q(iloc)=rho
         q(iloc+qskip)=rho*u
         q(iloc+2*qskip)=rho*v
         q(iloc+3*qskip)=rho*w
         q(iloc+4*qskip)=p/gm1+rho*v2
         !
      enddo
   enddo
enddo
!
return
end subroutine vortex_exact_cart
!
subroutine taylor_green(fsmach,gamma,q,x,js,je,ks,ke,ls,le,jmax,kmax,lmax,nvar,istor)
!========================================================================================================  
! This is the modified subroutine to initialize the lamb vortex (scully,vatistas formula)
!
implicit none
!
real*8, intent(in) :: fsmach,gamma
integer, intent(in) :: js,je,ks,ke,ls,le
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in)  :: x(3*jmax*kmax*lmax)
real*8, intent(out) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
! local variables
!
real*8 :: gm1,tinf,sinf,pi,cinf,rc,alpha,mu
real*8 :: rinf,uinf,vinf,winf,pinf,x0,y0,strnth0,z0
real*8 :: z00,strnth,afac,bfac,scale,sigma
real*8 :: xbar,ybar,zbar,rsq,ee,eex,eey,u,v,w,v2,tpr,t,p,rho
integer :: j,k,l,n
integer :: iloc,ilocx,iq,qskip,qmult
integer :: xskip,xmult,npts
!
! begin
!
rinf=1.0
uinf=fsmach
vinf=0.0
winf=0.0
pinf=1./gamma
x0=0.0
y0=0.0
z0=0.0
!
gm1=gamma-1.
!
npts=jmax*kmax*lmax
call getstride(qskip,qmult,istor,npts,nvar)
call getstride(xskip,xmult,istor,npts,3)
!
!
scale=1.0
do l=ls,le
   do k=ks,ke
      do j=js,je
         iq=(l-1)*jmax*kmax+(k-1)*jmax+(j-1)
         iloc=iq*qmult+1
         ilocx=iq*xmult+1
         xbar=x(ilocx)-x0
         ybar=x(ilocx+xskip)-y0
         zbar=x(ilocx+2*xskip)-z0
         u=uinf*sin(xbar)*cos(ybar)*cos(zbar)
         v=-uinf*cos(xbar)*sin(ybar)*cos(zbar)
         w=0d0
         p=pinf + (1d0/16d0)*rinf*uinf**2*(cos(2*xbar)+cos(2*ybar))*(cos(2*zbar)+2d0)
         v2=0.5*(u**2+v**2+w**2)
         rho=(gamma*p)
         q(iloc)=rho
         q(iloc+qskip)=rho*u
         q(iloc+2*qskip)=rho*v
         q(iloc+3*qskip)=rho*w
         q(iloc+4*qskip)=p/gm1+rho*v2
      enddo
   enddo
enddo
!
return
end subroutine taylor_green
!===============================================================================
! initialize metrics
!===============================================================================
subroutine init_metrics(timeMetric,dx,dy,dz,jmax,kmax,lmax,istor)
!
implicit none
integer, intent(in) :: jmax,kmax,lmax
real*8, intent(inout) :: timeMetric(3*jmax*kmax*lmax)
real*8, intent(in) :: dx,dy,dz
!
integer :: j,k,l
real*8 :: dxdy,dxdz,dydz,dxdydz
integer :: qskip,qmult,iq,iloc
integer :: tskip,tmult,npts
character*(*) :: istor
!
npts=jmax*kmax*lmax
call getstride(tskip,tmult,istor,npts,3)
!
iq=0
!
do l=1,lmax
   do k=1,kmax
      do j=1,jmax
         iloc=iq*qmult+1
         timeMetric(iloc)=0d0
         timeMetric(iloc+tskip)=0d0
         timeMetric(iloc+2*tskip)=0d0
      enddo
   enddo
enddo
return
end subroutine init_metrics
!>
!> get the storage order of the data
!>
subroutine getstride(qskip,qmult,istor,npts,nq)
implicit none
integer, intent(in) :: nq,npts
integer, intent(inout) :: qskip,qmult
character*(*) :: istor
if (istor=='row') then
   qskip=1
   qmult=nq
else
   qskip=npts
   qmult=1
endif
!
end subroutine getstride
