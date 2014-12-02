!>
!> This subroutine gives the 2nd order boundary conditions using
!> extrapolation from interior for five sides
!> and creates a Gaussian isentropic velocity perturbation on the inlet
!>
subroutine bc_per(q,x,fsmach,gamma,amp,freq,t,jmax,kmax,lmax)
!
implicit none
!
integer, intent(in) :: jmax,kmax,lmax
real*8, intent(inout) :: q(5,jmax,kmax,lmax)
real*8, intent(in)    :: x(3,jmax,kmax,lmax)
real*8, intent(in) :: fsmach,gamma,amp,freq,t
!
! local variables
!
integer :: j,k,l,j1
real*8 :: u,v,w,u2,p,TT,rho,p0,T0,xx,yy,zz,spread,gaussian
!
T0=(1+0.5*(gamma-1)*fsmach**2)
p0=T0**(gamma/(gamma-1))/gamma
spread=0.3
!
!
do j=1,2
 do k=1,kmax
   do l=1,lmax

      xx=x(1,j,k,l)
      yy=x(2,j,k,l)
      zz=x(3,j,k,l)
      gaussian = exp(-spread*(yy*yy+zz*zz))
      !gaussian=0.
      rho=1.0
      u=fsmach+amp*sin(freq*t)*gaussian
      v=0.
      w=0.
      u2=u*u+v*v+w*w
      j1=3
      p=(gamma-1)*(q(5,j1,k,l)-0.5*(q(2,j1,k,l)**2+q(3,j1,k,l)**2+q(4,j1,k,l)**2)/q(1,j1,k,l))
      !p=p0/(1+0.5*(gamma-1)*u2)**(gamma/(gamma-1))
      !TT=T0/(1+0.5*(gamma-1)*u2)
      !rho=gamma*p/TT

      q(1,j,k,l)=rho
      q(2,j,k,l)=rho*u
      q(3,j,k,l)=rho*v
      q(4,j,k,l)=rho*w
      q(5,j,k,l)=p/(gamma-1)+0.5*rho*u2
   enddo
 enddo
enddo
!
return
j=jmax
!
do k=1,kmax
   do l=1,lmax
      q(:,j,k,l)=q(:,j-1,k,l)
   enddo
enddo
!
k=1
!
do j=1,jmax
   do l=1,lmax
      q(:,j,k,l)=q(:,j,k+1,l)
   enddo
enddo
!
k=kmax
!
do j=1,jmax
   do l=1,lmax
      q(:,j,k,l)=q(:,j,k-1,l)
   enddo
enddo
!
l=1
do j=1,jmax
   do k=1,kmax
      q(:,j,k,l)=q(:,j,k,l+1)
   enddo
enddo
!
l=lmax
do j=1,jmax
   do k=1,kmax
      q(:,j,k,l)=q(:,j,k,l-1)
   enddo
enddo
!
return
end subroutine bc_per
!>
!> just a density perturbation
!>
subroutine bc_per1(q,x,fsmach,gamma,amp,freq,t,jmax,kmax,lmax,bctyp,fluxOrder,dissOrder)
!
implicit none
!
integer, intent(in) :: jmax,kmax,lmax,bctyp,fluxOrder,dissOrder
real*8, intent(inout) :: q(5,jmax,kmax,lmax)
real*8, intent(in)    :: x(3,jmax,kmax,lmax)
real*8, intent(in) :: fsmach,gamma,amp,freq,t
!
! local variables
!
integer :: j,k,l,j1,js,je
real*8 :: u,v,w,u2,p,TT,rho,p0,T0,xx,yy,zz,spread,gaussian,kernel
!
T0=(1+0.5*(gamma-1)*fsmach**2)
p0=T0**(gamma/(gamma-1))/gamma
spread=0.3

js = 1
je = 2

if (fluxOrder.eq.6) then
   je = 3
endif

!
do j=js,je
 do k=1,kmax
   do l=1,lmax

      xx=x(1,j,k,l)
      yy=x(2,j,k,l)
      zz=x(3,j,k,l)
      gaussian = exp(-spread*(yy*yy+zz*zz))
      if (bctyp==2) then
         kernel = 2./5.+4*cos(freq*t)-1.
      else
         kernel = sin(freq*t)
      endif

      rho=1.+amp*kernel*gaussian
      !gaussian=0.
      u=fsmach
      v=0.
      w=0.
      u2=u*u+v*v+w*w
      j1=3
      p=(gamma-1)*(q(5,j1,k,l)-0.5*(q(2,j1,k,l)**2+q(3,j1,k,l)**2+q(4,j1,k,l)**2)/q(1,j1,k,l))

      !p=p0/(1+0.5*(gamma-1)*u2)**(gamma/(gamma-1))
      !TT=T0/(1+0.5*(gamma-1)*u2)
      !rho=gamma*p/TT

      q(1,j,k,l)=rho
      q(2,j,k,l)=rho*u
      q(3,j,k,l)=rho*v
      q(4,j,k,l)=rho*w
      q(5,j,k,l)=p/(gamma-1)+0.5*rho*u2
   enddo
 enddo
enddo
!
!
do j=jmax-2,jmax
 do k=1,kmax
   do l=1,lmax
      q(1:4,j,k,l)=q(1:4,jmax-3,k,l)
      p=1./gamma
      q(5,j,k,l)=p/(gamma-1)+0.5*(q(2,j,k,l)**2+q(3,j,k,l)**2+q(4,j,k,l)**2)/q(1,j,k,l)
   enddo
 enddo
enddo
!
return
!
k=1
!
do j=1,jmax
   do l=1,lmax
      q(:,j,k,l)=q(:,j,k+1,l)
   enddo
enddo
!
k=kmax
!
do j=1,jmax
   do l=1,lmax
      q(:,j,k,l)=q(:,j,k-1,l)
   enddo
enddo
!
l=1
do j=1,jmax
   do k=1,kmax
      q(:,j,k,l)=q(:,j,k,l+1)
   enddo
enddo
!
l=lmax
do j=1,jmax
   do k=1,kmax
      q(:,j,k,l)=q(:,j,k,l-1)
   enddo
enddo
!
return
end subroutine bc_per1
