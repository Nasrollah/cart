program ifluxcheck
  real*8 :: q(5),q1(5)
  real*8 :: dq(5)
  real*8 :: flux(5),flux1(5)
  real*8 :: dflux(5),dflux1(5)
  real*8 :: gamma,gm1
  real*8 :: tm(3)
  real*8 :: eps=1e-5
  !
  gamma=1.4d0
  gm1=gamma-1d0
  tm=(/1d0,2d0,3d0/)
  q=(/1.d0,0.5d0,0.1d0,0.1d0,1d0/(gm1*gamma)+0.5d0*(0.1d0**2+0.1d0**2+0.5d0**2)/)
  dq=(/eps,eps,eps,eps,eps/)
  q1=q+dq
  !
  do idir=1,3
   call iflux(nq,idir,gm1,q,tm,flux)
   call iflux(nq,idir,gm1,q1,tm,flux1)
   dflux1=flux1-flux
   !
   call iflux_der(nq,idir,gm1,q,dq,tm,dflux)
   !
   write(6,"(5(1x,E15.7))") dflux-dflux1
  enddo
  !
end program ifluxcheck
!
subroutine iflux(nq,idir,gm1,q,tm,flux)
!
implicit none
integer, intent(in) :: nq
integer, intent(in) :: idir
real*8, intent(in) :: gm1
real*8, intent(in)  :: tm(3)
real*8, intent(in)  :: q(nq)
real*8, intent(out) :: flux(nq)
!
integer :: ip
real*8 :: vel,rvel,p
!
ip=idir+1
vel=q(ip)/q(1)
rvel=vel-tm(idir)
p=gm1*(q(5)-0.5d0*(q(2)*q(2)+q(3)*q(3)+q(4)*q(4))/q(1))
!
flux(1)=q(1)*rvel
flux(2)=q(2)*rvel
flux(3)=q(3)*rvel
flux(4)=q(4)*rvel
flux(5)=q(5)*rvel+p*vel
flux(ip)=flux(ip)+p
!
return
end subroutine iflux
!
subroutine iflux_der(nq,idir,gm1,q,dq,tm,dflux)
!
implicit none
integer, intent(in) :: nq
integer, intent(in) :: idir
real*8, intent(in) :: gm1
real*8, intent(in)  :: tm(3)
real*8, intent(in)  :: q(nq)
real*8, intent(in)  :: dq(nq)
real*8, intent(out) :: dflux(nq)
!
integer :: ip
real*8 :: drvel,dp,vel,p,rvel,u2
real*8 :: irho,irho2
!
ip=idir+1
irho=1d0/q(1)
irho2=irho*irho
vel=q(ip)*irho
rvel=vel-tm(idir)
u2=0.5d0*(q(2)*q(2)+q(3)*q(3)+q(4)*q(4))
p=gm1*(q(5)-u2*irho)
drvel=-q(ip)*irho2*dq(1)+dq(ip)*irho
dp=gm1*(dq(1)*u2*irho2&
     -dq(2)*q(2)*irho&
     -dq(3)*q(3)*irho&
     -dq(4)*q(4)*irho&
     +dq(5))
!
dflux(1)=dq(1)*rvel+q(1)*drvel
dflux(2)=dq(2)*rvel+q(2)*drvel
dflux(3)=dq(3)*rvel+q(3)*drvel
dflux(4)=dq(4)*rvel+q(4)*drvel
dflux(5)=dq(5)*rvel+q(5)*drvel+dp*vel+p*drvel
dflux(ip)=dflux(ip)+dp
!
return
end subroutine iflux_der
