subroutine fjdq(nq,idir,gm1,q,dq,tm,eps,dflux)
!
implicit none
integer, intent(in)  :: nq
integer, intent(in)  :: idir
real*8,  intent(in)  :: gm1
real*8,  intent(in)  :: tm(3)
real*8,  intent(in)  :: q(nq)
real*8,  intent(in)  :: dq(nq)
real*8,  intent(in)  :: eps
real*8,  intent(out) :: dflux(nq)
!
integer :: ip
real*8 :: drvel,vel,rvel,u2
real*8 :: irho,irho2,p,dp
!
ip=idir+1
irho=1d0/q(1)
irho2=irho*irho
vel=q(ip)*irho
rvel=vel-tm(idir)
u2=0.5d0*(q(2)*q(2)+q(3)*q(3)+q(4)*q(4))
p=gm1*(q(5)-u2*irho)
drvel=-q(ip)*irho2*dq(1)+dq(ip)*irho

! store dp(j,k,l) and send it in like p. if sending in p get rid of u2
dp=gm1*(dq(1)*u2*irho2&
     -dq(2)*q(2)*irho&
     -dq(3)*q(3)*irho&
     -dq(4)*q(4)*irho&
     +dq(5))
!
!dflux(1)=dq(1)*rvel+q(1)*drvel 
!dflux(2)=dq(2)*rvel+q(2)*drvel
!dflux(3)=dq(3)*rvel+q(3)*drvel
!dflux(4)=dq(4)*rvel+q(4)*drvel
!dflux(5)=dq(5)*rvel+q(5)*drvel+dp*vel+p*drvel
!dflux(ip)=dflux(ip)+dp

dflux(1)=dq(1)*(rvel+eps)+q(1)*drvel
dflux(2)=dq(2)*(rvel+eps)+q(2)*drvel
dflux(3)=dq(3)*(rvel+eps)+q(3)*drvel
dflux(4)=dq(4)*(rvel+eps)+q(4)*drvel
dflux(5)=dq(5)*(rvel+eps)+q(5)*drvel+dp*vel+p*drvel
dflux(ip)=dflux(ip)+dp
!
return
end subroutine fjdq
