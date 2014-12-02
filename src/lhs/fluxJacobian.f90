!===================================================================!
!> \brief
!! This subroutine computes the second-order central difference 
!! flux Jacobian for conservative variables on a Cartesian grid system. 
!!
!! Versions:\par
!!    - Leffell 09/22/2014
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include flux_jacobian.f90
!!
!====================================================================!
subroutine fluxJacobian(nq,gm1,q,vel,rvel,p,flux_order,&
                        diss_order,efac,idir,dfdq)
!
implicit none
!
integer, intent(in)    :: nq                     !< number of field variables stored
real*8,  intent(in)    :: gm1                  !< ratio of specific heats-1
real*8,  intent(in)    :: q(nq)                  !< solution variables
real*8,  intent(in)    :: vel                    !< fluid velocity
real*8,  intent(in)    :: rvel                   !< relative velocity
real*8,  intent(in)    :: p                      !< fluid pressure
integer, intent(in)    :: flux_order             !< order of physical flux
integer, intent(in)    :: diss_order             !< order of dissipation
real*8,  intent(in)    :: efac                   !< scaling for dissipation
integer, intent(in)    :: idir                   !< coordinate direction index
real*8,  intent(out)   :: dfdq(nq,nq)            !< Flux Jacobian
!
! local variables
!
real*8  :: dp(nq)
real*8  :: rho,irho,ep,uu,u,v,w,e,drvel_dq1,drvel_dqip
integer :: ip,jp,dir
!
ip         = idir + 1
rho        = q(1)
irho       = 1d0/rho
dfdq       = 0d0

! Nonzero components of drvel/dq (first and idir+1th)
drvel_dq1  = -q(ip)*irho**2
drvel_dqip = irho

u          = q(2)*irho
v          = q(3)*irho
w          = q(4)*irho
e          = q(5)
uu         = u**2 + v**2 + w**2

dp = (/ 0.5*gm1*uu,-gm1*u,-gm1*v,-gm1*w,gm1/)

!
! Mass: can hardcode to dfdq(1,1) = -xdot, dfdq(1,ip) = 1d0
!
dfdq(1,1)  = rho*drvel_dq1 + rvel
!dfdq(1,ip) = q(1)*drvel_dqip  
dfdq(1,ip) = 1.0d0


!
! Momentum
!
do dir = 1,3
   jp = dir+1
   dfdq(jp,1)  = q(jp)*drvel_dq1
   dfdq(jp,ip) = q(jp)*drvel_dqip
   dfdq(jp,jp) = dfdq(jp,jp) + rvel
enddo
dfdq(ip,:) = dfdq(ip,:) + dp

!
! Energy
!
! dvel = drvel
ep         = e + p
dfdq(5,:)  = vel*dp
dfdq(5,1)  = dfdq(5,1)  + (ep)*drvel_dq1
dfdq(5,ip) = dfdq(5,ip) + (ep)*drvel_dqip
dfdq(5,5)  = dfdq(5,5)  +rvel


return
end subroutine fluxJacobian
