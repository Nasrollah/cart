!===================================================================!
!> \brief
!! This subroutine computes turbulent kinetic energy, enstrophy and
!! dissipation rate for the taylor-green vortex problem
!! 
!! See case 3.3 documentation of the 3rd higher order workshop at
!! https://www.grc.nasa.gov/hiocfd/
!!
!! Versions:\par
!!    - Sitaraman 02/03/2015
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include compute_tke_params
!!
!====================================================================!
subroutine compute_tke_params(nf,t,fsmach,rey,nq,nvar,gamma,q,qstar,qwork,dx,dy,dz,jmax,kmax,lmax,&
                             flux_order,istor)
!
use mpi
use spatialCommunication, only : cartComm,myid_spatial,id,iprocs
implicit none
!
integer, intent(in) :: nf                        !< number of fringes
real*8, intent(in) :: t                          !< current time
real*8, intent(in) :: fsmach                     !< free stream Mach number
real*8, intent(in) ::  rey                       !< reynolds number (based on speed of sound)
integer, intent(in) :: nq                        !< number of field variables stored
integer, intent(in) :: nvar                      !< number of field variables to compute residuals for
real*8, intent(in) :: gamma                      !< ratio of specific heats
real*8,  intent(in) :: q(nq*jmax*kmax*lmax)      !< solution variables
real*8,  intent(inout) :: qwork(9*jmax*kmax*lmax) !< work array for storing velocity gradient
real*8,  intent(inout) :: qstar(nq*jmax*kmax*lmax) !< work array for storing velocity gradient
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
integer :: j,k,l,n,m,jj,kk,ierr
integer :: js,je,ks,ke,ls,le
real*8  :: gradu(3,3)
real*8 :: Lstar,omega,domega,V0,mu,rho0,pp,dv,fj,fk,fl
real*8, dimension(4) :: tke_params_local,tke_params,fac
real*8, dimension(3) :: vorticity
real*8 :: e1_local,e2_local,e3_local,enstrophy_local,tke_local
real*8 :: qprim(nvar)
integer :: jkmax,iq,iloc,qmult,qskip,ilocw
logical, save :: first_time=.true.
real*8 :: divv,bulk_strain
!
!
! begin
!
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
jkmax=jmax*kmax
!
ls=nf+1
le=lmax-nf
ks=nf+1
ke=kmax-nf
js=nf+1
je=jmax-nf
if (id(1)==iprocs(1)) je=je+1
if (id(2)==iprocs(2)) ke=ke+1
if (id(3)==iprocs(3)) le=le+1
!
! change variables to primitive and
! compute velocity gradients
!
call convertVariables(gamma,nq,jmax,kmax,lmax,q,0,istor)
call velocityGradient(dx,dy,dz,jmax,kmax,lmax,nq,q,qwork,flux_order,istor)
!
domega=dx*dy*dz
!
tke_local=0d0
enstrophy_local=0d0
e1_local=0d0
e2_local=0d0
e3_local=0d0
!
do l=ls,le
   fl=1d0
   if (l==ls .and. id(3)==1) fl=0.5
   if (l==le .and. id(3)==iprocs(3)) fl=0.5
   do k=ks,ke
      fk=1d0
      if (k==ks .and. id(2)==1) fk=0.5
      if (k==ke .and. id(2)==iprocs(2)) fk=0.5
      do j=js,je
         fj=1d0
         if (j==js .and. id(1)==1) fj=0.5
         if (j==je .and. id(1)==iprocs(1)) fj=0.5
         iq=(l-1)*jkmax+(k-1)*jmax+(j-1)
         iloc=iq*qmult+1
         !
         do n=1,nvar
            qprim(n)=q(iloc)
            iloc=iloc+qskip
         enddo
         !
         ilocw=iq*9
         m=0
         do kk=1,3
            do jj=1,3
               m=m+1
               gradu(jj,kk)=qwork(ilocw+m)
            enddo
         enddo
         !
         dv=fl*fk*fj*domega
         tke_local=tke_local+0.5d0*qprim(1)*(qprim(2)**2+qprim(3)**2+qprim(4)**2)*dv
         !tke_local=tke_local+dv
         !
         vorticity(1)=gradu(3,2)-gradu(2,3)
         vorticity(2)=gradu(1,3)-gradu(3,1)
         vorticity(3)=gradu(2,1)-gradu(1,2)
         !
         enstrophy_local=enstrophy_local+qprim(1)*0.5d0*dot_product(vorticity,vorticity)*dv
         !
         divv=gradu(1,1)+gradu(2,2)+gradu(3,3)
         bulk_strain=divv/3d0
         !
         ! norm of deviatoric part of the strain-rate tensor
         !
         e1_local=e1_local+(0.5d0*((gradu(1,2)+gradu(2,1))**2+&
                                  (gradu(1,3)+gradu(3,1))**2+&
                                  (gradu(2,3)+gradu(3,2))**2)+&
                                  (gradu(1,1)-bulk_strain)**2+&
                                  (gradu(2,2)-bulk_strain)**2+&
                                  (gradu(3,3)-bulk_strain)**2)*dv
                                      
         !
         ! bulk velocity contribution
         !
         e2_local=0d0
         !
         ! contribution from compressibility
         ! 
         pp=qprim(1)*qprim(5)/gamma
         e3_local=e3_local-pp*(gradu(1,1)+gradu(2,2)+gradu(3,3))*dv
         !
         iq=(l-1)*jkmax+(k-1)*jmax+(j-1)
         iloc=iq*qmult+1
         qstar(iloc)=0.5d0*qprim(1)*(qprim(2)**2+qprim(3)**2+qprim(4)**2) !-q1
         iloc=iloc+qskip
         qstar(iloc)=(0.5d0*((gradu(1,2)+gradu(2,1))**2+&
                                  (gradu(1,3)+gradu(3,1))**2+&
                                  (gradu(2,3)+gradu(3,2))**2)+&
                                  (gradu(1,1)-bulk_strain)**2+&
                                  (gradu(2,2)-bulk_strain)**2+&
                                  (gradu(3,3)-bulk_strain)**2)           !-q2
          iloc=iloc+qskip
          qstar(iloc)=pp*(gradu(1,1)+gradu(2,2)+gradu(3,3))              !-q3
          iloc=iloc+qskip
          qstar(iloc)=divv                                               !-q4
          iloc=iloc+qskip
          qstar(iloc)=qprim(1)*0.5d0*dot_product(vorticity,vorticity)    !-q5
      enddo
   enddo
enddo
!
!write(6,*) 'tke_local=',tke_local,js,je,ks,ke,ls,le
!
call convertVariables(gamma,nq,jmax,kmax,lmax,q,1,istor)
!
!tke_local=tke_local*domega
!enstrophy_local=enstrophy_local*domega
!e1_local=e1_local*domega
!e3_local=e3_local*domega
tke_params_local=(/tke_local,enstrophy_local,e1_local,e3_local/)
!
call mpi_reduce(tke_params_local,tke_params,4,MPI_DOUBLE_PRECISION,MPI_SUM,0,cartComm,ierr)
!
if (myid_spatial==0) then
   Lstar=1d0
   omega=(2d0*acos(-1d0)*Lstar)**3  !< volume of the whole simulation
   V0=fsmach
   mu=1d0/rey
   rho0=1d0
   fac=(/ V0*V0, V0*V0/Lstar/Lstar, V0**3/Lstar, V0**3/Lstar /)
   tke_params(1)=tke_params(1)/(rho0*omega)/fac(1)
   tke_params(2)=tke_params(2)/(rho0*omega)/fac(2)
   tke_params(3)=tke_params(3)*(2d0*mu/rho0/omega)/fac(3)
   tke_params(4)=tke_params(4)/(rho0*omega)/fac(4)
   if (first_time) then
      open(unit=1040,file='tke_params.dat',form='formatted')
      first_time=.false.
   endif
   write(1040,"(6(1x,E15.7))") t*V0/Lstar,tke_params,tke_params(3)+tke_params(4)
   call flush(1040)
endif
!
return
end subroutine compute_tke_params
