!>
!> RHS interface to samrai
!> same interface as ARC3D
!>
subroutine rhs(nq,jd,kd,ld,jk2d,fringe,                   &
   js,je,ks,ke,ls,le,q,s,iblank,                          &
   dxyz,gamma,dis2,dis4,ad_fso,ad_smo,                    &
   FSO,irhs,resid,x,y,z,xc,yc,zc,                         &
   ihover,omegaz,imaneuver,gridvelx,gridvely,gridvelz,    &
   refPrndl,refRey,refMach,ivis,iturb,qp,vort_mag,mut,    &
   tjac, turb_resid)
  !
  implicit none
  !
  integer, intent (in) :: nq, jd,kd,ld,irhs,js,je,ks,ke,ls,le,fringe
  integer, intent (in) :: ivis, iturb
  real*8, intent(in) :: dxyz,gamma,fso
  real*8,  INTENT (IN) :: DIS2,DIS4,ad_fso,ad_smo
  real*8, intent(in) :: refPrndl,refRey,refMach
  real*8, intent(inout) :: resid,turb_resid
  logical, intent(in) :: jk2D
  INTEGER, DIMENSION(JD,KD,LD), INTENT (IN) :: IBLANK
  !
  real*8, dimension(jd,kd,ld,nq), intent(in) :: q
  real*8, dimension(jd,kd,ld), intent(in) :: x,y,z
  real*8, dimension(jd,kd,ld,nq), intent(inout) :: s
  real*8,dimension(jd,kd,ld),intent(inout)        :: tjac
  real*8, dimension(jd,kd,ld,nq,7), intent(inout) :: qp
  real*8, dimension(jd,kd,ld),intent(inout)     :: vort_mag
  real*8, dimension(jd,kd,ld),intent(inout)       :: mut
  integer :: ihover,imaneuver
  real*8 :: xc,yc,zc,omegaz,gridvelx,gridvely,gridvelz
  !
  real*8 :: dx,dy,dz
  real*8 :: rey,pr,prtr
  real*8 :: tm(3)
  integer :: fluxOrder,dissOrder,viscOrder
  real*8 :: dissCoef
  integer :: ibval
  integer :: j,k,l,n
  character*10::istor
  integer :: nvar
  !
  nvar=5
  istor='column'
  dx=1d0
  dy=1d0
  dz=1d0
  tm=(/gridvelx,gridvely,gridvelz/)
  fluxOrder=6
  dissOrder=6
  viscOrder=4
  dissCoef=0.5d0
  rey=refRey/refMach
  pr=refPrndl
  prtr=refPrndl
  !
  call inviscidRHSunified(nq,nvar,gamma,q,s,qp,tm,dx,dy,dz,jd,kd,ld,&
           fluxOrder,dissOrder,dissCoef,istor)
  !
  if (ivis==1) then
     call viscousRHS(rey,pr,prtr,nq,nvar,gamma,q,qp,s,dx,dy,dz,jd,kd,ld,&
          min(4,viscOrder),istor)
  endif
  !
  do l=1,ld
   do k=1,kd
    do j=1,jd
     ibval=max(iblank(j,k,l),0)
     do n=1,5
      s(j,k,l,n)=-s(j,k,l,n)*ibval
     enddo
    enddo
   enddo
  enddo
  !
  call compute_norm(resid,s,jd,kd,ld,nq,istor)
  !
  return
end subroutine rhs
!>
!> invert the linear system using the appropriate 
!> method
!>
subroutine lhs(nq,jd,kd,ld,js,je,ks,ke,ls,le,q,s,iblank,dxyz,&
     gamma,gridvelx,gridvely,gridvelz,refprndl,refrey,&
     refmach,dt,iorder,ivis,iturb)
!
implicit none
!
integer, intent(in) :: nq
integer, intent(in) :: jd,kd,ld
integer, intent(in) :: js,je,ks,ke,ls,le
real*8, intent(inout) :: q(jd,kd,ld,nq)
real*8, intent(inout) :: s(jd,kd,ld,nq)
integer, intent(in) :: iblank(jd,kd,ld)
real*8, intent(in) :: dxyz
real*8, intent(in) :: gamma
real*8, intent(in) :: gridvelx,gridvely,gridvelz
real*8, intent(in) :: refprndl
real*8, intent(in) :: refrey
real*8, intent(in) :: refmach
real*8, intent(in) :: dt
integer, intent(in) :: iorder
integer, intent(in) :: iturb
integer, intent(in) :: ivis
!
integer :: nvar
real*8 :: rey,pr,prtr
real*8 :: tscal
real*8, allocatable :: tm(:,:,:,:)
real*8 :: ibval
real*8 :: dx,dy,dz,facearea,voli
integer :: j,k,l,n
!
nvar=5
dx=dxyz
dy=dxyz
dz=dxyz
faceArea=dx*dy
voli=1d0/(dx*dy*dz)
rey=refRey/refMach
pr=refPrndl
prtr=refPrndl
tscal=dt
allocate(tm(jd,kd,ld,3))
do l=1,ld
 do k=1,kd
   do j=1,jd
    tm(j,k,l,1)=gridvelx
    tm(j,k,l,2)=gridvely
    tm(j,k,l,3)=gridvelz
  enddo
 enddo
enddo
!
if (iorder==2) tscal=2d0/3d0*dt
!
!write(6,*) 'calling lusgs'
call lusgs_hyb_col(gamma,rey,q,s,tscal,tm,dx,dy,dz,nq,nvar,jd,kd,ld)
!
do l=ls,le
   do k=ks,ke
      do j=js,je
         ibval=max(iblank(j,k,l),0)*voli
         do n=1,5
	    s(j,k,l,n)=s(j,k,l,n)*ibval
            q(j,k,l,n)=q(j,k,l,n)+s(j,k,l,n)
         enddo
      enddo
   enddo
enddo
!
deallocate(tm)
!
return
end subroutine lhs
!
!> Add the bdf source terms for the appropriate
!> implicit method
!>
subroutine bdfsource(nq,nvar,pn,jd,kd,ld,js,je,ks,ke,ls,le,&
     dt,dxyz,iorder,q,qn,qnm1,s)
!
implicit none
!
integer, intent(in) :: nq,nvar,jd,kd,ld,js,je,ks,ke,ls,le,pn
real*8, intent(in)  :: dt
real*8, intent(in)  :: dxyz
integer, intent(in) :: iorder
real*8, intent(in) :: q(jd,kd,ld,nq)
real*8, intent(in) :: qn(jd,kd,ld,nq)
real*8, intent(in) :: qnm1(jd,kd,ld,nq)
real*8, intent(inout) :: s(jd,kd,ld,nq)
!
!>
!
real*8 :: facearea,dtvol
integer :: j,k,l,n
real*8 :: resid
!
facearea=dxyz*dxyz
!
!if (pn==0) then
! call compute_norm(resid,s,jd,kd,ld,nq,'column')
! write(6,*) pn,jd,kd,ld,resid,iorder
! write(6,*) q(27,38,15,:)
! write(6,*) qn(27,38,15,:)
!endif
!
if (iorder==0) then
 do l=ls,le
   do k=ks,ke
     do j=js,je
      do n=1,nvar
       s(j,k,l,n)=-facearea*s(j,k,l,n)
      enddo
     enddo
   enddo
  enddo
else if (iorder==1) then
   dtvol=dxyz*dxyz*dxyz/dt
   do l=ls,le
      do k=ks,ke
         do j=js,je
            do n=1,nvar
               s(j,k,l,n)=-facearea*s(j,k,l,n)-(q(j,k,l,n)-qn(j,k,l,n))*dtvol
            enddo
         enddo
      enddo
   enddo
else
   dtvol=0.5*dxyz*dxyz*dxyz/dt
   do l=ls,le
      do k=ks,ke
         do j=js,je
            do n=1,nvar
               s(j,k,l,n)=-facearea*s(j,k,l,n)&
                    -(3d0*q(j,k,l,n)-4d0*qn(j,k,l,n)+qnm1(j,k,l,n))*dtvol
            enddo
         enddo
      enddo
   enddo
endif
!
!if (pn==0) then
! call compute_inf_norm(resid,s,jd,kd,ld,nq,'column')
! write(6,*) pn,jd,kd,ld,resid
!endif
!
return
end subroutine bdfsource


