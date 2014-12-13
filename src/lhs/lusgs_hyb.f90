!>
!!
!! LU-SGS scheme 
!! Calculate the implicit inversion of the LHS
!! This involves two bidiagonal scalar inversions
!! 
subroutine lusgs_hyb(gamma,rey,q,s,tscale,timeMetric,dx,dy,dz,nq,nvar,jmax,kmax,lmax)
!
!
implicit none
!
! subroutine arguments
!
real*8, intent(in) :: gamma
real*8, intent(in) :: rey
integer, intent(in) :: jmax,kmax,lmax,nvar,nq
real*8, intent(in) :: q(nq,jmax,kmax,lmax)
real*8, intent(inout) :: s(nq,jmax,kmax,lmax)
real*8, intent(in) :: timeMetric(3,jmax,kmax,lmax)
real*8, intent(in) :: dx,dy,dz
real*8, intent(in) :: tscale
!
! local variables
!
integer, allocatable :: ms(:),me(:)
real*8, allocatable :: uv(:,:),vn(:,:,:),ge(:,:),qx(:,:),qz(:,:)
real*8, allocatable :: cx(:,:),cz(:,:),a(:,:),c(:,:)
real*8, allocatable :: b(:,:,:),d(:,:,:)
!
integer :: js,je,ks,ke,ls,le
integer :: i,j,k,l,m,n,j1,k1,mdim,kr,l1
real*8  :: ggm1,gm1
real*8  :: oat
real*8  :: eps2,epsv
real*8,save :: epse=0.01
real*8 :: s1,s2,s3,s4,s5
real*8 :: uu,vv,ww
real*8 :: uvw,cjkl
real*8 :: rr1,rr2,rr3
real*8 :: dxi,dyi,dzi
real*8 :: qq1,qq2,qq3,qqx,qqy,qqz
real*8 :: ug,vg,wg
real*8 :: qq,cc,sp1,sm2,chkx,chky,chkz,spec
real*8 :: a2,a3,a5,svt,er
real*8 :: rhoi,volumeinv
real*8 :: tsmax,tsmin
real*8 :: vscale
real*8 :: vnu
integer :: jd1,kd1,ld1,msweep
!
integer, parameter :: ixx=1
integer, parameter :: ixy=2
integer, parameter :: ixz=3
integer, parameter :: iyx=4
integer, parameter :: iyy=5
integer, parameter :: iyz=6
integer, parameter :: izx=7
integer, parameter :: izy=8
integer, parameter :: izz=9
integer, parameter :: ivol=10
!
! begin
!
eps2=epse*2.
epsv=1.+eps2 ! this is a fudge to account for linearization errors
gm1=gamma-1
ggm1=gamma*gm1
!
tsmax=-1e10
tsmin=1e10
!
ks=2
ke=kmax-1
js=2
je=jmax-1
ls=2
le=lmax-1
do l=1,lmax,lmax-1
  do k=1,kmax
   do j=1,jmax
    s(:,j,k,l)=0.
   enddo
  enddo
enddo
!
do k=1,kmax,kmax-1
 do j=1,jmax
  do l=1,lmax
   s(:,j,k,l)=0.
  enddo
 enddo
enddo
!
do j=1,jmax,jmax-1
  do k=1,kmax
   do l=1,lmax
    s(:,j,k,l)=0.
   enddo
  enddo
enddo
!
mdim=max(jmax,kmax)
mdim=max(mdim,lmax)
!
! allocate arrays
!
allocate(ms(je+le),me(je+le))
allocate(uv(jmax,lmax),ge(jmax,lmax),qx(jmax,lmax),qz(jmax,lmax))
allocate(cx(jmax,lmax),cz(jmax,lmax),d(jmax,kmax,lmax))
allocate(vn(3,jmax,lmax))
allocate(a(5,mdim),b(5,jmax,lmax),c(5,mdim))
!
! set up for hyper plane loop
!
do m=js+ls,je+le
   ms(m)=max(js,m-le)
   me(m)=min(je,m-js)
enddo
!
!
! set time accuracy
!
dxi=1d0/dx
dyi=1d0/dy
dzi=1d0/dz
!
! Setup D, diagonal term and store some arrays
!
do msweep=1,1
do k=ks,ke
   do l=ls-1,le+1
      do j=js-1,je+1
         rhoi=1./q(1,j,k,l)
         uu=q(2,j,k,l)*rhoi
         vv=q(3,j,k,l)*rhoi
         ww=q(4,j,k,l)*rhoi
         ug=timeMetric(1,j,k,l)
         vg=timeMetric(2,j,k,l)
         wg=timeMetric(3,j,k,l)
         uvw=0.5*(uu*uu+vv*vv + ww*ww)
         cjkl=sqrt(ggm1*(q(5,j,k,l)*rhoi-uvw))
         !
         qq1=uu*dxi
         qqx=abs((uu-ug)*dxi)
         rr1=dxi
         !
         qq2=vv*dyi
         qqy=abs((vv-vg)*dyi)
         rr2=dyi
         !
         qq3=ww*dzi
         qqz=abs((ww-wg)*dzi)
         rr3=dzi         
         !
         svt=tscale
         if (nq.gt.5) then
            vscale=(1d0+q(6,j,k,l))/(rey*q(1,j,k,l))
         else
            vscale=1d0/(rey*q(1,j,k,l))
         endif
         !
         ! diagonal contribution of the direct viscous terms
         ! no cross terms and off-diagonals for now
         !
         vn(1,j,l)=rr1*rr1*vscale
         vn(2,j,l)=rr2*rr2*vscale
         vn(3,j,l)=rr3*rr3*vscale
         vnu=2d0*(vn(1,j,l)+vn(2,j,l)+vn(3,j,l))
         !
         d(j,k,l)=1./(1+svt*(qqx+qqy+qqz+cjkl*(rr1+rr2+rr3))*epsv+vnu)
         uv(j,l)=uvw
         qx(j,l)=qq1
         qz(j,l)=qq3
         cx(j,l)=cjkl*rr1
         cz(j,l)=cjkl*rr3
         ge(j,l)=gamma*q(5,j,k,l)*rhoi-(gamma-1)*uvw
	 if (msweep.eq.1) s(1:nvar,j,k,l)=s(:,j,k,l)*svt
      enddo
   enddo
   kr=k-1
   do l=ls-1,le+1
      do j=js-1,je+1
         rhoi=q(1,j,kr,l)
         uu=q(2,j,kr,l)*rhoi
         vv=q(3,j,kr,l)*rhoi
         ww=q(4,j,kr,l)*rhoi         
         er=q(5,j,kr,l)*rhoi
         uvw=0.5*(uu*uu+vv*vv+ww*ww)
         !
         ug=timeMetric(1,j,kr,l)
         vg=timeMetric(2,j,kr,l)
         wg=timeMetric(3,j,kr,l)
         !
         rr2=dyi*dyi
         qq=vv*dyi
         cc=sqrt(ggm1*(er-uvw)*rr2)
         qqy=-dyi*vg+qq
         sp1=abs(qqy)+cc
         sm2=eps2*sp1
         chky=0.5+sign(0.5,qqy+cc)
         spec=chky*(qqy+sp1)+sm2+vn(2,j,l)
         !
         s1=s(1,j,kr,l)
         s2=s(2,j,kr,l)
         s3=s(3,j,kr,l)
         s4=s(4,j,kr,l)
         s5=s(5,j,kr,l)
         !
         a5=chky*(dyi*s3-qq*s1)
         a2=chky*gm1*(uvw*s1-(uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (B-sigma*I)*dQ(:,j,k-1,l)
         !
         b(1,j,l)=a5 + spec*s1
         b(2,j,l)=         uu*a5 + spec*s2
         b(3,j,l)=dyi*a2 + vv*a5 + spec*s3
         b(4,j,l)=         ww*a5 + spec*s4
         b(5,j,l)=qq*a2 + (gamma*er - gm1*uvw)*a5 + spec*s5
         !
      enddo
   enddo
   !
   ! hyper-plane forward sweep
   !
   do m=js+ls,je+le
      do j=ms(m),me(m)
         l=m-j
         j1=j-1
         rhoi=q(1,j1,k,l)
         uu=q(2,j1,k,l)*rhoi
         vv=q(3,j1,k,l)*rhoi
         ww=q(4,j1,k,l)*rhoi
         uvw=uv(j1,l)
         !
         qq=qx(j1,l)
         cc=cx(j1,l)
         ug=timeMetric(1,j1,k,l)
         vg=timeMetric(2,j1,k,l)
         wg=timeMetric(3,j1,k,l)
         qqx=-dxi*ug+qq
         sp1=abs(qqx)+cc
         sm2=eps2*sp1
         chkx=0.5 + sign(0.5,qqx+cc)
         spec=chkx*(qqx+sp1)+sm2+vn(1,j1,l)
         !
         s1=s(1,j1,k,l)
         s2=s(2,j1,k,l)
         s3=s(3,j1,k,l)
         s4=s(4,j1,k,l)
         s5=s(5,j1,k,l)
         !
         a5=chkx*(dxi*s2 - qq*s1)
         a2=chkx*gm1*(uvw*s1 - (uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (A-sigma*I)*dQ(:,j-1,k,l)
         !
         a(1,j) = a5 +spec*s1
         a(2,j) = dxi*a2 + uu*a5 + spec*s2
         a(3,j) =          vv*a5 + spec*s3
         a(4,j) =           ww*a5 + spec*s4
         a(5,j) = qq*a2 + ge(j1,l)*a5 + spec*s5
      enddo

      do j=ms(m),me(m)
         l=m-j
         l1=l-1
         !
         rhoi=q(1,j,k,l1)
         uu=q(2,j,k,l1)*rhoi
         vv=q(3,j,k,l1)*rhoi
         ww=q(4,j,k,l1)*rhoi
         uvw=uv(j,l1)
         !
         ug=timeMetric(1,j,k,l1)
         vg=timeMetric(2,j,k,l1)
         wg=timeMetric(3,j,k,l1)
         qq=qz(j,l1)
         cc=cz(j,l1)
         !
         qqz=-dzi*wg+qq
         sp1= abs(qqz)+ cc
         sm2= eps2*sp1
         chkz = 0.5 + sign(0.5, qqz+cc)
         spec=chkz*(qqz+sp1) + sm2 + vn(3,j,l1)
         !
         s1=s(1,j,k,l1)
         s2=s(2,j,k,l1)
         s3=s(3,j,k,l1)
         s4=s(4,j,k,l1)
         s5=s(5,j,k,l1)
         !
         a5=chkz*(dzi*s4 - qq*s1)
         a2=chkz*gm1*(uvw*s1 - (uu*s2+vv*s3+ww*s4) + s5)
         !
         ! (C-sigma*I)*dQ(:,j,k,l-1)
         !
         c(1,j) = a5 + spec*s1
         c(2,j) =          uu*a5 + spec*s2
         c(3,j) =          vv*a5 + spec*s3
         c(4,j) = dzi*a2 + ww*a5 + spec*s4
         c(5,j) = qq*a2 + ge(j,l1)*a5 + spec*s5
         !
      enddo
   !
      do j=ms(m),me(m)
       l=m-j
       svt=tscale
       !
       ! factor of half because its really (A+sigma*I)/2 and (C+sigma*I)/2 and (B+sigma*I)/2
       !
       s(1:nvar,j,k,l)=(s(1:nvar,j,k,l) + svt*0.5*(a(:,j)+b(:,j,l)+c(:,j)))*d(j,k,l)
       !
     enddo
   enddo
enddo
!
! backward sweep
!
do k=ke,ks,-1
   do l=ls-1,le+1
      do j=js-1,je+1
         rhoi=1./q(1,j,k,l)
         uu=q(2,j,k,l)*rhoi
         vv=q(3,j,k,l)*rhoi
         ww=q(4,j,k,l)*rhoi
         ug=timeMetric(1,j,k,l)
         vg=timeMetric(2,j,k,l)
         wg=timeMetric(3,j,k,l)
         uvw=0.5*(uu*uu+vv*vv + ww*ww)
         cjkl=sqrt(ggm1*(q(5,j,k,l)*rhoi-uvw))
         !
         qq1=dxi*uu
         qqx=abs(dxi*(uu-ug))
         rr1=dxi
         !
         qq3=dzi*ww
         qqz=abs(dzi*(ww-wg))
         rr3=dzi
         !
         uv(j,l)=uvw
         qx(j,l)=qq1
         qz(j,l)=qq3
         cx(j,l)=cjkl*rr1
         cz(j,l)=cjkl*rr3
         ge(j,l)=gamma*q(5,j,k,l)*rhoi-(gamma-1)*uvw
      enddo
   enddo
   !
   kr=k+1
   !
   do l=ls-1,le+1
      do j=js-1,je+1
         rhoi=q(1,j,kr,l)
         uu=q(2,j,kr,l)*rhoi
         vv=q(3,j,kr,l)*rhoi
         ww=q(4,j,kr,l)*rhoi         
         er=q(5,j,kr,l)*rhoi
         uvw=0.5*(uu*uu+vv*vv+ww*ww)
         !
         ug=timeMetric(1,j,kr,l)
         vg=timeMetric(2,j,kr,l)
         wg=timeMetric(3,j,kr,l)
         !
         rr2=dyi*dyi
         qq=dyi*vv
         cc=sqrt(ggm1*(er-uvw)*rr2)
         qqy=-dyi*vg+qq
         sp1=abs(qqy)+cc
         sm2=eps2*sp1
         chky=0.5-sign(0.5,qqy-cc)
         spec=chky*(qqy-sp1)-sm2 - vn(2,j,l)
         !
         s1=s(1,j,kr,l)
         s2=s(2,j,kr,l)
         s3=s(3,j,kr,l)
         s4=s(4,j,kr,l)
         s5=s(5,j,kr,l)
         !
         a5=chky*(dyi*s3-qq*s1)
         a2=chky*gm1*(uvw*s1-(uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (B+sigma*I)*dQ(:,j,k+1,l)
         !
         b(1,j,l)=a5 + spec*s1
         b(2,j,l)=         uu*a5 + spec*s2
         b(3,j,l)=dyi*a2 + vv*a5 + spec*s3
         b(4,j,l)=         ww*a5 + spec*s4
         b(5,j,l)=qq*a2 + (gamma*er - gm1*uvw)*a5 + spec*s5
         !
      enddo
   enddo
   !
   do m=je+le,js+ls,-1
      do j=ms(m),me(m)
         l=m-j
         j1=j+1
         rhoi=q(1,j1,k,l)
         uu=q(2,j1,k,l)*rhoi
         vv=q(3,j1,k,l)*rhoi
         ww=q(4,j1,k,l)*rhoi
         uvw=uv(j1,l)
         !
         qq=qx(j1,l)
         cc=cx(j1,l)
         ug=timeMetric(1,j1,k,l)
         vg=timeMetric(2,j1,k,l)
         wg=timeMetric(3,j1,k,l)
         qqx=-dxi*ug + qq
         sp1=abs(qqx)+cc
         sm2=eps2*sp1
         chkx=0.5 - sign(0.5,qqx-cc)
         spec=chkx*(qqx-sp1)-sm2 - vn(1,j1,l)
         !
         s1=s(1,j1,k,l)
         s2=s(2,j1,k,l)
         s3=s(3,j1,k,l)
         s4=s(4,j1,k,l)
         s5=s(5,j1,k,l)
         !
         a5=chkx*(dxi*s2 - qq*s1)
         a2=chkx*gm1*(uvw*s1 - (uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (A+sigma*I)*dQ(:,j+1,k,l)
         !
         a(1,j) = a5 +spec*s1
         a(2,j) = dxi*a2 + uu*a5 + spec*s2
         a(3,j) =          vv*a5 + spec*s3
         a(4,j) =          ww*a5 + spec*s4
         a(5,j) = qq*a2 + ge(j1,l)*a5 + spec*s5
      enddo

      do j=ms(m),me(m)
         l=m-j
         l1=l+1
         !
         rhoi=q(1,j,k,l1)
         uu=q(2,j,k,l1)*rhoi
         vv=q(3,j,k,l1)*rhoi
         ww=q(4,j,k,l1)*rhoi
         uvw=uv(j,l1)
         !
         ug=timeMetric(1,j,k,l1)
         vg=timeMetric(2,j,k,l1)
         wg=timeMetric(3,j,k,l1)
         qq=qz(j,l1)
         cc=cz(j,l1)
         !
         qqz=-dzi*wg+qq
         sp1= abs(qqz)+ cc
         sm2= eps2*sp1
         chkz = 0.5 - sign(0.5, qqz-cc)
         spec=chkz*(qqz-sp1) - sm2  - vn(3,j,l1)
         !
         s1=s(1,j,k,l1)
         s2=s(2,j,k,l1)
         s3=s(3,j,k,l1)
         s4=s(4,j,k,l1)
         s5=s(5,j,k,l1)
         !
         a5=chkz*(dzi*s4 - qq*s1)
         a2=chkz*gm1*(uvw*s1 - (uu*s2+vv*s3+ww*s4) + s5)
         !
         ! (C+sigma*I)*dQ(:,j,k,l+1)
         !
         c(1,j) = a5 + spec*s1
         c(2,j) =          uu*a5 + spec*s2
         c(3,j) =          vv*a5 + spec*s3
         c(4,j) = dzi*a2 + ww*a5 + spec*s4
         c(5,j) = qq*a2 + ge(j,l1)*a5 + spec*s5
         !
      enddo
      !
      do j=ms(m),me(m)
        l=m-j
        svt=tscale
        s(1:nvar,j,k,l)=s(1:nvar,j,k,l)-svt*0.5*(a(:,j)+b(:,j,l)+c(:,j))*d(j,k,l)
      enddo
      !
   enddo
   !
enddo
enddo
!
do l=1,lmax,lmax-1
  do k=1,kmax
   do j=1,jmax
    s(:,j,k,l)=0.
   enddo
  enddo
enddo
!
do k=1,kmax,kmax-1
 do j=1,jmax
  do l=1,lmax
   s(:,j,k,l)=0.
  enddo
 enddo
enddo
!
do j=1,jmax,jmax-1
  do k=1,kmax
   do l=1,lmax
    s(:,j,k,l)=0.
   enddo
  enddo
enddo
!
!
!write(6,*) 'tsmax/tsmin=',tsmax,tsmin,jd1,kd1,ld1
!
deallocate(ms,me,uv,vn,ge,qx,qz,cx,cz,d,a,b,c)
!
return
end subroutine lusgs_hyb

