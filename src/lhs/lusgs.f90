!>
!!
!! LU-SGS scheme 
!! Calculate the implicit inversion of the LHS
!! This involves two bidiagonal scalar inversions
!! 
subroutine lusgs(gamma,q,s,tscale,timeMetric,spaceMetric,nq,nvar,jmax,kmax,lmax)
!
!
implicit none
!
! subroutine arguments
!
real*8, intent(in) :: gamma
integer, intent(in) :: jmax,kmax,lmax,nvar,nq
real*8, intent(in) :: q(nvar,jmax,kmax,lmax)
real*8, intent(inout) :: s(nvar,jmax,kmax,lmax)
real*8, intent(in) :: timeMetric(3,jmax,kmax,lmax)
real*8, intent(in) :: spaceMetric(10,jmax,kmax,lmax)
real*8, intent(in) :: tscale(jmax,kmax,lmax)
!
! local variables
!
integer, allocatable :: ms(:),me(:)
real*8, allocatable :: uv(:,:),vn(:,:),ge(:,:),qx(:,:),qz(:,:)
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
real*8 :: ri1,ri2,ri3,rj1,rj2,rj3,rk1,rk2,rk3,rl1,rl2,rl3
real*8 :: qq1,qq2,qq3,qqx,qqy,qqz
real*8 :: ug,vg,wg
real*8 :: qq,cc,sp1,sm2,chkx,chky,chkz,spec
real*8 :: a2,a3,a5,svt,er
real*8 :: rhoi,volumeinv
real*8 :: tsmax,tsmin
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
allocate(uv(jmax,lmax),vn(jmax,lmax),ge(jmax,lmax),qx(jmax,lmax),qz(jmax,lmax))
allocate(cx(jmax,lmax),cz(jmax,lmax),d(jmax,kmax,lmax))
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
oat=1.
!
! Setup D, diagonal term and store some arrays
!
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
         volumeinv=1./spaceMetric(ivol,j,k,l)
         rj1=spaceMetric(ixx,j,k,l)*volumeinv
         rj2=spaceMetric(ixy,j,k,l)*volumeinv
         rj3=spaceMetric(ixz,j,k,l)*volumeinv
         qq1=rj1*uu+rj2*vv+rj3*ww
         qqx=abs(rj1*(uu-ug)+rj2*(vv-vg)+rj3*(ww-wg))
         rr1=sqrt(rj1**2+rj2**2+rj3**2)
         !
         rk1=spaceMetric(iyx,j,k,l)*volumeinv
         rk2=spaceMetric(iyy,j,k,l)*volumeinv
         rk3=spaceMetric(iyz,j,k,l)*volumeinv
         qq2=rk1*uu+rk2*vv+rk3*ww
         qqy=abs(rk1*(uu-ug)+rk2*(vv-vg)+rk3*(ww-wg))
         rr2=sqrt(rk1**2+rk2**2+rk3**2)
         !
         rl1=spaceMetric(izx,j,k,l)*volumeinv
         rl2=spaceMetric(izy,j,k,l)*volumeinv
         rl3=spaceMetric(izz,j,k,l)*volumeinv
         qq3=rl1*uu+rl2*vv+rl3*ww
         qqz=abs(rl1*(uu-ug)+rl2*(vv-vg)+rl3*(ww-wg))
         rr3=sqrt(rl1**2+rl2**2+rl3**2)
         !
         !tscal=CFL*spaceMetric(ivol,j,k,l)/specrad(j,k,l) ! change here for time accurate
         !tscale=CFL*(sqrt(vol(j,k))+0.002)/(sqrt(vol(j,k))+1.)
         !
         if (tscale(j,k,l) > tsmax) then
           tsmax=tscale(j,k,l)
           jd1=j
           kd1=k
           ld1=l
         endif
         tsmin=min(tsmin,tscale(j,k,l))
         svt=tscale(j,k,l)
         d(j,k,l)=1./(1+svt*(qqx+qqy+qqz+cjkl*(rr1+rr2+rr3))*epsv)
         uv(j,l)=uvw
         qx(j,l)=qq1
         qz(j,l)=qq3
         cx(j,l)=cjkl*rr1
         cz(j,l)=cjkl*rr3
         ge(j,l)=gamma*q(5,j,k,l)*rhoi-(gamma-1)*uvw
	 if (msweep.eq.1) s(1:nq,j,k,l)=s(:,j,k,l)*svt
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
         volumeinv=1./spaceMetric(ivol,j,kr,l)
         ri1=spaceMetric(iyx,j,kr,l)*volumeinv
         ri2=spaceMetric(iyy,j,kr,l)*volumeinv
         ri3=spaceMetric(iyz,j,kr,l)*volumeinv
         !
         ug=timeMetric(1,j,kr,l)
         vg=timeMetric(2,j,kr,l)
         wg=timeMetric(3,j,kr,l)
         !
         rr2=ri1*ri1+ri2*ri2+ri3*ri3
         qq=ri1*uu+ri2*vv+ri3*ww
         cc=sqrt(ggm1*(er-uvw)*rr2)
         qqy=-ri1*ug-ri2*vg-ri3*wg+qq
         sp1=abs(qqy)+cc
         sm2=eps2*sp1
         chky=0.5+sign(0.5,qqy+cc)
         spec=chky*(qqy+sp1)+sm2
         !
         s1=s(1,j,kr,l)
         s2=s(2,j,kr,l)
         s3=s(3,j,kr,l)
         s4=s(4,j,kr,l)
         s5=s(5,j,kr,l)
         !
         a5=chky*(ri1*s2+ri2*s3+ri3*s4-qq*s1)
         a2=chky*gm1*(uvw*s1-(uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (B+sigma*I)*dQ(:,j,k-1,l)
         !
         b(1,j,l)=a5 + spec*s1
         b(2,j,l)=ri1*a2 + uu*a5 + spec*s2
         b(3,j,l)=ri2*a2 + vv*a5 + spec*s3
         b(4,j,l)=ri3*a2 + ww*a5 + spec*s4
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
         volumeinv=1./spaceMetric(ivol,j1,k,l)
         ri1=spaceMetric(ixx,j1,k,l)*volumeinv
         ri2=spaceMetric(ixy,j1,k,l)*volumeinv
         ri3=spaceMetric(ixz,j1,k,l)*volumeinv
         qq=qx(j1,l)
         cc=cx(j1,l)
         ug=timeMetric(1,j1,k,l)
         vg=timeMetric(2,j1,k,l)
         wg=timeMetric(3,j1,k,l)
         qqx=-ri1*ug-ri2*vg-ri3*wg + qq
         sp1=abs(qqx)+cc
         sm2=eps2*sp1
         chkx=0.5 + sign(0.5,qqx+cc)
         spec=chkx*(qqx+sp1)+sm2
         !
         s1=s(1,j1,k,l)
         s2=s(2,j1,k,l)
         s3=s(3,j1,k,l)
         s4=s(4,j1,k,l)
         s5=s(5,j1,k,l)
         !
         a5=chkx*(ri1*s2+ri2*s3+ri3*s4 - qq*s1)
         a2=chkx*gm1*(uvw*s1 - (uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (A+sigma*I)*dQ(:,j-1,k,l)
         !
         a(1,j) = a5 +spec*s1
         a(2,j) = ri1*a2 + uu*a5 + spec*s2
         a(3,j) = ri2*a2 + vv*a5 + spec*s3
         a(4,j) = ri3*a2 + ww*a5 + spec*s4
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
         volumeinv=1./spaceMetric(ivol,j,k,l1)
         ri1=spaceMetric(izx,j,k,l1)*volumeinv
         ri2=spaceMetric(izy,j,k,l1)*volumeinv
         ri3=spaceMetric(izz,j,k,l1)*volumeinv
         ug=timeMetric(1,j,k,l1)
         vg=timeMetric(2,j,k,l1)
         wg=timeMetric(3,j,k,l1)
         qq=qz(j,l1)
         cc=cz(j,l1)
         !
         qqz=-ri1*ug-ri2*vg-ri3*wg+qq
         sp1= abs(qqz)+ cc
         sm2= eps2*sp1
         chkz = 0.5 + sign(0.5, qqz+cc)
         spec=chkz*(qqz+sp1) + sm2
         !
         s1=s(1,j,k,l1)
         s2=s(2,j,k,l1)
         s3=s(3,j,k,l1)
         s4=s(4,j,k,l1)
         s5=s(5,j,k,l1)
         !
         a5=chkz*(ri1*s2+ri2*s3+ri3*s4 - qq*s1)
         a2=chkz*gm1*(uvw*s1 - (uu*s2+vv*s3+ww*s4) + s5)
         !
         ! (C+sigma*I)*dQ(:,j,k,l-1)
         !
         c(1,j) = a5 + spec*s1
         c(2,j) = ri1*a2 + uu*a5 + spec*s2
         c(3,j) = ri2*a2 + vv*a5 + spec*s3
         c(4,j) = ri3*a2 + ww*a5 + spec*s4
         c(5,j) = qq*a2 + ge(j,l1)*a5 + spec*s5
         !
      enddo
   !
      do j=ms(m),me(m)
       l=m-j
       svt=tscale(j,k,l)
       !
       ! factor of half because its really (A+sigma*I)/2 and (C+sigma*I)/2 and (B+sigma*I)/2
       !
       s(1:nq,j,k,l)=(s(1:nq,j,k,l) + svt*0.5*(a(:,j)+b(:,j,l)+c(:,j)))*d(j,k,l)
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
         volumeinv=1./spaceMetric(ivol,j,k,l)
         rj1=spaceMetric(ixx,j,k,l)*volumeinv
         rj2=spaceMetric(ixy,j,k,l)*volumeinv
         rj3=spaceMetric(ixz,j,k,l)*volumeinv
         qq1=rj1*uu+rj2*vv+rj3*ww
         qqx=abs(rj1*(uu-ug)+rj2*(vv-vg)+rj3*(ww-wg))
         rr1=sqrt(rj1**2+rj2**2+rj3**2)
         !
         rl1=spaceMetric(izx,j,k,l)*volumeinv
         rl2=spaceMetric(izy,j,k,l)*volumeinv
         rl3=spaceMetric(izz,j,k,l)*volumeinv
         qq3=rl1*uu+rl2*vv+rl3*ww
         qqz=abs(rl1*(uu-ug)+rl2*(vv-vg)+rl3*(ww-wg))
         rr3=sqrt(rl1**2+rl2**2+rl3**2)
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
         volumeinv=1./spaceMetric(ivol,j,kr,l)
         ri1=spaceMetric(iyx,j,kr,l)*volumeinv
         ri2=spaceMetric(iyy,j,kr,l)*volumeinv
         ri3=spaceMetric(iyz,j,kr,l)*volumeinv
         !
         ug=timeMetric(1,j,kr,l)
         vg=timeMetric(2,j,kr,l)
         wg=timeMetric(3,j,kr,l)
         !
         rr2=ri1*ri1+ri2*ri2+ri3*ri3
         qq=ri1*uu+ri2*vv+ri3*ww
         cc=sqrt(ggm1*(er-uvw)*rr2)
         qqy=-ri1*ug-ri2*vg-ri3*wg+qq
         sp1=abs(qqy)+cc
         sm2=eps2*sp1
         chky=0.5-sign(0.5,qqy-cc)
         spec=chky*(qqy-sp1)-sm2
         !
         s1=s(1,j,kr,l)
         s2=s(2,j,kr,l)
         s3=s(3,j,kr,l)
         s4=s(4,j,kr,l)
         s5=s(5,j,kr,l)
         !
         a5=chky*(ri1*s2+ri2*s3+ri3*s4-qq*s1)
         a2=chky*gm1*(uvw*s1-(uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (B+sigma*I)*dQ(:,j,k+1,l)
         !
         b(1,j,l)=a5 + spec*s1
         b(2,j,l)=ri1*a2 + uu*a5 + spec*s2
         b(3,j,l)=ri2*a2 + vv*a5 + spec*s3
         b(4,j,l)=ri3*a3 + ww*a5 + spec*s4
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
         volumeinv=1./spaceMetric(ivol,j1,k,l)
         ri1=spaceMetric(ixx,j1,k,l)*volumeinv
         ri2=spaceMetric(ixy,j1,k,l)*volumeinv
         ri3=spaceMetric(ixz,j1,k,l)*volumeinv
         qq=qx(j1,l)
         cc=cx(j1,l)
         ug=timeMetric(1,j1,k,l)
         vg=timeMetric(2,j1,k,l)
         wg=timeMetric(3,j1,k,l)
         qqx=-ri1*ug-ri2*vg-ri3*wg + qq
         sp1=abs(qqx)+cc
         sm2=eps2*sp1
         chkx=0.5 - sign(0.5,qqx-cc)
         spec=chkx*(qqx-sp1)-sm2
         !
         s1=s(1,j1,k,l)
         s2=s(2,j1,k,l)
         s3=s(3,j1,k,l)
         s4=s(4,j1,k,l)
         s5=s(5,j1,k,l)
         !
         a5=chkx*(ri1*s2+ri2*s3+ri3*s4 - qq*s1)
         a2=chkx*gm1*(uvw*s1 - (uu*s2 + vv*s3 + ww*s4) + s5)
         !
         ! (A+sigma*I)*dQ(:,j+1,k,l)
         !
         a(1,j) = a5 +spec*s1
         a(2,j) = ri1*a2 + uu*a5 + spec*s2
         a(3,j) = ri2*a2 + vv*a5 + spec*s3
         a(4,j) = ri3*a2 + ww*a5 + spec*s4
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
         volumeinv=1./spaceMetric(ivol,j,k,l1)
         ri1=spaceMetric(izx,j,k,l1)*volumeinv
         ri2=spaceMetric(izy,j,k,l1)*volumeinv
         ri3=spaceMetric(izz,j,k,l1)*volumeinv
         ug=timeMetric(1,j,k,l1)
         vg=timeMetric(2,j,k,l1)
         wg=timeMetric(3,j,k,l1)
         qq=qz(j,l1)
         cc=cz(j,l1)
         !
         qqz=-ri1*ug-ri2*vg-ri3*wg+qq
         sp1= abs(qqz)+ cc
         sm2= eps2*sp1
         chkz = 0.5 - sign(0.5, qqz-cc)
         spec=chkz*(qqz-sp1) - sm2
         !
         s1=s(1,j,k,l1)
         s2=s(2,j,k,l1)
         s3=s(3,j,k,l1)
         s4=s(4,j,k,l1)
         s5=s(5,j,k,l1)
         !
         a5=chkz*(ri1*s2+ri2*s3+ri3*s4 - qq*s1)
         a2=chkz*gm1*(uvw*s1 - (uu*s2+vv*s3+ww*s4) + s5)
         !
         ! (C+sigma*I)*dQ(:,j,k,l+1)
         !
         c(1,j) = a5 + spec*s1
         c(2,j) = ri1*a2 + uu*a5 + spec*s2
         c(3,j) = ri2*a2 + vv*a5 + spec*s3
         c(4,j) = ri3*a2 + ww*a5 + spec*s4
         c(5,j) = qq*a2 + ge(j,l1)*a5 + spec*s5
         !
      enddo
      !
      do j=ms(m),me(m)
        l=m-j
        svt=tscale(j,k,l)
        s(1:nq,j,k,l)=s(1:nq,j,k,l)-svt*0.5*(a(:,j)+b(:,j,l)+c(:,j))*d(j,k,l)
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
end subroutine lusgs

