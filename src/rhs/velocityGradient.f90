!===================================================================!
!> \brief
!! This subroutine computes velocity gradient tensor 
!! on a Cartesian grid topology 
!! 2nd, 4th and 6th order gradients are implemented and
!! can easily be extended to even high order. See matlab/octave
!! code below to get the coefficients of a given scheme
!!
!! Versions:\par
!!    - Sitaraman 08/22/2014
!!
!! Uses:\par
!!    None
!!
!! Source code:\par
!!   \include inviscidRHSupwind.f90
!!
!====================================================================!
subroutine velocityGradient(dx,dy,dz,jmax,kmax,lmax,nq,q,qwork,order,istor)
!
implicit none
!
real*8,  intent(in)    :: dx                     !< coordinate 1 spacing
real*8,  intent(in)    :: dy                     !< coordinate 2 spacing
real*8,  intent(in)    :: dz                     !< coordinate 3 spacing
integer, intent(in) :: jmax                      !< coordinate 1 dimension
integer, intent(in) :: kmax                      !< coordinate 2 dimension
integer, intent(in) :: lmax                      !< coordinate 3 dimension
integer, intent(in) :: nq                        !< number of q-variables stored
real*8,  intent(in)    :: q(nq*jmax*kmax*lmax)   !< Array of primitive variables ([\rho,u,v,w,T,...])
real*8,  intent(inout) :: qwork(9*jmax*kmax*lmax) !< work array that returns velocity gradients
integer, intent(in)    :: order                   !< order of the gradients computed
character*(*), intent(in) :: istor                !< storage type 'row' or 'column'
                                                  !< row=> q(:,j,k,l), column=>q(j,k,l,:) 
!
integer :: j,k,l,n,m,mdim
real*8 :: ds(3),dqds(3),dsinv
integer :: js,je,ks,ke,ls,le,nf,joff,joff2,joff3,iqoffset,offset
integer :: iscount,nfm2,nfp1,fwidth
integer :: idim0(2,3),stride(3),ldim(3),kdim(3)
integer :: dim0(3),dim2,dim1
integer :: idir,jkmax,ip,ipp,jstride,iq,iloc,qmult,qstride,ilocw,iqstore
!
real*8, allocatable :: scof(:,:),cof(:,:)
real*8, allocatable :: qv(:,:)
!
if (order==4) then
   allocate(scof(5,5))
   allocate(cof(5,5))
   scof(:,1)=(/-25d0,48d0,-36d0,16d0,-3d0/)
   scof(:,2)=(/-3d0,-10d0,18d0,-6d0,1d0/)
   scof(:,3)=(/1d0,-8d0,0d0,8d0,-1d0/)
   scof(:,4)=(/-1d0,6d0,-18d0,10d0,3d0/)
   scof(:,5)=(/3d0,-16d0,36d0,-48d0,25d0/)
   scof=scof*1.0d0/12.0d0
   nf=2
else if (order==6) then
  allocate(scof(7,7))
  allocate(cof(7,7))
   scof(:,1)=(/-147d0,360d0,-450d0,400d0,-225d0,72d0,-10d0/)
   scof(:,2)=(/-10d0,-77d0,150d0,-100d0,50d0,-15d0,2d0/)
   scof(:,3)=(/2d0,-24d0,-35d0,80d0,-30d0,8d0,-1d0/)
   scof(:,4)=(/-1d0,9d0,-45d0,0d0,45d0,-9d0,1d0/)
   scof(:,5)=(/1d0,-8d0,30d0,-80d0,35d0,24d0,-2d0/)
   scof(:,6)=(/-2d0,15d0,-50d0,100d0,-150d0,77d0,10d0/)
   scof(:,7)=(/10d0,-72d0,225d0,-400d0,450d0,-360d0,147d0/)
   scof=scof*1.0d0/60.0d0
   nf=3
else
   allocate(scof(3,3))
   allocate(cof(3,3))
   scof(:,1)=(/-3d0,4d0,-1d0/)
   scof(:,2)=(/-1d0,0d0,1d0/)
   scof(:,3)=(/1d0,-4d0,3d0/)
   scof=scof*0.5d0
   nf=1
endif
!
nfm2=2*nf+1
nfp1=nf+1
mdim=max(max(jmax,kmax),lmax)
allocate(qv(3,mdim))
!
! use storage string to set the variable
! location order
!
if (istor=='row') then
   qstride=1
   qmult=nq
else
   qstride=jmax*kmax*lmax
   qmult=1
endif
!
idim0(1,1)=1
idim0(2,1)=jmax
idim0(1,2)=1
idim0(2,2)=kmax
idim0(1,3)=1
idim0(2,3)=lmax
!
jkmax=jmax*kmax
!
dim0=(/jmax,kmax,lmax/)
ds=(/dx,dy,dz/)
stride=(/1,jmax,jkmax/)
ldim=(/jkmax,jkmax,jmax/)
kdim=(/jmax,1,1/)
!
! compute gradients of velocity and temperature
! and store it in "qwork"
!
do idir=1,3
   if (idir == 2) then
      ipp=mod(idir,3)+1
      ip=mod(idir+1,3)+1
   else
      ip=mod(idir,3)+1
      ipp=mod(idir+1,3)+1
   endif
   !
   js=idim0(1,idir)
   je=idim0(2,idir)
   ks=idim0(1,ip)
   ke=idim0(2,ip)
   ls=idim0(1,ipp)
   le=idim0(2,ipp)
   !
   jstride=stride(idir)   
   joff=jstride*qmult
   iqoffset=(idir-1)*3
   dim2=ldim(idir)
   dim1=kdim(idir)
   dsinv=1d0/ds(idir)
   cof=scof*dsinv
   !
   do l=ls,le
      do k=ks,ke
         !
         iq=(l-1)*dim2+(k-1)*dim1
         iqstore=iq
         !
         ! collect velocity vector along
         ! the current line in "qv"
         !
         do j=js,je
            iloc=iq*qmult+1                           
            iloc=iloc+qstride
            do n=1,3
               qv(n,j)=q(iloc)
               iloc=iloc+qstride
            enddo
            iq=iq+jstride
         enddo
         !
         iq=iqstore
         !
         iscount=1
         do j=js,js+nf-1
            dqds=0d0
            do m=1,nfm2
               dqds(:)=dqds(:)+cof(m,iscount)*qv(:,m)
            enddo
            ilocw=iq*9+iqoffset
            qwork(ilocw+1:ilocw+3)=dqds(:)
            iscount=iscount+1
            iq=iq+jstride
         enddo
         !
         ! interior
         !
         do j=js+nf,je-nf
            dqds=0d0
            do m=1,nf
               dqds(:)=dqds(:)+cof(nfp1-m,iscount)*qv(:,j-m)+&
                    cof(nfp1+m,iscount)*qv(:,j+m)
            enddo
            ilocw=iq*9+iqoffset
            qwork(ilocw+1:ilocw+3)=dqds(:)
            iq=iq+jstride
         enddo
         !
         ! do the end boundary
         !
         iscount=iscount+1
         do j=je-nf+1,je
            dqds=0d0
            do m=1,nfm2
               dqds(:)=dqds(:)+cof(m,iscount)*qv(:,je-nfm2+m)
            enddo
            ilocw=iq*9+iqoffset
            qwork(ilocw+1:ilocw+3)=dqds(:)
            iscount=iscount+1
            iq=iq+jstride
         enddo
         !
      enddo
   enddo
enddo
!
deallocate(scof,qv)
!
end subroutine velocityGradient
 ! % Octave code to 
 ! % find coefficients for nth order scheme for first derivative
 ! % run as $ octave ttable.m
 ! % Jay Sitaraman 08/22/14
 ! order=6;
 ! n=order+1;
 ! s=[];
 ! d=[];
 ! taylor_error=[];
 ! for p=0:-1:-(n-1)
 ! x=[0:n-1]+p;
 ! for i=1:n
 !  fac=1;
 !  for j=1:n
 !    m(j,i)=x(i)^(j-1)/fac;
 !    fac=fac*j;
 !  end
 ! end
 ! %
 ! b(1:n)=0.;
 ! b(2)=1.0;
 ! %
 ! m1=inv(m)*b';
 ! [N,D]=rat(m1);
 ! ld=D(1);
 ! for i=2:n
 !  ld=lcm(ld,D(i));
 ! end
 ! N=N*ld./D;
 ! du=zeros(1,order+2);
 ! for i=1:n
 !  du=du+N(i)*taylor(x(i),order+2)/ld;
 ! end
 ! taylor_error=[taylor_error;du];
 ! s=[s;N'];
 ! d=[d ld];
 ! end
 ! printf("\n---- Scheme of %d order including boundaries --- \n\n",order); 
 ! taylor_error
 ! for j=1:n
 !  printf("scof(:,%d)=(/",j);
 !  for i=1:n
 !   if (i==n)
 !   printf("%dd0",s(j,i));
 !   else
 !   printf("%dd0,",s(j,i));
 !   end
 !  end
 !  printf("/)\n");
 ! end
 ! d


