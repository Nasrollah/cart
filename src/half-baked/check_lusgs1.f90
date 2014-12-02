!
! solve \nabla u = 1
!
program check_lusgs
  implicit none
  real*8, allocatable :: u(:,:),uf(:,:),ub(:,:),deltau(:,:)
  real*8, allocatable :: f(:,:)
  integer :: jmax,kmax
  integer :: j,k,m
  real*8 :: du
  integer :: msweep
  !
  jmax=32
  kmax=32
  msweep=2000
  !
  allocate(u(jmax,kmax))
  allocate(f(jmax,kmax))
  allocate(uf(jmax,kmax))
  allocate(ub(jmax,kmax))
  allocate(deltau(jmax,kmax))
  !
  u=0d0
  !
  do k=1,kmax
     do j=1,jmax
        u(j,k)=0d0
        f(j,k)=1d0
        ub(j,k)=0d0
        uf(j,k)=0d0
        deltau(j,k)=-1d0
     enddo
  enddo
  !  
  do m=1,msweep
     du=0d0
     do k=2,kmax-1
      do j=2,jmax-1
           !u(j,k)=-f(j,k)+(u(j-1,k)+u(j,k-1))*0.25d0
           u(j,k)=-f(j,k)+(u(j-1,k)+u(j,k-1)+u(j+1,k)+u(j,k+1))*0.25d0
           deltau(j,k)=-f(j,k)+(deltau(j-1,k)+deltau(j,k-1)+deltau(j+1,k)+deltau(j,k+1))*0.25d0
        enddo
     enddo
     !do k=kmax-1,2,-1
     ! do j=jmax-1,2,-1
     !      u(j,k)=u(j,k)+(u(j+1,k)+u(j,k+1))*0.25d0
     !      !u(j,k)=-f(j,k)+(u(j-1,k)+u(j,k-1)+u(j+1,k)+u(j,k+1))*0.25d0
     !      !deltau(j,k)=deltau(j,k)+(deltau(j+1,k)+deltau(j,k+1))*0.25d0
     !   enddo
     !enddo
     do j=2,jmax-1
       do k=2,kmax -1
        u(j,k)=u(j,k)+deltau(j,k)
        f(j,k)=(1d0 - (u(j,k-1)+u(j,k+1)+u(j+1,k)+u(j-1,k)-4d0*u(j,k))*0.25d0)
        du = du+f(j,k)**2
        deltau(j,k)=0d0
      enddo
     enddo
     write(6,*) m,sqrt(du/jmax/kmax)
  enddo
end program check_lusgs
