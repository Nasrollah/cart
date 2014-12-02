!
! solve \nabla u = 1
!
program check_lusgs
  implicit none
  real*8, allocatable :: u(:,:),uf(:,:),ub(:,:)
  integer :: jmax,kmax
  integer :: j,k,m
  real*8 :: du
  integer :: msweep
  !
  jmax=32
  kmax=32
  msweep=4000
  !
  allocate(u(jmax,kmax))
  allocate(uf(jmax,kmax))
  allocate(ub(jmax,kmax))
  !
  u=0d0
  !
  do k=2,kmax-1
     do j=2,jmax-1
        u(j,k)=-1d0
        ub(j,k)=0d0
        uf(j,k)=0d0
     enddo
  enddo
  !  
  do m=1,msweep
     du=0d0
     do j=2,jmax-1
        do k=2,kmax-1
           u(j,k)=-1d0+(u(j-1,k)+u(j+1,k)+u(j,k+1)+u(j,k-1))*0.25d0
        enddo
     enddo
     do j=2,jmax-1
       do k=2,kmax -1
        du = du+(1d0 - (u(j,k-1)+u(j,k+1)+u(j+1,k)+u(j-1,k)-4d0*u(j,k))*0.25d0)**2
      enddo
     enddo
     write(6,*) m,sqrt(du/jmax/kmax)
  enddo
end program check_lusgs
