!>
!> Store plot3d data
!> 
subroutine storep3d(x,q,istep,fsmach,alpha,rey,totime,jmax,kmax,lmax,nvar,nf,istor)
!
implicit none
integer, intent(in) :: istep
integer, intent(in) :: nf
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in) :: fsmach,alpha,rey,totime
real*8, intent(in) :: x(3*jmax*kmax*lmax)
real*8, intent(in) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
integer :: j,k,l,n
integer :: jkmax,npts
character*16 :: integer_string,filename
integer :: xskip,xmult,qskip,qmult
!
! write output plot3d files
!
jkmax=jmax*kmax
npts=jkmax*lmax
call getstride(xskip,xmult,istor,npts,3)
call getstride(qskip,qmult,istor,npts,nvar)
!
write(integer_string,"(I7)") 1000000+istep
filename='x'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax-2*nf,kmax-2*nf,lmax-2*nf
write(1) ((((x(((l-1)*jkmax+(k-1)*jmax+j-1)*xmult+1+(n-1)*xskip),&
     j=nf+1,jmax-nf),k=nf+1,kmax-nf),l=nf+1,lmax-nf),n=1,3)
close(1)
!
filename='q'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax-2*nf,kmax-2*nf,lmax-2*nf
!
write(1) fsmach,alpha,rey,totime
write(1) ((((q(((l-1)*jkmax+(k-1)*jmax+j-1)*qmult+1+(n-1)*qskip),&
     j=nf+1,jmax-nf),k=nf+1,kmax-nf),l=nf+1,lmax-nf),n=1,5)
close(1)

return
end subroutine storep3d

!>
!> Store plot3d data
!> 
subroutine storep3di(x,q,istep,fsmach,alpha,rey,totime,jmax,kmax,lmax,nvar,nf,istor)
!
implicit none
integer, intent(in) :: istep
integer, intent(in) :: nf
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in) :: fsmach,alpha,rey,totime
real*8, intent(in) :: x(3*jmax*kmax*lmax)
real*8, intent(in) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
integer :: j,k,l,n
integer :: jkmax,npts
character*16 :: integer_string,filename
integer :: xskip,xmult,qskip,qmult
integer, allocatable :: iblank(:,:,:)
!
! write output plot3d files
!
allocate(iblank(jmax,kmax,lmax))
iblank=1
!
do l=1,lmax
   do k=1,kmax
      do j=1,nf
         iblank(j,k,l)=-1
      end do
      do j=jmax-nf+1,jmax
         iblank(j,k,l)=-1
      enddo
   enddo
enddo
!
do l=1,lmax
   do j=1,jmax
      do k=1,nf
         iblank(j,k,l)=-1
      end do
      do k=kmax-nf+1,kmax
         iblank(j,k,l)=-1
      enddo
   enddo
enddo
!
do j=1,jmax
   do k=1,kmax
      do l=1,nf
         iblank(j,k,l)=-1
      end do
      do l=lmax-nf+1,lmax
         iblank(j,k,l)=-1
      enddo
   enddo
enddo
!
jkmax=jmax*kmax
npts=jkmax*lmax
call getstride(xskip,xmult,istor,npts,3)
call getstride(qskip,qmult,istor,npts,nvar)
!
write(integer_string,"(I7)") 1000000+istep
filename='x'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax,kmax,lmax
write(1) ((((x(((l-1)*jkmax+(k-1)*jmax+j-1)*xmult+1+(n-1)*xskip),&
     j=1,jmax),k=1,kmax),l=1,lmax),n=1,3),&
     (((iblank(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
close(1)
!
filename='q'//trim(adjustl(integer_string(2:)))//'.p3d'
open(unit=1,file=filename,form='unformatted')
write(1) jmax,kmax,lmax
!
write(1) fsmach,alpha,rey,totime
write(1) ((((q(((l-1)*jkmax+(k-1)*jmax+j-1)*qmult+1+(n-1)*qskip),&
     j=1,jmax),k=1,kmax),l=1,lmax),n=1,5)
close(1)
!
deallocate(iblank)
return
end subroutine storep3di
!>
!> Store plot3d data using parallel i/o
!> this implementation uses seek, have to change to more
!> efficient window type method
!> 
subroutine storep3d_parallel(x,q,istep,fsmach,alpha,rey,totime,jmax,kmax,lmax,nvar,nf,istor)
!
use spatialCommunication, only : id,cartComm,numprocs_spatial,myid_spatial
implicit none
include 'mpif.h'
integer, intent(in) :: istep
integer, intent(in) :: nf
integer, intent(in) :: jmax,kmax,lmax,nvar
real*8, intent(in) :: fsmach,alpha,rey,totime
real*8, intent(in) :: x(3*jmax*kmax*lmax)
real*8, intent(in) :: q(nvar*jmax*kmax*lmax)
character*(*) :: istor
!
integer :: jkmax,npts
character*16 :: integer_string,xfilename,qfilename
integer :: xskip,xmult,qskip,qmult
!
integer :: ierr
integer :: itmp(3)
integer :: xsize,qsize,hsize
integer :: xsize0(1),qsize0(1)
integer(kind=MPI_OFFSET_KIND) :: disp,disp0,dispx,dispq
integer :: i,j,k,l,n,m
integer :: xfile,qfile
integer, allocatable :: jm(:),km(:),lm(:),header(:)
real*8, allocatable :: xbuf(:),qbuf(:)
integer :: js,je,ks,ke,ls,le,jdim,kdim,ldim
real*8, dimension(4) :: freestream_vec
integer :: fvecsize(1)
!
! write output plot3d files
!
allocate(jm(0:numprocs_spatial-1),km(0:numprocs_spatial-1),lm(0:numprocs_spatial-1))
!
js=nf
je=jmax-nf
ks=nf
ke=kmax-nf
ls=nf
le=lmax-nf
if (id(1)==1) js=js+1
if (id(2)==1) ks=ks+1
if (id(3)==1) ls=ls+1
!
jdim=je-js+1
kdim=ke-ks+1
ldim=le-ls+1
!
jkmax=jmax*kmax
npts=jkmax*lmax
call getstride(xskip,xmult,istor,npts,3)
call getstride(qskip,qmult,istor,npts,nvar)
allocate(xbuf(3*jdim*kdim*ldim))
allocate(qbuf(5*jdim*kdim*ldim))
!
m=0
do n=1,3
   do l=ls,le
      do k=ks,ke
         do j=js,je
            m=m+1
            xbuf(m)=x(((l-1)*jkmax+(k-1)*jmax+j-1)*xmult+1+(n-1)*xskip)
         enddo
      enddo
   enddo
enddo
!
m=0
do n=1,3
   do l=ls,le
      do k=ks,ke
         do j=js,je
            m=m+1
            qbuf(m)=q(((l-1)*jkmax+(k-1)*jmax+j-1)*qmult+1+(n-1)*qskip)
         enddo
      enddo
   enddo
enddo
!
write(integer_string,"(I7)") 1000000+istep
xfilename='x'//trim(adjustl(integer_string(2:)))//'.p3d'
qfilename='q'//trim(adjustl(integer_string(2:)))//'.p3d'
!    
! collect the dimension information at the root
! can probably use a MPI_Gather to scale better (change later)
!        
if (myid_spatial==0) then
   jm(myid_spatial)=jdim
   km(myid_spatial)=kdim
   lm(myid_spatial)=ldim
   do i=1,numprocs_spatial-1
      call mpi_recv(itmp,3,MPI_integer,i,0,cartComm,MPI_STATUS_IGNORE,ierr)
      jm(i)=itmp(1)
      km(i)=itmp(2)
      lm(i)=itmp(3)
   enddo
else
   itmp(1)=jdim
   itmp(2)=kdim
   itmp(3)=ldim
   call mpi_send(itmp,3,MPI_integer,0,0,cartComm,MPI_STATUS_IGNORE,ierr)
endif
!
call mpi_bcast(jm,numprocs_spatial,MPI_INTEGER,0,cartComm,ierr)
call mpi_bcast(km,numprocs_spatial,MPI_INTEGER,0,cartComm,ierr)
call mpi_bcast(lm,numprocs_spatial,MPI_INTEGER,0,cartComm,ierr)
!
! delete the files if they exist
!
call MPI_FILE_DELETE(xfilename,MPI_INFO_NULL,ierr)
call MPI_FILE_DELETE(qfilename,MPI_INFO_NULL,ierr)
!
! create x and q files for parallel i/o
!
call MPI_FILE_OPEN(cartComm,xfilename, &
     MPI_MODE_WRONLY + MPI_MODE_CREATE, &
     MPI_INFO_NULL,xfile,ierr)
!
call MPI_FILE_OPEN(cartComm,qfilename, &
     MPI_MODE_WRONLY + MPI_MODE_CREATE, &
     MPI_INFO_NULL,qfile,ierr)
!
! set the offsets and byte counts 
! for the unformatted i/o format
!
hsize=3*numprocs_spatial+1+4
xsize=3*(jdim*kdim*ldim)
qsize=5*(jdim*kdim*ldim)
freestream_vec=(/fsmach,alpha,rey,totime/)
fvecsize(1)=32
dispx=hsize*4 
dispq=hsize*4
do i=0,myid_spatial-1
 dispx=dispx+(jm(i)*km(i)*lm(i))*24+8
 dispq=dispq+(jm(i)*km(i)*lm(i))*40+8+40
enddo
!
! write the header in processor 0
!
if (myid_spatial==0) then
   disp0=0
   call MPI_FILE_SEEK(xfile,disp0,MPI_SEEK_SET,ierr)
   call MPI_FILE_SEEK(qfile,disp0,MPI_SEEK_SET,ierr)
   allocate(header(hsize))
   header(1)=4
   header(2)=numprocs_spatial
   header(3)=4
   header(4)=numprocs_spatial*3*4
   m=4
   do i=0,numprocs_spatial-1
      header(m+1)=jm(i)
      header(m+2)=km(i)
      header(m+3)=lm(i)
      m=m+3
   end do
   header(m+1)=numprocs_spatial*3*4
   call MPI_FILE_WRITE(xfile,header,hsize,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
   call MPI_FILE_WRITE(qfile,header,hsize,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
endif
!
! write the data in all processors
!
call MPI_FILE_SEEK(xfile,dispx,MPI_SEEK_SET,ierr)
xsize0(1)=xsize*8
call MPI_FILE_WRITE(xfile,xsize0,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_WRITE(xfile,xbuf,xsize,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_WRITE(xfile,xsize0,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(xfile,ierr)
!
call MPI_FILE_SEEK(qfile,dispq,MPI_SEEK_SET,ierr)
call  MPI_FILE_WRITE(qfile,fvecsize,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_WRITE(qfile,freestream_vec,4,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call  MPI_FILE_WRITE(qfile,fvecsize,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
qsize0(1)=qsize*8
call  MPI_FILE_WRITE(qfile,qsize0,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_WRITE(qfile,qbuf,qsize,MPI_REAL8,MPI_STATUS_IGNORE,ierr)
call  MPI_FILE_WRITE(qfile,qsize0,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
call MPI_FILE_CLOSE(xfile,ierr)
!
call mpi_barrier(cartComm,ierr)
!
deallocate(xbuf,qbuf)
!
return
end subroutine storep3d_parallel
