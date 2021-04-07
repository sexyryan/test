program r2r_scatter

   use, intrinsic :: iso_c_binding
   USE r2r_myfft
   USE r2r_dts
   implicit none
   ! MPI Preparation
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
   call fftw_mpi_init()

CALL r2r_setup

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

  if(myid.eq.0) then
    do l=1,nprocs
       write(*,"('count',I4)")rcounts(l)
       write(*,"('dis',I4)") displs(l)
    end do
  end if

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
   !Input Initialization
   ! input data
  if (myid==0) then
     do ix = 1, NC
       do iy = 1, NR
       if (ix>=1 .AND. ix <= NC-1 ) rin_all(iy,ix) = 2
       if (iy>=1 .AND. iy <= NR-1 ) rin_all(iy,ix) = 2
      ! if (ix>=5 .AND. ix <= 8 ) rin_all(iy,ix) = (2)*10 + ix - 4
      ! if (ix>=9 .AND. ix <= 11) rin_all(iy,ix) = (3)*10 + ix - 8
       if (ix==NC .OR. iy == NR )rin_all(iy,ix) = 1
       end do
     end do
  end if


  CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

  if(myid==0)then
    write(*,*)"========input========="
    do i=1,NR
      do j=1,NC
        write(*,'(F9.2)',advance='no') rin_all(i,j)
        enddo
      write(*,*)
    enddo
  endif

rout_all(:,:) = dst(rin_all(:,:))

! Print array on root
if(myid==0)then
  write(*,*)"========forward========="
  do i=1,NR
    do j=1,NC
      write(*,'(F9.2)',advance='no') rout_all(i,j)
      enddo
    write(*,*)
  enddo
endif
!Finish forward
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

CALL MPI_SCATTERV( rout_all, rcounts_n, displs_n, MPI_REAL8, &
 & rout, local_M_n, MPI_REAL8, &
 & 0, MPI_COMM_WORLD, IERR)

 CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

if(myid==0)then
write(*,*)"========proc 0 ========="
do i=1,NR
  do j=1,local_M
    write(*,'(F9.2)',advance='no') rout(i,j)
    enddo
  write(*,*)
enddo
endif
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

if(myid==1)then
write(*,*)"========proc 1 ========="
do i=1,NR
  do j=1,local_M
    write(*,'(F9.2)',advance='no') rout(i,j)
    enddo
  write(*,*)
enddo
endif
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

if(myid==2)then
write(*,*)"========proc 2 ========="
do i=1,NR
  do j=1,local_M
    write(*,'(F9.2)',advance='no') rout(i,j)
    enddo
  write(*,*)
enddo
endif


 CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
 call fftw_mpi_execute_r2r(planb, rout, rin)
 CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

 call MPI_gatherv (rin,local_M_n, MPI_REAL8,&
 & rin_all, rcounts_n, displs_n ,MPI_REAL8,&
 & 0, MPI_COMM_WORLD, ierr)

 CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
 if(myid==0)then
   write(*,*)"========backward========="
   do i=1,NR
     do j=1,NC
       write(*,'(F9.2)',advance='no') rin_all(i,j)/(6*12*2*2)
       enddo
     write(*,*)
   enddo
 endif
   ! Memory Deallocation
call r2r_deallocate
call MPI_Finalize(ierr)

end program r2r_scatter

!mpif90 -I/usr/local/include -L/usr/local/lib -o r2r r2r_scatter.f90 -lfftw3_mpi -lfftw3 -lm
