module r2r_myfft
  !use to allocate plans and do the transformation
    use, intrinsic :: iso_c_binding
    USE r2r_dts
    implicit none
    include 'mpif.h'
    include 'fftw3-mpi.f03'
 !  Variable Definition
    integer(C_INTPTR_T) :: i, j, l, ix, iy, alloc_local, local_M, local_M_n, local_j_offset
 !   real(C_DOUBLE), allocatable :: rin(:), rin_all(:)
 !   real(C_DOUBLE), allocatable :: rout(:), rout_all(:)
    type(C_PTR) :: planf, planb
    type(C_PTR) :: cdataf, cdatab

    integer,allocatable:: rcounts(:) ! array of local_M's (for mpi_gatrherv)
    integer,allocatable:: displs(:)  ! array of local_j_offset (for mpi_gatherv)
    integer,allocatable:: rcounts_n(:) ! array of local_M's (integer size for mpi_gatrherv)
    integer,allocatable:: displs_n(:)  ! array of local_j_offset (integer size for mpi_gatherv)
    integer :: myid, nprocs, ierr

contains

subroutine r2r_setup

     ! Memory Allocation
     allocate(rcounts(nprocs), displs(nprocs))
     allocate(rin_all(NR, NC), rout_all(NR, NC))

     displs(1) = 0

     alloc_local = fftw_mpi_local_size_2d(NC, NR, &
                 & MPI_COMM_WORLD, local_M, local_j_offset)


     cdataf = fftw_alloc_real(alloc_local)
     call c_f_pointer(cdataf, rin, [NR, local_M])

     cdatab = fftw_alloc_real(alloc_local)
     call c_f_pointer(cdatab, rout, [NR, local_M])


    ! Plan Creation
    planf = fftw_mpi_plan_r2r_2d &
       & (NC, NR, rin, rout, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE)
    planb = fftw_mpi_plan_r2r_2d &
       & (NC, NR, rout, rin, MPI_COMM_WORLD, FFTW_RODFT00, FFTW_REDFT00, FFTW_ESTIMATE)
    if ((.not. c_associated(planf)) .or. (.not. c_associated(planb))) then
       write(*,*) "plan creation error!!"
       stop
    end if

    CALL MPI_ALLGATHER (local_M, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER,&
     & MPI_COMM_WORLD, IERR)

     displs(1) = 0
     do j=1,nprocs
       if((j-1).ne.0)displs(j) = displs(j-1) + rcounts(j-1)
     enddo

         local_M_n = NR * local_M
         rcounts_n = NR * rcounts
         displs_n  = NR * displs

end subroutine r2r_setup


subroutine r2r_deallocate
  USE r2r_dts
  deallocate(rin, rout, rin_all, rout_all)
  call fftw_destroy_plan(planf)
  call fftw_destroy_plan(planb)

end subroutine r2r_deallocate

function dst(psi)

  REAL(C_DOUBLE) :: psi(:,:)
  REAL(C_DOUBLE) :: dst(1:NR,1:NC)
  rin_all = psi
!scatter
  CALL MPI_SCATTERV( rin_all, rcounts_n, displs_n, MPI_REAL8, &
   & rin, local_M_n, MPI_REAL8, &
   & 0, MPI_COMM_WORLD, IERR)


   if(myid==0)then
   write(*,*)"========proc 0 ========="
   do i=1,NR
     do j=1,local_M
       write(*,'(F9.2)',advance='no') rin(i,j)
       enddo
     write(*,*)
   enddo
   endif
   CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

   if(myid==1)then
   write(*,*)"========proc 1 ========="
   do i=1,NR
     do j=1,local_M
       write(*,'(F9.2)',advance='no') rin(i,j)
       enddo
     write(*,*)
   enddo
   endif
   CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

   if(myid==2)then
   write(*,*)"========proc 2 ========="
   do i=1,NR
     do j=1,local_M
       write(*,'(F9.2)',advance='no') rin(i,j)
       enddo
     write(*,*)
   enddo
   endif

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
call fftw_mpi_execute_r2r(planf, rin, rout)
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

  call MPI_gatherv (rout,local_M_n, MPI_REAL8,&
  & rout_all, rcounts_n, displs_n,MPI_REAL8,&
  & 0, MPI_COMM_WORLD, ierr)

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

  dst = rout_all

end FUNCTION

end module
