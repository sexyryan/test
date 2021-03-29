program example_mpi_real_to_real_2d

   use, intrinsic :: iso_c_binding
   implicit none

   include 'mpif.h'
   include 'fftw3-mpi.f03'
   !include 'aslfftw3-mpi.f03'

   ! Parameter Definition
   integer(C_INTPTR_T), parameter :: NR = 5 !number of rows

   integer(C_INTPTR_T), parameter :: NC = 11 !number of columns

   ! Variable Definition
   integer(C_INTPTR_T) :: i, j, l, ix, iy, alloc_local, local_M, local_j_offset
!   real(C_DOUBLE), allocatable :: rin(:), rin_all(:)
!   real(C_DOUBLE), allocatable :: rout(:), rout_all(:)

   real(C_DOUBLE), POINTER :: rin(:,:), rin_all(:,:)
   real(C_DOUBLE), POINTER :: rout(:,:), rout_all(:,:)

   type(C_PTR) :: planf, planb
   type(C_PTR) :: cdataf, cdatab


   integer,allocatable:: rcounts(:) ! array of local_M's (for mpi_gatrherv)
   integer,allocatable:: displs(:)  ! array of local_j_offset (for mpi_gatherv)

   integer :: myid, nprocs, ierr


   ! MPI Preparation
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
   call fftw_mpi_init()


   ! Memory Allocation
   allocate(rcounts(nprocs), displs(nprocs))
   allocate(rin_all(NR ,NC), rout_all(NR , NC))

  displs(1) = 0

  alloc_local = fftw_mpi_local_size_2d(NC, NR, &
              & MPI_COMM_WORLD, local_M, local_j_offset)


  cdataf = fftw_alloc_real(alloc_local)
  call c_f_pointer(cdataf, rin, [NR, local_M])

  cdatab = fftw_alloc_real(alloc_local)
  call c_f_pointer(cdatab, rout, [NR, local_M])


   ! Plan Creation
   planf = fftw_mpi_plan_r2r_2d &
      & (NC, NR, rin, rout, MPI_COMM_WORLD, FFTW_RODFT10, FFTW_REDFT10, FFTW_ESTIMATE)
   planb = fftw_mpi_plan_r2r_2d &
      & (NC, NR, rout, rin, MPI_COMM_WORLD, FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE)
   if ((.not. c_associated(planf)) .or. (.not. c_associated(planb))) then
      write(*,*) "plan creation error!!"
      stop
   end if
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)


CALL MPI_ALLGATHER (local_M, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER,&
 & MPI_COMM_WORLD, IERR)

 displs(1) = 0
 do j=1,nprocs
   if((j-1).ne.0)displs(j) = displs(j-1) + rcounts(j-1)
 enddo



CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

if(myid.eq.0) then
  do l=1,nprocs
     write(*,"('count',I4)")rcounts(l)
     write(*,"('dis',I4)") displs(l)
  end do
end if

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

 !Input Initialization
do ix = 1, local_M
  do iy = 1, NR
  rin(iy,ix) = (myid+1)*10 + ix
  end do
end do

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)
if(myid==0)then
  write(*,*)"========proc 0 ========="
  do i=1,NR
    do j=1,local_M
      write(*,'(F6.2)',advance='no') rin(i,j)
      enddo
    write(*,*)
  enddo
endif
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

if(myid==1)then
  write(*,*)"========proc 1 ========="
  do i=1,NR
    do j=1,local_M
      write(*,'(F6.2)',advance='no') rin(i,j)
      enddo
    write(*,*)
  enddo
endif
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

if(myid==2)then
  write(*,*)"========proc 2 ========="
  do i=1,NR
    do j=1,local_M
      write(*,'(F6.2)',advance='no') rin(i,j)
      enddo
    write(*,*)
  enddo
endif

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

local_M = NR * local_M
rcounts = NR * rcounts
displs  = NR * displs

call MPI_gatherv (rin,local_M, MPI_REAL8,&
& rin_all, rcounts, displs,MPI_REAL8,&
& 0, MPI_COMM_WORLD, ierr)


if(myid==0)then
  write(*,*)"========input========="
  do i=1,NR
    do j=1,NC
      write(*,'(F9.2)',advance='no') rin_all(i,j)
      enddo
    write(*,*)
  enddo
endif

! FFT Execution (forward)
CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)

call fftw_mpi_execute_r2r(planf, rin, rout)

CALL MPI_BARRIER (MPI_COMM_WORLD, IERR)


  call MPI_gatherv (rout,local_M, MPI_REAL8,&
  & rout_all, rcounts, displs,MPI_REAL8,&
  & 0, MPI_COMM_WORLD, ierr)

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

 call fftw_mpi_execute_r2r(planb, rout, rin)

 call MPI_gatherv (rin,local_M, MPI_REAL8,&
 & rin_all, rcounts, displs,MPI_REAL8,&
 & 0, MPI_COMM_WORLD, ierr)

 if(myid==0)then
   write(*,*)"========backward========="
   do i=1,NR
     do j=1,NC
       write(*,'(F9.2)',advance='no') rin_all(i,j)
       enddo
     write(*,*)
   enddo
 endif


   ! Memory Deallocation
 deallocate(rin, rout, rin_all, rout_all)

 call MPI_Finalize(ierr)
end program example_mpi_real_to_real_2d


!mpif90 -I/usr/local/include -L/usr/local/lib -o r2r example_mpi_real_to_real_2d.f90 -lfftw3_mpi -lfftw3 -lm
