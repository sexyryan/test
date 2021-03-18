program example_mpi_real_to_real_2d
   use, intrinsic :: iso_c_binding
   implicit none
   include 'mpif.h'
   include 'fftw3-mpi.f03'
   !include 'aslfftw3-mpi.f03'

   ! Parameter Definition
   integer(C_INTPTR_T), parameter :: NX = 4, NY = 4

   ! Variable Definition
   integer(C_INTPTR_T) :: i, ix, iy, alloc_local, local_M, local_j_offset
!   real(C_DOUBLE), allocatable :: rin(:), rin_all(:)
!   real(C_DOUBLE), allocatable :: rout(:), rout_all(:)

   real(C_DOUBLE), pointer :: rin(:,:), rin_all(:,:)
   real(C_DOUBLE), pointer :: rout(:,:), rout_all(:,:)
   real(C_DOUBLE) :: fout
   type(C_PTR) :: planf, planb
   type(C_PTR) :: cdataf, cdatab


   real(C_DOUBLE),dimension(4,4)::ones
   real(C_DOUBLE),dimension(4,4)::a
   integer(C_INTPTR_T) ::r,c


   integer :: my_rank, num_procs, ierr
   integer, allocatable :: recv_counts(:), recv_offsets(:)


   ! MPI Preparation
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
   call fftw_mpi_init()

   do r = 1 ,4
     do c= 1,4
       ones(r,c)=r
     end do
  end do

  a = ones

   ! Memory Allocation
   allocate(recv_counts(num_procs), recv_offsets(num_procs))
   alloc_local = fftw_mpi_local_size_2d(NY, NX, MPI_COMM_WORLD, local_M, local_j_offset)
  ! allocate(rin(alloc_local), rout(alloc_local))
  !allocate(rin_all(NX * NY), rout_all(NX * NY))

  cdataf = fftw_alloc_real(alloc_local)
  call c_f_pointer(cdataf, rin, [NX, local_M])

  cdatab = fftw_alloc_real(alloc_local)
  call c_f_pointer(cdatab, rout, [NX, local_M])

  allocate(rin_all(NX ,NY), rout_all(NX , NY))

   ! Plan Creation
   planf = fftw_mpi_plan_r2r_2d &
      & (NY, NX, rin, rout, MPI_COMM_WORLD, FFTW_RODFT10, FFTW_REDFT10, FFTW_ESTIMATE)
   planb = fftw_mpi_plan_r2r_2d &
      & (NY, NX, rout, rin, MPI_COMM_WORLD, FFTW_RODFT01, FFTW_REDFT01, FFTW_ESTIMATE)
   if ((.not. c_associated(planf)) .or. (.not. c_associated(planb))) then
      write(*,*) "plan creation error!!"
      stop
   end if

  ! (just preparation for print of whole input/result data)
   call MPI_Allgather &
     & (int(local_M), 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
   do i = 1, num_procs
      recv_counts(i) = recv_counts(i) * int(NX)
   end do
   recv_offsets(1) = 0
   do i = 2, num_procs
      recv_offsets(i) = recv_offsets(i - 1) + recv_counts(i - 1)
   end do

   ! Input Initialization
   do iy = 1, local_M
        do ix = 1, NX
        !  call initial(ix, (iy + local_j_offset), NX, NY, fout)
        !  rin(ix, iy) = fout
      !rin(ix , iy+local_j_offset ) = 1.0d0
      rin(ix, iy) = a( ix,iy)
        end do
   end do

   ! (print of whole input data)
   call MPI_Allgatherv &
     & (rin, recv_counts(my_rank + 1), MPI_DOUBLE_PRECISION, &
     & rin_all, recv_counts, recv_offsets, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
   if (my_rank == 0) then
      write(*,*) "===INPUT==="
      call print_real(rin_all, NX, NY, .true., 1.0d0)
   end if

   ! FFT Execution (forward)
   call fftw_mpi_execute_r2r(planf, rin, rout)

   ! (print of whole result data)
   call MPI_Allgatherv &
      & (rout, recv_counts(my_rank + 1), MPI_DOUBLE_PRECISION,&
      & rout_all, recv_counts, recv_offsets, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
   if (my_rank == 0) then
      write(*,*) "===FORWARD==="
      call print_real(rout_all, NX, NY, .false., 1.0d0)
   end if

   ! FFT Execution (backward)
   call fftw_mpi_execute_r2r(planb, rout, rin)

   ! (print of whole result data)
   call MPI_Allgatherv &
      & (rin, recv_counts(my_rank + 1), MPI_DOUBLE_PRECISION, &
      & rin_all, recv_counts, recv_offsets, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
   if (my_rank == 0) then
      write(*,*) "===BACKWARD==="
      call print_real(rin_all, NX, NY, .true., dble(4 * NX * NY))
   end if

   ! Plan Destruction
   call fftw_destroy_plan(planf)
   call fftw_destroy_plan(planb)

   ! Memory Deallocation
   deallocate(rin, rout, rin_all, rout_all)

   call MPI_Finalize(ierr)
end program example_mpi_real_to_real_2d

!
! === Print of Real Data ===
!
subroutine print_real(r, nx, ny, io, dscale)
   use, intrinsic :: iso_c_binding
   implicit none

   integer(C_INTPTR_T), intent(in) :: nx, ny
   real(C_DOUBLE), intent(in) :: r(nx * ny)
   logical, intent(in) :: io
   real(C_DOUBLE), intent(in) :: dscale

   character(14), parameter :: c1 = '##### rin ####'
   character(14), parameter :: c2 = '#### rout ####'
   character(2), parameter :: c3 = 'ix', c4 = 'iy'
   character(6), parameter :: c6 = '-real-'
   integer(C_INTPTR_T) :: ix, iy, i

   if (io) then
      write(*,'(A29)') c1
   else
      write(*,'(A29)') c2
   end if
   write(*,'(A5,A5,A15)') c3, c4, c6
   do iy = 0, ny - 1
   do ix = 0, nx - 1
      i = ix + nx * iy + 1
      write(*,'(2I5,6X,E12.3)') ix, iy, r(i) / dscale
   end do
   end do
end subroutine print_real

subroutine initial(ix, iy, NX, NY, fout)
use, intrinsic :: iso_c_binding
integer(C_INTPTR_T), intent(in) :: ix, iy, NX, NY
real(C_DOUBLE), intent(out) :: fout

fout = ix + iy -2

return
end
