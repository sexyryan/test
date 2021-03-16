module myfft

      use, intrinsic :: iso_c_binding
      implicit none

      include 'mpif.h'
      include 'fftw3-mpi.f03'

contains

subroutine myfft_setup(nproc, myid, a)
       ! Parameter Definition

       integer(C_INTPTR_T), parameter :: L = 1024
       integer(C_INTPTR_T), parameter :: M = 1024

       !Variable Definition
       integer(C_INTPTR_T) :: alloc_local, local_M, local_j_offset
       integer(C_INTPTR_T) :: i, j

       complex(C_DOUBLE_COMPLEX), pointer :: fdata(:,:)
       complex(C_DOUBLE_COMPLEX) :: fout


       type(C_PTR) :: plan, cdata
       integer :: ierr, a
       integer, intent(in) :: myid, nproc

       real(C_DOUBLE) :: t1, t2, t3, t4, tplan, tmid, texec

!
!   get local data size and allocate (note dimension reversal)
   alloc_local = fftw_mpi_local_size_2d(M, L, &
    &                  MPI_COMM_WORLD, local_M, local_j_offset)
   cdata = fftw_alloc_complex(alloc_local)
   call c_f_pointer(cdata, fdata, [L,local_M])
!

!   create MPI plan for in-place forward DFT (note dimension reversal)
         t1 = MPI_wtime()
   plan = fftw_mpi_plan_dft_2d(M, L, fdata, fdata, &
    &         MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE)
         t2 = MPI_wtime()
!
! initialize data to some function my_function(i,j)
   do j = 1, local_M
     do i = 1, L
       call initial(i, (j + local_j_offset), L, M, fout)
       fdata(i, j) = fout
     end do
   end do
!
! compute transform (as many times as desired)
         t3 = MPI_wtime()
         call fftw_mpi_execute_dft(plan, fdata, fdata)
         t4 = MPI_wtime()
!
! print solutinos
         tplan = t2 - t1
         tmid = t3 - t2
         texec = t4 - t3
         a = a + 1
         write(*,*) 'in module a = : ', a, 'in', myid, 'of ', nproc

         if (myid.eq.0) print*,'Tplan=',tplan,'   Tmid=',tmid,'   Texec=',texec
!
! deallocate and destroy plans
      call fftw_destroy_plan(plan)
      call fftw_mpi_cleanup()
      call fftw_free(cdata)

 end subroutine
!
subroutine initial(i, j, L, M, fout)
  use, intrinsic :: iso_c_binding
  integer(C_INTPTR_T) :: i, j, L, M
  complex(C_DOUBLE_COMPLEX) :: fout
  real(C_DOUBLE), parameter :: amp = 0.25
  real(C_DOUBLE) :: xx, yy, LL, MM, r1

    xx = real(i, C_DOUBLE) - real((L+1)/2, C_DOUBLE)
    yy = real(j, C_DOUBLE) - real((M+1)/2, C_DOUBLE)
    LL = real(L, C_DOUBLE)
    MM = real(M, C_DOUBLE)

    r1 = sqrt(((xx/LL)**2.) + ((yy/MM)**2.))
    if (r1 .le. amp) then
      fout = CMPLX(1., 0. , C_DOUBLE_COMPLEX)
          else
            fout = CMPLX(0., 1. , C_DOUBLE_COMPLEX)
          endif
          return
end subroutine

subroutine add2(a)
  implicit none
  integer :: a

  a = a + 2
  return
end subroutine

function add3(a)
  implicit none
  integer :: a
  integer :: add3
  a = a + 3
  return
end function

end module

!to compile
!mpif90 -I/usr/local/include -L/usr/local/lib -o ex2_f ex2_f.f90 -lfftw3_mpi -lfftw3 -lm
!to compile seperately
!mpif90 -I/usr/local/include -L/usr/local/lib -c myfft.f90 -lfftw3_mpi -lfftw3 -lm
!mpif90 -I/usr/local/include -L/usr/local/lib -c main.f90 -lfftw3_mpi -lfftw3 -lm
!mpif90 -I/usr/local/include -L/usr/local/lib myfft.o main.o -lfftw3_mpi -lfftw3 -lm
!to execute
!mpirun -n #nprocs ./main
