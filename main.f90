program main

  use, intrinsic :: iso_c_binding
  use myfft

  implicit none

  ! Initialize
  integer :: ierr, myid, nproc
  real(C_DOUBLE) :: texec

  integer :: a = 0
  integer :: b = 0

  call mpi_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call fftw_mpi_init()


  if(myid .eq. 0) then
    write(*,*) 'Nya-Hallo~'
    write(*,*) a
    call add2(a)
    write(*,*) a
    b = add3(a)
    write(*,*) a,'and' ,b
  end if

  write(*,*) 'In main program  a = : ', a, 'in', myid, 'of ', nproc

    write(*,*) 'Hello World from process: ', myid, 'of ', nproc


    CALL myfft_setup(nproc, myid, a)

  call mpi_finalize(ierr)

  end






 ! **********************************************************************
