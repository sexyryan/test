module myfft
  !use to allocate plans and do the transformation
  use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'
    include 'fftw3-mpi.f03'

contains

subroutine r2r_setup


end subroutine r2r_setup


subroutine r2r_execute


end subroutine r2r_execute

end module
