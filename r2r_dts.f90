module r2r_dts
  use, intrinsic :: iso_c_binding
  implicit none
  !used to declare variables
  ! Parameter Definition
  integer(C_INTPTR_T), parameter :: NR = 5 !number of rows

  integer(C_INTPTR_T), parameter :: NC = 11 !number of columns

  integer(C_INTPTR_T),parameter :: NRR = 5 !number of rows

  integer(C_INTPTR_T),parameter :: NCC = 11 !number of columns

  real(C_DOUBLE), POINTER :: rin(:,:), rin_all(:,:)
  real(C_DOUBLE), POINTER :: rout(:,:), rout_all(:,:)


  REAL(C_DOUBLE) ,dimension(2 :NRR, 2:NCC):: test


end module
