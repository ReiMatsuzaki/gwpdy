#include "macros.fpp"

module UTest_gwp
  use Mod_GWP
  implicit none
contains
  subroutine UTest_GWP_run
    type(Obj_GWP) :: gwp
    integer, parameter :: dim = 1
    integer, parameter :: num = 3
    integer ierr
    call GWP_new(gwp, dim, num, 'c', ierr); check_err(ierr)
    call GWP_setup(gwp, ierr);              check_err(ierr)
    call GWP_delete(gwp, ierr);             check_err(ierr)
  end subroutine UTest_GWP_run
end module UTest_gwp

program main
  use UTest_GWP
  call UTest_GWP_run
end program main
