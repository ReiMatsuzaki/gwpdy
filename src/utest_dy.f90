#include "macros.fpp"

module Mod_UTestDy
  implicit none
contains
  subroutine UTestDy_run
    call UTestDy_mono_1dim
  end subroutine UTestDy_run
  subroutine UTestDy_mono_1dim
    use Mod_DyMono
    write(*,*) "hello"
  end subroutine UTestDy_mono_1dim
end module Mod_UTestDy

program main
  use Mod_UTestDy
  call UTestDy_run
end program main
