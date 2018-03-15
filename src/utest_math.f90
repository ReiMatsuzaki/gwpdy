#include "macros.fpp"
#include "macros_utest.fpp"
module Mod_UTestMath
  use Mod_UTest
  implicit none
contains
  subroutine UTestMath_run
    use Mod_Timer
    type(Obj_Timer) :: timer
    integer :: ierr

    call Timer_new(timer, "UTestDy", .true., ierr)
    write(*,*)
    write(*,*) "UTestMath begin"
    write(*,*)

    call Timer_begin(timer, "lapack_zggev", ierr)
    call test_zggev(ierr)
    call Timer_end(timer, "lapack_zggev", ierr)
    
    write(*,*)
    write(*,*) "UTestMath end"
    write(*,*) 

    call Timer_delete(timer, ierr)
    
  end subroutine UTestMath_run
  subroutine test_zggev(ierr)
    use Mod_math, only : lapack_zggev
    integer, intent(out) :: ierr
    integer, parameter :: n = 2
    complex(kind(0d0)) :: H(n,n), S(n,n)
    complex(kind(0d0)) :: w(n), UL(n,n), UR(n,n)

    ierr = 0

    H(1,1) = 1
    H(1,2) = 2
    H(2,1) = 2
    H(2,2) = 3

    S(1,1) = 1
    S(1,2) = 0
    S(2,1) = 0
    S(2,2) = 1

    call lapack_zggev(n, H, S, w, UL, UR, ierr)
    
  end subroutine test_zggev
end module Mod_UTestMath

program main
  use Mod_UTestMath
  call UTestMath_run
end program main
