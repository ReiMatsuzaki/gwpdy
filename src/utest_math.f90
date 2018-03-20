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
    use Mod_math, only : lapack_zggev, lapack_zggev_shift
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

    H = 0
    S = 0
    H(1,1) = 2
    H(2,2) = 2
    S(1,1) = 1
    S(2,2) = 1
    call lapack_zggev(n, H, S, w, UL, UR, ierr)

    H = 0
    S = 0
    H(1,1) = (  2.5000000000000022E-004,  0.0000000000000000     )
    H(1,2) = ( -4.5952824878625620E-015, -4.6383036291635780E-017)
    H(2,1) = (  4.5952824878625620E-015,  4.6383036291635780E-017)
    !    H(2,2) = (  2.4999999999999675E-004,  8.1315162936412833E-019)
    H(2,2) = (  2.4999999999999675E-004,  8.1315162936412833E-019)
    S(1,1) = 1
    S(1,2) = (  0.0000000000000000     , -1.8379742171768383E-013)
    S(2,1) = (  0.0000000000000000     ,  1.8379742171768383E-013)
    S(2,2) = 1
    call lapack_zggev_shift(n, H, S, H(1,1), w, UL, UR, ierr); CHK_ERR(ierr)
    EXPECT_EQ_D(0.0d0, sum(abs(matmul(H, UR(:,1))-w(1)*matmul(S, UR(:,1)))), ierr)
    EXPECT_EQ_D(0.0d0, sum(abs(matmul(H, UR(:,2))-w(2)*matmul(S, UR(:,2)))), ierr)
    EXPECT_EQ_D(1.0d0, abs(dot_product(UL(:,1), matmul(S, UR(:,1)))), ierr)
    EXPECT_EQ_D(1.0d0, abs(dot_product(UL(:,2), matmul(S, UR(:,2)))), ierr)
    EXPECT_EQ_D(0.0d0, abs(dot_product(UL(:,1), matmul(S, UR(:,2)))), ierr)

    !H:
    !(  2.5000000000000022E-004,  0.0000000000000000     ) ( -4.5952824878625620E-015, -4.6383036291635780E-017)
    !    (  4.5952824878625620E-015,  4.6383036291635780E-017) (  2.4999999999999675E-004,  8.1315162936412833E-019)
    !    S:
    !(  1.0000000000000000     ,  0.0000000000000000     ) (  0.0000000000000000     , -1.8379742171768383E-013)
    ! (  0.0000000000000000     ,  1.8379742171768383E-013) (  1.0000000000000002     ,  0.0000000000000000     )
     
  end subroutine test_zggev
end module Mod_UTestMath

program main
  use Mod_UTestMath
  call UTestMath_run
end program main
