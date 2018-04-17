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

    call Timer_begin(timer, "dfact", ierr)
    call test_dfact(ierr)
    call Timer_end(timer, "dfact", ierr)

    call Timer_begin(timer, "gamma", ierr)
    call test_gamma(ierr)
    call Timer_end(timer, "gamma", ierr)

    call Timer_begin(timer, "gtoint_shift", ierr)
    call test_gtoint_shift(ierr)
    call Timer_end(timer, "gtoint_shift", ierr)
    
    
    write(*,*)
    write(*,*) "UTestMath end"
    write(*,*) 

    call Timer_delete(timer, ierr)
    
  end subroutine UTestMath_run
  subroutine test_zggev(ierr)
    use Mod_math, only : lapack_zggev, lapack_zggev_shift, sort_zggev
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
    H(1,1) = (0.025d0, 0.0d0)
    H(1,2) = (0.0d0,   0.005d0)
    H(2,1) = (0.0d0,  -0.005d0)
    H(2,2) = (0.035d0, 0.0d0)
    S(1,1) = 1
    S(1,2) = (  0.2d0, 0.3d0)
    S(2,1) = (  0.2d0,-0.3d0)
    S(2,2) = 1
    call lapack_zggev(n, H, S, w, UL, UR, ierr); CHK_ERR(ierr)
    call sort_zggev(n, w, UL, UR, ierr)
    EXPECT_EQ_C(w(1),  dot_product(UL(:,1), matmul(H, UR(:,1))), ierr)
    EXPECT_EQ_C(w(2),  dot_product(UL(:,2), matmul(H, UR(:,2))), ierr)
    EXPECT_EQ_D(0.0d0, abs(dot_product(UL(:,1), matmul(H, UR(:,2)))), ierr)
    EXPECT_EQ_D(0.0d0, abs(dot_product(UL(:,2), matmul(H, UR(:,1)))), ierr)
    
    EXPECT_EQ_D(1.0d0, abs(dot_product(UL(:,1), matmul(S, UR(:,1)))), ierr)
    EXPECT_EQ_D(1.0d0, abs(dot_product(UL(:,2), matmul(S, UR(:,2)))), ierr)
    EXPECT_EQ_D(0.0d0, abs(dot_product(UL(:,1), matmul(S, UR(:,2)))), ierr)
    EXPECT_EQ_D(0.0d0, abs(dot_product(UL(:,2), matmul(S, UR(:,1)))), ierr)
    
    !call lapack_zggev_shift(n, H, S, H(1,1), w, UL, UR, ierr); CHK_ERR(ierr)
    !EXPECT_EQ_D(0.0d0, sum(abs(matmul(H, UR(:,1))-w(1)*matmul(S, UR(:,1)))), ierr)
    !EXPECT_EQ_D(0.0d0, sum(abs(matmul(H, UR(:,2))-w(2)*matmul(S, UR(:,2)))), ierr)
    !EXPECT_EQ_D(1.0d0, abs(dot_product(UL(:,1), matmul(S, UR(:,1)))), ierr)
    !EXPECT_EQ_D(1.0d0, abs(dot_product(UL(:,2), matmul(S, UR(:,2)))), ierr)
    !EXPECT_EQ_D(0.0d0, abs(dot_product(UL(:,1), matmul(S, UR(:,2)))), ierr)

    !H:
    !(  2.5000000000000022E-004,  0.0000000000000000     ) ( -4.5952824878625620E-015, -4.6383036291635780E-017)
    !    (  4.5952824878625620E-015,  4.6383036291635780E-017) (  2.4999999999999675E-004,  8.1315162936412833E-019)
    !    S:
    !(  1.0000000000000000     ,  0.0000000000000000     ) (  0.0000000000000000     , -1.8379742171768383E-013)
    ! (  0.0000000000000000     ,  1.8379742171768383E-013) (  1.0000000000000002     ,  0.0000000000000000     )
     
  end subroutine test_zggev
  subroutine test_dfact(ierr)
    use Mod_math, only : dfact
    integer, intent(out) :: ierr
    integer, parameter :: maxn = 7
    integer :: npp(0:maxn)
    call dfact(maxn, npp(:), ierr); CHK_ERR(ierr)
    EXPECT_EQ_I(1, npp(0), ierr)
    EXPECT_EQ_I(1, npp(1), ierr)
    EXPECT_EQ_I(1*2, npp(2), ierr)
    EXPECT_EQ_I(1*3, npp(3), ierr)
    EXPECT_EQ_I(1*2*4, npp(4), ierr)
    EXPECT_EQ_I(1*3*5, npp(5), ierr)
    EXPECT_EQ_I(1*2*4*6, npp(6), ierr)
    EXPECT_EQ_I(1*3*5*7, npp(7), ierr)
  end subroutine test_dfact
  subroutine test_gamma(ierr)
    use Mod_const, only : PI
    use Mod_math, only : gamma_half_int
    integer, intent(out) :: ierr
    integer, parameter :: nmax=5
    double precision  :: g(0:nmax)
    call gamma_half_int(nmax, g, ierr); CHK_ERR(ierr)
    EXPECT_EQ_D(sqrt(PI), g(0), ierr)
    EXPECT_EQ_D(sqrt(PI)/2, g(1), ierr)
    EXPECT_EQ_D(sqrt(PI)*3.0d0/4, g(2), ierr)
    EXPECT_EQ_D(sqrt(PI)*15.0d0/8, g(3), ierr)
  end subroutine test_gamma
  subroutine test_gtoint_shift(ierr)
    ! See 2018/4/16/qpbranch for detail
    use Mod_math, only : gtoint_shift
    integer, intent(out) :: ierr
    complex(kind(0d0)) a
    double precision q0, b
    integer n
    double precision, parameter :: tol = 1.0d-13
    integer, parameter :: maxn = 6
    double precision refs(0:maxn)
    complex(kind(0d0)) :: calc(0:maxn)
    a = 1.1d0
    b = 1.2d0
    q0 = 0.3d0
    
    refs(0:maxn) = (/1.10988699911633d0, &
         -0.173721443339947d0, &
         0.268470964852411d0, &
         -0.117552604646441d0, &
         0.193489297804928d0, &
         -0.132504937609851d0, &
         0.231054357413855d0 /)

    call gtoint_shift(maxn, a, b, q0, calc, ierr)

    do n = 0, maxn
       EXPECT_NEAR_D(refs(n), real(calc(n)), tol, ierr)
       if(ierr.ne.0) then
          write(0,*) "n:", n
       end if
    end do
   
  end subroutine test_gtoint_shift
end module Mod_UTestMath

program main
  use Mod_UTestMath
  call UTestMath_run
end program main
