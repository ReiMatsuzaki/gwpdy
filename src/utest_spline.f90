#include "macros_utest.fpp"
module Mod_UTestSpline
  use Mod_Utest
  implicit none
contains
  subroutine UTestSpline_run
    call test_spline
    call test_spline2
  end subroutine UTestSpline_run
  subroutine test_spline
    use Mod_spline
    type(Obj_Spline) spl
    integer, parameter :: n = 20
    double precision :: x(n), y(n), x0
    double precision, parameter :: a = 1.2d0
    integer i, ierr
    double precision ref, res

    x = (/(  (i-1)*0.1d0,  i = 1, n  )/)
    y = (/(  sin(a*x(i)),    i = 1, n  )/)

    call Spline_new(spl, n, x, y, "const0", ierr)

    x0 = 0.32d0
    
    ref = sin(a*x0)
    call Spline_at(spl, x0, res, ierr)
    EXPECT_NEAR_D(ref, res, 1.0d-5, ierr)
    
    ref = sin(a*x0)
    call Spline_deriv_at(spl, 0, x0, res, ierr)
    EXPECT_NEAR_D(ref, res, 1.0d-5, ierr)

    ref = a*cos(a*x0)
    call Spline_deriv_at(spl, 1, x0, res, ierr)
    EXPECT_NEAR_D(ref, res, 1.0d-5, ierr)

    ref = -a*a*sin(a*x0)
    call Spline_deriv_at(spl, 2, x0, res, ierr)
    EXPECT_NEAR_D(ref, res, 1.0d-4, ierr)
    
    call Spline_delete(spl)
    
  end subroutine test_spline
  subroutine test_spline2
    use Mod_spline
    type(Obj_Spline) spl
    double precision res
    integer ierr
    
    call Spline_new_file(spl, "data/x12.csv", "const0", ierr)

    call Spline_at(spl, -3.8431372549019613d0, res, ierr)    
    EXPECT_EQ_D(0.0638358366961411d0, res, ierr)

    call Spline_at(spl, -3.6862745098039227d0, res, ierr)
    EXPECT_EQ_D(0.0694612051533508d0, res, ierr)

    call Spline_at(spl, -3.529411764705884d0, res, ierr)
    EXPECT_EQ_D(0.0765542993986583d0, res, ierr)
    
    call Spline_at(spl, -3.3725490196078418d0, res, ierr)
    EXPECT_EQ_D(0.0857038924140142d0, res, ierr)

    call Spline_at(spl, -3.2156862745098032d0, res, ierr)
    EXPECT_EQ_D(0.0977228373804536d0, res, ierr)
    
    call Spline_delete(spl)
    
  end subroutine test_spline2
end module Mod_UTestSpline

program main
  use Mod_UTestSpline
  call UTestSpline_run
end program main
