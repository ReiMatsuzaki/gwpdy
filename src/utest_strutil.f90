#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_UtestStrUtil
  use Mod_Utest
  use Mod_StrUtil
  implicit none
contains
  subroutine TestStrUtil_run

    call test_str_split
    call test_basic
    call test_str_split

  end subroutine TestStrUtil_run
  subroutine test_basic
    integer :: i
    double precision :: d
    integer :: ierr
    call str2i("2", i, ierr); CHK_ERR(ierr)
    EXPECT_EQ_I(2, i, ierr)
    call str2d("2.1", d, ierr); CHK_ERR(ierr)
    EXPECT_EQ_D(2.1d0, d, ierr)
  end subroutine test_basic
  subroutine test_str_split
    character(20) :: line
    character(20) :: lines(10)
    integer n
    integer ierr
    line = "x1,x2,3.14"
    call str_split(line, ",", n, lines(:), ierr) ; CHK_ERR(ierr)
    EXPECT_EQ_I(3, n, ierr)
    EXPECT_EQ_S("x1", lines(1), ierr)
    EXPECT_EQ_S("x2", lines(2), ierr)
    EXPECT_EQ_S("3.14", lines(3), ierr)

    line = "momo"
    call str_split(line, "n", n,lines(:), ierr); CHK_ERR(ierr)
    EXPECT_EQ_I(n, 1, ierr)
    EXPECT_EQ_S("momo", lines(1), ierr)
    
  end subroutine test_str_split
  subroutine test_str2vec
    character(100) :: line
    double precision :: xs(5)
    integer n, ierr
    line = "linspace:1.0:3.0:3"
    call str2vec(line, n, xs, ierr); CHK_ERR(ierr)
    EXPECT_EQ_I(n,3, ierr)
    EXPECT_EQ_D(1.0d0, xs(1), ierr)
    EXPECT_EQ_D(2.0d0, xs(2), ierr)
    EXPECT_EQ_D(3.0d0, xs(3), ierr)

    line = "scalar:1.12"
    call str2vec(line, n, xs, ierr); CHK_ERR(ierr)
    EXPECT_EQ_I(n,1,ierr)
    EXPECT_EQ_D(1.12d0, xs(1), ierr)
    
  end subroutine test_str2vec
end module Mod_UtestStrUtil

program main
  use Mod_UtestStrUtil

  call TestStrUtil_run
  
end program main
