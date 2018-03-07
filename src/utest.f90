#include "macros.fpp"

module mod_Utest
  implicit none
contains
  subroutine expect_near_d(a,b,eps,ierr)
    double precision, intent(in) :: a, b, eps
    integer, intent(out) :: ierr
    ierr = 0
    if(abs(a-b) > eps) then
       begin_err("not equal error")
       ierr = 1
       write(0, '("a:  ", E20.10)') a
       write(0, '("b:  ", E20.10)') b
       write(0, '("eps: ", E20.5)') eps
       return
    end if
  end subroutine expect_near_d
  subroutine expect_near_c(a,b,eps,ierr)
    complex(kind(0d0)), intent(in) :: a, b
    double precision, intent(in) :: eps
    integer, intent(out) :: ierr
    ierr = 0
    if(abs(a-b) > eps) then
       begin_err("not equal error")
       ierr = 1
       write(0, *) a
       write(0, *) b
       write(0, '("eps: ", E20.5)') eps
       return
    end if
  end subroutine expect_near_c
  subroutine expect_eq_i(a,b,ierr)
    integer, intent(in) :: a, b
    integer, intent(out) :: ierr
    ierr = 0
    if(a.ne.b) then
       begin_err("not equal error")
       ierr = 1
       write(0,*) "a:", a
       write(0,*) "b:", b
       return
    end if
  end subroutine expect_eq_i
  subroutine expect_eq_s(a,b,ierr)
    character(*), intent(in) :: a, b
    integer, intent(out) :: ierr
    ierr = 0
    if(a.ne.b) then
       begin_err("not equal error")
       ierr = 1
       write(0,*) "a:", a
       write(0,*) "b:", b
       return
    end if
  end subroutine expect_eq_s
end module mod_Utest
