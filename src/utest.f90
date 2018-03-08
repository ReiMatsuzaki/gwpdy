#include "macros.fpp"

module mod_Utest
  implicit none
contains
  subroutine utest_near_d(a,b,eps,file,line,ierr)
    double precision, intent(in) :: a, b, eps
    character(*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(out) :: ierr
    ierr = 0
    if(abs(a-b) > eps) then
       write(0,'(A, ":", I0, ": ", A)') file, line, "a and b are not near"
       ierr = 1
       write(0, '("a:  ", E20.10)') a
       write(0, '("b:  ", E20.10)') b
       write(0, '("eps: ", E20.5)') eps
       return
    end if
  end subroutine utest_near_d
  subroutine utest_near_c(a,b,eps,file,line,ierr)
    complex(kind(0d0)), intent(in) :: a, b
    double precision, intent(in) :: eps
    character(*), intent(in) :: file
    integer, intent(in) :: line    
    integer, intent(out) :: ierr
    ierr = 0
    if(abs(a-b) > eps) then
       write(0,'(A, ":", I0, ": ", A)') file, line, "a and b are not near"
       ierr = 1
       write(0, *) a
       write(0, *) b
       write(0, '("eps: ", E20.5)') eps
       return
    end if
  end subroutine utest_near_c
  subroutine utest_eq_i(a,b,file,line,ierr)
    integer, intent(in) :: a, b
    character(*), intent(in) :: file
    integer, intent(in) :: line
    integer, intent(out) :: ierr
    ierr = 0
    if(a.ne.b) then
       write(0,'(A, ":", I0, ": ", A)') file, line, "a and b are not near"
       ierr = 1
       write(0,'("a: ", I0)') a
       write(0,'("b: ", I0)') b
       return
    end if
  end subroutine utest_eq_i
  subroutine utest_eq_s(a,b,file,line,ierr)
    character(*), intent(in) :: a, b
    character(*), intent(in) :: file
    integer, intent(in) :: line    
    integer, intent(out) :: ierr
    ierr = 0
    if(a.ne.b) then
       write(0,'(A, ":", I0, ": ", A)') file, line, "a and b are equal"
       ierr = 1
       write(0,*) "a:", a
       write(0,*) "b:", b
       return
    end if
  end subroutine utest_eq_s
end module mod_Utest
