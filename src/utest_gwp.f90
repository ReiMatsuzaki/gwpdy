#include "macros.fpp"
#include "macros_utest.fpp"

module UTestGwp
  use Mod_UTest
  implicit none
contains
  subroutine UTestGWP_run
    use Mod_Timer
    type(Obj_Timer) :: timer
    integer :: ierr
    
    call Timer_new(timer, "UTestGWP", .true., ierr); CHK_ERR(ierr)
    write(*,*) 
    write(*,*) "UTestGWP begin"
    write(*,*)

    call Timer_begin(timer, "run", ierr)
    call test_run
    call Timer_end(timer, "run", ierr)

    call Timer_begin(timer, "overlap", ierr)
    call test_overlap(ierr); CHK_ERR(ierr)
    call Timer_end(timer, "overlap", ierr)

    call Timer_begin(timer, "p2", ierr)
    call test_p2(ierr); CHK_ERR(ierr)
    call Timer_end(timer, "p2", ierr)
    
  end subroutine UTestGWP_run
  subroutine test_run
    use Mod_GWP  
    type(Obj_GWP) :: gwp
    integer, parameter :: dim = 1
    integer, parameter :: num = 3
    integer ierr
    call GWP_new(gwp, dim, num, 'c', ierr); CHK_ERR(ierr)
    call GWP_setup(gwp, ierr);              CHK_ERR(ierr)
    call GWP_delete(gwp, ierr);             CHK_ERR(ierr)
  end subroutine test_run
  subroutine test_overlap(ierr)
    use Mod_GWP
    use Mod_math, only : gtoint
    type(Obj_GWP) :: gwp
    integer, parameter :: dim = 1
    integer, parameter :: num = 3
    complex(kind(0d0)) :: S(num,num)
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: ref
    ierr = 0
    call GWP_new(gwp, dim, num, 'c', ierr); CHK_ERR(ierr)
    gwp%g(:,1,1) = 0.5d0
    gwp%R(1,1) = 0.0d0; gwp%P(1,1) = 0.0d0;
    gwp%R(2,1) = 0.0d0; gwp%P(2,1) = 0.4d0;
    gwp%R(3,1) = 0.3d0; gwp%P(3,1) = 0.4d0;
    call GWP_setup(gwp, ierr); CHK_ERR(ierr)

    call GWP_overlap(gwp, S, ierr); CHK_ERR(ierr)    

    EXPECT_EQ_C((1.0d0, 0.0d0), S(1,1), ierr); CHK_ERR(ierr)
    EXPECT_EQ_C((1.0d0, 0.0d0), S(2,2), ierr); CHK_ERR(ierr)
    EXPECT_EQ_C((1.0d0, 0.0d0), S(3,3), ierr); CHK_ERR(ierr)

    ! - see ./support/int_gwp.py
    ref = (0.9377226265225954d0, -0.05633097098542215d0)
    EXPECT_EQ_C(ref, S(1,3), ierr); CHK_ERR(ierr)
    
  end subroutine test_overlap
  subroutine test_p2(ierr)
    use Mod_GWP  
    type(Obj_GWP) :: gwp
    integer, parameter :: dim = 1
    integer, parameter :: num = 3
    complex(kind(0d0)) :: S(dim, num,num)
    integer, intent(out) :: ierr
    ierr = 0
    call GWP_new(gwp, dim, num, 'c', ierr); CHK_ERR(ierr)
    gwp%g(:,1,1) = 0.5d0
    gwp%R(1,1) = 0.0d0; gwp%P(1,1) = 0.0d0;
    gwp%R(2,1) = 0.0d0; gwp%P(2,1) = 0.4d0;
    gwp%R(3,1) = 0.3d0; gwp%P(3,1) = 0.4d0;
    call GWP_setup(gwp, ierr); CHK_ERR(ierr)

    call GWP_p2(gwp, S, ierr); CHK_ERR(ierr)
    
  end subroutine test_p2
end module UTestGwp

program main
  use UTestGWP
  call UTestGWP_run
end program main
