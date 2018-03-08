#include "macros.fpp"
module Mod_Harm
  implicit none
  double precision k_, m_, r0_
contains
  subroutine Harm_new
    k_ = 1
    m_ = 1
    r0_ = 0
  end subroutine Harm_new
  subroutine Harm_H_X(Q, HeIJ, XkIJ, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
    integer, intent(out) :: ierr
    ierr = 0
    HeIJ(1,1) = k_/2*(Q(1)-r0_)**2
    XkIJ(1,1,1) = 0.0d0
  end subroutine Harm_H_X
  subroutine Harm_ctraj(R0, P0, t, R, P, ierr)
    double precision, intent(in) :: R0(:), P0(:)
    double precision, intent(in) :: t
    double precision, intent(out):: R(:), P(:)
    integer, intent(out) :: ierr
    double precision :: w, A, B
    ierr = 0
    w    = sqrt(k_/m_)
    A    = R0(1)
    B    = P0(1)*m_/w
    R(1) = A*cos(w*t) + B*sin(w*t)
    P(1) = w/m_*(-A*sin(w*t) + B*cos(w*t))
    
  end subroutine Harm_ctraj
end module Mod_harm

module Mod_Uncoupled2e
  implicit none
contains
  subroutine Uncoupled2e_H_X(Q, HeIJ, XkIJ, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
    integer, intent(out) :: ierr
    ierr = 0
    HeIJ(:,:) = 0
    HeIJ(1,1) = Q(1)**2
    HeIJ(2,2) = Q(1)**2+1
    XkIJ(:,:,:) = 0
  end subroutine Uncoupled2e_H_X
end module Mod_Uncoupled2e

module Mod_UTestDy
  use Mod_Utest
  implicit none
contains
  subroutine UTestDy_run
    use Mod_Timer    
    type(Obj_Timer) :: timer
    integer ierr
    call Timer_new(timer, "UTestDy", .true., ierr)
    write(*,*)
    write(*,*) "UTestDy begin"
    write(*,*)
    call Timer_begin(timer, "1f1e", ierr)
    call test_1f1e
    call Timer_end(timer, "1f1e", ierr)
    call Timer_begin(timer, "1f2e", ierr)
    call test_1f2e_uncouple
    call Timer_end(timer, "1f2e", ierr)
    write(*,*)
    write(*,*) "UTestDy end"
    write(*,*)    
  end subroutine UTestDy_run
  subroutine test_1f1e
    use Mod_DyMono, only : R_, P_, c_, dt_, nt_, n1t_, &
         DyMono_new, DyMono_setup, DyMono_update, DyMono_delete
    use Mod_Harm, only : Harm_new, Harm_H_X, Harm_ctraj, k_, m_, r0_
    integer, parameter :: nf = 1, ne=1
    double precision, parameter :: eps = 1.0d-10
    integer ierr, it, iit
    double precision :: R(nf), P(nf)

    ! -- Initialize --
    call DyMono_new(nf, ne, ierr); check_err(ierr)
    R_(1)   = 1
    P_(1)   = 1
    nt_ = 100
    dt_ = 0.01d0
    n1t_ = 2
    c_(1) = 1
    call DyMono_setup(ierr)
    call Harm_new
    k_ = 1
    m_ = 1
    r0_ = 0

    ! -- Calculate --
    call Harm_ctraj(R_(:), P_(:), dt_*nt_, R(:), P(:), ierr)
    do it = 1, nt_
       do iit = 1, n1t_
          call DyMono_update(Harm_H_X, ierr)
       end do
    end do
    call expect_near_d(R_(1), R(1), 1.0d-2, ierr)
    call expect_near_d(P_(1), P(1), 1.0d-2, ierr)
    
    ! -- Finalize --
    call DyMono_delete(ierr)
    
  end subroutine test_1f1e
  subroutine test_1f2e_uncouple
    use Mod_DyMono
    use Mod_Uncoupled2e
    integer it, iit, ierr
    complex(kind(0d0)) H(2,2)
    ! -- Initialize --
    call DyMono_new(1, 2, ierr); check_err(ierr)
    R_(1) = 1
    P_(1) = 0
    nt_ = 1
    n1t_ = 1
    call DyMono_setup(ierr)

    call HIJ(Uncoupled2e_H_X, H, ierr)
    call expect_near_c( (0.0d0,0.0d0), H(1,2), 1.0d-10, ierr); check_err(ierr)
    call expect_near_c( (0.0d0,0.0d0), H(2,1), 1.0d-10, ierr); check_err(ierr)

    ! -- Calculate --
    do it = 1, nt_
       do iit = 1, n1t_
          call DyMono_update(Uncoupled2e_H_X, ierr)
       end do
    end do

    call expect_near_d(1.0d0, abs(c_(1)), 1.0d-10, ierr)
    call expect_near_d(0.0d0, abs(c_(2)), 1.0d-10, ierr)
    
    ! -- Finalize --
    call DyMono_delete(ierr)
    
  end subroutine test_1f2e_uncouple
end module Mod_UTestDy

program main
  use Mod_UTestDy
  call UTestDy_run
end program main
