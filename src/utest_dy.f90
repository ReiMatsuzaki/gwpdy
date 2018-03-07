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
  subroutine Harm_HeIJ(Q, res, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    ierr = 0
    res(1,1) = k_/2*(Q(1)-r0_)**2
  end subroutine Harm_HeIJ
  subroutine Harm_XkIJ(Q, res, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: res(:,:,:)
    integer, intent(out) :: ierr
    ierr = 0
    res = 0*Q(1)
  end subroutine Harm_XkIJ
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

module Mod_UTestDy
  use Mod_Utest
  implicit none
contains
  subroutine UTestDy_run
    call UTestDy_mono_1f1e
  end subroutine UTestDy_run
  subroutine UTestDy_mono_1f1e
    use Mod_DyMono, only : gwp_, c_, dt_, &
         DyMono_new, DyMono_setup, DyMono_update, DyMono_delete
    use Mod_Harm, only : Harm_new, Harm_HeIJ, Harm_XkIJ, Harm_ctraj, k_, m_, r0_
    integer, parameter :: nf = 1, nebasis=1, nt=100
    double precision, parameter :: eps = 1.0d-10
    integer ierr, it
    double precision :: R(nf), P(nf)

    ! -- Initialize --
    call DyMono_new(nf, nebasis, ierr); check_err(ierr)
    gwp_%g(1,1,1) = 1
    gwp_%R(1,1)   = 1
    gwp_%P(1,1)   = 1
    gwp_%c(1)     = 1
    dt_ = 0.01d0
    c_(1) = 1
    call DyMono_setup(ierr)
    call Harm_new
    k_ = 1
    m_ = 1
    r0_ = 0

    ! -- Calculate --
    call Harm_ctraj(gwp_%R(1,:), gwp_%P(1,:), dt_*nt, R(:), P(:), ierr)
    do it = 1, nt
       call DyMono_update(Harm_HeIJ, Harm_XkIJ, ierr)
    end do
    call expect_near_d(gwp_%R(1,1), R(1), 1.0d-2, ierr)
    call expect_near_d(gwp_%P(1,1), P(1), 1.0d-2, ierr)
    
    ! -- Finalize --
    call DyMono_delete(ierr)
    
  end subroutine UTestDy_mono_1f1e
end module Mod_UTestDy

program main
  use Mod_UTestDy
  call UTestDy_run
end program main
