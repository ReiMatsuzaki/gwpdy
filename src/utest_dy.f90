#include "macros.fpp"
#include "macros_utest.fpp"

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

module Mod_Tully1
contains
  subroutine tully1_calc_H_X(Q, HeIJ, XkIJ, ierr)
    implicit none
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
    integer, intent(out) :: ierr
    double precision A,B,C,D,x
    x = Q(1)
    A = 0.01d0
    B = 1.6d0
    C = 0.005d0
    D = 1.0d0
    
    ierr = 0
    HeIJ(:,:) = 0
    XkIJ(:,:,:) = 0
    
    if(x>0) then
       HeIJ(1,1) = A*(1-exp(-B*x))
    else
       HeIJ(1,1) = -A*(1-exp(B*x))
    end if
    HeIJ(1,2) = C*exp(-D*x**2)
    HeIJ(2,1) = HeIJ(1,2)
    HeIJ(2,2) = -HeIJ(1,1)
    
  end subroutine tully1_calc_H_X
end module Mod_Tully1

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

    call Timer_begin(timer, "1f2e_fit", ierr)
    call test_1f2e_fit
    call Timer_end(timer, "1f2e_fit", ierr)

    call Timer_begin(timer, "branch_1f2e", ierr)
    call test_branch_1f2e
    call Timer_end(timer, "branch_1f2e", ierr)
    
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
    call DyMono_new(nf, ne, ierr); CHK_ERR(ierr)
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
    EXPECT_NEAR_D(R_(1), R(1), 1.0d-2, ierr)
    EXPECT_NEAR_D(P_(1), P(1), 1.0d-2, ierr)
    
    ! -- Finalize --
    call DyMono_delete(ierr)
    
  end subroutine test_1f1e
  subroutine test_1f2e_uncouple
    use Mod_DyMono
    use Mod_Uncoupled2e
    integer it, iit, ierr
    complex(kind(0d0)) H(2,2)
    ! -- Initialize --
    call DyMono_new(1, 2, ierr); CHK_ERR(ierr)
    R_(1) = 1
    P_(1) = 0
    nt_ = 1
    n1t_ = 1
    call DyMono_setup(ierr)

    call HIJ(Uncoupled2e_H_X, H, ierr)
    EXPECT_NEAR_C((0.0d0,0.0d0), H(1,2), 1.0d-10, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_C((0.0d0,0.0d0), H(2,1), 1.0d-10, ierr); CHK_ERR(ierr)

    ! -- Calculate --
    do it = 1, nt_
       do iit = 1, n1t_
          call DyMono_update(Uncoupled2e_H_X, ierr)
       end do
    end do

    EXPECT_NEAR_D(1.0d0, abs(c_(1)), 1.0d-10, ierr)
    EXPECT_NEAR_D(0.0d0, abs(c_(2)), 1.0d-10, ierr)
    
    ! -- Finalize --
    call DyMono_delete(ierr)
    
  end subroutine test_1f2e_uncouple
  subroutine test_1f2e_fit
    use Mod_DyMono
    use Mod_MoleFit
    use Mod_Tully1
    use Mod_sys, only : open_w
    integer, parameter :: ifile_xs = 29831, ifile_h = 29832, ifile_xk = 29833
    integer ierr, ix, iit, it
    integer, parameter :: nx = 100
    double precision :: x(1)
    complex(kind(0d0)) :: HeIJ(2,2), XkIJ(1,2,2)
    double precision :: R0(1), P0(1), R1(1), P1(1)
    complex(kind(0d0)) :: c0(2), c1(2)

    ! -- write data --
    call open_w(ifile_xs, "data/mole_tully1_xs.csv", ierr)
    call open_w(ifile_h,  "data/mole_tully1_heij.csv", ierr)
    call open_w(ifile_xk, "data/mole_tully1_xkij.csv", ierr)
    write(ifile_xs, '("val")')
    write(ifile_h , '("i,j,k,val")')
    write(ifile_xk, '("i,j,k,l,val")')
    do ix = 1, nx
       x(1) = (ix-nx/2) * 0.2d0
       call Tully1_calc_H_X(x(:), HeIJ, XkIJ, ierr)
       write(ifile_xs,   '(f20.10)') x(1)
       write(ifile_h , '(i0,",",i0,",",i0,",",f20.10)') ix,1,1,real(HeIJ(1,1))
       write(ifile_h , '(i0,",",i0,",",i0,",",f20.10)') ix,1,2,real(HeIJ(1,2))
       write(ifile_h , '(i0,",",i0,",",i0,",",f20.10)') ix,2,2,real(HeIJ(2,2))
       write(ifile_xk, '(i0,",",i0,",",i0,",",i0,",",f20.10)') ix,1,1,2,real(XkIJ(1,1,2))
    end do
    close(ifile_xs)
    close(ifile_h)
    close(ifile_xk)

    ! -- build from file --
    call MoleFit_new_file("data/mole_tully1_xs.csv", &
         "data/mole_tully1_heij.csv", "data/mole_tully1_xkij.csv", &
         ierr); CHK_ERR(ierr)
    
    ! -- fitting calc --
    call DyMono_new(1, 2, ierr); CHK_ERR(ierr)
    R_(1) = -2
    P_(1) = 5
    nt_ = 10
    n1t_ = 2
    call DyMono_setup(ierr)
    do it = 1, nt_
       do iit = 1, n1t_
          call DyMono_update(MoleFit_H_X, ierr); CHK_ERR(ierr)
       end do
    end do
    R0=R_; P0=P_; c0=c_
    call DyMono_delete(ierr)

    ! -- original calc --
    call DyMono_new(1, 2, ierr); CHK_ERR(ierr)
    R_(1) = -2
    P_(1) = 5
    nt_ = 10
    n1t_ = 2
    call DyMono_setup(ierr)
    do it = 1, nt_
       do iit = 1, n1t_
          call DyMono_update(Tully1_calc_H_X, ierr) ; CHK_ERR(ierr)
       end do
    end do
    R1=R_; P1=P_; c1=c_
    call DyMono_delete(ierr)

    EXPECT_NEAR_D(R0(1), R1(1), 1.0d-4, ierr)
    EXPECT_NEAR_D(P0(1), P1(1), 1.0d-4, ierr)
    EXPECT_NEAR_C(c0(1), c1(1), 2.0d-4, ierr)
    EXPECT_NEAR_C(c0(2), c1(2), 1.0d-4, ierr)
    
  end subroutine test_1f2e_fit
  subroutine test_branch_1f2e
    use Mod_DyBranch
    use Mod_Tully1
    integer it, iit, ierr
    ! -- Initialize --
    call DyBranch_new(1, 2, 1, ierr); CHK_ERR(ierr)
    R_(1,1) = 1
    P_(1,1) = 0
    nt_ = 10
    n1t_ = 1
    call DyBranch_setup(ierr)

    ! -- Calculate --
    do it = 1, nt_
       !       if(it.eq.5) then
       !          call branch(Tully1_calc_H_X, ierr); CHK_ERR(ierr)
       !       end if
       do iit = 1, n1t_
          !      call DyBranch_update(Tully1_calc_H_X, ierr)
          call DyBranch_update(Tully1_calc_H_X, ierr); CHK_ERR(ierr)
       end do
    end do

    EXPECT_NEAR_D(1.0d0, abs(cc_(1)), 1.0d-10, ierr)
    EXPECT_NEAR_D(0.0d0, abs(cc_(2)), 1.0d-10, ierr)
    
    ! -- Finalize --
    call DyBranch_delete(ierr)
    
  end subroutine test_branch_1f2e
  subroutine test_branch_1f2e_branch
    use Mod_DyBranch
    use Mod_Tully1
    integer it, iit, ierr
    ! -- Initialize --
    call DyBranch_new(1, 2, 1, ierr); CHK_ERR(ierr)
    R_(1,1) = -0.1d0
    P_(1,1) = +10.0d0
    nt_ = 10
    n1t_ = 1
    call DyBranch_setup(ierr)

    call branch(Tully1_calc_H_X, ierr); CHK_ERR(ierr)

    ! -- Calculate --
    do it = 1, nt_
       do iit = 1, n1t_
          !      call DyBranch_update(Tully1_calc_H_X, ierr)
          call DyBranch_update(Tully1_calc_H_X, ierr); CHK_ERR(ierr)
       end do
    end do

    EXPECT_NEAR_D(1.0d0, abs(cc_(1)), 1.0d-10, ierr)
    EXPECT_NEAR_D(0.0d0, abs(cc_(2)), 1.0d-10, ierr)
    
    ! -- Finalize --
    call DyBranch_delete(ierr)
        
  end subroutine test_branch_1f2e_branch
end module Mod_UTestDy

program main
  use Mod_UTestDy
  call UTestDy_run
end program main
