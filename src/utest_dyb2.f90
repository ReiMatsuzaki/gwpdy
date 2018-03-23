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
    write(*,*) "UTestDyB2 begin"
    write(*,*)

    call Timer_begin(timer, "uncouple", ierr)
    call test_uncouple
    call Timer_end(timer, "uncouple", ierr)

    call Timer_begin(timer, "couple", ierr)
    call test_couple
    call Timer_end(timer, "couple", ierr)    
    
    write(*,*)
    write(*,*) "UTestDyB2 end"
    write(*,*)    
  end subroutine UTestDy_run
  subroutine test_uncouple
    use Mod_DyBranch2
    use Mod_Uncoupled2e
    integer ierr
    double precision :: probe(2), norm2, probe0(2)
    ! -- Initialize --
    call DyBranch_new(1, 2, 2, ierr); CHK_ERR(ierr)
    R_(1,1) = -7; R_(2,1) = -6
    P_(1,1) = +3; P_(2,1) = +2
    nt_ = 1
    n1t_ = 1
    dt_ = 2
    cc_(1) = sqrt(1/3.0d0); cc_(2) = sqrt(2/3.0d0)
    c_(1,1)= 1;     c_(2,1)=0
    c_(1,2)= 0;     c_(2,2)=1
    inte_RP_ = "RK4"
    call DyBranch_setup(ierr)

    ! -- Calculate --
    call calc_probe(probe0(:), ierr); CHK_ERR(ierr)
    call DyBranch_update(Uncoupled2e_H_X, ierr); CHK_ERR(ierr)
    call calc_probe(probe(:), ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(probe0(1), probe(1), 1.d-10, ierr)
    EXPECT_NEAR_D(probe0(2), probe(2), 1.d-10, ierr)
    norm2 = sum(probe)
    EXPECT_NEAR_D(1.0d0, norm2, 1.d-10, ierr)
    
    ! -- Finalize --
    call DyBranch_delete(ierr)
    
  end subroutine test_uncouple
  subroutine test_couple
    use Mod_DyBranch2
    use Mod_PWGTO
    use Mod_Tully1
    use Mod_Math, only : vmv
    integer ierr
    integer, parameter :: nf=1, ne=2, np=2
    complex(kind(0d0)) :: HeIJ(np,ne,ne), XkIJ(np,nf,ne,ne), S(np,np), P2(nf,np,np)
    complex(kind(0d0)) :: dR(np,np), dP(np,np)
    double precision :: dotR(np,nf), dotP(np,nf)
    complex(kind(0d0)) :: HH(np,np), HK(np,ne,ne), SS(np,np), TTe(np,np), TTn(np,np)
    type(Obj_PWGTO) :: gwp
    double precision :: norm2
    integer KK
    
    ! -- Initialize --
    call DyBranch_new(nf, ne, np, ierr); CHK_ERR(ierr)
    R_(1,1) = +0.3d0;
    P_(1,1) = +1.0d0
    c_(1,:) = (/(1.0d0,0.0d0), (0.0d0,0.0d0)/)
    R_(2,1) = +0.2d0;
    P_(2,1) = -1.2d0
    c_(2,:) = (/(0.2d0,0.0d0), (1.0d0,0.0d0)/)
    dt_ = 1.0d0 !=> norm=0.1002185667E+01
    m_  = 2000.0d0
    call DyBranch_setup(ierr)

    ! -- nuclear part --
    call make_GWP("0RP", gwp, ierr); CHK_ERR(ierr)
    call PWGTO_overlap(gwp,   1, 1, S, ierr)
    call PWGTO_overlap(gwp,   1, 2, dR(:,:),ierr); CHK_ERR(ierr)
    call PWGTO_overlap(gwp,   1, 3, dP(:,:),ierr); CHK_ERR(ierr)
    call PWGTO_kineticp2(gwp, 1, 1, P2(1,:,:), ierr)
    call PWGTO_delete(gwp, ierr)

    write(*,*) S

    ! -- electron part --
    do KK = 1, npath_
       call dot_RP(Tully1_calc_H_X, KK, dotR(KK,:), dotP(KK,:), ierr)
       call Tully1_calc_H_X(R_(KK,:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr)
       call local_HIJ(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P_(KK,:), HK(KK,:,:), ierr)
    end do

    ! -- global Hamiltonian --
    call global_HIJ(HeIJ(:,:,:), XkIJ(:,:,:,:), S(:,:), p2(:,:,:), &
         SS(:,:), HH(:,:), ierr)
    call global_T_el(HK, S, TTe, ierr)
    call global_T_nu(dotR, dotP, dR, dP, TTn, ierr)

    ! -- check Hamiltonian --    
    EXPECT_EQ_D(0.0d0, aimag(HH(2,2)), ierr)
    EXPECT_EQ_D(0.0d0, aimag(HH(1,1)), ierr)
    EXPECT_EQ_C(HH(1,2), conjg(HH(2,1)), ierr)

    !EXPECT_EQ_D(0.0d0,   real(TT(1,1)),  ierr)
    !EXPECT_EQ_D(0.0d0,   real(TT(2,2)),  ierr)
    !EXPECT_EQ_C(TT(1,2), -conjg(TT(2,1)), ierr)

    ! -- update --
    call Update1st(Tully1_calc_H_X, ierr); CHK_ERR(ierr)
    call calc_norm2(norm2, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, norm2, 3.0d-3, ierr)
    
    ! -- Finalize --
    call DyBranch_delete(ierr)
        
  end subroutine test_couple
end module Mod_UTestDy

program main
  use Mod_UTestDy
  call UTestDy_run
end program main
