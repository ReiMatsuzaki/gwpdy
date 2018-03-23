#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_Uncoupled2e
contains
  subroutine Uncoupled2e_H_X(Q, HeIJ, XkIJ, ierr)
    implicit none
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

module Mod_UTestCl
  use Mod_Utest
  implicit none
contains
  subroutine UTestCl_run
    use Mod_Timer    
    type(Obj_Timer) :: timer
    integer ierr
    call Timer_new(timer, "UTestDy", .true., ierr)
    write(*,*)
    write(*,*) "UTestCl begin"
    write(*,*)

    call Timer_begin(timer, "uncouple", ierr)
    call test_uncouple
    call Timer_end(timer, "uncouple", ierr)

    call Timer_begin(timer, "couple", ierr)
    call test_couple
    call Timer_end(timer, "couple", ierr)    

    call Timer_begin(timer, "psa", ierr)
    call test_psa
    call Timer_end(timer, "psa", ierr)    
    
    write(*,*)
    write(*,*) "UTestCl end"
    write(*,*)    
  end subroutine UTestCl_run
  subroutine test_uncouple
    use Mod_ClPSANB
    use Mod_Uncoupled2e
    integer ierr
    double precision :: probe(2), norm2, probe0(2)
    ! -- Initialize --
    call ClPSANB_new(1, 2, 2, ierr); CHK_ERR(ierr)
    R_(1,1) = -7; R_(2,1) = -6
    P_(1,1) = +3; P_(2,1) = +2
    nt_ = 1
    n1t_ = 1
    dt_ = 2
    CC_(1) = 1;     CC_(2) = 2
    c_(1,1)= 1;     c_(2,1)= 3
    c_(1,2)= 0;     c_(2,2)= 1
    inte_RP_ = "RK4"
    call ClPSANB_setup(ierr)

    ! -- Calculate --
    call calc_probe(probe0(:), ierr); CHK_ERR(ierr)
    call update_set(Uncoupled2e_H_X, ierr); CHK_ERR(ierr)
    call calc_probe(probe(:), ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(probe0(1), probe(1), 1.d-10, ierr)
    EXPECT_NEAR_D(probe0(2), probe(2), 1.d-10, ierr)
    norm2 = sum(probe)
    EXPECT_NEAR_D(1.0d0, norm2, 1.d-10, ierr)
    
    ! -- Finalize --
    call ClPSANB_delete(ierr)
    
  end subroutine test_uncouple
  subroutine test_couple
    use Mod_ClPSANB
    use Mod_Tully1
    integer ierr
    integer, parameter :: nf=1, ne=2, np=2
    double precision :: probe(ne)
    
    ! -- Initialize --
    call ClPSANB_new(nf, ne, np, ierr); CHK_ERR(ierr)
    R_(1,1) = +0.3d0;
    P_(1,1) = +1.0d0
    c_(1,:) = (/(1.0d0,0.0d0), (0.0d0,0.0d0)/)
    CC_(1)  = 1
    R_(2,1) = +0.2d0;
    P_(2,1) = -1.2d0
    c_(2,:) = (/(0.2d0,0.0d0), (1.0d0,0.0d0)/)
    CC_(2)  = 2
    dt_ = 1.0d0 !=> norm=0.1002185667E+01
    m_  = 2000.0d0
    call ClPSANB_setup(ierr)

    ! -- update --
    call update_set(tully1_calc_H_X, ierr); CHK_ERR(ierr)
    call calc_probe(probe(:), ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, sum(probe(:)), 1.0d-13, ierr)
    
    ! -- Finalize --
    call ClPSANB_delete(ierr)
        
  end subroutine test_couple
  subroutine test_psa
    use Mod_ClPSANB
    use Mod_Tully1
    integer ierr
    integer, parameter :: nf=1, ne=2
    integer it
    double precision :: probe(ne)
    
    ! -- Initialization --
    call ClPSANB_new(nf, ne, 1, ierr); CHK_ERR(ierr)
    R_(1,1) = -5.0d0
    P_(1,1) = +10.0d0
    c_(1,:) = (/(1.0d0, 0.0d0), (0.0d0, 0.0d0)/)
    dt_ = 40.0d0
    m_ = 2000
    call ClPSANB_setup(ierr); CHK_ERR(ierr)

    ! -- Update -
    call calc_probe(probe, ierr); CHK_ERR(ierr)
    write(*,*) "prob:", probe
    write(*,*) "R:", R_(:2,:)
    write(*,*) "P:", P_(:2,:)

    write(*,*) 
    write(*,*) "SET"
    write(*,*) 
    do it = 1, 20
       call update_set(tully1_calc_H_X, ierr); CHK_ERR(ierr)
    end do

    call calc_probe(probe, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, sum(probe), 1.0d-10, ierr)
    write(*,*) "prob:", probe
    write(*,*) "R:", R_(:2,:)
    write(*,*) "P:", P_(:2,:)

    write(*,*) 
    write(*,*) "branch->PSA->finish"
    write(*,*) 
    call branch(ierr); CHK_ERR(ierr)
    do it = 1, 10
       call update_psa(tully1_calc_H_X, ierr); CHK_ERR(ierr)
    end do
    call finish_psa(ierr)
    
    call calc_probe(probe, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, sum(probe), 1.0d-10, ierr)
    write(*,*) "prob:", probe
    write(*,*) "R:", R_(:2,:)
    write(*,*) "P:", P_(:2,:)

    write(*,*) 
    write(*,*) "SET"
    write(*,*) 
    do it = 1, 20
       call update_set(tully1_calc_H_X, ierr); CHK_ERR(ierr)
    end do
    
    call calc_probe(probe, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, sum(probe), 1.0d-10, ierr)
    write(*,*) "prob:", probe
    write(*,*) "R:", R_(:2,:)
    write(*,*) "P:", P_(:2,:)

    ! -- Finalize --
    call ClPSANB_delete(ierr)
    
  end subroutine test_psa
end module Mod_UTestCl

program main
  use Mod_UTestCl
  call UTestCl_run
end program main


