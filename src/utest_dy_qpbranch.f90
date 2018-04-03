#include "macros.fpp"
#include "macros_utest.fpp"


module Mod_UTestDyQPBranch
  use Mod_Utest
  use Mod_DyQPBranch
  implicit none
contains
  subroutine run
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

    call Timer_begin(timer, "psa", ierr)
    call test_psa
    call Timer_end(timer, "psa", ierr)
    
    
    write(*,*)
    write(*,*) "UTestDyB2 end"
    write(*,*)    
  end subroutine run
  subroutine test_uncouple
    use Mod_Uncoupled2e
    integer ierr
    double precision :: probe(2), norm2, probe0(2)
    ! -- Initialize --
    call DyQPBranch_new(1, 2, 1, 2, ierr); CHK_ERR(ierr)
    np_ = 2
    R_(1,1,1) = -7; R_(2,1,1) = -6
    P_(1,1,1) = +3; P_(2,1,1) = +2
    nt_ = 1
    n1t_ = 1
    dt_ = 2
    Ctot_(1) = sqrt(1/3.0d0);
    Ctot_(2) = sqrt(2/3.0d0)
    c_(1,1,:) = (/1,0/)
    c_(2,1,:) = (/0,1/)
    nnuc_(:) = 1
    inte_RP_ = "RK4"
    call DyQPBranch_setup(ierr); CHK_ERR(ierr)

    ! -- Calculate --
    call calc_probe(probe0(:), ierr); CHK_ERR(ierr)
    call update_set(Uncoupled2e_H_X, 1, ierr); CHK_ERR(ierr)
    call update_set(Uncoupled2e_H_X, 2, ierr); CHK_ERR(ierr)
    call calc_probe(probe(:), ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(probe0(1), probe(1), 1.d-10, ierr)
    EXPECT_NEAR_D(probe0(2), probe(2), 1.d-10, ierr)
    norm2 = sum(probe)
    EXPECT_NEAR_D(1.0d0, norm2, 1.d-10, ierr)
    
    ! -- Finalize --
    call DyQPBranch_delete(ierr)
    
  end subroutine test_uncouple
  subroutine test_couple
    use Mod_Tully1
    use Mod_Math, only : vmv
    integer ierr
    integer, parameter :: nf=1, ne=2, np=2
    double precision :: norm2
    
    ! -- Initialize --
    call DyQPBranch_new(nf, ne, 1, np, ierr); CHK_ERR(ierr)
    R_(1,1,1) = +0.3d0;
    P_(1,1,1) = +1.0d0
    Ctot_(1) = 1
    c_(1,1,:) = (/(1.0d0,0.0d0), (0.0d0,0.0d0)/)
    Ctot_(1) = 2
    
    R_(2,1,1) = +0.2d0;
    P_(2,1,1) = -1.2d0
    c_(2,1,:) = (/(0.2d0,0.0d0), (1.0d0,0.0d0)/)
    Ctot_(2) = 1
    
    dt_ = 1.0d0 !=> norm=0.1002185667E+01
    m_  = 2000.0d0
    call DyQPBranch_setup(ierr)

    ! -- update --
    call update_set(Tully1_calc_H_X, 1, ierr); CHK_ERR(ierr)
    call update_set(Tully1_calc_H_X, 1, ierr); CHK_ERR(ierr)
    call calc_norm2(norm2, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, norm2, 3.0d-3, ierr)
    
    ! -- Finalize --
    call DyQPBranch_delete(ierr)
        
  end subroutine test_couple
  subroutine test_psa
    use Mod_DyQPBranchPSA
    use Mod_Tully1
    use Mod_Math, only : vmv
    integer ierr
    integer, parameter :: nf=1, ne=2, np=2
    double precision :: norm2
    
    ! -- Initialize --
    call DyQPBranch_new(nf, ne, 1, np, ierr); CHK_ERR(ierr)
    call DyQPBRanchPSA_new(2, ierr); CHK_ERR(ierr)
    np_ = 1    
    R_(1,1,1) = +0.3d0;
    P_(1,1,1) = +1.0d0
    Ctot_(1) = 1
    c_(1,1,:) = (/(1.0d0,0.0d0), (0.0d0,0.0d0)/)
    Ctot_(1) = 1
    gauss_mode_ = "thawed"
    
    dt_ = 1.0d0 !=> norm=0.1002185667E+01
    m_  = 2000.0d0
    call DyQPBranch_setup(ierr)
    call DyQPBranchPSA_setup(ierr)

    ! -- update --
    call update_set(Tully1_calc_H_X, 1, ierr); CHK_ERR(ierr)
    call begin_clpsa(1, ierr); CHK_ERR(ierr)
    call update_clpsa(Tully1_calc_H_X, 1, ierr)
    call end_clpsa(1, ierr); CHK_ERR(ierr)
    !    call DyQPBranch_dump(6, ierr)
    call update_set(Tully1_calc_H_X, 1, ierr); CHK_ERR(ierr)
    !    call DyQPBranch_dump(6, ierr)
    call calc_norm2(norm2, ierr); CHK_ERR(ierr)
    EXPECT_NEAR_D(1.0d0, norm2, 3.0d-5, ierr)
    
    ! -- Finalize --
    call DyQPBranch_delete(ierr)
  end subroutine test_psa
end module Mod_UTestDyQPBranch

program main
  use Mod_UTestDyQPBranch
  call run
end program main
