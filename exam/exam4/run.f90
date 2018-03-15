#include "macros.fpp"

subroutine run
  use Mod_DyBranch
  use Mod_MoleFit
  integer ierr

  write(*,*) "exam4"

  call MoleFit_new_file("xs.csv", "heij.csv", "xkij.csv", ierr)
  
  call DyBranch_new(1, 2, 1, ierr)
  R_(1, 1) = -7
  P_(1, 1) = 15
  m_    = 2000
  dt_   = 20.0d0
  nt_   = 100
  n1t_  = 1
  cc_(1) = 1
  c_(1,1) = 1
  c_(1,2) = 0
  inte_RP_ = "RK4"
  call DyBranch_setup(ierr)

  call write_input(ierr)
  call print_input

  do it = 0, nt_
     call write_status(it)
     call print_status(it)

     if(it.eq.nt_/2) then
        call branch(MoleFit_H_X, ierr); CHK_ERR(ierr)
     end if
     
     do i1t = 1, n1t_
        call DyBranch_update(MoleFit_H_X, ierr); CHK_ERR(ierr)
        if(ierr.ne.0) then
           write(0,*) "Error"
           write(0,*) "(it,i1t):", it, i1t
           write(0,*) "R_:", R_(:npath_,:)
           write(0,*) "P_:", P_(:npath_,:)
           write(0,*) "c_:", c_(:npath_,:)
           stop
        end if
     end do
  end do
  
  !  call DyBranch_run(MoleFit_H_X)

  call DyBranch_delete(ierr)
end subroutine run
program main
  call run
end program main
