
program main
  use Mod_DyMono
  use Mod_MoleFit
  integer ierr

  call MoleFit_new_file("xs.csv", "heij.csv", "xkij.csv", ierr)
  
  call DyMono_new(1, 2, ierr)
  R_(1) = -7
  P_(1) = 20
  m_    = 2000
  dt_   = 20.0d0
  nt_   = 100
  n1t_  = 1
  c_(1) = 1
  inte_RP_ = "RK4"
  call DyMono_setup(ierr)

  call DyMono_run(MoleFit_H_X)

  call DyMono_delete(ierr)
  
end program main
