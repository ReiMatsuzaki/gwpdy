
program main
  use Mod_DyMono
  use Mod_MoleFit
  integer ierr
  !  double precision :: Q(1)=(/0.3d0/)
  !  complex(kind(0d0)) :: HeIJ(2,2)
  !  complex(kind(0d0)) :: XkIJ(1,2,2)

  call MoleFit_new_file("xs.csv", "heij.csv", "xkij.csv", ierr)
  !  call MoleFit_H_X(Q, HeIJ, XkIJ, ierr)
  
  call DyMono_new(1, 2, ierr)
  R_(1) = -7
  P_(1) = 15
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
