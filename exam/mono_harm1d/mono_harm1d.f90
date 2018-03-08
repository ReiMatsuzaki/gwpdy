module Mod_Harm
  implicit none
contains
  subroutine Harm_H_X(Q, HeIJ, XkIJ, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
    integer, intent(out) :: ierr
    double precision k, r0

    k = 1
    r0 = 0
    
    ierr = 0
    HeIJ(1,1) = k/2*(Q(1)-r0)**2
    XkIJ(1,1,1) = 0.0d0
  end subroutine Harm_H_X
end module Mod_harm

program main
  use Mod_DyMono
  use Mod_Harm
  integer ierr
  
  call DyMono_new(1, 1, ierr)
  R_(1)   = 1
  P_(1)   = 1
  c_(1)   = 1
  nt_ = 100
  dt_ = 0.3d0  
  n1t_ = 10
  inte_RP_ = "RK4"
  !  inte_RP_ = "euler"
  call DyMono_setup(ierr)

  call DyMono_run(Harm_H_X)

  call DyMono_delete(ierr)
  
end program main
