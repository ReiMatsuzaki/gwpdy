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
end module Mod_harm

program main
  use Mod_DyMono
  use Mod_Harm
  integer ierr
  
  call DyMono_new(1, 1, ierr)
  gwp_%g(1,1,1) = 1
  gwp_%R(1,1)   = 1
  gwp_%P(1,1)   = 1
  gwp_%c(1)     = 1
  dt_ = 0.01d0
  nt_ = 30*4
  n1t_ = 5*10
  c_(1) = 1
  inte_RP_ = "RK4"
  call DyMono_setup(ierr)

  call Harm_new
  k_ = 1
  m_ = 1
  r0_ = 0

  call DyMono_run(Harm_H_X)

  call DyMono_delete(ierr)
  
end program main
