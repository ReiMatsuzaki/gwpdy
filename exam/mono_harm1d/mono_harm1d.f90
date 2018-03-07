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
  dt_ = 0.1d0
  nt_ = 30
  n1t_ = 20
  c_(1) = 1
  call DyMono_setup(ierr)

  call Harm_new
  k_ = 1
  m_ = 1
  r0_ = 0

  call DyMono_run(Harm_HeIJ, Harm_XkIJ)

  call DyMono_delete(ierr)
  
end program main
