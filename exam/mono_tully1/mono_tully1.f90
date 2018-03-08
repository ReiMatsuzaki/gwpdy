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

program main
  use Mod_DyMono
  use Mod_Tully1
  integer ierr
  
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

  call DyMono_run(tully1_calc_H_X)

  call DyMono_delete(ierr)
  
end program main
