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
