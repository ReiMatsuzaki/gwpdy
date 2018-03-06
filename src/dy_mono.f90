#include "macros.fpp"

module Mod_DyMono
  ! molecular wave function \Psi are expanded as 
  !     \Psi(r,Q,t) = G(Q,t) \sum_I c_I(t) \Phi_I(r;Q) .
  ! G is auss wave packet
  use Mod_GWP
  implicit none
  integer :: nf_        ! number of freedom
  integer :: nebasis_   ! number of electron Hilbert space
  type(Obj_GWP) :: gwp_
  complex(kind(0d0)), allocatable :: c_(:)
  integer :: nd_        ! order of numerical differate
  double precision :: dt_, dR_, dP_ ! finite difference for t R and P
contains
  subroutine DyMono_run
    
  end subroutine DyMono_run
  subroutine DyMono_new(nf, nebasis, ierr)
    integer, intent(in) :: nf, nebasis
    integer, intent(out) :: ierr
    ierr = 0
    nf_  = nf
    nebasis_ = nebasis
    call GWP_new(gwp_, nf, 1, 'c', ierr) ; check_err(ierr)
    allocate(c_(nebasis))
    c_ = 0

    nd_ = 2
    dt_ = 0.1d0
    dR_ = 0.01d0
    dP_ = 0.01d0
    
  end subroutine DyMono_new
  subroutine DyMono_setup(ierr)
    integer, intent(out) :: ierr
    double precision :: norm2
    double precision, parameter :: tol = 1.0d-15
    ierr = 0
    call GWP_setup(gwp_, ierr); check_err(ierr)
    norm2 = sum(abs(c_(:))**2)
    if(norm2<tol) then
       begin_err("norm is too small")
       ierr = 1
       write(0,*) "norm2:", norm2
       return
    end if
    c_(:) = c_(:) / sqrt(norm2)
  end subroutine DyMono_setup
  subroutine DyMono_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    call GWP_delete(gwp_, ierr); check_err(ierr)
    deallocate(c_)
  end subroutine DyMono_delete
  subroutine DyMono_update(HIJ, ierr)
    use Mod_const, only : II
    use Mod_math, only : lapack_zheev
    interface
       subroutine HIJ(R,P,res, ierr)
         double precision, intent(in) :: R(:),P(:)
         complex(kind(0d0)), intent(out) :: res(:,:)
         integer, intent(out) :: ierr
       end subroutine HIJ
    end interface
    integer, intent(out) :: ierr
    complex(kind(0d0)), allocatable :: dH_dR(:,:), dH_dP(:,:), H1(:,:), H2(:,:)
    integer k
    double precision, allocatable :: R(:), P(:), dotR(:), dotP(:)
    double precision :: w(nebasis_)
    complex(kind(0d0)) :: U(nebasis_, nebasis_), UH(nebasis_, nebasis_)

    ierr = 0

    allocate(dH_dR(nebasis_,nebasis_), dH_dP(nebasis_,nebasis_))
    allocate(H1(nebasis_,nebasis_),    H2(nebasis_,nebasis_))
    allocate(R(nf_), P(nf_), dotR(nf_), dotP(nf_))

    if(nd_.eq.2) then
       do k = 1, nf_
          R(:) = gwp_%R(1,:); R(k) = R(k) + dR_
          P(:) = gwp_%P(1,:)
          call HIJ(R,P,H1(:,:), ierr)
          R(:) = gwp_%R(1,:); R(k) = R(k) - dR_ 
          call HIJ(R,P,H2(:,:), ierr)
          dH_dR(:,:) = (H1(:,:) - H2(:,:)) / (2*dR_)
          dotR(k) = real(dot_product(c_(:), matmul(dH_dR(:,:), c_(:))))

          R(:) = gwp_%R(1,:)
          P(:) = gwp_%P(1,:); P(k) = P(k) + dP_
          call HIJ(R,P,H1(:,:), ierr)
          P(:) = gwp_%P(1,:); P(k) = P(k) - dP_ 
          call HIJ(R,P,H2(:,:), ierr)
          dH_dP(:,:) = (H1(:,:) - H2(:,:)) / (2*dP_)
          dotP(k) = real(dot_product(c_(:), matmul(dH_dP(:,:), c_(:))))
       end do
    else
       begin_err("invalid nd_")
       ierr = 1
    end if

    call HIJ(gwp_%R(1,:), gwp_%P(1,:), H1(:,:), ierr); check_err(ierr)
    call lapack_zheev(nebasis_, H1(:,:), w(:), U(:,:), ierr); check_err(ierr)
    UH(:,:) = conjg(transpose(U(:,:)))
    c_(:) = matmul(UH(:,:), c_(:))
    c_(:) = exp(-II*w(:)*dt_) * c_(:)
    c_(:) = matmul(U(:,:), c_(:))

    gwp_%R(1,:) = gwp_%R(1,:) + dotR(:) * dt_
    gwp_%P(1,:) = gwp_%P(1,:) + dotP(:) * dt_

    deallocate(dH_dR, dH_dP, H1, H2)
    deallocate(R, P, dotR, dotP)
    
  end subroutine DyMono_update
end module Mod_DyMono
