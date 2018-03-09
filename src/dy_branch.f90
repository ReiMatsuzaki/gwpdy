#include "macros.fpp"
module Mod_DyBranch
  ! molecular wave function are expressed as
  !     Psi(r,Q,t) = \sum_K C_K(t) G[Q,x_K(t)] \phi_K[r;Q,c_K(t)]
  !
  ! variables
  !   K : index of path
  !   I : index of electron basis
  !   i : dimension
  implicit none
  integer nf_    ! # of freedom
  integer ne_    ! # of electron basis
  integer npath_ ! # of path
  double precision m_ ! mass
  double precision, allocatable :: g_(:), R_(:,:), P_(:,:), theta_(:) ! GWP
  complex(kind(0d0)), allocatable :: CC_(:), c_(:,:) ! coefficient
  integer :: nd_                    ! order of numerical differate
  character(8) :: inte_RP_          ! euler or RK4
  double precision :: dydt_, dR_, dP_ ! finite difference for t, R and P
  double precision :: dt_
  integer :: nt_, n1t_              ! number of time for print, number for each step
  integer, parameter :: MAX_PATH=100
  private :: dydt_  
contains
  ! ==== IO ====
  ! ==== Main ====
  subroutine DyBranch_new(nf, ne, npath, ierr)
    integer, intent(in) :: nf, ne, npath
    integer, intent(out) :: ierr
    ierr = 0
    
    nf_ = nf
    ne_ = ne
    npath_ = npath

    allocate(g_(MAX_PATH), R_(MAX_PATH,nf), P_(MAX_PATH,nf), theta_(MAX_PATH))
    g_(:) = 1
    R_(:,:) = 0
    P_(:,:) = 0
    theta_(:) = 0

    allocate(CC_(MAX_PATH), c_(MAX_PATH, ne))
    CC_(:) = 0
    c_(:,:) = 0
    CC_(1) = 1
    c_(:,1) = 1

    m_ = 1

    nd_ = 2
    inte_RP_ = "euler"
    dt_ = 1.0d0
    nt_ = 10
    n1t_ = 1
    dR_ = 0.01d0
    dP_ = 0.01d0
    
  end subroutine DyBranch_new
  subroutine DyBranch_setup(ierr)
    integer, intent(out) :: ierr
    integer :: K
    double precision :: norm2
    double precision, parameter :: tol = 1.0d-10
    ierr = 0
    do K = 1, npath_
       norm2 = sum(abs(c_(K,:))**2)
       if(norm2 < tol) then
          MSG_ERR("norm is too small")
          ierr = 1
          write(0,*) "K:", K
          return
       end if
       c_(K,:) = c_(K,:) / sqrt(norm2)
    end do

    if(inte_RP_.ne."euler" .and. inte_RP_.ne."RK4") then
       MSG_ERR("unsupported inte_RP_")
       ierr = 1
       write(0,*) "inte_RP_:", inte_RP_
       return
    end if

    dydt_ = dt_/n1t_    
    
  end subroutine DyBranch_setup
  subroutine DyBranch_update(calc_H_X, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    double precision :: dR(npath_,nf_), dP(npath_,nf_)

    ierr = 0

    if(inte_RP_.eq."euler") then
       call inte_nuc_euler(calc_H_X, dR(:,:), dP(:,:), ierr); CHK_ERR(ierr)
    else if(inte_RP_.eq."RK4") then
       call inte_nuc_RK4(  calc_H_X, dR(:,:), dP(:,:), ierr); CHK_ERR(ierr)
    end if

    call inte_ele_diag(calc_H_X, ierr)

    R_(1:npath_,:) = R_(1:npath_,:) + dR(:,:)
    P_(1:npath_,:) = P_(1:npath_,:) + dP(:,:)
    
  end subroutine DyBranch_update
  subroutine DyBranch_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(g_, R_, P_, theta_)
    deallocate(CC_, c_)
  end subroutine DyBranch_delete
  ! ==== Calc ====
  subroutine inte_nuc_euler(calc_H_X, dR, dP, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    double precision, intent(out) :: dR(:,:), dP(:,:)
    double precision :: dotR(nf_), dotP(nf_)
    integer, intent(out) :: ierr
    integer KK
    do KK = 1, npath_
       call dot_RP(calc_H_X, KK, dotR, dotP, ierr)
       dR(KK,:) = dotR(:)*dydt_
       dP(KK,:) = dotP(:)*dydt_
    end do
  end subroutine inte_nuc_euler
  subroutine inte_nuc_RK4(calc_H_X, dR, dP, ierr)
    interface
        subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    double precision, intent(out) :: dR(:,:), dP(:,:)
    integer, intent(out) :: ierr
    double precision :: kR(4,nf_), kP(4,nf_)
    integer KK
    ierr = 0
    do KK = 1, npath_
       call dot_RP(calc_H_X, KK, kR(1,:), kP(1,:), ierr)
    
       R_(KK,:) = R_(KK,:) + kR(1,:) * dydt_/2
       P_(KK,:) = P_(KK,:) + kP(1,:) * dydt_/2
       call dot_RP(calc_H_X, KK, kR(2,:), kP(2,:), ierr)
       R_(KK,:) = R_(KK,:) - kR(1,:) * dydt_/2
       P_(KK,:) = P_(KK,:) - kP(1,:) * dydt_/2
       
       R_(KK,:) = R_(KK,:) + kR(2,:) * dydt_/2
       P_(KK,:) = P_(KK,:) + kP(2,:) * dydt_/2
       call dot_RP(calc_H_X, KK, kR(3,:), kP(3,:), ierr)
       R_(KK,:) = R_(KK,:) - kR(2,:) * dydt_/2
       P_(KK,:) = P_(KK,:) - kP(2,:) * dydt_/2
       
       R_(KK,:) = R_(KK,:) + kR(3,:) * dydt_/2
       P_(KK,:) = P_(KK,:) + kP(3,:) * dydt_/2
       call dot_RP(calc_H_X, KK, kR(3,:), kP(4,:), ierr)
       R_(KK,:) = R_(KK,:) - kR(3,:) * dydt_/2
       P_(KK,:) = P_(KK,:) - kP(3,:) * dydt_/2    
       
       dR(KK,:) =  (kR(1,:) + 2*kR(2,:) + 2*kR(3,:) + kR(4,:)) * dydt_/6
       dP(KK,:) =  (kP(1,:) + 2*kP(2,:) + 2*kP(3,:) + kP(4,:)) * dydt_/6    
    end do
    
  end subroutine inte_nuc_RK4
  subroutine branch(calc_H_X, ierr)
    use Mod_math, only : lapack_zheev
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    integer K, lam
    double precision ::  w(ne_)
    complex(kind(0d0)) :: H(ne_, ne_), U(ne_, ne_), c0(ne_)
    complex(kind(0d0)) :: CC0
    
!    K = 1
    !    do i = 1, nf_
!       R_(i) = R_(i) + dR_
!       call HIJ(calc_H_X, K, HK, ierr); CHK_ERR(ierr)
!       call lapack_zheev(ne_, H(:,:), w1(:), U(:,:), ierr); CHK_ERR(ierr)
!       R_(i) = R_(i) - 2*dR_
!       call HIJ(calc_H_X, K, HK, ierr); CHK_ERR(ierr)
!       call lapack_zheev(ne_, H(:,:), w2(:), U(:,:), ierr); CHK_ERR(ierr)
!       R_(i) = R_(i) + dR_
!       dH_dR(i,:) = (w1(:)-w2(:))/(2*dR_)
!       
!       P_(i) = P_(i) + dP_
!       call HIJ(calc_H_X, K, HK, ierr); CHK_ERR(ierr)
!       call lapack_zheev(ne_, H(:,:), w1(:), U(:,:), ierr); CHK_ERR(ierr)
!       P_(i) = P_(i) - 2*dP_
!       call HIJ(calc_H_X, K, HK, ierr); CHK_ERR(ierr)
!       call lapack_zheev(ne_, H(:,:), w2(:), U(:,:), ierr); CHK_ERR(ierr)
!       P_(i) = P_(i) + dP_
!       dH_dP(i,:) = (w1(:)-w2(:))/(2*dP_)
!    end do
    ierr = 0
    K = 1
    call HIJ(calc_H_X, K, H, ierr); CHK_ERR(ierr)
    call lapack_zheev(ne_, H(:,:), w(:), U(:,:), ierr); CHK_ERR(ierr)
    !    UH(:,:) = conjg(transpose(U(:,:)))
    
    npath_ = ne_
    c0(:) = c_(K,:)
    CC0   = CC_(K)
    !    R0(:) = R_(K,:)
    !    P0(:) = P_(K,:)

    do lam = 1, ne_
       c_(lam,:) = U(:,lam)
       CC_(lam) = CC0 * c0(lam)
       !       R_(lam,:)   = R0(:) + dH_dP(:) * dydt_
       !       P_(lam,:)   = P0(:) - dH_dR(:) * dydt_
    end do
    
  end subroutine branch
  subroutine inte_ele_diag(calc_H_X, ierr)
    use Mod_const, only : II
    use Mod_math, only : lapack_zheev
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: U(ne_, ne_), UH(ne_, ne_), H(ne_,ne_)
    double precision :: w(ne_)
    integer K
    ierr = 0

    do K = 1, npath_
       call HIJ(calc_H_X, K, H(:,:), ierr); CHK_ERR(ierr)
       call lapack_zheev(ne_, H(:,:), w(:), U(:,:), ierr); CHK_ERR(ierr)
       UH(:,:) = conjg(transpose(U(:,:)))
       c_(K,:) = matmul(UH(:,:), c_(K,:))
       c_(K,:) = exp(-II*w(:)*dydt_) * c_(K,:)
       c_(K,:) = matmul(U(:,:), c_(K,:))
    end do
    
  end subroutine inte_ele_diag
  subroutine dot_RP(calc_H_X, K, dotR, dotP, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in) :: K
    double precision, intent(out) :: dotR(:), dotP(:)
    integer, intent(out) :: ierr
    integer i
    complex(kind(0d0)) :: H1(ne_,ne_), H2(ne_,ne_)
    complex(kind(0d0)) :: dH_dR(ne_,ne_), dH_dP(ne_,ne_)
    if(nd_.eq.2) then
       do i = 1, nf_
          R_(K,i) = R_(K,i) + dR_
          call HIJ(calc_H_X, K, H1(:,:), ierr)
          R_(K,i) = R_(K,i) - 2*dR_
          call HIJ(calc_H_X, K, H2(:,:), ierr)
          R_(K,i) = R_(K,i) + dR_
          dH_dR(:,:) = (H1(:,:) - H2(:,:)) / (2*dR_)
          dotP(i) = -real(dot_product(c_(K,:), matmul(dH_dR(:,:), c_(K,:))))

          P_(K,i) = P_(K,i) + dP_
          call HIJ(calc_H_X, K, H1(:,:), ierr)
          P_(K,i) = P_(K,i) -2*dP_
          call HIJ(calc_H_X, K, H2(:,:), ierr)
          P_(K,i) = P_(K,i) + dP_
          dH_dP(:,:) = (H1(:,:) - H2(:,:)) / (2*dP_)
          dotR(i) = real(dot_product(c_(K,:), matmul(dH_dP(:,:), c_(K,:))))
       end do
    else
       MSG_ERR("invalid nd_")
       ierr = 1
    end if        
  end subroutine dot_RP
  subroutine HIJ(calc_H_X, K, res, ierr)
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in) :: K
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer I
    integer ierr
    ierr = 0
    do I = 1, ne_
       res(:,I) = 0
       res(I,I) = 1
       call Hc(calc_H_X, K, res(:,I), ierr)
    end do
    
  end subroutine HIJ
  subroutine Hc(calc_H_X, K, c, ierr)
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    complex(kind(0d0)), intent(inout) :: c(:)
    integer, intent(in) :: K
    complex(kind(0d0)) :: HeIJ(ne_, ne_), XkIJ(nf_, ne_, ne_)
    integer ierr, i
    complex(kind(0d0)) :: c0(ne_), Xc0(ne_)
    ierr = 0
    call calc_H_X(R_(K,:), HeIJ(:,:), XkIJ(:,:,:), ierr)
    c0(:) = c(:)

    c(:) = matmul(HeIJ(:,:), c0(:))
    do i = 1, nf_
       Xc0(:) = matmul(XkIJ(i,:,:), c0(:))
       c(:) = c(:) + P_(K,i)*P_(K,i)/(2*m_) * c0(:)
       c(:) = c(:) - II * P_(K,i) * Xc0(:)
       c(:) = c(:) - 0.5d0 * matmul(XkIJ(i,:,:), Xc0(:))
    end do
        
  end subroutine Hc
end module Mod_DyBranch
