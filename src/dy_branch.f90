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
  subroutine write_input(ierr)
    use Mod_sys, only : open_w, mkdirp_if_not
    integer, parameter :: ifile = 921
    integer ierr
    call mkdirp_if_not("out")
    call open_w(ifile, "out/common.json", ierr); CHK_ERR(ierr)
    write(ifile, '("{")')
    write(ifile,*) '"nf":', nf_, ","
    write(ifile,*) '"ne":', ne_, ","
    write(ifile,*) '"dt":', dt_, ","
    write(ifile,*) '"nt":', nt_, ","
    write(ifile,*) '"n1t":', n1t_,","
    write(ifile,*) '"dR":', dR_, ","
    write(ifile,*) '"dP":', dP_
    write(ifile, '("}")')
    close(ifile)
  end subroutine write_input
  subroutine print_input
    write(*,*) "-- basics --"
    write(*,*) "nf:", nf_
    write(*,*) "ne:", ne_
    write(*,*) "-- init conds --"
    write(*,*) "g:", g_(1:npath_)
    write(*,*) "R:", R_(1:npath_,:)
    write(*,*) "P:", P_(1:npath_,:)
    write(*,*) "C:", cc_(1:npath_)
    write(*,*) "c:", c_(1:npath_,:)
    write(*,*) "-- time --"
    write(*,*) "dt:", dt_
    write(*,*) "nt:", nt_
    write(*,*) "n1t:", n1t_
    write(*,*) "-- calc method --"
    write(*,*) "nd:", nd_
    write(*,*) "dR:", dR_
    write(*,*) "dP:", dP_
    write(*,*)    
  end subroutine print_input
  subroutine write_status(it)
    use Mod_sys, only : open_w, mkdirp_if_not
    integer it, k, KK, I
    character(50) fn
    integer ierr, ifile
    ifile = 1242

    write(fn, '("out/", I0)') it
    call mkdirp_if_not(fn)
    
    write(fn, '("out/", I0, "/", A)') it, "g.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,val")')
    do k = 1, nf_
       write(ifile,'(I0,",",F20.10)') k,g_(k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "r.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,val")')
    do KK = 1, npath_
       do k = 1, nf_
          write(ifile,'(I0,",",I0,",",F20.10)') KK,k,R_(KK,k)
       end do
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "p.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,val")')
    do KK = 1, npath_
       do k = 1, nf_
          write(ifile,'(I0,",",I0,",",F20.10)') KK,k,P_(KK,k)
       end do
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "c.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,re,im")')
    do KK = 1, npath_
       do I = 1, ne_
          write(ifile, '(I0,",",I0,",",F20.10,",",F20.10)') KK, I, real(c_(KK,I)), aimag(c_(KK,I))
       end do
    end do

    write(fn, '("out/", I0, "/", A)') it, "cc.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,re,im")')
    do KK = 1, npath_
       write(ifile,'(I0,",",F20.10,",",F20.10)') KK, real(cc_(KK)),  aimag(cc_(KK))
    end do
    close(ifile)
    ifile = ifile + 1    
    
  end subroutine write_status
  subroutine print_status(it)
    integer, intent(in) :: it
    integer KK, k, I
    write(*,'("t: ", F20.10)') it*n1t_*dt_
    do KK = 1, npath_
       write(*,'("C(",I0,"):", F20.10, F20.10)') KK, real(CC_(KK)), aimag(CC_(KK))
       do k = 1, nf_
          write(*,'("R(",I0,",",I0,"): ", F20.10)') KK, k, R_(KK, k)
          write(*,'("P(",I0,",",I0,"): ", F20.10)') KK, k, P_(KK, k)
       end do
       do I = 1, ne_
          write(*,'("c(",I0,",",I0,"):", F20.10, F20.10)') KK, I, c_(KK,I)
       end do
    end do
    write(*,*) "---------------------"
  end subroutine print_status
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
    ! update with diagonal method for coefficient and Eulaer for (R,P).
    ! H_{KL} = <G_K c_{I,K} Phi_I | H | G_L c_{J,L} Phi_J >
    !        = c_{I,K}c_{I,L} <G_K|T|G_L> +
    !          c_{I,K}c_{J,L} (P_K+P_L)/2 (X_{IJ}(Q_K) + X_{IJ}(Q_L))/2 (K,L)
    use Mod_GWP
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    double precision :: dotR(npath_,nf_), dotP(npath_,nf_)
    integer KK, LL, k
    complex(kind(0d0)) :: HK(npath_,ne_,ne_), HH(npath_,npath_)
    complex(kind(0d0)) :: c_KL(npath_,napth_), tmp
    complex(kind(0d0)) :: P2(nf_, npath_,npath_), S(path_,npath_)
    type(Obj_GWP) :: gwp

    ierr = 0

    do KK = 1, npath_
       call dot_RP(calc_H_X, KK, dotR(KK,:), dotP(KK,:), ierr); CHK_ERR(ierr)
       call HIJ(calc_H_X, KK, HK(KK,:,:), ierr); CHK_ERR(ierr)
    end do
    call make_GWP(gwp, ierr); CHK_ERR(ierr)
    call gwp_overlap(gwp, S(:,:), ierr); CHK_ERR(ierr)
    call gwp_p2(gwp, P2(:,:,:), ierr); CHK_ERR(ierr)
    call GWP_delete(gwp, ierr); CHK_ERR(ierr)
    c_KL(:,:) = 0
    do KK = 1, npath_
       do LL = 1, npath_
          c_KL(KK,LL) = dot_product(c_(:,KK), c_(:,LL))
       end do
    end do    
    HH(:,:) = 0
    do KK = 1, npath_
       do LL = 1, npath_
          tmp = 0
          do k = 1, nf_
             tmp = tmp + c_KL(KK,LL)*P2(k, KK, LL)/(2*m_)
             tmp = tmp - ii* (P_(KK,k)+P_(LL,k))/2 * (X_())
          end do
          HH(KK,LL) = tmp
       end do
    end do
    
    call inte_ele_diag(calc_H_X, ierr); CHK_ERR(ierr)

    R_(1:npath_,:) = R_(1:npath_,:) + dR(:,:)
    P_(1:npath_,:) = P_(1:npath_,:) + dP(:,:)
    
  end subroutine DyBranch_update
  subroutine DyBranch_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(g_, R_, P_, theta_)
    deallocate(CC_, c_)
  end subroutine DyBranch_delete
  ! ==== 1st order ====  
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
    write(*,*) "in RK4"
    do KK = 1, npath_
       call dot_RP(calc_H_X, KK, kR(1,:), kP(1,:), ierr); CHK_ERR(ierr)
    
       R_(KK,:) = R_(KK,:) + kR(1,:) * dydt_/2
       P_(KK,:) = P_(KK,:) + kP(1,:) * dydt_/2
       call dot_RP(calc_H_X, KK, kR(2,:), kP(2,:), ierr); CHK_ERR(ierr)
       R_(KK,:) = R_(KK,:) - kR(1,:) * dydt_/2
       P_(KK,:) = P_(KK,:) - kP(1,:) * dydt_/2
       
       R_(KK,:) = R_(KK,:) + kR(2,:) * dydt_/2
       P_(KK,:) = P_(KK,:) + kP(2,:) * dydt_/2
       call dot_RP(calc_H_X, KK, kR(3,:), kP(3,:), ierr); CHK_ERR(ierr)
       R_(KK,:) = R_(KK,:) - kR(2,:) * dydt_/2
       P_(KK,:) = P_(KK,:) - kP(2,:) * dydt_/2
       
       R_(KK,:) = R_(KK,:) + kR(3,:) * dydt_
       P_(KK,:) = P_(KK,:) + kP(3,:) * dydt_
       call dot_RP(calc_H_X, KK, kR(4,:), kP(4,:), ierr); CHK_ERR(ierr)
       R_(KK,:) = R_(KK,:) - kR(3,:) * dydt_
       P_(KK,:) = P_(KK,:) - kP(3,:) * dydt_
       
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
    integer K
    double precision ::  w(ne_)
    complex(kind(0d0)) :: H(ne_, ne_), U(ne_, ne_)
    
    ierr = 0
    K = 1
    call HIJ(calc_H_X, K, H, ierr); CHK_ERR(ierr)
    call lapack_zheev(ne_, H(:,:), w(:), U(:,:), ierr); CHK_ERR(ierr)
    
    npath_ = ne_

    do K = 2, ne_
       R_(K,:) = R_(1,:)
       P_(K,:) = P_(1,:)       
    end do

    CC_(1:npath_)  = matmul( conjg(transpose(U(:,:))), c_(1,:))
    c_(1:npath_,:) = transpose(U(:,:))
    
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
    integer KK
    ierr = 0

    do KK = 1, npath_
       call HIJ(calc_H_X, KK, H(:,:), ierr); CHK_ERR(ierr)
       call lapack_zheev(ne_, H(:,:), w(:), U(:,:), ierr); CHK_ERR(ierr)
       UH(:,:) = conjg(transpose(U(:,:)))
       c_(KK,:) = matmul(UH(:,:), c_(KK,:))
       c_(KK,:) = exp(-II*w(:)*dydt_) * c_(KK,:)
       c_(KK,:) = matmul(U(:,:), c_(KK,:))
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
          call HIJ(calc_H_X, K, H1(:,:), ierr); CHK_ERR(ierr)
          R_(K,i) = R_(K,i) - 2*dR_
          call HIJ(calc_H_X, K, H2(:,:), ierr); CHK_ERR(ierr)
          R_(K,i) = R_(K,i) + dR_
          dH_dR(:,:) = (H1(:,:) - H2(:,:)) / (2*dR_); CHK_ERR(ierr)
          dotP(i) = -real(dot_product(c_(K,:), matmul(dH_dR(:,:), c_(K,:))))

          P_(K,i) = P_(K,i) + dP_
          call HIJ(calc_H_X, K, H1(:,:), ierr); CHK_ERR(ierr)
          P_(K,i) = P_(K,i) -2*dP_
          call HIJ(calc_H_X, K, H2(:,:), ierr); CHK_ERR(ierr)
          P_(K,i) = P_(K,i) + dP_
          dH_dP(:,:) = (H1(:,:) - H2(:,:)) / (2*dP_)
          dotR(i) = real(dot_product(c_(K,:), matmul(dH_dP(:,:), c_(K,:))))
       end do
    else
       MSG_ERR("invalid nd_")
       ierr = 1
    end if        
  end subroutine dot_RP
  subroutine make_gwp(gwp, ierr)
    use Mod_GWP
    type(Obj_GWP) :: gwp
    integer, intent(out) :: ierr
    integer K, KK
    call GWP_new(gwp, nf_, npath_, "c", ierr); CHK_ERR(ierr)
    do KK = 1, npath_
       do k = 1, nf_
          gwp%g(KK,k,k) = g_(KK)
          gwp%R(KK,k)   = R_(KK,k)
          gwp%P(KK,k)   = P_(KK,k)
          gwp%c(KK)     = theta_(KK)
       end do
    end do
    call GWP_setup(gwp, ierr); CHK_ERR(ierr)
    
  end subroutine make_gwp
  subroutine HIJ(calc_H_X, KK, res, ierr)
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in) :: KK
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer,intent(out) ::  ierr
    !    type(Obj_GWP) :: gwp
    complex(kind(0d0)) :: H(ne_,ne_), X(nf_,ne_,ne_)
    !    complex(kind(0d0)) :: s(npath_,npath_), p2(nf_,npath_,npath_)
    integer k
    

    call calc_H_X(R_(KK,:), H(:,:), X(:,:,:), ierr); CHK_ERR(ierr)
    !    call make_gwp(gwp, ierr); CHK_ERR(ierr)
    !    call gwp_overlap(gwp, s, ierr); CHK_ERR(ierr)
    !    call gwp_p2(gwp, p2, ierr); CHK_ERR(ierr)
    
    ierr = 0
    res(:,:) = H(:,:)
    do k = 1, nf_
       res(:,:) = res(:,:) + 1/(2*m_)*P_(KK,k)**2
       res(:,:) = res(:,:) - II * X(k,:,:)*P_(KK,k)/m_
       res(:,:) = res(:,:) - 1/(2*m_) * matmul(X(k,:,:), X(k,:,:))
    end do
    
  end subroutine HIJ
  subroutine global_HIJ(calc_H_X, res, ierr)
    use Mod_GWP
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    type(Obj_GWP) :: gwp
    complex(kind(0d0)) :: H(ne_,ne_), X(nf_,ne_,ne_)
    complex(kind(0d0)) :: HH(npath_), XX(nf_,npath_)
    complex(kind(0d0)) :: s(npath_,npath_), p2(nf_,npath_,npath_)
    integer :: KK, LL, k
    complex(kind(0d0)) :: tmp
    ierr = 0
    
    call make_gwp(gwp, ierr); CHK_ERR(ierr)
    call gwp_overlap(gwp, s, ierr); CHK_ERR(ierr)
    call gwp_p2(gwp, p2, ierr); CHK_ERR(ierr)

    do KK = 1, npath_
       call calc_H_X(R_(KK,:), H(:,:), X(:,:,:), ierr); CHK_ERR(ierr)
       HH(KK) = dot_product(c_(KK,:), matmul(H(:,:), c_(KK,:)))
       do k = 1, nf_
          XX(k,KK) = dot_product(c_(KK,:), matmul(X(k,:,:), c_(KK,:)))
       end do
    end do

    res(:,:) = 0
    do KK = 1, npath_
       do LL = 1, npath_
          tmp = H(KK,LL)
          do k = 1, nf_
             tmp = tmp + p2(k,KK,LL)/(2*m_)
             tmp = tmp + 1/m_*(P_(KK,k)+P_(LL,k))/2*(XX(k,KK)+XX(k,LL))/2 * s(KK,LL)
          end do
          res(KK,LL) = res(KK,LL) + tmp
       end do
    end do
  end subroutine global_HIJ
end module Mod_DyBranch

module Mod_Updater1st
  implicit none
  use Mod_DyBranch
  use Mod_GWP
  complex(kind(0d0)), allocatable  :: HeIJ_(:,:,:), XkIJ_(:,:,:,:)
  complex(kind(0d0)), allocatable  :: s_(:,:), p2_(:,:,:)
  type(Obj_GWP) :: gwp_
contains
  subroutine updater1st_new(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    call GWP_new(gwp_, nf_, npath_, "d", ierr); CHK_ERR(ierr)    
  end subroutine updater1st_new
  subroutine updaate_1st(calc_H_X, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    double precision :: dotR(npath_,nf_), dotP(npath_,nf_)
    integer KK, LL, k
    complex(kind(0d0)) :: HK(npath_,ne_,ne_), HH(npath_,npath_)
    
    ierr  = 0
    do KK = 1, npath_
       call dot_RP(calc_H_X, KK, dotR, dotP, ierr); CHK_ERR(ierr)
    end do

    call update_emat(calc_H_X, ierr); CHK_ERR(ierr)
    call update_nmat(ierr); CHK_ERR(ierr)

    call calc_HK(HK, ierr); CHK_ERR(ierr)
    call calc_global_H(dotR, HH, ierr); CHK_ERR(ierr)
    
  end subroutine updaate_1st
  subroutine update_emat(calc_H_X, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    integer KK
    do KK = 1, npath_
       call calc_H_X(R_(KK,:), HeIJ, XkIJ, ierr); CHK_ERR(ierr)
    end do
  end subroutine update_emat
  subroutine update_nmat(ierr)
    integer, intent(out) :: ierr
    integer KK, k
    do KK = 1, npath_
       do k = 1, nf_
          gwp%g(KK,k,k) = g_(KK)
          gwp%R(KK,k)   = R_(KK,k)
          gwp%P(KK,k)   = P_(KK,k)
          gwp%c(KK)     = theta_(KK)
       end do
    end do
    call GWP_setup(gwp, ierr); CHK_ERR(ierr)
    
  end subroutine update_nmat
  subroutine calc_HK(res, ierr)    
    complex(kind(0d0)), intent(out) :: res(:,:,:)
    integer, intent(out) :: ierr
    integer KK, k
    ierr = 0
    
    do KK = 1, npath_
       res(KK,:,:) = HeIJ_(KK,:,:)
       do k = 1, nf_
          res(KK,:,:) = res(KK,:,:) + 1/(2*m_)*P_(KK,k)**2
          res(KK,:,:) = res(KK,:,:) - II * XkIJ_(KK,k,:,:)*P_(KK,k)/m_
          res(KK,:,:) = res(KK,:,:) - 1/(2*m_) * matmul(XkIJ_(KK,k,:,:), XkIJ_(KK,k,:,:))
       end do
    end do

  end subroutine calc_HK
  subroutine calc_global_t(dotR, dotc, res, ierr)
    use Mod_const, only : II
    double precision, intent(in) :: dotR(:,:)
    complex(kind(0d0)), intent(in) :: dotc(:,:)
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    integer KK, LL
    do KK = 1, npath_
       do LL = 1, npath_
          res(KK, LL) = -dot_product(dotR(LL,:), P_(LL,:)) * S_(KK,LL) &
               * dot_product(c_(KK,:), c_(LL,:)) &
               -II* sum(conjg(c(KK,:))*dotc(LL,:)) * S_(KK,LL)
       end do
    end do
  end subroutine calc_global_t
  subroutine calc_global_H(res, ierr)
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    integer KK, LL, I, J
    complex(kind(0d0)) :: H(npath_,npath_,ne_,ne_), tmp, x, dr
    ierr = 0
    do KK = 1, npath_
       do LL = 1, npath_
          H(KK,LL,:,:) = (HeIJ(KK,:,:)+HeIJ(LL,:,:))/2
          do k = 1, nf_
             do I = 1, ne_             
                H(KK,LL,I,I) = H(KK,LL,I,I) + p2(k,KK,LL)/(2*m_)
             end do
          end do             
          dr = (P(KK,k)+P(LL,k))/(2*m_)
          x  = 
          H(KK,LL,I,J) = H(KK,LL,I,J) &
               -II*dr*(XkIJ(KK,:,:)+XkIJ(LL,:,:))/2*s(KK,LL)               
       end do
    end do

    do KK = 1, nf_
       do LL = 1, nf_
          res(KK,LL) = dot_product(c_(KK,:), matmul(H(KK,LL,:,:)*c_(LL,:))) 
       end do
    end do
    
  end subroutine calc_global_H
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
          call HIJ(calc_H_X, K, H1(:,:), ierr); CHK_ERR(ierr)
          R_(K,i) = R_(K,i) - 2*dR_
          call HIJ(calc_H_X, K, H2(:,:), ierr); CHK_ERR(ierr)
          R_(K,i) = R_(K,i) + dR_
          dH_dR(:,:) = (H1(:,:) - H2(:,:)) / (2*dR_); CHK_ERR(ierr)
          dotP(i) = -real(dot_product(c_(K,:), matmul(dH_dR(:,:), c_(K,:))))

          P_(K,i) = P_(K,i) + dP_
          call HIJ(calc_H_X, K, H1(:,:), ierr); CHK_ERR(ierr)
          P_(K,i) = P_(K,i) -2*dP_
          call HIJ(calc_H_X, K, H2(:,:), ierr); CHK_ERR(ierr)
          P_(K,i) = P_(K,i) + dP_
          dH_dP(:,:) = (H1(:,:) - H2(:,:)) / (2*dP_)
          dotR(i) = real(dot_product(c_(K,:), matmul(dH_dP(:,:), c_(K,:))))
       end do
    else
       MSG_ERR("invalid nd_")
       ierr = 1
    end if        
  end subroutine dot_RP
end module Mod_Updater1st
