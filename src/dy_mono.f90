#include "macros.fpp"



module Mod_DyMono
  ! molecular wave function \Psi are expanded as 
  !     \Psi(r,Q,t) = G(Q,t) \sum_I c_I(t) \Phi_I(r;Q) .
  ! G is auss wave packet
  implicit none
  integer :: nf_         ! number of freedom
  integer :: ne_         ! number of electron Hilbert space
  double precision :: m_ ! mass
  double precision :: g_ ! exponent of GWP
  double precision, allocatable :: R_(:), P_(:) ! phase space location of GWP
  complex(kind(0d0)), allocatable :: c_(:)      ! expansion coefficient of e-state
  integer :: nd_                    ! order of numerical differate
  character(8) :: inte_RP_          ! euler or RK4
  double precision :: dydt_, dR_, dP_ ! finite difference for t, R and P
  double precision :: dt_
  integer :: nt_, n1t_              ! number of time for print, number for each step
  private :: dydt_
contains
  ! ==== Driver =====
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
    write(*,*) "g:", g_
    write(*,*) "R:", R_(:)
    write(*,*) "P:", P_(:)
    write(*,*) "c:", c_(:)
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
  subroutine write_res(it)
    use Mod_sys, only : open_w, mkdirp_if_not
    integer it, k, i
    character(50) fn
    integer ierr, ifile
    ifile = 1242

    write(fn, '("out/", I0)') it
    call mkdirp_if_not(fn)
    
    write(fn, '("out/", I0, "/", A)') it, "g.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("val")')
    do k = 1, nf_
       write(ifile,*) g_
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "r.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("val")')
    do k = 1, nf_
       write(ifile,*) R_(k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "p.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("val")')
    do k = 1, nf_
       write(ifile,*) P_(k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "c.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("re,im")')
    do i = 1, ne_
       write(ifile,*) real(c_(i)), ",", aimag(c_(i))
    end do
    close(ifile)
    ifile = ifile + 1    
    
  end subroutine write_res
  subroutine print_res(it)
    integer, intent(in) :: it
    integer k, i
    write(*,'("t: ", F20.10)') it*n1t_*dt_
    do k = 1, nf_
       write(*,'("R",I0,": ", F20.10)') k, R_(k)
       write(*,'("P",I0,": ", F20.10)') k, P_(k)
    end do
    do i = 1, ne_
       write(*,'("re_c",I0,": ", F20.10)') i, real( c_(i))
       write(*,'("im_c",I0,": ", F20.10)') i, aimag(c_(i))
    end do
    write(*,*) "---------------------"
  end subroutine print_res
  function dir_exists(it) result(res)
    integer, intent(in) :: it
    logical res
    character(100) :: fn

    write(fn, '("out/", I0)') it
    
    if(0.eq.access(fn, " ")) then
       res = .false.
    else
       res = .true.
    end if
        
  end function dir_exists
  subroutine read_res(it)
    use Mod_sys, only : mkdirp_if_not, open_w
    integer, intent(in) :: it
    character(100) :: fn
    integer ifile, k, i
    double precision re, im

    integer ierr
    ifile = 23422
    write(fn, '("out/", I0)') it
    call mkdirp_if_not(fn)
    
    write(fn, '("out/", I0, "/", A)') it, "g.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    read(ifile)
    read(ifile,*) g_
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "r.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    read(ifile)
    do k = 1, nf_
       read(ifile, *) R_(k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "p.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    read(ifile)
    do k = 1, nf_
       write(ifile,*) P_(k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "c.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    read(ifile)
    do i = 1, ne_
       read(ifile, *) re, im
       c_(i) = dcmplx(re, im)
    end do
    close(ifile)
    ifile = ifile + 1
    
  end subroutine read_res
  subroutine DyMono_run(calc_H_X)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer it, i1t, ierr, it0

    write(*,*) ""
    write(*,*) "  GWPDY_Mono BEGIN"
    write(*,*) ""
    write(*,*) "==== Inputs ===="

    call write_input(ierr); CHK_ERR(ierr)
    call print_input

    write(*, *) "==== read if prev exists ===="
    it0 = 0
    do it = 0, nt_
       if(dir_exists(it) .and. .not. dir_exists(it+1)) then
          call read_res(it)
          it0 = it
       end if
    end do
    
    write(*,*) "==== Calc ===="
    do it = it0, nt_       
       call write_res(it)
       call print_res(it)
       do i1t = 1, n1t_
          call DyMono_update(calc_H_X, ierr); CHK_ERR(ierr)
       end do
    end do

    write(*,*) ""
    write(*,*) "  GWPDY_Mono END"
    write(*,*) ""    
    
  end subroutine DyMono_run
  ! ==== Main ====
  subroutine DyMono_new(nf, ne, ierr)
    integer, intent(in) :: nf, ne
    integer, intent(out) :: ierr
    ierr = 0
    nf_  = nf
    ne_ = ne

    g_ = 1
    allocate(R_(nf), P_(nf), c_(ne))
    R_(:) = 0
    P_(:) = 0
    c_(:) = 0
    c_(1) = 1
    m_ = 1
    
    nd_ = 2
    inte_RP_ = "euler"
    dt_ = 0.1d0
    nt_  = 10
    n1t_ = 1
    dR_ = 0.01d0
    dP_ = 0.01d0
    
  end subroutine DyMono_new
  subroutine DyMono_setup(ierr)
    integer, intent(out) :: ierr
    double precision :: norm2
    double precision, parameter :: tol = 1.0d-15
    ierr = 0
    norm2 = sum(abs(c_(:))**2)
    if(norm2<tol) then
       MSG_ERR("norm is too small")
       ierr = 1
       write(0,*) "norm2:", norm2
       return
    end if
    c_(:) = c_(:) / sqrt(norm2)

    if(inte_RP_.ne."euler" .and. inte_RP_.ne."RK4") then
       MSG_ERR("unsupported inte_RP_")
       ierr = 1
       write(0,*) "inte_RP_:", inte_RP_
       return
    end if

    dydt_ = dt_/n1t_
  end subroutine DyMono_setup
  subroutine DyMono_update(calc_H_X, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    double precision :: dR(nf_), dP(nf_)

    ierr = 0

    if(inte_RP_.eq."euler") then
       call inte_nuc_euler(calc_H_X, dR(:), dP(:), ierr)
    else if(inte_RP_.eq."RK4") then
       call inte_nuc_RK4(calc_H_X, dR(:), dP(:), ierr)
    end if

    call inte_ele_diag(calc_H_X, ierr)
    
    R_(:) = R_(:) + dR(:)
    P_(:) = P_(:) + dP(:)

  end subroutine DyMono_update
  subroutine DyMono_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(R_,P_,c_)
  end subroutine DyMono_delete
  ! ==== Private ====
  subroutine inte_nuc_euler(calc_H_X, dR, dP, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    double precision, intent(out) :: dR(:), dP(:)
    double precision :: dotR(nf_), dotP(nf_)
    integer, intent(out) :: ierr
    call dot_RP(calc_H_X, dotR, dotP, ierr)
    dR(:) = dotR(:)*dydt_
    dP(:) = dotP(:)*dydt_
  end subroutine inte_nuc_euler
  subroutine inte_nuc_RK4(calc_H_X, dR, dP, ierr)
    interface
        subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    double precision, intent(out) :: dR(:), dP(:)
    integer, intent(out) :: ierr
    double precision :: kR(4,nf_), kP(4,nf_)
    
    call dot_RP(calc_H_X, kR(1,:), kP(1,:), ierr)
    
    R_(:) = R_(:) + kR(1,:) * dydt_/2
    P_(:) = P_(:) + kP(1,:) * dydt_/2
    call dot_RP(calc_H_X, kR(2,:), kP(2,:), ierr)
    R_(:) = R_(:) - kR(1,:) * dydt_/2
    P_(:) = P_(:) - kP(1,:) * dydt_/2

    R_(:) = R_(:) + kR(2,:) * dydt_/2
    P_(:) = P_(:) + kP(2,:) * dydt_/2
    call dot_RP(calc_H_X, kR(3,:), kP(3,:), ierr)
    R_(:) = R_(:) - kR(2,:) * dydt_/2
    P_(:) = P_(:) - kP(2,:) * dydt_/2

    R_(:) = R_(:) + kR(3,:) * dydt_
    P_(:) = P_(:) + kP(3,:) * dydt_
    call dot_RP(calc_H_X, kR(4,:), kP(4,:), ierr)
    R_(:) = R_(:) - kR(3,:) * dydt_
    P_(:) = P_(:) - kP(3,:) * dydt_

    dR(:) =  (kR(1,:) + 2*kR(2,:) + 2*kR(3,:) + kR(4,:)) * dydt_/6
    dP(:) =  (kP(1,:) + 2*kP(2,:) + 2*kP(3,:) + kP(4,:)) * dydt_/6    
    
  end subroutine inte_nuc_RK4
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
    complex(kind(0d0)) :: U(ne_, ne_), UH(ne_, ne_), H1(ne_,ne_)
    double precision :: w(ne_)
    ierr = 0

    call HIJ(calc_H_X, H1(:,:), ierr); CHK_ERR(ierr)
    call lapack_zheev(ne_, H1(:,:), w(:), U(:,:), ierr); CHK_ERR(ierr)
    UH(:,:) = conjg(transpose(U(:,:)))
    c_(:) = matmul(UH(:,:), c_(:))
    c_(:) = exp(-II*w(:)*dydt_) * c_(:)
    c_(:) = matmul(U(:,:), c_(:))
    
  end subroutine inte_ele_diag
  subroutine dot_RP(calc_H_X, dotR, dotP, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    double precision, intent(out) :: dotR(:), dotP(:)
    integer, intent(out) :: ierr
    integer k
    complex(kind(0d0)) :: H1(ne_,ne_), H2(ne_,ne_)
    complex(kind(0d0)) :: dH_dR(ne_,ne_), dH_dP(ne_,ne_)
    if(nd_.eq.2) then
       do k = 1, nf_
          R_(k) = R_(k) + dR_
          call HIJ(calc_H_X, H1(:,:), ierr)
          R_(k) = R_(k) - 2*dR_
          call HIJ(calc_H_X, H2(:,:), ierr)
          R_(k) = R_(k) + dR_
          dH_dR(:,:) = (H1(:,:) - H2(:,:)) / (2*dR_)
          dotP(k) = -real(dot_product(c_(:), matmul(dH_dR(:,:), c_(:))))

          P_(k) = P_(k) + dP_
          call HIJ(calc_H_X, H1(:,:), ierr)
          P_(k) = P_(k) -2*dP_
          call HIJ(calc_H_X, H2(:,:), ierr)
          P_(k) = P_(k) + dP_
          dH_dP(:,:) = (H1(:,:) - H2(:,:)) / (2*dP_)
          dotR(k) = real(dot_product(c_(:), matmul(dH_dP(:,:), c_(:))))
       end do
    else
       MSG_ERR("invalid nd_")
       ierr = 1
    end if        
  end subroutine dot_RP
  subroutine HIJ(calc_H_X, res, ierr)
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    complex(kind(0d0)), intent(out) :: res(:,:)
    !    complex(kind(0d0)) :: HeIJ0(ne_, ne_), XkIJ0(nf_, ne_, ne_)
    !   integer k, I
    integer I
    integer ierr
    ierr = 0
    do I = 1, ne_
       res(:,I) = 0
       res(I,I) = 1
       call Hc(calc_H_X, res(:,I), ierr)
    end do
    
!    call HeIJ(gwp_%R(1,:), HeIJ0(:,:),   ierr)
!    call XkIJ(gwp_%R(1,:), XkIJ0(:,:,:), ierr)
!
!    res(:,:) = HeIJ0(:,:)
!    do k = 1, nf_
!       res(:,:) = res(:,:) + gwp_%P(1,k)**2 * 0.5d0
!       res(:,:) = res(:,:) - II * XkIJ0(k,:,:)*gwp_%P(1,k)
!       res(:,:) = res(:,:) - 0.5d0 * matmul(XkIJ0(k,:,:), XkIJ0(k,:,:))
!    end do
    
  end subroutine HIJ
  subroutine Hc(calc_H_X, c, ierr)
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    complex(kind(0d0)), intent(inout) :: c(:)
    complex(kind(0d0)) :: HeIJ(ne_, ne_), XkIJ(nf_, ne_, ne_)
    integer k
    integer ierr
    complex(kind(0d0)) :: c0(ne_), Xc0(ne_)
    ierr = 0
    call calc_H_X(R_(:), HeIJ(:,:), XkIJ(:,:,:), ierr)
    c0(:) = c(:)

    c(:) = matmul(HeIJ(:,:), c0(:))
    do k = 1, nf_
       Xc0(:) = matmul(XkIJ(k,:,:), c0(:))
       c(:) = c(:) + P_(k)*P_(k)/(2*m_) * c0(:)
       c(:) = c(:) - II * P_(k)/m_ * Xc0(:)
       c(:) = c(:) - 1/(2*m_) * matmul(XkIJ(k,:,:), Xc0(:))
    end do
        
  end subroutine Hc
end module Mod_DyMono
