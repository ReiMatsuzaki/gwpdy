#include "macros.fpp"

module Mod_DyMono
  ! molecular wave function \Psi are expanded as 
  !     \Psi(r,Q,t) = G(Q,t) \sum_I c_I(t) \Phi_I(r;Q) .
  ! G is auss wave packet
  use Mod_GWP
  implicit none
  integer :: nf_   ! number of freedom
  integer :: ne_   ! number of electron Hilbert space
  double precision :: m_ ! mass
  type(Obj_GWP) :: gwp_
  complex(kind(0d0)), allocatable :: c_(:)
  integer :: nd_        ! order of numerical differate
  character(8) :: inte_RP_    ! euler or RK4
  double precision :: dt_, dR_, dP_ ! finite difference for t R and P
  integer :: nt_, n1t_  ! number of time for print, number for each step
contains
  ! ==== Driver =====
  subroutine write_input(ierr)
    use Mod_sys, only : open_w, mkdirp_if_not
    integer, parameter :: ifile = 921
    integer ierr
    call mkdirp_if_not("out")
    call open_w(ifile, "out/common.json", ierr); check_err(ierr)
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
    write(*,*) "gwp%g:", gwp_%g(1,1,1)
    write(*,*) "gwp%R:", gwp_%R(1,:)
    write(*,*) "gwp%P:", gwp_%P(1,:)
    write(*,*) "gwp%c:", gwp_%c(1)
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
    
    write(fn, '("out/", I0, "/", A)') it, "gwp_g.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    write(ifile, '("val")')
    do k = 1, nf_
       write(ifile,*) gwp_%g(1,1,1)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "gwp_r.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    write(ifile, '("val")')
    do k = 1, nf_
       write(ifile,*) gwp_%R(1,k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "gwp_p.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    write(ifile, '("val")')
    do k = 1, nf_
       write(ifile,*) gwp_%P(1,k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "gwp_c.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    write(ifile, '("re,im")')
    do i = 1, ne_
       write(ifile,'(E20.10,",",E20.10)') real(c_(i)), aimag(c_(i))
    end do
    close(ifile)
    ifile = ifile + 1            
    
  end subroutine write_res
  subroutine print_res(it)
    integer, intent(in) :: it
    integer k, i
    write(*,'("t: ", F20.10)') (it-1)*n1t_*dt_
    do k = 1, nf_
       write(*,'("R",I0,": ", F20.10)') k, gwp_%R(1,k)
       write(*,'("P",I0,": ", F20.10)') k, gwp_%P(1,k)
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
    double precision re, im
    integer i, ifile, k
    integer ierr
    ifile = 23422
    write(fn, '("out/", I0)') it
    call mkdirp_if_not(fn)
    
    write(fn, '("out/", I0, "/", A)') it, "gwp_g.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    read(ifile)
    gwp_%g(:,:,:) = 0
    do k = 1, nf_
       read(ifile,*) gwp_%g(1,k,k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "gwp_r.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    read(ifile)
    do k = 1, nf_
       read(ifile, *) gwp_%R(1,k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "gwp_p.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
    read(ifile)
    do k = 1, nf_
       write(ifile,*) gwp_%P(1,k)
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '("out/", I0, "/", A)') it, "gwp_c.csv"
    call open_w(ifile, fn, ierr); check_err(ierr)
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

    call write_input(ierr); check_err(ierr)
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
          call DyMono_update(calc_H_X, ierr); check_err(ierr)
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
    
    call GWP_new(gwp_, nf, 1, 'c', ierr) ; check_err(ierr)
    allocate(c_(ne))
    c_ = 0
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
    call GWP_setup(gwp_, ierr); check_err(ierr)
    norm2 = sum(abs(c_(:))**2)
    if(norm2<tol) then
       begin_err("norm is too small")
       ierr = 1
       write(0,*) "norm2:", norm2
       return
    end if
    c_(:) = c_(:) / sqrt(norm2)

    if(inte_RP_.ne."euler" .and. inte_RP_.ne."RK4") then
       begin_err("unsupported inte_RP_")
       ierr = 1
       write(0,*) "inte_RP_:", inte_RP_
       return
    end if
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
    
    gwp_%R(1,:) = gwp_%R(1,:) + dR(:)
    gwp_%P(1,:) = gwp_%P(1,:) + dP(:)

  end subroutine DyMono_update
  subroutine DyMono_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    call GWP_delete(gwp_, ierr); check_err(ierr)
    deallocate(c_)
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
    dR(:) = dotR(:)*dt_
    dP(:) = dotP(:)*dt_
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
    
    gwp_%R(1,:) = gwp_%R(1,:) + kR(1,:) * dt_/2
    gwp_%P(1,:) = gwp_%P(1,:) + kP(1,:) * dt_/2
    call dot_RP(calc_H_X, kR(2,:), kP(2,:), ierr)
    gwp_%R(1,:) = gwp_%R(1,:) - kR(1,:) * dt_/2
    gwp_%P(1,:) = gwp_%P(1,:) - kP(1,:) * dt_/2

    gwp_%R(1,:) = gwp_%R(1,:) + kR(2,:) * dt_/2
    gwp_%P(1,:) = gwp_%P(1,:) + kP(2,:) * dt_/2
    call dot_RP(calc_H_X, kR(3,:), kP(3,:), ierr)
    gwp_%R(1,:) = gwp_%R(1,:) - kR(2,:) * dt_/2
    gwp_%P(1,:) = gwp_%P(1,:) - kP(2,:) * dt_/2

    gwp_%R(1,:) = gwp_%R(1,:) + kR(3,:) * dt_/2
    gwp_%P(1,:) = gwp_%P(1,:) + kP(3,:) * dt_/2
    call dot_RP(calc_H_X, kR(3,:), kP(4,:), ierr)
    gwp_%R(1,:) = gwp_%R(1,:) - kR(3,:) * dt_/2
    gwp_%P(1,:) = gwp_%P(1,:) - kP(3,:) * dt_/2    

    dR(:) =  (kR(1,:) + 2*kR(2,:) + 2*kR(3,:) + kR(4,:)) * dt_/6
    dP(:) =  (kP(1,:) + 2*kP(2,:) + 2*kP(3,:) + kP(4,:)) * dt_/6
    
    
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

    call HIJ(calc_H_X, H1(:,:), ierr); check_err(ierr)
    call lapack_zheev(ne_, H1(:,:), w(:), U(:,:), ierr); check_err(ierr)
    UH(:,:) = conjg(transpose(U(:,:)))
    c_(:) = matmul(UH(:,:), c_(:))
    c_(:) = exp(-II*w(:)*dt_) * c_(:)
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
          gwp_%R(1,k) = gwp_%R(1,k) + dR_
          call HIJ(calc_H_X, H1(:,:), ierr)
          gwp_%R(1,k) = gwp_%R(1,k) - 2*dR_
          call HIJ(calc_H_X, H2(:,:), ierr)
          gwp_%R(1,k) = gwp_%R(1,k) + dR_
          dH_dR(:,:) = (H1(:,:) - H2(:,:)) / (2*dR_)
          dotP(k) = -real(dot_product(c_(:), matmul(dH_dR(:,:), c_(:))))

          gwp_%P(1,k) = gwp_%P(1,k) + dP_
          call HIJ(calc_H_X, H1(:,:), ierr)
          gwp_%P(1,k) = gwp_%P(1,k) -2*dP_
          call HIJ(calc_H_X, H2(:,:), ierr)
          gwp_%P(1,k) = gwp_%P(1,k) + dP_
          dH_dP(:,:) = (H1(:,:) - H2(:,:)) / (2*dP_)
          dotR(k) = real(dot_product(c_(:), matmul(dH_dP(:,:), c_(:))))
       end do
    else
       begin_err("invalid nd_")
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
       res(:,I) = 1
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
    call calc_H_X(gwp_%R(1,:), HeIJ(:,:), XkIJ(:,:,:), ierr)
    c0(:) = c(:)

    c(:) = matmul(HeIJ(:,:), c0(:))
    do k = 1, nf_
       Xc0(:) = matmul(XkIJ(k,:,:), c0(:))
       c(:) = c(:) + gwp_%P(1,k)**2 * 0.5d0 * c0(:)
       c(:) = c(:) - II * gwp_%P(1,k) * Xc0(:)
       c(:) = c(:) - 0.5d0 * matmul(XkIJ(k,:,:), Xc0(:))
    end do
        
  end subroutine Hc
end module Mod_DyMono
