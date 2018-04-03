#include "macros.fpp"
! molecular dynamics with quantum path branching representation.
!    \Psi(r,Q,t)       = \sum_K C_K(t) \Psi^{(K)}(r,Q,t)
!    \Psi^{(K)}(r,Q,t) = \sum_{AI} d^{K}_{AI}(t) \chi_A(Q,x(t))\Phi_I(r;Q)
module Mod_DyQPBranch
  implicit none
  ! -- data size --
  integer nf_, ne_, maxnuc_, maxnp_
  integer, allocatable :: nnuc_(:)       ! nnuc_(K)
  ! -- variable --
  double precision, allocatable :: R_(:,:,:), P_(:,:,:)  ! x(K,A)
  complex(kind(0d0)), allocatable :: g_(:,:)
  complex(kind(0d0)), allocatable :: Ctot_(:), C_(:,:,:) ! C_(K,A,I)
  integer, allocatable :: mode_(:)        ! mode_(K) : 0:NON, 1:SET or 2:PSA
  ! -- time constant --
  double precision :: m_ ! mass  
  double precision :: dt_, dR_, dP_
  integer :: nt_, n1t_
  character(8) :: inte_RP_     ! euler or RK4  
  character(6) :: gauss_mode_ ! frozen, thawed
  ! -- intermediate --  
  double precision :: dydt_
  integer :: np_
contains
  ! ==== Main ====
  subroutine DyQPBranch_new(nf, ne, maxnuc, maxnp, ierr)
    integer, intent(in) :: nf, ne, maxnuc, maxnp
    integer, intent(out) :: ierr

    ierr = 0
    nf_ = nf
    ne_ = ne
    maxnuc_ = maxnuc
    maxnp_ = maxnp
    np_ = 1

    allocate(nnuc_(maxnp))
    allocate(R_(maxnp,maxnuc,nf), P_(maxnp,maxnuc,nf_))
    allocate(Ctot_(maxnp), C_(maxnp, maxnuc, ne))    
    allocate(g_(maxnp,maxnuc))
    allocate(mode_(maxnp))
    nnuc_(:) = 0
    nnuc_(1) = 1
    g_ = 1
    R_ = 0
    P_ = 0
    C_ = 0
    Ctot_ = 0
    C_(1,1,1) = 1
    Ctot_(1) = 1
    mode_(:) = 0
    mode_(1) = 1
    
    m_ = 1
    dt_ = 1
    dR_ = 0.01d0
    dP_ = 0.01d0
    nt_ = 10
    n1t_ = 1

    gauss_mode_ = "frozen"
    inte_RP_ = "euler"
    
  end subroutine DyQPBranch_new
  subroutine DyQPBranch_setup(ierr)
    integer, intent(out) :: ierr
    integer KK
    double precision :: norm2
    ierr = 0
    do KK = 1, np_
       call calc_norm2_path(KK, norm2, ierr); CHK_ERR(ierr)
       if(norm2<1.0d-10) then
          MSG_ERR("norm is too small")
          write(0,*) "norm: ", norm2
          ierr = 1; return
       end if
       C_(KK,:,:) = C_(KK,:,:)/sqrt(norm2)
    end do

    call calc_norm2(norm2, ierr); CHK_ERR(ierr)
    if(norm2<1.0d-10) then
       MSG_ERR("global norm is too small")
       write(*,*) "KK, Ctot, C"
       do KK = 1, np_
          write(*,*) KK, Ctot_(KK), C_(KK,:,:)
       end do
       ierr = 1; return
    end if
    Ctot_(:) = Ctot_(:)/sqrt(norm2)

    if(inte_RP_.ne."euler" .and. inte_RP_.ne."RK4") then
       MSG_ERR("unsupported inte_RP_")
       ierr = 1; return
    end if

    if(gauss_mode_.ne."frozen" .and. gauss_mode_.ne."thawed") then
       MSG_ERR("unsupported gauss_mode")
       ierr = 1; return
    end if

    dydt_ = dt_/n1t_
    
  end subroutine DyQPBranch_setup
  subroutine DyQPBranch_dump(ifile, ierr)
    integer, intent(in) :: ifile
    integer, intent(out) :: ierr
    integer KK
    ierr = 0
    write(ifile,*) "==== DyQPBRanch ===="
    write(ifile,*) "nf:", nf_
    write(ifile,*) "ne:", ne_
    write(ifile,*) "maxnuc:", maxnuc_
    write(ifile,*) "maxnp:", maxnp_
    write(ifile,*) "np:", np_
    write(ifile,*) "inte_RP:", inte_RP_
    write(ifile,*) "gauss_mode_:", gauss_mode_
    write(ifile,*) "paths:"
    do KK = 1, np_
       write(ifile,'(A,I0,A)') "-- path(K = ", KK, ") --"
       write(ifile,*) "Ctot:", Ctot_(KK)
       write(ifile,*) "C:", C_(KK,:,:)
       write(ifile,*) "g:", g_(KK,:nnuc_(KK))
       write(ifile,*) "R:", R_(KK,:nnuc_(KK), :)
       write(ifile,*) "P:", P_(KK,:nnuc_(KK), :) 
    end do
    write(ifile,*) "==================="
  end subroutine DyQPBranch_dump
  subroutine DyQPBranch_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(nnuc_)
    deallocate(R_, P_, Ctot_, C_)
    deallocate(g_, mode_)
  end subroutine DyQPBranch_delete
  ! ==== calc ====
  subroutine update_set(calc_H_X, KK, ierr)
    use Mod_PWGTO2
    use Mod_const, only : II
    use Mod_Math, only : intet_diag
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in) :: KK
    integer, intent(out) :: ierr
    double precision  :: dotR(nf_), dotP(nf_), dotx, doty
    integer :: k
    complex(kind(0d0)) :: HeIJ(np_,ne_,ne_), XkIJ(np_,nf_,ne_,ne_)
    complex(kind(0d0)) :: H0(ne_,ne_), H1(ne_,ne_), H2(ne_,ne_), C(ne_)
    double precision :: R(nf_), P(nf_), dg

    if(nnuc_(KK).ne.1) then
       MSG_ERR("nnuc_(KK)!=1")
       ierr = 1; return
    end if
    
    R(:) = R_(KK,1,:); P(:) = P_(KK,1,:); C(:) = C_(KK,1,:)
    do k = 1, nf_
       R(k) = R(k) + dR_
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H1, ierr)
       R(k) = R(k) -2*dR_
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H2, ierr)
       R(k) = R(k) + dR_
       dotP(:) = -real(dot_product(C, matmul(H1-H2, C)))/(2*dR_)
    end do
        
    call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
    call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P(:), H0(:,:), ierr)
    do k = 1, nf_
       P(k) = P(k) + dP_
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H1, ierr); CHK_ERR(ierr)
       P(k) = P(k) -2*dP_
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H2, ierr); CHK_ERR(ierr)
       P(k) = P(k) + dP_
       dotR(:) = real(dot_product(C, matmul(H1-H2, C)))/(2*dP_)
    end do

    select case(gauss_mode_)
    case("frozen")
       dotx=0; doty=0
    case("thawed")
       dg = 2*real(g_(KK,1))*sum(R)*dR_/sum(R*R)  !dg.(R,R) = g.2(R,dR) => dg = 2g.(R,dR)/(R,R)
       g_(Kk,1) = g_(kK,1) + dg
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H1, ierr)
       g_(Kk,1) = g_(kK,1) - 2*dg
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H2, ierr)
       g_(Kk,1) = g_(kK,1) + dg
       doty = -4*real(g_(KK,1))**2*real(dot_product(C, matmul(H1-H2, C)))/(2*dg)

       g_(Kk,1) = g_(kK,1) + II*dg
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H1, ierr)
       g_(Kk,1) = g_(kK,1) - 2*II*dg
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H2, ierr)
       g_(Kk,1) = g_(kK,1) + II*dg
       dotx = +4*real(g_(KK,1))**2*real(dot_product(C, matmul(H1-H2, C)))/(2*dg)
       
    case default
       MSG_ERR("unsupported gauss_mode_")
       ierr = 1; return
    end select
    
    ! -- update --       
    call intet_diag(ne_, H0(:,:), dydt_, C_(KK,1,:), ierr); CHK_ERR(ierr)
    g_(KK,1)   = g_(KK,1)   + dydt_*(dotx + II*doty)
    R_(KK,1,:) = R_(KK,1,:) + dydt_*dotR(:)
    P_(KK,1,:) = P_(KK,1,:) + dydt_*dotP(:)
    
  end subroutine update_set
  subroutine update_vp(calc_H_X, KK, ierr)
    use Mod_PWGTO2
    use Mod_const, only : II
    use Mod_Math, only : intet_diag
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in) :: KK
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: HeIJ(ne_,ne_), XkIJ(1,ne_,ne_)
    integer :: numvar, i, j, Ie, Je
    integer, allocatable :: ivars(:), iop0
    double precision :: X(:,:), y(:)
    complex(kind(0d0)) :: tmp, S(nnuc_(KK), nnuc_(KK)), Si0(nnuc_(KK), nnuc_(KK)), S0j(nnuc_(KK), nnuc_(KK)), Sij(nnuc_(KK), nnuc_(KK)), P2(nnuc_(KK), nnuc_(KK)), P1(nnuc_(KK), nnuc_(KK)), H00(nnuc_(KK), nnuc_(KK)), Hi0(nnuc_(KK), nnuc_(KK))
    complex(kind(0d0)), allocatable :: v1(:), v2(:)
    type(Obj_PWGTO) :: nbasis
    double precision HeIJs(nnuc_(KK),ne_,ne_), XkIJs(nnuc_(KK),nf_,ne_,ne_)
    
    ierr = 0
    if(gauss_mode_.eq."thawed") then
       MSG_ERR("not impl"); ierr = 1; return
    end if
    if(nf_.ne.1) then
       MSG_ERR("not impl"); ierr = 1; return
    end if

    numvar = 2
    allocate(iops(numvar), X(numvar,numvar), y(numvar), v(nnuc_(KK)))
    iop0 = 3
    iopP2= 4
    iopP1 = 5
    call make_nbasis((/"dR ", "dP ", "0  ", "P2 ", "P1 "/, 1, KK, mbasis, ierr)
    CHK_ERR(ierr)

    do A = 1, nnuc_(KK)
       call calc_H_X(R_(KK,:), HeIJs(KK,:,:), XkIJs(KK,:,:,:), ierr)
       CHK_ERR(ierr)
    end do

    call PWGTO_overlap(nbasis, iop0, iop0,  S00, ierr);CHK_ERR(ierr)
    call PWGTO_overlap(nbasis, iop0, iopP1, P1,  ierr);CHK_ERR(ierr)
    call PWGTO_overlap(nbasis, iop0, iopP2, P2,  ierr);CHK_ERR(ierr)
    
    do i = 1, numvar
       call PWGTO_overlap(nbasis, i, iop0, Si0, ierr); CHK_ERR(ierr)
       do j = 1, numvar                    
          call PWGTO_overlap(nbasis, iop0, j, S0j, ierr); CHK_ERR(ierr)
          call PWGTO_overlap(nbasis, i,    j, Sij, ierr); CHK_ERR(ierr)
          tmp = 0
          do Je = 1, ne_
             v(:) = matmul(S0j, cs_(KK,:,Je))
             call lapack_zgesv_1(nnuc_(KK), S00(:,:), v(:), ierr)
             v(:) = matmul(Si0, v(:))
             tmp = tmp + dot_product(cs_(KK,:,Je), v(:))
             tmp = tmp + dot_product(cs_(KK,:,Je), matmul(Sij(:,:), cs_(KK,:,Je))
          end do
          X(i,j) = -aimag(tmp)
       end do

       do Ie = 1, ne_
          do Je = 1, ne_
             
             do A = 1, nnuc_(KK)
                do B = 1, nnuc_(KK)
                   H00(A,B) = (1/(2*m_))*matmul(P2(A,B) &
                        -II/m_* (XkIJs(,1,Ie,Je))/2 *matmul(P1(:,:), cs_(KK,:,Je))
                   v(:) = cs_(KK,:,Je))
                end do
             end do
          end do
       end do
    end do
    
  end subroutine update_vp
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
    write(ifile,*) '"maxnp":', maxnp_, ","
    write(ifile,*) '"dt":', dt_, ","
    write(ifile,*) '"nt":', nt_, ","
    write(ifile,*) '"n1t":', n1t_,","
    write(ifile,*) '"dR":', dR_, ","
    write(ifile,*) '"dP":', dP_
    write(ifile, '("}")')
    close(ifile)
  end subroutine write_input
  subroutine write_status(it, ierr)
    use Mod_Math, only : dten2csv, cvec2csv, dvec2csv, cten2csv, cmat2csv, ivec2csv
    use Mod_sys, only : open_w, mkdirp_if_not
    integer, intent(in) :: it
    integer, intent(out) :: ierr
    character(50) :: fn, dir_it
    integer :: ifile
    double precision :: probe(ne_)
    ifile = 1242

    write(dir_it,'("out/", I0)') it
    call mkdirp_if_not(dir_it)

    write(fn, '(A, "/", A)') trim(dir_it), "param.json"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, *) "{"
    write(ifile, '(A,I0,A)') '"np": ', np_
    write(ifile, *) "}"
    close(ifile)
    ifile = ifile + 1

    write(fn,'(A, "/", A)') trim(dir_it), "mode.csv"
    call ivec2csv(mode_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "g.csv"
    call cmat2csv(g_, fn, ierr); CHK_ERR(ierr)
    
    write(fn, '(A, "/", A)') trim(dir_it), "r.csv"
    call dten2csv(R_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "p.csv"
    call dten2csv(P_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "ctot.csv"
    call cvec2csv(Ctot_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "c.csv"
    call cten2csv(c_, fn, ierr); CHK_ERR(ierr)

    call calc_probe(probe, ierr); CHK_ERR(ierr)
    write(fn,'("out/", I0, "/", A)') it, "probe.csv"
    call dvec2csv(probe, fn, ierr); CHK_ERR(ierr)
    
  end subroutine write_status
  subroutine print_status(it, ierr)
    integer, intent(in) :: it
    integer, intent(out) :: ierr
    integer KK
    ierr = 0
    write(*,'("t: ", F20.10)') it*n1t_*dt_    
    do KK = 1, np_
       write(*,*) "-- K=", KK, " --"
       write(*,*) "mode:", mode_(KK)
       write(*,'("Ctot:", F20.10, F20.10)') real(Ctot_(KK)), aimag(Ctot_(KK))       
       write(*,*) "R:", R_(KK, :, :)
       write(*,*) "P:", P_(KK, :, :)
       write(*,*) "C:", C_(KK,:,:)
    end do
    write(*,*) "========================="
  end subroutine print_status
  ! ==== utils ====
  subroutine make_nbasis(ops_typ, maxnd, KK, res, ierr)
    use Mod_PWGTO2
    character(3), intent(in) :: ops_typ(:)
    integer, intent(in) :: maxnd
    type(Obj_PWGTO) :: res
    integer, intent(out) :: ierr
    integer KK, A

    ierr = 0
    if(nnuc_(KK).eq.0) then
       MSG_ERR("nnuc_(KK)==0")
       write(0,*) "KK:", KK
       ierr = 1; return
    end if
    
    call PWGTO_new(res, nnuc_(KK), size(ops_typ), maxnd, ierr)
    CHK_ERR(ierr)
    do A = 1, nnuc_(KK)
       res%gs(A)     = g_(KK,A)
       res%Rs(A)     = R_(KK,A,1)
       res%Ps(A)     = P_(KK,A,1)
       res%ops_typ(A)= ops_typ(A)
       res%thetas(A) = 0
    end do
    call PWGTO_setup(res, ierr); CHK_ERR(ierr)
    
  end subroutine make_nbasis
  subroutine hamiltonian(HeIJ, XkIJ, P, res, ierr)
    use Mod_Const, only : II
    complex(kind(0d0)), intent(in) :: HeIJ(:,:), XkIJ(:,:,:)
    double precision, intent(in) :: P(:)
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    integer k, id(ne_,ne_), I

    ierr = 0

    id(:,:) = 0
    do I = 1, ne_
       id(I,I) = 1
    end do

    res(:,:) = HeIJ(:,:)
    do k = 1, nf_
       res(:,:) = res(:,:) + id(:,:)*1/(2*m_)*P(k)**2
       res(:,:) = res(:,:) - II*XkIJ(k,:,:)*P(k)/m_
       res(:,:) = res(:,:) - 1/(2*m_) * matmul(XkIJ(k,:,:), XkIJ(k,:,:))
    end do
    
  end subroutine hamiltonian
  subroutine calc_norm2_path(KK, res, ierr)
    use Mod_PWGTO2
    integer, intent(in) :: KK
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    type(Obj_PWGTO) :: nbasis
    complex(kind(0d0)) :: S(nnuc_(KK), nnuc_(KK))
    integer I, nn
    ierr = 0
    nn = nnuc_(KK)    
    call make_nbasis((/"0  "/), 1, KK, nbasis, ierr); CHK_ERR(ierr)
    call PWGTO_overlap(nbasis, 1, 1, S(:,:), ierr); CHK_ERR(ierr)
    call PWGTO_delete(nbasis, ierr)
    res = 0
    do I = 1, ne_
       res = res + real(dot_product(C_(KK,:nn,I), matmul(S, C_(KK,:nn,I))))
    end do
    
  end subroutine calc_norm2_path
  subroutine calc_probe(res, ierr)
    use Mod_PWGTO2
    double precision, intent(out) :: res(:)
    integer, intent(out) :: ierr
    integer KK
    ierr = 0
    res = 0
    do KK = 1, np_
       if(nnuc_(KK).ne.1) then
          MSG_ERR("nnuc(KK).ne.1")
          write(0,*) "KK:", KK
          write(0,*) "nnuc(KK):", nnuc_(KK)
          ierr = 1; return
       end if
       
       res(:) = res(:) + abs(Ctot_(KK)*C_(KK,1,:))**2
    end do
    
    !type(Obj_PWGTO) :: nbasis
    !integer KK, LL, I
    !complex(kind(0d0)) :: cumsum
    !complex(kind(0d0)) :: S(np_, maxnuc_, maxnuc_)
    !
    !res = 0
    !ierr = 0
!
!    do KK = 1, np_
!       call make_nbasis("0", KK, nbasis, ierr); CHK_ERR(ierr)
!       call PWGTO_overlap(nbasis, 1, 1, S(KK,:nnuc_(KK),:nnuc_(KK)), ierr); CHK_ERR(ierr)
!       call PWGTO_delete(nbasis, ierr); CHK_ERR(ierr)
!    end do
!    
!    do I = 1, ne_
!       cumsum = 0
!       do KK = 1, np_
!          do LL = 1, np_
!             do A = 1, nnuc_(KK)
!                do B = 1, nnuc_(KK)
!             cumsum = cumsum + conjg(CC_(KK)*d_(KK,I))*S(,)*CC_(LL)*d_(LL,I)
!          end do
!       end do
!       res(I) = real(cumsum)
!    end do
!
!    call PWGTO_delete(nbasis, ierr)
    
  end subroutine calc_probe
  subroutine calc_norm2(res, ierr)
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    double precision :: probe(ne_)
    ierr = 0
    call calc_probe(probe, ierr); CHK_ERR(ierr)
    res = sum(probe)
  end subroutine calc_norm2
end module Mod_DyQPBranch
  
module Mod_DyQPBranchPSA
  use Mod_DyQPBranch
  implicit none
  integer nlam_ ! size of eigen space.
  double precision, allocatable   :: Rlam_(:,:,:), Plam_(:,:,:)  ! R(KK,lam,:), P_(KK,lam,:) : eigen phase
  complex(kind(0d0)), allocatable :: glam_(:,:) ! g(KK,lam)   ! eigen width
  complex(kind(0d0)), allocatable :: U_(:,:,:)  ! U(KK,lam,I) ! eigen vector
  complex(kind(0d0)), allocatable :: Clam_(:,:) ! Clam(KK,I)  
contains
  subroutine DyQPBranchPSA_new(nlam, ierr)
    integer, intent(in) :: nlam
    integer, intent(out) :: ierr
    ierr = 0
    nlam_ = nlam
    if(nlam>ne_) then
       MSG_ERR("nlam > ne_")
       ierr = 1; return
    end if
    allocate(Rlam_(maxnp_,nlam,nf_))
    allocate(Plam_(maxnp_,nlam,nf_))
    allocate(glam_(maxnp_,nlam))
    allocate(U_(maxnp_,nlam,ne_))
    allocate(Clam_(maxnp_,nlam))
  end subroutine DyQPBranchPSA_New
  subroutine DyQPBranchPSA_setup(ierr)
    integer, intent(out) :: ierr
    ierr = 0       
  end subroutine DyQPBranchPSA_setup
  subroutine DyQPBranchPSA_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(Rlam_, Plam_, glam_, U_, Clam_)
  end subroutine DyQPBranchPSA_delete
  subroutine begin_clpsa(KK, ierr)
    integer, intent(in) :: KK
    integer, intent(out) :: ierr    
    integer :: lam
    
    ierr = 0
    if(nnuc_(KK).ne.1) then
       MSG_ERR("nnuc_(KK) != 1")
       ierr = 1; return
    end if
    if(mode_(KK).ne.1) then
       MSG_ERR("mode_(KK).ne.1")
       ierr = 1; return
    end if
    if(maxnp_<np_+nlam_-1) then
       MSG_ERR("too small maxnp")
       ierr = 1
       write(0,*) "maxnp:", maxnp_
       write(0,*) "np:", np_
       write(0,*) "nlam:", nlam_
       return 
    end if

    mode_(KK) = 2 ! to PSA mode

    do lam = 1, nlam_
       Rlam_(KK,lam,:) = R_(KK,1,:)
       Plam_(KK,lam,:) = P_(KK,1,:)
       glam_(KK,lam)   = g_(KK,1)
    end do
    
  end subroutine begin_clpsa
  subroutine update_clpsa(calc_H_X, KK, ierr)
    use Mod_math
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in)  :: KK
    integer, intent(out) :: ierr
    double precision :: w0(ne_), w1(ne_), w2(ne_)
    double precision :: dotR(nlam_,nf_), dotP(nlam_,nf_), dotx(nlam_), doty(nlam_)
    double precision :: R(nf_), P(nf_)
    complex(kind(0d0)) :: U(ne_,ne_), H(ne_,ne_), U1(ne_,ne_)
    complex(kind(0d0)) :: HeIJ(ne_,ne_), XkIJ(nf_,ne_,ne_)
    integer :: k, lam
    double precision :: dg, cc

    if(mode_(KK).ne.2) then
       MSG_ERR("mode_(KK).ne.2") 
       ierr = 1; return
    end if

    ! -- solve Hamilton equation on the averaged path --    
    R(:) = R_(KK,1,:)
    P(:) = P_(KK,1,:)
    call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)    
    do k = 1, nf_
       
       !P(k) = P(k) + dP_
       !call hamiltonian(HeIJ, XkIJ, P, H, ierr)
       !call lapack_zheev(ne_, H, w1, U, ierr); CHK_ERR(ierr)
       !P(k) = P(k) -2*dP_
       !call hamiltonian(HeIJ, XkIJ, P, H, ierr)
       !call lapack_zheev(ne_, H, w2, U, ierr); CHK_ERR(ierr)
       !P(k) = P(k) + dP_
       !do KK = 1, ne_
       !   dotR(KK,:) = (w1(KK)-w2(KK))/(2*dP_)
       !end do
       dotR(:,k) = Plam_(KK,:,k)/m_

       R(k) = R(k) + dR_
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr)
       call lapack_zheev(ne_, H, w1, U, ierr); CHK_ERR(ierr)          
       R(k) = R(k) -2*dR_
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr)
       call lapack_zheev(ne_, H, w2, U1, ierr); CHK_ERR(ierr)
!       if(abs(dot_product(U(:,1), U1(:,1))) < abs(dot_product(U(:,1), U1(:,2)))) then
!          tmp = w2(1)
!          w2(1) = w2(2)
!          w2(2) = tmp
!       end if
       R(k) = R(k) + dR_
       do lam = 1, nlam_ ! assume nnuc_(KK) == ne_
          dotP(lam,k) = -(w1(lam)-w2(lam))/(2*dR_)
       end do
    end do
    
    if(gauss_mode_.eq."thawed") then
       dg = 2*real(g_(KK,1))*sum(R)*dR_/sum(R*R)  !dg.(R,R) = g.2(R,dR) => dg = 2g.(R,dR)/(R,R)
       g_(KK,1) = g_(KK,1) + dg
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
       call lapack_zheev(ne_, H, w1, U, ierr); CHK_ERR(ierr)
       g_(KK,1) = g_(KK,1) - 2*dg
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
       call lapack_zheev(ne_, H, w2, U, ierr); CHK_ERR(ierr)
       g_(KK,1) = g_(KK,1) + dg
       doty(:) = -4*real(g_(KK,1))**2*(w1(:)-w2(:))/(2*dg)

       g_(KK,1) = g_(KK,1) + II*dg
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
       call lapack_zheev(ne_, H, w1, U, ierr); CHK_ERR(ierr)
       g_(KK,1) = g_(KK,1) - 2*II*dg
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
       call lapack_zheev(ne_, H, w2, U, ierr); CHK_ERR(ierr)
       g_(KK,1) = g_(KK,1) + II*dg
       dotx(:) = +4*real(g_(KK,1))**2*(w1(:)-w2(:))/(2*dg)
    else if(gauss_mode_.eq."frozen") then
       dotx(:) = 0; doty(:) = 0
    else
       MSG_ERR("unsupported gauss_mode_"); ierr=1; return
    end if
    
    ! -- diag hamiltonian at averaged point --    
    R(:) = R_(KK,1,:)
    P(:) = P_(KK,1,:)
    call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
    call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
    call lapack_zheev(ne_, H, w0, U_(KK,:,:), ierr); CHK_ERR(ierr)

    ! -- re-expand electronic wave packet with eigen state.
    !   \chi(Q;<x(t)>) \Theta(r;Q) = \chi(Q;<x(t)>) \sum_I CC(KK,I) \Phi_I(r;Q)
    !                              = \chi(Q;<x(t)>) \sum_Lam C_Lam \Phi_Lam(r;Q,t)
    ! To determine D_Lam, take inner product with Phi_Lam and get
    !   C_Lam = \sum_I C_I <\chi(<x>)Phi_Lam|\chi(<x>)\Phi_I>
    !         = \sum_I U_{I,lam}^* C_I
    Clam_(KK,:) = matmul(transpose(conjg(U_(KK,:,:))), C_(KK,1,:))

    ! -- update R_lam,P_lam and <R>,<P> --
    do lam = 1, nlam_
       glam_(KK,lam)   = glam_(KK,lam)   + (dotx(lam) + II*doty(lam))
       Rlam_(KK,lam,:) = Rlam_(KK,lam,:) + dotR(lam,:) * dydt_
       Plam_(KK,lam,:) = Plam_(KK,lam,:) + dotP(lam,:) * dydt_
       cc = abs(Clam_(KK,lam))**2
       g_(KK,1)   = g_(KK,1)   + cc * dydt_*(dotx(lam)+II*doty(lam))
       R_(KK,1,:) = R_(KK,1,:) + cc * dydt_*dotR(lam,:)
       P_(KK,1,:) = P_(KK,1,:) + cc * dydt_*dotP(lam,:)
    end do
    
    ! -- update C(t) by Hamiltonian at averaged point --
    call intet_diag(ne_, H(:,:), dydt_, C_(KK,1,:), ierr); CHK_ERR(ierr)
    
  end subroutine update_clpsa
  subroutine end_clpsa(KK, ierr)
    integer, intent(in) :: KK
    integer, intent(out) :: ierr    
    integer K, lam

    ierr = 0

    if(mode_(KK).ne.2) then
       MSG_ERR("mode_(KK).ne.2")
       ierr = 1; return
    end if

    mode_(KK) = 1 ! SET
    nnuc_(KK) = 1
    R_(KK,1,:) = Rlam_(KK,1,:)
    P_(KK,1,:) = Plam_(KK,1,:)
    g_(KK,1)   = glam_(KK,1)    
    C_(KK,1,:) = U_(KK,1,:)
    Ctot_(KK) = Clam_(KK,1)

    do lam = 2, nlam_
       K = np_+lam-1
       mode_(K) = 1 ! SET
       nnuc_(K) = 1
       R_(K,1,:) = Rlam_(KK,lam,:)
       P_(K,1,:) = Plam_(KK,lam,:)
       g_(K,1)   = glam_(KK,lam)       
       C_(K,1,:) = U_(KK,lam,:)
       Ctot_(K)  = Clam_(KK,lam)
    end do
    
    np_ = np_ + nlam_-1
    
  end subroutine end_clpsa
end module Mod_DyQPBranchPSA
