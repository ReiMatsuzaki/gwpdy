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
  complex(kind(0d0)), allocatable :: c_(:), d_(:,:,:) ! c_(K), d_(K,A,I)
  ! -- time constant --
  double precision :: m_ ! mass  
  double precision :: dt_, dR_, dP_
  integer :: nt_, n1t_
  character(8) :: inte_RP_     ! euler or RK4
  character(3) :: mode_        ! SET or PSA
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
    maxnp_ = maxnp
    np_ = 1

    allocate(nnuc_(maxnp))
    allocate(R_(maxnp,maxnuc,nf), P_(maxnp,maxnuc,nf_), c_(maxnp), d_(maxnp, maxnp, ne))    
    allocate(g_(maxnp,maxnuc))
    nnuc_(:) = 0
    nnuc_(1) = 1
    g_ = 1
    R_ = 0
    P_ = 0
    d_ = 0
    C_ = 0
    d_(1,1,1) = 1
    C_(1) = 1
    
    m_ = 1
    dt_ = 1
    dR_ = 0.01d0
    dP_ = 0.01d0
    nt_ = 10
    n1t_ = 1

    inte_RP_ = "euler"
    mode_ = "SET"
    
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
          ierr = 1; return
       end if
       d_(KK,:,:) = d_(KK,:,:)/sqrt(norm2)
    end do

    call calc_norm2(norm2, ierr); CHK_ERR(ierr)
    if(norm2<1.0d-10) then
       MSG_ERR("global norm is too small")
       ierr = 1; return
    end if
    C_(:) = C_(:)/sqrt(norm2)

    if(inte_RP_.ne."euler" .and. inte_RP_.ne."RK4") then
       MSG_ERR("unsupported inte_RP_")
       ierr = 1; return
    end if

    dydt_ = dt_/n1t_
    
  end subroutine DyQPBranch_setup
  subroutine DyQPBranch_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(nnuc_)
    deallocate(R_, P_, d_, C_)
    deallocate(g_)
  end subroutine DyQPBranch_delete
  ! ==== calc ====
  subroutine update_set(calc_H_X, KK, ierr)
    use Mod_PWGTO
    use Mod_const, only : II
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(in) :: KK
    integer, intent(out) :: ierr
    double precision  :: dotR(nf_), dotP(nf_)
    integer :: k
    complex(kind(0d0)) :: HeIJ(np_,ne_,ne_), XkIJ(np_,nf_,ne_,ne_)
    complex(kind(0d0)) :: H0(ne_,ne_), H1(ne_,ne_), H2(ne_,ne_), d(ne_)
    double precision :: R(nf_), P(nf_)
    type(Obj_PWGTO) :: nbasis

    if(nnuc_(KK).ne.1) then
       MSG_ERR("nnuc_(KK)!=1")
       ierr = 1; return
    end if
    
    R(:) = R_(KK,1,:); P(:) = P_(KK,1,:); d(:) = d_(KK,1,:)
    do k = 1, nf_
       R(k) = R(k) + dR_
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H1, ierr)
       R(k) = R(k) -2*dR_
       call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H2, ierr)
       R(k) = R(k) + dR_
       dotP(:) = -real(dot_product(d, matmul(H1-H2, d)))/(2*dR_)
    end do
        
    call calc_H_X(R(:), HeIJ(KK,:,:), XkIJ(KK,:,:,:), ierr); CHK_ERR(ierr)
    call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P(:), H0(:,:), ierr)
    do k = 1, nf_
       P(k) = P(k) + dP_
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H1, ierr); CHK_ERR(ierr)
       P(k) = P(k) -2*dP_
       call hamiltonian(HeIJ(KK,:,:), XkIJ(KK,:,:,:), P, H2, ierr); CHK_ERR(ierr)
       P(k) = P(k) + dP_
       dotR(:) = real(dot_product(d, matmul(H1-H2, d)))/(2*dP_)
    end do
    
    ! -- update --       
    call intet_diag(ne_, H0(:,:), d_(KK,1,:), ierr); CHK_ERR(ierr)
    R_(KK,1,:) = R_(KK,1,:) + dotR(:)*dydt_
    P_(KK,1,:) = P_(KK,1,:) + dotP(:)*dydt_
    
  end subroutine update_set
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
    use Mod_sys, only : open_w, mkdirp_if_not
    integer, intent(in) :: it
    integer, intent(out) :: ierr
    character(50) :: fn, dir_it
    integer :: ifile
    integer :: KK, I, k
    double precision :: probe(ne_)
    ifile = 1242

    write(dir_it,'("out/", I0)') it
    call mkdirp_if_not(dir_it)

    write(fn, '(A, "/", A)') trim(dir_it), "param.json"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, *) "{"
    write(ifile, '(A,A,A,A,I0,A)') '"', "np",   '"', ":", np_, ","
    write(ifile, '(A,A,A,A,A,A,A)') '"',   "mode", '"', ":", '"', mode_, '"'
    write(ifile, *) "}"
    close(ifile)
    ifile = ifile + 1

    write(fn, '(A, "/", A)') trim(dir_it), "r.csv"
    call dten2csv(R_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "p.csv"
    call dten2csv(P_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "c.csv"
    call cvec2csv(C_, fn, ierr); CHK_ERR(ierr)

    write(fn, '(A, "/", A)') trim(dir_it), "d.csv"
    call cten2csv(d_, fn, ierr); CHK_ERR(ierr)

    call calc_probe(probe, ierr); CHK_ERR(ierr)
    write(fn,'("out/", I0, "/", A)') it, "probe.csv"
    call dvec2csv(probe, fn, ierr); CHK_ERR(ierr)
    
  end subroutine write_status
  ! ==== utils ====
  subroutine calc_norm2_path(KK, res, ierr)
    use Mod_PWGTO
    integer, intent(in) :: KK
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    type(Obj_PWGTO) :: nbasis
    complex(kind(0d0)) :: S(nnuc_(KK), nnuc_(KK))
    integer I
    ierr = 0
    
    call make_nbasis("0", KK, nbasis, ierr); CHK_ERR(ierr)
    call PWGTO_overlap(nbasis, 1, 1, S(:,:), ierr); CHK_ERR(ierr)
    call PWGTO_delete(nbasis, ierr)
    res = 0
    do I = 1, ne_
       res = res + real(dot_product(d_(KK,:,I), matmul(S, d_(KK,:,I))))
    end do
    
  end subroutine calc_norm2_path
  subroutine make_nbasis(mode, KK, res, ierr)
    use Mod_PWGTO
    character(*), intent(in) :: mode
    type(Obj_PWGTO) :: res
    integer, intent(out) :: ierr
    integer, parameter :: maxnd=2
    integer numNCs
    integer KK, A

    ierr = 0
    if(mode.eq."0") then
       numNCs = 1
    else if(mode.eq."0RP") then
       numNCs = 3
    else
       MSG_ERR("unsupported")
       ierr = 1; return
    end if

    if(KK.eq.0) then
       MSG_ERR("KK==0 is not impl")
       ierr = 1 ; return
    end if
    
    call PWGTO_new(res, nnuc_(KK), numNCs, maxnd, ierr); CHK_ERR(ierr)
    do A = 1, nnuc_(KK)
       res%gs(A)     = g_(KK,A)
       res%Rs(A)     = R_(KK,A,1)
       res%Ps(A)     = P_(KK,A,1)
       res%thetas(KK) = 0
    end do
    if(mode.eq."0RP") then
       res%typ(2) = "dR"
       res%typ(3) = "dP"
    end if
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
  subroutine intet_diag(n, H, c, ierr)
    use Mod_const, only : II
    use Mod_math, only  : lapack_zheev
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: H(n,n)
    complex(kind(0d0)), intent(inout) :: c(n)
    integer, intent(out) :: ierr
    double precision :: w(n)
    complex(kind(0d0)) :: U(n,n), UH(n,n)

    call lapack_zheev(n, H, w, U, ierr); CHK_ERR(ierr)
    UH(:,:) = conjg(transpose(U(:,:)))
    c(:)   = matmul(UH(:,:),       c(:))
    c(:)   = exp(-II*w(:)*dydt_) * c(:)
    c(:)   = matmul(U(:,:),        c(:))
    
  end subroutine intet_diag
  subroutine dten2csv(x, fn, ierr)
    double precision, intent(in) :: x(:,:,:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i, j, k
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,j,k,val")')
    do i = 1, size(x, 1)
       do j = 1, size(x, 2)
          do k = 1, size(x, 3)
             write(ifile,'(I0,",",I0,",",I0,",",F20.10)') i,j,k,x(i,j,k)
          end do
       end do
    end do
    close(ifile)
    
  end subroutine dten2csv
  subroutine cten2csv(x, fn, ierr)
    complex(kind(0d0)), intent(in) :: x(:,:,:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i, j, k
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,j,k,re,im")')
    do i = 1, size(x, 1)
       do j = 1, size(x, 2)
          do k = 1, size(x, 3)
             write(ifile,'(I0,",",I0,",",I0,",",F20.10,F20.10)') i,j,k,real(x(i,j,k)),aimag(x(i,j,k))
          end do
       end do
    end do
    close(ifile)
    
  end subroutine cten2csv
  subroutine dvec2csv(x, fn, ierr)
    double precision, intent(in) :: x(:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,val")')
    do i = 1, size(x, 1)
       write(ifile,'(I0,",",F20.10)') i,x(i)
    end do
    close(ifile)
    
  end subroutine dvec2csv
  subroutine cvec2csv(x, fn, ierr)
    complex(kind(0d0)), intent(in) :: x(:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,j,k,re,im")')
    do i = 1, size(x, 1)
       write(ifile,'(I0,",",F20.10,F20.10)') i,real(x(i)),aimag(x(i))
    end do
    close(ifile)
    
  end subroutine cvec2csv
  ! ==== utils(global) ====
  subroutine calc_probe(res, ierr)
    use Mod_PWGTO
    double precision, intent(out) :: res(:)
    integer, intent(out) :: ierr

    MSG_ERR("not impl")
    ierr = 1
    res = 0
    
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
  double precision, allocatable   :: RR_(:,:), PP_(:,:)  ! RR_(KK,:), PP_(KK,:) : average path
  complex(kind(0d0)), allocatable :: gg_(:), dd_(:,:)    ! gg_(KK), dd_(KK,I)
contains
  subroutine DyQPBranchPSA_new(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    allocate(RR_(maxnp_,nf_), PP_(maxnp_,nf_), gg_(maxnp_))
  end subroutine DyQPBranchPSA_New
  subroutine DyQPBranchPSA_delete(ierr)
    integer, intent(out) :: ierr
    deallocate(RR_, PP_, gg_)
  end subroutine DyQPBranchPSA_delete
  subroutine begin_clpsa(KK, ierr)
    integer, intent(in) :: KK
    integer, intent(out) :: ierr    
    integer :: A
    
    ierr = 0
    if(nnuc_(KK).ne.1) then
       MSG_ERR("nnuc_(KK) != 1")
       ierr = 1; return
    end if
    if(mode_.ne."SET") then
       MSG_ERR("mode must be SET")
       ierr = 1; return
    end if

    mode_ = "PSA"     
    nnuc_(1) = ne_
    do A = 1, nnuc_(1)
       R_(1, A, :) = R_(1,1,:)
    end do
    
    RR_(KK,:) = R_(KK,1,:)
    PP_(KK,:) = R_(KK,1,:)
    dd_(KK,:) = d_(KK,1,:)
    do A = 2, nnuc_(KK)
       R_(KK,A,:) = R_(KK,A,:)
       P_(KK,A,:) = P_(KK,:)
       d_(KK,A,:) = d_(KK,A,:)
    end do
    
  end subroutine begin_clpsa
  subroutine update_clpsa(calc_H_X, KK, ierr)
    use Mod_math
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
    double precision :: dotR(nnuc_(KK),nf_), dotP(nnuc_(KK),nf_), R(nf_), P(nf_)  ! ne .eq. nnuc_(KK)
    complex(kind(0d0)) :: U(ne_,ne_), H(ne_,ne_), U1(ne_,ne_)
    complex(kind(0d0)) :: HeIJ(ne_,ne_), XkIJ(nf_,ne_,ne_)
    complex(kind(0d0)) :: D(ne_)
    integer :: k, A
    double precision :: tmp

    ! -- solve Hamilton equation on the averaged path --    
    R(:) = RR_(KK,:)
    P(:) = PP_(KK,:)    
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
       dotR(:,k) = P_(KK,:,k)/m_

       R(k) = R(k) + dR_
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr)
       call lapack_zheev(ne_, H, w1, U, ierr); CHK_ERR(ierr)          
       R(k) = R(k) -2*dR_
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H, ierr)
       call lapack_zheev(ne_, H, w2, U1, ierr); CHK_ERR(ierr)
       if(abs(dot_product(U(:,1), U1(:,1))) < abs(dot_product(U(:,1), U1(:,2)))) then
          write(0,*) "swap!"
          tmp = w2(1)
          w2(1) = w2(2)
          w2(2) = tmp
       end if
       R(k) = R(k) + dR_
       do A = 1, nnuc_(KK) ! assume nnuc_(KK) == ne_
          dotP(A,k) = -(w1(A)-w2(A))/(2*dR_)
       end do
    end do

    ! -- diag hamiltonian at averaged point --    
    R(:) = RR_(KK,:)
    P(:) = PP_(KK,:)
    call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
    call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
    call lapack_zheev(ne_, H, w0, U, ierr); CHK_ERR(ierr)

    ! -- re-expand electronic wave packet with eigen state.
    ! \chi(Q;<x(t)>) \Theta(r;Q) = \chi(Q;<x(t)>) \sum_I dd(KK,I) \Phi_I(r;Q)
    !                            = \chi(Q;<x(t)>) \sum_Lambda D_Lam \Theta_Lam
    ! To determine D_Lam, take inner product and get
    !      dd(KK,I) = \sum_\Lambda D_\Lambda <\Phi_I|Theta_\Lambda>  
    !               = \sum_\Lambda D_\Lambda U_{I\Lambda}         
    call lapack_zgesv_1(ne_, U, dd_(KK,:), D(:), ierr); CHK_ERR(ierr)

    ! -- update R,P and <R>,<P> --
    do A = 1, nnuc_(KK)
       R_(KK,A,:) = R_(KK,A,:) + dotR(KK,A,:)*dydt_
       P_(KK,A,:) = P_(KK,A,:) + dotP(KK,A,:)*dydt_
       RR_(KK,:)  = RR_(KK,:)  + abs(D(A))**2 * dotR(KK,A,:)*dydt_
       PP_(KK,:)  = PP_(KK,:)  + abs(D(A))**2 * dotP(KK,A,:)*dydt_
    end do
    
    ! -- update d_ave(t) by Hamiltonian at averaged point --
    call intet_diag(ne_, H(:,:), dd_(KK,:), ierr); CHK_ERR(ierr)

    ! -- update d(t) by expansion of eigen state --
    ! \chi(Q;<x>) \sum_I dd_I\Phi_I(r;Q) = \sum_{AI} d_{AI}\chi_A\Phi_I
    ! => \chi(Q;<x>) dd_I = sum_A d_{AI} \chi_A
    ! => dd_I = sum_A <\chi(<x>)|\chi_A> d_{AI}
    ! see Eq.(*)
    do A = 1, ne_
       dd_(KK,A,:) = matmul(U(:,:), D(:))
    end do
    
  end subroutine update_clpsa
  subroutine end_clpsa(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    mode_ = "SET"    
  end subroutine end_clpsa
end module Mod_DyQPBranchPSA
