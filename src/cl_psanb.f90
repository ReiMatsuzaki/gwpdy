#include "macros.fpp"
module Mod_ClPSANB
  implicit none
  ! -- data size --
  integer nf_, ne_, np_
  ! -- variable --
  double precision, allocatable :: R_(:,:), P_(:,:), Rave_(:), Pave_(:)
  complex(kind(0d0)), allocatable :: c_(:,:), cave_(:), cc_(:)
  ! -- time constant --
  double precision :: m_ ! mass  
  double precision :: dt_, dR_, dP_
  integer :: nt_, n1t_
  integer :: maxnp_
  character(8) :: inte_RP_     ! euler or RK4
  character(3) :: mode_        ! SET or PSA
  ! -- intermediate --
  double precision :: dydt_
  !  integer :: status_(:)        ! 0: active, 1:death
contains
  ! ==== Main ====
  subroutine ClPSANB_new(nf, ne, maxnp, ierr)
    integer, intent(in) :: nf, ne, maxnp
    integer, intent(out) :: ierr

    ierr = 0
    nf_ = nf
    ne_ = ne
    maxnp_ = maxnp
    np_ = 1
    
    allocate(R_(maxnp,nf), P_(maxnp,nf_), c_(maxnp,ne_), cave_(ne_), cc_(maxnp))
    allocate(Rave_(nf), Pave_(nf))
    !    allocate(status_(maxnp))    
    R_ = 0
    P_ = 0
    c_ = 0
    CC_ = 0
    c_(1,1) = 1
    CC_(1) = 1
    Rave_ = 0
    Pave_ = 0
    !    status_(:) = 1
    !    status_(1) = 0 ! active only first path

    m_ = 1
    dt_ = 1
    dR_ = 0.01d0
    dP_ = 0.01d0
    nt_ = 10
    n1t_ = 1

    inte_RP_ = "euler"
    mode_ = "SET"
    
  end subroutine ClPSANB_new
  subroutine ClPSANB_setup(ierr)
    integer, intent(out) :: ierr
    integer KK
    double precision :: norm2
    ierr = 0
    do KK = 1, np_
       norm2 = real(dot_product(c_(KK,:), c_(KK,:)))
       if(norm2<1.0d-10) then
          MSG_ERR("norm is too small")
          ierr = 1; return
       end if
       c_(KK,:) = c_(KK,:)/sqrt(norm2)
    end do

    norm2 = 0
    do KK = 1, np_
       norm2 = norm2 + abs(CC_(KK))**2
    end do
    if(norm2<1.0d-10) then
       MSG_ERR("|CC| is too small")
       ierr = 1; return
    end if
    CC_(:) = CC_(:)/sqrt(norm2)

    if(inte_RP_.ne."euler" .and. inte_RP_.ne."RK4") then
       MSG_ERR("unsupported inte_RP_")
       ierr = 1; return
    end if

    dydt_ = dt_/n1t_
    
  end subroutine ClPSANB_setup
  subroutine ClPSANB_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(R_, P_, Rave_, Pave_, c_, cave_, cc_)
  end subroutine ClPSANB_delete
  ! ==== calc ====
  subroutine update(calc_H_X, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    
    if(mode_.eq."SET") then
       call update_set(calc_H_X, ierr); CHK_ERR(ierr)
    else if(mode_.eq."PSA") then
       call update_psa(calc_H_X, ierr); CHK_ERR(ierr)
    end if
    
  end subroutine update
  subroutine update_set(calc_H_X, ierr)
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    double precision  :: dotR(nf_), dotP(nf_)
    integer :: KK, k
    complex(kind(0d0)) :: HeIJ(ne_,ne_), XkIJ(nf_,ne_,ne_)
    complex(kind(0d0)) :: H0(ne_,ne_), H1(ne_,ne_), H2(ne_,ne_), c(ne_)
    double precision :: R(nf_), P(nf_)
        
    do KK = 1, np_
       R(:) = R_(KK,:)
       P(:) = P_(KK,:)
       c(:) = c_(KK,:)
       call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
       call hamiltonian(HeIJ, XkIJ, P, H0, ierr)
       do k = 1, nf_
          P(k) = P(k) + dP_
          call hamiltonian(HeIJ, XkIJ, P, H1, ierr)
          P(k) = P(k) -2*dP_
          call hamiltonian(HeIJ, XkIJ, P, H2, ierr)
          P(k) = P(k) + dP_
          dotR(:) = real(dot_product(c, matmul(H1-H2, c)))/(2*dP_)
       end do
       do k = 1, nf_
          R(k) = R(k) + dR_
          call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
          call hamiltonian(HeIJ, XkIJ, P, H1, ierr)
          R(k) = R(k) -2*dR_
          call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
          call hamiltonian(HeIJ, XkIJ, P, H2, ierr)
          R(k) = R(k) + dR_
          dotP(:) = -real(dot_product(c, matmul(H1-H2, c)))/(2*dR_)
       end do
       
       if(.not. dotR(1).eq.dotR(1)) then
          MSG_ERR("dotR is NaN")
          ierr = 1; return
       end if
       if(.not. dotP(1).eq.dotP(1)) then
          MSG_ERR("dotP is NaN")
          ierr = 1; return
       end if

       ! -- update --
       call intet_diag(ne_, H0(:,:), c_(KK,:), ierr); CHK_ERR(ierr)
       R_(KK,:) = R_(KK,:) + dotR(:)*dydt_
       P_(KK,:) = P_(KK,:) + dotP(:)*dydt_
    end do

  end subroutine update_set
  subroutine begin_psa(ierr)
    integer, intent(out) :: ierr
    integer :: KK
    
    ierr = 0
    if(np_.ne.1) then
       MSG_ERR("np_ must be 1")
       ierr = 1; return
    end if

    mode_ = "PSA"     
    np_ = ne_
    !    status_(1:ne_) = 0
    Rave_(:) = R_(1,:)
    Pave_(:) = P_(1,:)
    cave_(:)  = c_(1,:)
    do KK = 2, ne_
       R_(KK,:) = R_(1,:)
       P_(KK,:) = P_(1,:)
    end do
    
  end subroutine begin_psa
  subroutine update_psa(calc_H_X, ierr)
    use Mod_math
    interface
       subroutine calc_H_X(Q, HeIJ, XkIJ, ierr)
         double precision, intent(in) :: Q(:)
         complex(kind(0d0)), intent(out) :: HeIJ(:,:), XkIJ(:,:,:)
         integer, intent(out) :: ierr
       end subroutine calc_H_X
    end interface
    integer, intent(out) :: ierr
    double precision :: w0(ne_), w1(ne_), w2(ne_)
    double precision :: dotR(ne_,nf_), dotP(ne_,nf_), R(nf_), P(nf_)
    complex(kind(0d0)) :: U(ne_,ne_), H(ne_,ne_), U1(ne_,ne_)
    complex(kind(0d0)) :: HeIJ(ne_,ne_), XkIJ(nf_,ne_,ne_)
    complex(kind(0d0)) :: D(ne_)
    integer :: KK, k
    double precision :: tmp

    ! -- solve Hamilton equation on the averaged path --    
    R(:) = Rave_(:)
    P(:) = Pave_(:)
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
       dotR(:,k) = P_(:np_,k)/m_

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
       do KK = 1, ne_
          dotP(KK,k) = -(w1(KK)-w2(KK))/(2*dR_)
       end do
    end do

    ! -- diag coupled hamiltonian --    
    R(:) = Rave_(:)
    P(:) = Pave_(:)
    call calc_H_X(R(:), HeIJ(:,:), XkIJ(:,:,:), ierr); CHK_ERR(ierr)
    call hamiltonian(HeIJ, XkIJ, P, H, ierr); CHK_ERR(ierr)
    call lapack_zheev(ne_, H, w0, U, ierr); CHK_ERR(ierr)
    call lapack_zgesv_1(ne_, U, cave_(:), D(:), ierr); CHK_ERR(ierr)

    ! -- update R,P --
    do KK = 1, np_
       R_(KK,:) = R_(KK,:) + dotR(KK,:)*dydt_
       P_(KK,:) = P_(KK,:) + dotP(KK,:)*dydt_
       !       write(*,*) KK, dotP(KK,:)*dydt_
       Rave_(:) = Rave_(:) + abs(D(KK))**2 * dotR(KK,:)*dydt_
       Pave_(:) = Pave_(:) + abs(D(KK))**2 * dotP(KK,:)*dydt_
    end do

    ! -- update C_ave(t), C^K_I(t) --
    call intet_diag(ne_, H(:,:), cave_(:), ierr); CHK_ERR(ierr)
    do KK = 1, ne_
       c_(KK,:) = U(:,KK)
       cc_(KK)  = D(KK)
    end do
    
  end subroutine update_psa
  subroutine end_psa(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    mode_ = "SET"    
  end subroutine end_psa
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
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,val")')
    do KK = 1, maxnp_
       do k = 1, nf_
          write(ifile,'(I0,",",I0,",",F20.10)') KK,k,R_(KK,k)
       end do
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '(A, "/", A)') trim(dir_it), "p.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,val")')
    do KK = 1, maxnp_
       do k = 1, nf_
          write(ifile,'(I0,",",I0,",",F20.10)') KK,k,P_(KK,k)
       end do
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '(A, "/", A)') trim(dir_it), "c.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,re,im")')
    do KK = 1, maxnp_
       do I = 1, ne_
          write(ifile,'(I0,",",I0,",",F20.10,",",F20.10)') KK,I,real(C_(KK,I)),aimag(C_(KK,I))
       end do
    end do
    close(ifile)
    ifile = ifile + 1

    write(fn, '(A, "/", A)') trim(dir_it), "cc.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile, '("i,j,re,im")')
    do KK = 1, maxnp_
       write(ifile,'(I0,",",I0,",",F20.10,",",F20.10)') KK,I,real(CC_(KK)),aimag(CC_(KK))
    end do
    close(ifile)
    ifile = ifile + 1

    call calc_probe(probe, ierr); CHK_ERR(ierr)
    write(fn,'("out/", I0, "/", A)') it, "probe.csv"
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)
    write(ifile,'("i,val")')
    do I = 1, ne_
       write(ifile,'(I0,",",F20.10)') I, probe(I)
    end do
    close(ifile)
    ifile = ifile + 1
    
  end subroutine write_status
  ! ==== utils ====
  subroutine calc_probe(res, ierr)
    double precision, intent(out) :: res(:)
    integer, intent(out) :: ierr
    integer KK, I
    complex(kind(0d0)) :: cumsum
    
    res = 0
    ierr = 0
    do I = 1, ne_
       cumsum = 0
       do KK = 1, np_
          cumsum = cumsum + abs(c_(KK,I)*cc_(KK))**2
       end do
       res(I) = real(cumsum)
    end do
    
  end subroutine calc_probe
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
end module Mod_ClPSANB
  
