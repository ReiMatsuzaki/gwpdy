#include "macros.fpp"
module Mod_PWGTO
  ! Plane Wave Gauss Type Orbital
  !   G^k_i(Q) = sum_j c_kji(Q-Ri)^{nkji} exp[-g(Q-Ri)^2 + iPi(Q-Ri)]
  implicit none
  type Obj_PWGTO
     ! - data size -
     integer :: num, nf, maxnd
     ! - variable -
     integer, allocatable :: ns(:)
     complex(kind(0d0)), allocatable :: gs(:) 
     double precision,   allocatable :: Rs(:), Ps(:), thetas(:)
     character(3), allocatable :: typ(:)
     ! - intermediate -
     integer, allocatable :: nlc(:,:)     ! nlc(n,A) number of nth polynomial
     integer, allocatable :: nns(:,:,:) ! nns(n,A,:) powers of nth polynomial
     complex(kind(0d0)), allocatable :: cs(:,:,:) ! cs(n,A,:) : coeff of nthpoly
     complex(kind(0d0)), allocatable :: gP(:,:), RP(:,:), eP(:,:)
     complex(kind(0d0)), allocatable :: d(:,:,:,:,:) ! d(A,B,0:maxnd,0:maxnd,0:maxnd*2)
  end type Obj_PWGTO
contains
  ! -- main --
  subroutine PWGTO_new(this, num, numNCs, maxnd, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(in) :: num, numNCs, maxnd
    integer, intent(out) :: ierr

    ierr = 0
    this%num   = num
    this%maxnd = maxnd
    allocate(this%ns(num))
    allocate(this%gs(num), this%Rs(num), this%Ps(num), this%thetas(num))
    allocate(this%gP(num,num), this%RP(num,num), this%eP(num,num))
    allocate(this%d(num,num,0:maxnd,0:maxnd,0:maxnd*2))
    this%ns = 0
    this%gs = (1.0d0, 0.0d0)
    this%Rs = 0.0d0
    this%Ps = 0.0d0
    this%thetas = 0.0d0

    allocate(this%typ(numNCs))
    this%typ(1) = "0"
    allocate(this%nlc(numNCs, num))
    allocate(this%nns(numNCs, num, 3))
    allocate(this%cs(numNCs, num, 3))
    
  end subroutine PWGTO_new
  subroutine PWGTO_setup(this, ierr)
    use Mod_const, only : II
    use Mod_math, only : gtoint
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    integer A, B, nA, n, inc
    complex(kind(0d0)) :: gA, c0
    complex(kind(0d0)), allocatable :: gg(:)
    
    ierr = 0
    allocate(gg(0:maxval(this%ns(:))*2+2))

    ! -- normalization term / derivative basis --
    do A = 1, this % num
       nA = this%ns(A)
       gA = this%gs(A)
       call gtoint(2*nA, conjg(gA)+gA, gg, ierr); CHK_ERR(ierr)
       c0 = 1.0d0/sqrt(gg(2*nA))

       do inc = 1, size(this%typ)
          select case(this%typ(inc))
          case("0")
             this%nlc(inc,A) = 1
             this%cs(inc,A, 1) = c0
             this%nns(inc,A, 1) = this % ns(A)
          case("dR")
             if(nA .eq. 0) then
                this%nlc(inc,A) = 2
             else
                this%nlc(inc,A) = 3
             end if
             this%cs(inc,A, 1) = c0*2*gA
             this%nns(inc,A, 1) = nA+1
             this%cs(inc,A, 2) = -c0*ii*this%Ps(A)
             this%nns(inc,A, 2) = nA
             if(nA .ne. 0) then
                this%cs(inc,A, 3) = nA*c0
                this%nns(inc,A, 3) = nA-1
             end if
          case("dP")
             this%nlc(inc,A) = 1
             this%cs(inc,A, 1) = c0*ii
             this%nns(inc,A, 1) = nA+1
          case default
             MSG_ERR("unsupported typ")
             ierr = 1; return
          end select
       end do
    end do

    ! -- combination of GTO --
    do B = 1, this % num
       do A = 1, this % num
          call prod_gauss(&
               this%gs(A), this%Rs(A), this%Ps(A), this%thetas(A), &
               this%gs(B), this%Rs(B), this%Ps(B), this%thetas(B), &
               this%gP(A,B), this%RP(A,B), this%eP(A,B))
          n = this % maxnd
          call hermite_coef_d(this%gP(A,B), this%RP(A,B), &
               this%Rs(A), this%Rs(B), &
               n,n, this%d(A,B,:,:,:))
       end do
    end do
    
  end subroutine PWGTO_setup
  subroutine PWGTO_dump(this, ifile, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(in) :: ifile
    integer, intent(out) :: ierr
    ierr = 0
    write(ifile,*) "==== PWGTO ===="
    write(ifile,*) "ns:", this%ns
    write(ifile,*) "g:", this%gs
    write(ifile,*) "R:", this%Rs
    write(ifile,*) "P:", this%Ps
    write(ifile,*) "==============="
  end subroutine PWGTO_dump
  subroutine PWGTO_delete(this, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(this%gs, this%Rs, this%Ps, this%thetas)
    deallocate(this%nlc, this%nns, this%cs)
    deallocate(this%gP, this%RP, this%eP, this%d)
  end subroutine PWGTO_delete
  ! -- calc  --
  subroutine PWGTO_overlap(this, ibra, iket, X, ierr)
    type(Obj_PWGTO), intent(in) :: this
    integer, intent(in) :: ibra, iket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    complex(kind(0d0)) acc, tmp, hint, cAcB

    ierr = 0
    call check_matrix(this, X, ierr); CHK_ERR(ierr)

    do A = 1, this % num
       do B = 1, this % num
          call hermite_1dint(this%gP(A,B), hint)
          acc = 0.0d0
          do i = 1, this%nlc(ibra,A)
             do j = 1, this%nlc(iket,B)
                nA = this%nns(ibra,A, i)
                nB = this%nns(iket,B, j)
                cAcB = conjg(this%cs(ibra,A,i)) * this%cs(iket,B,j)
                tmp = cAcB * this%eP(A,B) * this%d(A,B, nA,nB,0) * hint
                acc = acc + tmp
             end do
          end do
          X(A,B) = acc
       end do
    end do
    
  end subroutine PWGTO_overlap
  subroutine PWGTO_multipole(this, ibra, m, iket, X, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(in) :: ibra, iket
    integer ,intent(in) :: m
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    
    integer :: A, B, i, j, nA, nB
    complex(kind(0d0)) :: acc, tmp, hint(0:m, 0:m), cAcB
    integer mm
    
    ierr = 0
    call check_matrix(this, X, ierr); CHK_ERR(ierr)

    do A = 1, this % num
       do B = 1, this % num
          hint(:, :) = 0.0d0
          call hermite_1drm(this%gP(A,B), this%RP(A,B), m ,hint)
          acc = 0.0d0
          do i = 1, this%nlc(ibra,A)
             do j = 1, this%nlc(iket,B)
                nA = this%nns(ibra,A, i)
                nB = this%nns(iket,B, j)
                cAcB = conjg(this%cs(ibra,A,i)) * this%cs(iket,B,j)
                mm = min(nA+nB, m)
                tmp = dot_product(this%d(A,B, nA,nB,0:mm), hint(m,0:mm))
                acc = acc + cAcB * tmp 
             end do
          end do
          X(A,B) = this%eP(A,B) * acc
       end do
    end do
    
  end subroutine PWGTO_multipole
  subroutine PWGTO_kineticP2(this, ibra, iket, X, ierr)
    use Mod_const, only : II
    type(Obj_PWGTO) :: this
    integer, intent(in) :: ibra, iket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    double precision :: PA, PB 
    complex(kind(0d0)) acc, tmp, hint, cAcB, gA, gB
    
    ierr = 0
    call check_matrix(this, X, ierr); CHK_ERR(ierr)
    do A = 1, this % num
       do B = 1, this % num
          gA = conjg(this % gs(A))
          gB = this % gs(B)
          PA = this % Ps(A)
          PB = this % Ps(B)          
          call hermite_1dint(this%gP(A,B), hint)
          acc = 0.0d0
          do i = 1, this%nlc(ibra,A)
             do j = 1, this%nlc(iket,B)
                nA = this%nns(ibra,A, i)
                nB = this%nns(iket,B, j)
                cAcB = conjg(this%cs(ibra,A,i)) * this%cs(iket,B,j)

                tmp = 4*gA*gB * this%d(A,B, nA+1,nB+1,0) &
                     -2*ii*gA*PB * this%d(A,B, nA+1,nB,0) &
                     +2*ii*PA*gB * this%d(A,B, nA,nB+1,0) &
                     +(PB*PA) * this%d(A,B, nA,nB,0)
                if(nA.ne.0) then
                   tmp = tmp +nA*(ii*PB) * this%d(A,B, nA-1,nB,0)
                   tmp = tmp -2*nA*gB * this%d(A,B, nA-1,nB+1,0)
                end if
                if(nB.ne.0) then
                   tmp = tmp +(-ii*PA)*nB * this%d(A,B, nA,nB-1,0)
                   tmp = tmp -2*gA*nB * this%d(A,B, nA+1,nB-1,0)
                end if
                if(nB.ne.0.and.nA.ne.0) then
                   tmp = tmp +(nA*nB) * this%d(A,B, nA-1,nB-1,0)
                end if
                acc = acc + cAcB * tmp * this%eP(A,B)*hint
             end do
          end do
          X(A,B) = acc
       end do
    end do
    
  end subroutine PWGTO_kineticP2
  function PWGTO_nterm(this, A) result(res)
    ! gives normalization term
    type(Obj_PWGTO) :: this
    integer :: A
    complex(kind(0d0)) :: res
    if(this%typ(1).ne."0") then
       MSG_ERR("invalid condition")
       write(0,*) "stop program..."
       stop
    end if
    res = this%cs(1,A,1)
  end function PWGTO_nterm
  ! -- utils --
  subroutine check_matrix(this, M, ierr)
    type(Obj_PWGTO), intent(in) :: this
    complex(kind(0d0)), intent(in) :: M(:,:)
    integer, intent(out) :: ierr

    ierr = 0
    if(size(M,1) .ne. this%num .or. size(M,2) .ne. this%num) then
       MSG_ERR("invalid size")
       ierr = 1; return
    end if
    
  end subroutine check_matrix
  subroutine prod_gauss(gA, RA, PA, tA, gB, RB, PB, tB, gP, RP, eP)    
    ! compute gauss parameter of product conjg(GA).GB
    use Mod_const, only : II
    complex(kind(0d0)), intent(in) :: gA, gB
    double precision, intent(in)   :: RA, PA, tA, RB, PB, tB
    complex(kind(0d0)) :: gP, RP, eP
    complex(kind(0d0)) :: cgA
    cgA =conjg(gA) 
    gP = cgA + gB
    RP = (2*cgA*RA+2*gB*RB-II*PA+ii*PB) / (2*gP)
    eP = exp(-cgA*RA*RA - gB*RB*RB &
         +ii*RA*PA  - ii*PB*RB &
         -ii*tA     + ii*tB &
         +gP*RP**2)
  end subroutine prod_gauss
  recursive subroutine hermite_coef_d_0(gP,wPk,RAk,RBk,nAk,nBk,Nk, res)
    ! compute coefficient of Hermitian GTO exapanded for cartesian GTO
    complex(kind(0d0)), intent(in)  :: gP, wPk
    double precision, intent(in)    :: RAk, RBk
    integer, intent(in) :: nAk, nBk, Nk
    complex(kind(0d0)) res, res0, res1, res2

    if(nAk .eq. 0 .and. nBk .eq. 0 .and. Nk .eq. 0) then
       res = 1.0d0; return       
    end if

    if(Nk .lt. 0 .or. Nk .gt. nAk + nBk) then
       res = 0.0d0; return
    end if

    if(nAk > 0) then
       call hermite_coef_d_0(gP, wPk, RAk, RBk, nAk-1, nBk, Nk-1, res0)
       call hermite_coef_d_0(gP, wPk, RAk, RBk, nAk-1, nBk, Nk,   res1)
       call hermite_coef_d_0(gP, wPk, RAk, RBk, nAk-1, nBk, Nk+1, res2)
       res = & 
            +1/(2*gP) *res0&
            +(wPk-RAk)*res1 &
            +(Nk+1)   *res2
    else
       call hermite_coef_d_0(gP, wPk, RAk, RBk, nAk, nBk-1, Nk-1, res0)
       call hermite_coef_d_0(gP, wPk, RAk, RBk, nAk, nBk-1, Nk,   res1)
       call hermite_coef_d_0(gP, wPk, RAk, RBk, nAk, nBk-1, Nk+1, res2)
       res = &
            +1/(2*gP) * res0 &
            +(wPk-RBk)* res1 &
            +(Nk+1)   * res2
    end if
    
  end subroutine hermite_coef_d_0
  subroutine hermite_coef_d(gP, wP,RA,RB, maxnA,maxnB, res)
    ! batch calculation of hermitian coefficient
    complex(kind(0d0)), intent(in) :: gP, wP
    double precision, intent(in)   :: RA, RB
    integer, intent(in) :: maxnA, maxnB
    complex(kind(0d0)), intent(out) :: res(0:,0:,0:)
    integer :: maxNk
    complex(kind(0d0))  v
    
    integer nA, nB, Nk, nAB

    maxNk = maxnA + maxnB
    
    res(0,0,0) = 1.0d0
    
    do nA = 0, maxnA
       do nB = 0, maxnB
          do Nk = nA+nB+1, maxNk
             res(nA,nB,Nk) = 0.0d0
          end do
       end do
    end do

    do nAB = 1, maxnA+maxnB
       do nA = 0, min(maxnA, nAB)
          nB = nAB-nA
          if(nB .le. maxnB) then
             do Nk = 0, min(nAB, maxNk)
                v = 0.0d0
                if(nA > 0) then
                   v = v + (wP-RA) * res(nA-1,nB,Nk)
                   if(Nk-1 >= 0) then
                      v = v + 1/(2*gP) * res(nA-1,nB,Nk-1)
                   end if
                   if(Nk+1 <= maxNk) then
                      v = v + (Nk+1)  * res(nA-1,nB,Nk+1)
                   end if
                else if(nB > 0) then
                   v = v + (wP-RB) * res(nA,nB-1,Nk)
                   if(Nk-1 >= 0) then
                      v = v + 1/(2*gP) * res(nA,nB-1,Nk-1)
                   end if
                   if(Nk+1 <= maxNk) then
                      v = v + (Nk+1)  * res(nA,nB-1,Nk+1)
                   end if
                end if
                res(nA, nB, Nk) = v
             end do
          end if
       end do
    end do
    
   ! do nA = 0, maxnA
   !    do nB = 0, maxnB
   !       do Nk = 0, maxNk
   !          d(nA,nB,Nk) = coef_d(zP, wP,RA,RB, nA,nB,Nk)
   !       end do
   !    end do
   ! end do
    
  end subroutine hermite_coef_d
  subroutine hermite_1dint(gP, res)
    use Mod_const, only : pi
    complex(kind(0d0)), intent(in) :: gP
    complex(kind(0d0)), intent(out) :: res
    res = sqrt(pi/gP)
  end subroutine hermite_1dint
  recursive subroutine hermite_1drm_0(gP, wP, m, n, res, ierr)
    ! x H_n(x) = (2w)^{-1}H_{n+1}(x) + nH_{n-1}(x)
    ! S(m,n)   = x^m H_n
    !          = (2w)^{-1} x^{m-1}H_{n+1} + nH_x^{m-1}{n-1}
    !          = (2w)^{-1} S(m-1,n+1) + nS(m-1,n-1)
    complex(kind(0d0)), intent(in) :: gP, wP
    integer, intent(in) :: m, n    
    complex(kind(0d0)), intent(out) :: res
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: res0
    
    ierr = 0
    
    if(n < 0) then
       MSG_ERR("n is negative")
       ierr = 1
       return
    end if
    if(m < 0) then
       MSG_ERR("m is negative")
       ierr = 1
       return 
    end if

    if(n>m) then
       res = 0.0d0
       return
    end if

    if(n.eq.0 .and. m.eq.0) then
       call hermite_1dint(gP, res)
       return
    end if

    call hermite_1drm_0(gP, wP, m-1, n+1, res0, ierr); CHK_ERR(ierr)
    res = 1/(2*gP) * res0
    call hermite_1drm_0(gP, wP, m-1, n  , res0, ierr); CHK_ERR(ierr)
    res = res + wP * res0
    if(n.ne.0) then
       call hermite_1drm_0(gP, wP, m-1, n-1, res0, ierr); CHK_ERR(ierr)
       res = res + n * res0
    end if
    
  end subroutine hermite_1drm_0
  subroutine hermite_1drm(gP, wP, mx, res)
    ! batch calculation of hermite 1drm
    complex(kind(0d0)), intent(in) :: gP, wP
    integer, intent(in) :: mx
    complex(kind(0d0)), intent(out) :: res(0:,0:)
    integer :: m, Nk, nk1

   do Nk = 0, mx
      do m = 0, mx
         if(Nk > m) then
            res(m, Nk) = 0.0d0
         end if
      end do
   end do
   call hermite_1dint(gP, res(0, 0))
   do m = 1, mx
      nk1 = min(m, mx)
      res(m, 0:nk1) = wP * res(m-1,0:nk1)
      do Nk = 0, nk1-1
         res(m, Nk) = res(m, Nk) + 1/(2*gP)*res(m-1,Nk+1)
      end do
      do Nk = 1, nk1
         res(m, Nk) = res(m, Nk) + Nk*res(m-1,Nk-1)
      end do
   end do

!    do m = 0, mx
!       do Nk = 0, mx
!          res(m,Nk) = calc_hermite_1drm(zP, wP, m, Nk)
!       end do
!    end do
    
  end subroutine hermite_1drm
end module Mod_PWGTO
