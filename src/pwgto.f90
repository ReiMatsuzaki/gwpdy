#include "macros.fpp"
module Mod_PWGTO
  ! Plane Wave Gauss Type Orbital
  !   G^k_i(Q) = sum_j c_kji(Q-Ri)^{nkji} exp[-g(Q-Ri)^2 + iPi(Q-Ri)]
  implicit none
  type Obj_NCs
     ! represent sum_j cj x^{nj}
     character(3) :: typ    ! 0: normal, dR: d.of R,  dP: d. of P
     integer, allocatable :: nlc(:)  ! nlc(A)  : number of j for GA
     integer, allocatable :: ns(:,:) ! ns(A,:) : {nj} for GA
     complex(kind(0d0)), allocatable :: cs(:,:) ! cs(A,:) : {cj} for GA
  end type Obj_NCs
  type Obj_PWGTO
     ! - size -
     integer :: num, maxnd
     ! - variable -
     integer, allocatable :: ns(:)
     complex(kind(0d0)), allocatable :: zs(:) 
     double precision,   allocatable :: Rs(:), Ps(:), gs(:)
     ! - intermediate -
     type(Obj_NCs), allocatable :: ncs(:) ! ncs(1): non, ncs(2):dR, ncs(3):dP
     complex(kind(0d0)), allocatable :: zP(:,:), RP(:,:), eP(:,:)
     complex(kind(0d0)), allocatable :: d(:,:,:,:,:) ! d(A,B,0:maxnd,0:maxnd,0:maxnd*2)
  end type Obj_PWGTO
contains
  ! -- NCs --
  subroutine NCs_new(this, num, numNC, ierr)
    type(Obj_NCs) :: this
    integer, intent(in) :: num, numNC
    integer, intent(out) :: ierr
    ierr = 0
    this%typ = ""
    allocate(this%nlc(num))
    allocate(this%ns(num, numNC))
    allocate(this%cs(num, numNC))
  end subroutine NCs_new
  ! -- main --
  subroutine PWGTO_new(this, num, numNCs, maxnd, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(in) :: num, numNCs, maxnd
    integer, intent(out) :: ierr
    integer icn

    ierr = 0
    this%num   = num
    this%maxnd = maxnd
    allocate(this%ns(num))
    allocate(this%zs(num), this%Rs(num), this%Ps(num), this%gs(num))
    allocate(this%zP(num,num), this%RP(num,num), this%eP(num,num))
    allocate(this%d(num,num,0:maxnd,0:maxnd,0:maxnd*2))
    this%ns = 0
    this%zs = (0.0d0, 0.0d0)
    this%Rs = 0.0d0
    this%Ps = 0.0d0
    this%gs = 0.0d0

    allocate(this%ncs(numNCs))
    do icn = 1, numNCs
       call NCs_new(this%ncs(icn), num, 3, ierr); CHK_ERR(ierr)
    end do
    this%ncs(1) %typ = "0"
    
  end subroutine PWGTO_new
  subroutine PWGTO_setup(this, ierr)
    use Mod_const, only : II
    use Mod_math, only : gtoint
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    integer A, B, nA, n, inc
    complex(kind(0d0)) :: zA, c0
    complex(kind(0d0)), allocatable :: gg(:)
    
    ierr = 0
    allocate(gg(0:maxval(this%ns(:))*2+2))

    ! -- normalization term / derivative basis --
    do A = 1, this % num
       nA = this%ns(A)
       zA = this%zs(A)
       call gtoint(2*nA, conjg(zA)+zA, gg, ierr); CHK_ERR(ierr)
       c0 = 1.0d0/sqrt(gg(2*nA))

       do inc = 1, size(this%ncs)
          select case(this%ncs(inc)%typ)
          case("0")
             this%ncs(inc)%nlc(A) = 1
             this%ncs(inc)%cs(A, 1) = c0
             this%ncs(inc)%ns(A, 1) = this % ns(A)
          case("dR")
             if(nA .eq. 0) then
                this%ncs(inc)%nlc(A) = 2
             else
                this%ncs(inc)%nlc(A) = 3
             end if
             this%ncs(inc)%cs(A, 1) = c0*2*zA
             this%ncs(inc)%ns(A, 1) = nA+1
             this%ncs(inc)%cs(A, 2) = -c0*ii*this%Ps(A)
             this%ncs(inc)%ns(A, 2) = nA
             if(nA .ne. 0) then
                this%ncs(inc) % cs(A, 3) = nA*c0
                this%ncs(inc) % ns(A, 3) = nA-1
             end if
          case("dP")
             this%ncs(inc)%nlc(A) = 1
             this%ncs(inc)%cs(A, 1) = c0*ii
             this%ncs(inc)%ns(A, 1) = nA+1
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
               this%zs(A), this%Rs(A), this%Ps(A), this%gs(A), &
               this%zs(B), this%Rs(B), this%Ps(B), this%gs(B), &
               this%zP(A,B), this%RP(A,B), this%eP(A,B))
          n = this % maxnd
          call hermite_coef_d(this%zP(A,B), this%RP(A,B), &
               this%Rs(A), this%Rs(B), &
               n,n, this%d(A,B,:,:,:))
       end do
    end do
    
  end subroutine PWGTO_setup
  subroutine PWGTO_delete(this, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(this%zs, this%Rs, this%Ps, this%gs)
    deallocate(this%ncs, this%zP, this%RP, this%eP, this%d)
  end subroutine PWGTO_delete
  ! -- calc (private) --
  subroutine overlap(this, bra, ket, X, ierr)
    type(Obj_PWGTO), intent(in) :: this
    type(Obj_NCs), intent(in) ::   bra, ket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    complex(kind(0d0)) acc, tmp, hint

    ierr = 0
    do A = 1, this % num
       do B = 1, this % num
          call hermite_1dint(this%zP(A,B), hint)
          acc = 0.0d0
          do i = 1, bra%nlc(A)
             do j = 1, ket % nlc(B)
                nA = bra % ns(A, i)
                nB = ket % ns(B, j)
                tmp = this%eP(A,B) * this%d(A,B, nA,nB,0) * hint
                acc = acc + conjg(bra%cs(A,i)) * ket%cs(B,j) * tmp
             end do
          end do
          X(A,B) = acc
       end do
    end do
    
  end subroutine overlap
  subroutine multipole(this, bra, m, ket, X, ierr)
    type(Obj_PWGTO), intent(in) :: this
    type(Obj_NCs), intent(in) ::   bra, ket
    integer, intent(in) :: m
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    complex(kind(0d0)) :: acc, tmp, hint(0:m, 0:m)
    integer mm
    ierr = 0
    do A = 1, this % num
       do B = 1, this % num
          hint(:, :) = 0.0d0
          call hermite_1drm(this%zP(A,B), this%RP(A,B), m ,hint)
          acc = 0.0d0
          do i = 1, bra%nlc(A)
             do j = 1, ket % nlc(B)
                nA = bra % ns(A, i)
                nB = ket % ns(B, j)
                mm = min(nA+nB, m)
                tmp = dot_product(this%d(A,B, nA,nB,0:mm), hint(m,0:mm))
                acc = acc + conjg(bra%cs(A,i)) * ket%cs(B,j) * tmp 
             end do
          end do
          X(A,B) = this%eP(A,B) * acc
       end do
    end do
    
  end subroutine multipole
  subroutine kineticP2(this, bra, ket, X, ierr)
    ! d/dr x^n Exp[-zx^2 + iPx]
    !    = {nx^{n-1} +iPx^n -2zx^{n+1}} Exp[]
    ! <gA|T|gB> = nAnB <gA|x^{-2}|gB>
    !            +(i.nA.PB - i.nB.PA) <x^{-1}>
    !            +(-2(nA.zB+nB.zA) <1>
    !            +i(-PA+PB)<1>
    !            +i(-2PA.zB + 2PB.zA) <x>
    !            +4.zA.zB <x^2>
    use Mod_const, only : ii
    type(Obj_PWGTO), intent(in) :: this
    type(Obj_NCs), intent(in) :: bra, ket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    double precision :: PA, PB 
    complex(kind(0d0)) acc, tmp, hint, cA, cB, zA, zB
    ierr = 0
    do A = 1, this % num
       do B = 1, this % num
          zA = conjg(this % zs(A))
          zB = this % zs(B)
          PA = this % Ps(A)
          PB = this % Ps(B)          
          call hermite_1dint(this%zP(A,B), hint)
          acc = 0.0d0
          do i = 1, bra % nlc(A)
             do j = 1, ket % nlc(B)
                cA = conjg(bra % cs(A, i))
                cB = ket % cs(B, j)
                nA = bra % ns(A, i)
                nB = ket % ns(B, j)
                tmp = 4*zA*zB * this%d(A,B, nA+1,nB+1,0) &
                     -2*ii*zA*PB * this%d(A,B, nA+1,nB,0) &
                     +2*ii*PA*zB * this%d(A,B, nA,nB+1,0) &
                     +(PB*PA) * this%d(A,B, nA,nB,0)
                if(nA.ne.0) then
                   tmp = tmp +nA*(ii*PB) * this%d(A,B, nA-1,nB,0)
                   tmp = tmp -2*nA*zB * this%d(A,B, nA-1,nB+1,0)
                end if
                if(nB.ne.0) then
                   tmp = tmp +(-ii*PA)*nB * this%d(A,B, nA,nB-1,0)
                   tmp = tmp -2*zA*nB * this%d(A,B, nA+1,nB-1,0)
                end if
                if(nB.ne.0.and.nA.ne.0) then
                   tmp = tmp +(nA*nB) * this%d(A,B, nA-1,nB-1,0)
                end if
                acc = acc + tmp * this%eP(A,B)*hint*cA*cB
             end do
          end do
          X(A,B) = acc
       end do
    end do
    
  end subroutine kineticP2
  ! -- calc (public) --
  subroutine PWGTO_overlap(this, ibra, iket, X, ierr)
    type(Obj_PWGTO), intent(in) :: this
    integer, intent(in) :: ibra, iket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    complex(kind(0d0)) acc, tmp, hint

    ierr = 0
    call check_matrix(this, X, ierr); CHK_ERR(ierr)

    do A = 1, this % num
       do B = 1, this % num
          call hermite_1dint(this%zP(A,B), hint)
          acc = 0.0d0
          do i = 1, this%ncs(ibra)%nlc(A)
             do j = 1, this%ncs(iket)%nlc(B)
                nA = this%ncs(ibra) % ns(A, i)
                nB = this%ncs(iket) % ns(B, j)
                tmp = this%eP(A,B) * this%d(A,B, nA,nB,0) * hint
                acc = acc + conjg(this%ncs(ibra)%cs(A,I)) * this%ncs(iket)%cs(B,J) * tmp
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
    complex(kind(0d0)) :: acc, tmp, hint(0:m, 0:m)
    integer mm
    
    ierr = 0
    call check_matrix(this, X, ierr); CHK_ERR(ierr)

    do A = 1, this % num
       do B = 1, this % num
          hint(:, :) = 0.0d0
          call hermite_1drm(this%zP(A,B), this%RP(A,B), m ,hint)
          acc = 0.0d0
          do i = 1, this%ncs(ibra)%nlc(A)
             do j = 1, this%ncs(iket)%nlc(B)
                nA = this%ncs(ibra)%ns(A,i) !bra % ns(A, i)
                nB = this%ncs(iket)%ns(B,j)
                mm = min(nA+nB, m)
                tmp = dot_product(this%d(A,B, nA,nB,0:mm), hint(m,0:mm))
                acc = acc + conjg(this%ncs(ibra)%cs(A,i)) * this%ncs(iket)%cs(B,j) * tmp 
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
    complex(kind(0d0)) acc, tmp, hint, cA, cB, zA, zB
    
    ierr = 0
    call check_matrix(this, X, ierr); CHK_ERR(ierr)
    do A = 1, this % num
       do B = 1, this % num
          zA = conjg(this % zs(A))
          zB = this % zs(B)
          PA = this % Ps(A)
          PB = this % Ps(B)          
          call hermite_1dint(this%zP(A,B), hint)
          acc = 0.0d0
          do i = 1, this%ncs(ibra)%nlc(A) 
             do j = 1, this%ncs(iket)%nlc(B)
                nA = this%ncs(ibra)%ns(A,i) 
                nB = this%ncs(iket)%ns(B,j)
                cA = conjg(this%ncs(ibra)%cs(A,i))
                cB = this%ncs(iket)%cs(B,j)
                tmp = 4*zA*zB * this%d(A,B, nA+1,nB+1,0) &
                     -2*ii*zA*PB * this%d(A,B, nA+1,nB,0) &
                     +2*ii*PA*zB * this%d(A,B, nA,nB+1,0) &
                     +(PB*PA) * this%d(A,B, nA,nB,0)
                if(nA.ne.0) then
                   tmp = tmp +nA*(ii*PB) * this%d(A,B, nA-1,nB,0)
                   tmp = tmp -2*nA*zB * this%d(A,B, nA-1,nB+1,0)
                end if
                if(nB.ne.0) then
                   tmp = tmp +(-ii*PA)*nB * this%d(A,B, nA,nB-1,0)
                   tmp = tmp -2*zA*nB * this%d(A,B, nA+1,nB-1,0)
                end if
                if(nB.ne.0.and.nA.ne.0) then
                   tmp = tmp +(nA*nB) * this%d(A,B, nA-1,nB-1,0)
                end if
                acc = acc + tmp * this%eP(A,B)*hint*cA*cB
             end do
          end do
          X(A,B) = acc
       end do
    end do
    
  end subroutine PWGTO_kineticP2
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
  function PWGTO_nterm(this, A) result(res)
    ! gives normalization term
    type(Obj_PWGTO) :: this
    integer :: A
    complex(kind(0d0)) :: res
    res = this%ncs(1)%cs(A,1)
  end function PWGTO_nterm
  subroutine get_iNC(dchar, res, ierr)
    ! get index Obj_PWGTO.ncs which represent dchar
    character, intent(in) :: dchar
    integer, intent(out) :: res, ierr
    ierr = 0
    if('0'.eq.dchar) then
       res = 1
    else if('R'.eq.dchar) then
       res = 2
    else if('P'.eq.dchar) then
       res = 3
    else
       MSG_ERR("not impl")
       ierr = 1
       return
    end if
  end subroutine get_iNC
  subroutine prod_gauss(zA, RA, PA, gA, zB, RB, PB, gB, zP, RP, eP)    
    ! compute gauss parameter of product conjg(GA).GB
    use Mod_const, only : II
    complex(kind(0d0)), intent(in) :: zA, zB
    double precision, intent(in)   :: RA, PA, gA, RB, PB, gB
    complex(kind(0d0)) :: zP, RP, eP
    complex(kind(0d0)) :: czA
    czA =conjg(zA) 
    zP = czA + zB
    RP = (2*czA*RA+2*zB*RB-II*PA+ii*PB) / (2*zP)
    eP = exp(-czA*RA*RA - zB*RB*RB &
         +ii*RA*PA  - ii*PB*RB &
         -ii*gA     + ii*gB &
         +zP*RP**2)
  end subroutine prod_gauss
  recursive subroutine hermite_coef_d_0(zP,wPk,RAk,RBk,nAk,nBk,Nk, res)
    ! compute coefficient of Hermitian GTO exapanded for cartesian GTO
    complex(kind(0d0)), intent(in)  :: zP, wPk
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
       call hermite_coef_d_0(zP, wPk, RAk, RBk, nAk-1, nBk, Nk-1, res0)
       call hermite_coef_d_0(zP, wPk, RAk, RBk, nAk-1, nBk, Nk,   res1)
       call hermite_coef_d_0(zP, wPk, RAk, RBk, nAk-1, nBk, Nk+1, res2)
       res = & 
            +1/(2*zP) *res0&
            +(wPk-RAk)*res1 &
            +(Nk+1)   *res2
    else
       call hermite_coef_d_0(zP, wPk, RAk, RBk, nAk, nBk-1, Nk-1, res0)
       call hermite_coef_d_0(zP, wPk, RAk, RBk, nAk, nBk-1, Nk,   res1)
       call hermite_coef_d_0(zP, wPk, RAk, RBk, nAk, nBk-1, Nk+1, res2)
       res = &
            +1/(2*zP) * res0 &
            +(wPk-RBk)* res1 &
            +(Nk+1)   * res2
    end if
    
  end subroutine hermite_coef_d_0
  subroutine hermite_coef_d(zP, wP,RA,RB, maxnA,maxnB, res)
    ! batch calculation of hermitian coefficient
    complex(kind(0d0)), intent(in) :: zP, wP
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
                      v = v + 1/(2*zP) * res(nA-1,nB,Nk-1)
                   end if
                   if(Nk+1 <= maxNk) then
                      v = v + (Nk+1)  * res(nA-1,nB,Nk+1)
                   end if
                else if(nB > 0) then
                   v = v + (wP-RB) * res(nA,nB-1,Nk)
                   if(Nk-1 >= 0) then
                      v = v + 1/(2*zP) * res(nA,nB-1,Nk-1)
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
  subroutine hermite_1dint(zP, res)
    use Mod_const, only : pi
    complex(kind(0d0)), intent(in) :: zP
    complex(kind(0d0)), intent(out) :: res
    res = sqrt(pi/zP)
  end subroutine hermite_1dint
  recursive subroutine hermite_1drm_0(zP, wP, m, n, res, ierr)
    ! x H_n(x) = (2w)^{-1}H_{n+1}(x) + nH_{n-1}(x)
    ! S(m,n)   = x^m H_n
    !          = (2w)^{-1} x^{m-1}H_{n+1} + nH_x^{m-1}{n-1}
    !          = (2w)^{-1} S(m-1,n+1) + nS(m-1,n-1)
    complex(kind(0d0)), intent(in) :: zP, wP
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
       call hermite_1dint(zP, res)
       return
    end if

    call hermite_1drm_0(zP, wP, m-1, n+1, res0, ierr); CHK_ERR(ierr)
    res = 1/(2*zP) * res0
    call hermite_1drm_0(zP, wP, m-1, n  , res0, ierr); CHK_ERR(ierr)
    res = res + wP * res0
    if(n.ne.0) then
       call hermite_1drm_0(zP, wP, m-1, n-1, res0, ierr); CHK_ERR(ierr)
       res = res + n * res0
    end if
    
  end subroutine hermite_1drm_0
  subroutine hermite_1drm(zP, wP, mx, res)
    ! batch calculation of hermite 1drm
    complex(kind(0d0)), intent(in) :: zP, wP
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
   call hermite_1dint(zP, res(0, 0))
   do m = 1, mx
      nk1 = min(m, mx)
      res(m, 0:nk1) = wP * res(m-1,0:nk1)
      do Nk = 0, nk1-1
         res(m, Nk) = res(m, Nk) + 1/(2*zP)*res(m-1,Nk+1)
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
