#include "macros.fpp"
module Mod_PWGTO
  ! Plane Wave Gauss Type Orbital
  !   G_i(Q) = Ni (Q-Ri)^ni exp[-g(Q-Ri)^2 + iPi(Q-Ri)]
  implicit none
  type Obj_ICs
     integer, allocatable :: nlc(:)
     integer, allocatable :: ns(:,:)
     complex(kind(0d0)), allocatable :: cs(:,:)
  end type Obj_ICs
  type Obj_PWGTO
     ! - size -
     integer :: num, maxnd
     ! - variable -
     integer, allocatable :: ns(:)
     complex(kind(0d0)), allocatable :: zs(:) 
     double precision,   allocatable :: Rs(:), Ps(:), gs(:)
     ! - intermediate -
     type(Obj_ICs) :: nc0, ncR, ncP
     complex(kind(0d0)), allocatable :: zP(:,:), RP(:,:), eP(:,:)
     complex(kind(0d0)), allocatable :: d(:,:,:,:,:)
  end type Obj_PWGTO
contains
  ! -- ICs --
  subroutine ICs_new(this, num, n, ierr)
    type(Obj_ICs) :: this
    integer, intent(in) :: num, n
    integer, intent(out) :: ierr
    ierr = 0
    allocate(this%nlc(num))
    allocate(this%ns(num, n))
    allocate(this%cs(num, n))
  end subroutine ICs_new
  ! -- main --
  subroutine PWGTO_new(this, num, maxnd, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(in) :: num, maxnd
    integer, intent(out) :: ierr
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
    
    call ICs_new(this%nc0, num, 1, ierr)
    call ICs_new(this%ncR, num, 3, ierr)
    call ICs_new(this%ncP, num, 1, ierr)
    
  end subroutine PWGTO_new
  subroutine PWGTO_setup(this, ierr)
    use Mod_const, only : II
    use Mod_math, only : gtoint
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    integer A, B, nA, n
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

       this % nc0 % nlc(A) = 1
       this % nc0 % cs(A, 1) = c0
       this % nc0 % ns(A, 1) = this % ns(A)

       if(nA .eq. 0) then
          this % ncR % nlc(A) = 2
       else
          this % ncR % nlc(A) = 3
       end if
       this % ncR % cs(A, 1) = c0*2*zA
       this % ncR % ns(A, 1) = nA+1
       this % ncR % cs(A, 2) = -c0*ii*this%Ps(A)
       this % ncR % ns(A, 2) = nA
       if(nA .ne. 0) then
          this % ncR % cs(A, 3) = nA*c0
          this % ncR % ns(A, 3) = nA-1
       end if

       this % ncP % nlc(A) = 1
       this % ncP % cs(A, 1) = c0*ii
       this % ncP % ns(A, 1) = nA+1
    end do

    ! -- combination of GTO --
    do B = 1, this % num
       do A = 1, this % num
          call prod_gauss(&
               this%zs(A), this%Rs(A), this%Ps(A), this%gs(A), &
               this%zs(B), this%Rs(B), this%Ps(B), this%gs(B), &
               this%zP(A,B), this%RP(A,B), this%eP(A,B))
          n = this % maxnd
          call hermite_coef(this%zP(A,B), this%RP(A,B), &
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
  end subroutine PWGTO_delete
  ! -- calc (private) --
  subroutine overlap(this, bra, ket, X, ierr)
    type(Obj_PWGTO), intent(in) :: this
    type(Obj_ICs), intent(in) ::   bra, ket
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
    type(Obj_ICs), intent(in) ::   bra, ket
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
  subroutine kinetic(this, bra, ket, X, ierr)
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
    type(Obj_ICs), intent(in) :: bra, ket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nA, nB
    double precision :: PA, PB 
    complex(kind(0d0)) acc, tmp, hint, cA, cB, zA, zB
    
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
          X(A,B) = 0.5d0*acc
       end do
    end do
    
  end subroutine kinetic
  ! -- calc (public) --
  subroutine PWGTO_overlap(this, dbra, dket, X, ierr)
    type(Obj_PWGTO), intent(in) :: this
    character, intent(in) :: dbra, dket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr

    ierr = 0
    if(size(X,1) .ne. this%num .or. size(X,2) .ne. this%num) then
       MSG_ERR("invalid size")
       ierr = 1; return
    end if

    if(dbra.eq.'0' .and. dket.eq.'0') then
       call overlap(this, this%nc0, this%nc0, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'0' .and. dket.eq.'R') then
       call overlap(this, this%nc0, this%ncR, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'0' .and. dket.eq.'P') then
       call overlap(this, this%nc0, this%ncP, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'R' .and. dket.eq.'R') then
       call overlap(this, this%ncR, this%ncR, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'R' .and. dket.eq.'P') then
       call overlap(this, this%ncR, this%ncP, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'P' .and. dket.eq.'P') then
       call overlap(this, this%ncP, this%ncP, X, ierr); CHK_ERR(ierr)
    else
       MSG_ERR("not implemented")
       ierr = 1; return
    end if
    
  end subroutine PWGTO_overlap
  subroutine PWGTO_multipole(this, dbra, m, dket, X, ierr)
    type(Obj_PWGTO) :: this
    character, intent(in) :: dbra, dket
    integer ,intent(in) :: m
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr
    ierr = 0
    if(size(X,1) .ne. this%num .or. size(X,2) .ne. this%num) then
       MSG_ERR("invalid size of matrix")
       ierr = 1; return
    end if

    if(dbra.eq.'0' .and. dket.eq.'0') then
       call multipole(this, this%nc0, m, this%nc0, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'0' .and. dket.eq.'R') then
       call multipole(this, this%nc0, m, this%ncR, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'0' .and. dket.eq.'P') then
       call multipole(this, this%nc0, m, this%ncP, X, ierr); CHK_ERR(ierr)
    else
       MSG_ERR("not impl")
       ierr = 1; return
    end if

  end subroutine PWGTO_multipole
  subroutine PWGTO_kinetic(this, dbra, dket, X, ierr)
    type(Obj_PWGTO) :: this
    character, intent(in) :: dbra, dket
    complex(kind(0d0)), intent(out) :: X(:,:)
    integer, intent(out) :: ierr

    if(size(X,1) .ne. this%num .or. size(X,2) .ne. this%num) then
       MSG_ERR("invalid size of matrix")
       ierr = 1; return 
    end if

    if(dbra.eq.'0' .and. dket.eq.'0') then
       call kinetic(this, this%nc0, this%nc0, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'0' .and. dket.eq.'R') then
       call kinetic(this, this%nc0, this%ncR, X, ierr); CHK_ERR(ierr)
    else if(dbra.eq.'0' .and. dket.eq.'P') then
       call kinetic(this, this%nc0, this%ncP, X, ierr); CHK_ERR(ierr)
    else
       MSG_ERR("not impl")
       ierr = 1; return
    end if
    
  end subroutine PWGTO_kinetic
  ! -- utils --
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
  recursive subroutine calc_coef_d(zP,wPk,RAk,RBk,nAk,nBk,Nk, res)
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
       call calc_coef_d(zP, wPk, RAk, RBk, nAk-1, nBk, Nk-1, res0)
       call calc_coef_d(zP, wPk, RAk, RBk, nAk-1, nBk, Nk,   res1)
       call calc_coef_d(zP, wPk, RAk, RBk, nAk-1, nBk, Nk+1, res2)
       res = & 
            +1/(2*zP) *res0&
            +(wPk-RAk)*res1 &
            +(Nk+1)   *res2
    else
       call calc_coef_d(zP, wPk, RAk, RBk, nAk, nBk-1, Nk-1, res0)
       call calc_coef_d(zP, wPk, RAk, RBk, nAk, nBk-1, Nk,   res1)
       call calc_coef_d(zP, wPk, RAk, RBk, nAk, nBk-1, Nk+1, res2)
       res = &
            +1/(2*zP) * res0 &
            +(wPk-RBk)* res1 &
            +(Nk+1)   * res2
    end if
    
  end subroutine calc_coef_d
  subroutine hermite_coef(zP, wP,RA,RB, maxnA,maxnB, res)
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
    
  end subroutine hermite_coef
  subroutine hermite_1dint(zP, res)
    use Mod_const, only : pi
    complex(kind(0d0)), intent(in) :: zP
    complex(kind(0d0)), intent(out) :: res
    res = sqrt(pi/zP)
  end subroutine hermite_1dint
  recursive subroutine calc_hermite_1drm(zP, wP, m, n, res, ierr)
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

    call calc_hermite_1drm(zP, wP, m-1, n+1, res0, ierr); CHK_ERR(ierr)
    res = 1/(2*zP) * res0
    call calc_hermite_1drm(zP, wP, m-1, n  , res0, ierr); CHK_ERR(ierr)
    res = res + wP * res0
    if(n.ne.0) then
       call calc_hermite_1drm(zP, wP, m-1, n-1, res0, ierr); CHK_ERR(ierr)
       res = res + n * res0
    end if
    
  end subroutine calc_hermite_1drm
  subroutine hermite_1drm(zP, wP, mx, res)
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
