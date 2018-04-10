! gwpdy/src/pwgto2.f90
!
! Plane Wave Gauss Type Orbital
!   GA(Q) = CA (Q-RA)^nA Exp[-gA(Q-RA)^2 + iPA(Q-RA)].
! The module also express the operators OPi operating G_A.
!   OPi.GA = sum_j CjiA (Q-RA)^{njiA} Exp[-gA(Q-RA)^2 + i PA(Q-RA)]
! Following matrix can be calcualted
!   { <OPi.GA|OPj.GB> | A,B=1,num}.
! Matrix elements are avaluated using ??? method. see Tachikawa(2008).
! Hermitian polynomial \Lambda(Q-R;gP) is defined as    
!   \Lambda_{A,B,N}(Q-R) = Exp[gP(Q-R)^2] (d/dRAB)^N Exp[-gp(Q-RAB)^2]
!   (Q-RA)^nA (Q-RB)^nB = sum_{N=0}^{nA+nB} d(nA,nB,N) \Lambda_{A,B,N}(Q-RAB)
! where
!   RAB = (gA.RA+gB.RB)/(gA+gB)
! matrix elements of f are evaluated as
!      <gA|f|gB> = Int[ eAB f(Q) (Q-RA)^nA (Q-RB)^nB Exp[-gAB (Q-RAB)^2]
!                = eAB \sum_{N=0}^{nA+nB} d_N^{nA,nB} Int[f(Q)\Lambda_{A,B,N}(Q_RAB)Exp[-gAB (Q-RAB)^2]
!                = eAB \sum_{N=0}^{nA+nB} [f]_N
!      [f]_N     = d_N^{nA,nB} (d/dRAB)^N Int[f(Q)Exp[-gAB (Q-RAB)^2]]
! hermite_coef_d : gives d(nA,nB,N)
! hermite_1drm   : gives [r^m]_N
!
! see 2018/3/30/qpbranch
#include "macros.fpp"  
module Mod_PWGTO2
  implicit none
  type Obj_PWGTO
     ! - data size -
     integer :: num, nf, maxnd, numops
     ! - variable -
     integer, allocatable :: ns(:)
     complex(kind(0d0)), allocatable :: gs(:) 
     double precision,   allocatable :: Rs(:), Ps(:), thetas(:)
     character(3), allocatable :: ops_typ(:)  ! dR, dP, P, P2
     ! - operator data -
     integer, allocatable :: ops_num(:,:)  ! op_num(opi,A)
     integer, allocatable :: ops_ns(:,:,:) ! op_ns( opi,A,:)
     complex(kind(0d0)), allocatable :: ops_cs(:,:,:) ! cs(opi,A,:)
     ! - intermediate -
     complex(kind(0d0)), allocatable :: gAB(:,:), RAB(:,:), eAB(:,:), hAB(:,:)
     complex(kind(0d0)), allocatable :: d(:,:,:,:,:) ! d(A,B,nA,nB,N)
  end type Obj_PWGTO
contains
  ! -- main --
  subroutine PWGTO_new(this, num, nf, maxnd, numops, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(in) :: num, nf, numops, maxnd
    integer, intent(out) :: ierr

    ierr = 0    
    if(nf.ne.1) then
       MSG_ERR("Unsupported")
       ierr = 1; return
    end if
    
    this%num= num
    this%nf = nf
    this%maxnd  = maxnd
    this%numops = numops
    
    allocate(this%ns(num))
    allocate(this%gs(num), this%Rs(num), this%Ps(num), this%thetas(num))
    this%ns = 0
    this%gs = (1.0d0, 0.0d0)
    this%Rs = 0.0d0
    this%Ps = 0.0d0
    this%thetas = 0.0d0
    
    allocate(this%ops_typ(numops))
    allocate(this%ops_num(numops, num))
    allocate(this%ops_ns(numops,  num, 5))
    allocate(this%ops_cs(numops,  num, 5))
    this%ops_typ = "0"

    allocate(this%gAB(num,num), this%RAB(num,num), this%eAB(num,num))
    allocate(this%hAB(num,num))
    allocate(this%d(num,num,0:maxnd,0:maxnd,0:2*maxnd))
    
  end subroutine PWGTO_new
  subroutine PWGTO_setup(this, ierr)
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    ierr = 0
    call setup_op(this, ierr); CHK_ERR(ierr)
    call setup_comb(this, ierr); CHK_ERR(ierr)
  end subroutine PWGTO_setup
  subroutine setup_op(this, ierr)
    use Mod_const, only : II
    use Mod_math, only : gtoint
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    integer :: A, nA, iop
    complex(kind(0d0)) :: gA, nt
    complex(kind(0d0)), allocatable :: gg(:)
    double precision tmp
    double precision pA

    ierr = 0
    
    ! -- normalization term and operation --
    allocate(gg(0:maxval(this%ns(:))*2+2))
    do A = 1, this%num
       nA = this%ns(A)
       gA = this%gs(A)
       pA = this%ps(A)
       call gtoint(2*nA, conjg(gA)+gA, gg, ierr); CHK_ERR(ierr)
       nt = 1.0d0/sqrt(gg(2*nA))
       do iop = 1, this%numops
          select case(this%ops_typ(iop))
          case("0")
             this%ops_num(iop,A)  = 1
             this%ops_ns(iop,A,1) = nA
             this%ops_cs(iop,A,1) = nt
          case("1")
             this%ops_num(iop,A) = 1
             this%ops_ns(iop,A,1) = nA+1
             this%ops_cs(iop,A,1) = nt
          case("2")
             this%ops_num(iop,A) = 1
             this%ops_ns(iop,A,1) = nA+2
             this%ops_cs(iop,A,1) = nt
          case("P1")
             if(nA .eq. 0) then
                this%ops_num(iop,A) = 2
             else
                this%ops_num(iop,A) = 3
             end if
             this%ops_ns(iop,A, 1) = nA+1
             this%ops_cs(iop,A, 1)  = II*2*gA*nt
             this%ops_ns(iop,A, 2) = nA
             this%ops_cs(iop,A, 2) = this%Ps(A)*nt
             if(nA .ne. 0) then                
                this%ops_ns(iop,A, 3) = nA-1
                this%ops_cs(iop,A, 3) = II*nA*nt
             end if
          case("P2")
             if(nA .eq. 0) then
                this%ops_num(iop,A) = 3
             else if(nA .eq. 1) then
                this%ops_num(iop,A) = 4
             else
                this%ops_num(iop,A) = 5
             end if
             this%ops_ns(iop,A, 1) = nA+2
             this%ops_cs(iop,A, 1) = -4*gA*gA * nt
             this%ops_ns(iop,A, 2) = nA+1
             this%ops_cs(iop,A, 2) = -(-4*II*gA*pA) * nt
             this%ops_ns(iop,A, 3) = nA
             this%ops_cs(iop,A, 3) = -(-2*gA*(nA+1) -2*gA*nA -pA*pA) * nt
             if(nA > 0) then                
                this%ops_ns(iop,A, 4) = nA-1
                this%ops_cs(iop,A, 4) = -(2*II*nA*pA) * nt
             end if
             if(nA > 1) then
                this%ops_ns(iop,A, 5) = nA-2
                this%ops_cs(iop,A, 5) = -nA*(nA-1) * nt
             end if
          case("dR")
             if(nA .eq. 0) then
                this%ops_num(iop,A) = 2
             else
                this%ops_num(iop,A) = 3
             end if
             this%ops_ns(iop,A, 1) = nA+1
             this%ops_cs(iop,A, 1)  = 2*gA*nt
             this%ops_ns(iop,A, 2) = nA
             this%ops_cs(iop,A, 2) = -II*this%Ps(A)*nt
             if(nA .ne. 0) then
                this%ops_cs(iop,A, 3) = -nA*nt
                this%ops_ns(iop,A, 3) = nA-1
             end if
          case("dP")
             this%ops_num(iop,A) = 1
             this%ops_cs(iop,A, 1)  = -II*nt
             this%ops_ns(iop,A, 1) = nA+1
          case("dgr")
             this%ops_num(iop,A)  = 2
             this%ops_cs(iop,A,1) = -nt
             this%ops_ns(iop,A,1) = nA+2
             call calc_nterm(1, nA, real(gA), tmp, ierr)
             this%ops_cs(iop,A,2) = tmp
             this%ops_ns(iop,A,2) = nA
          case("dgi")
             this%ops_num(iop,A)  = 1
             this%ops_cs(iop,A,1) = -II*nt
             this%ops_ns(iop,A,1) = nA+2             
          case default
             MSG_ERR("unsupported typ")
             ierr = 1; return
          end select
       end do
    end do
  end subroutine setup_op
  subroutine setup_comb(this, ierr)
    use Mod_const, only : II
    use Mod_math, only : gtoint
    type(Obj_PWGTO) :: this
    integer, intent(out) :: ierr
    integer :: A, B, n

    ierr = 0
    do B = 1, this % num
       do A = 1, this % num
          call prod_gauss(&
               this%gs(A), this%Rs(A), this%Ps(A), this%thetas(A), &
               this%gs(B), this%Rs(B), this%Ps(B), this%thetas(B), &
               this%gAB(A,B), this%RAB(A,B), this%eAB(A,B))
          n = this % maxnd
          call hermite_coef_d(this%gAB(A,B), this%RAB(A,B), &
               this%Rs(A), this%Rs(B), &
               n,n, this%d(A,B,:,:,:))
          call hermite_1dint(this%gAB(A,B), this%hAB(A,B))
       end do
    end do
    
  end subroutine setup_comb
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
    deallocate(this%ns, this%gs, this%Rs, this%Ps, this%thetas)
    deallocate(this%ops_typ)
    deallocate(this%ops_num, this%ops_ns, this%ops_cs)
    deallocate(this%gAB, this%RAB, this%eAB, this%hAB, this%d)
  end subroutine PWGTO_delete
  ! -- calc  --
  subroutine calc_nterm(ndiff, n, gr, res, ierr)
    ! give difference of normalization term
    !   (d/d(Re[g]))^(nd)N(t)
    ! where
    !   N(t) = 1/sqrt(S)
    !   S = Int_{-oo}^{+oo} x^{2n} exp[-2Re[g]x^2] dx.
    ! overlap S can be expressed by Gamma function
    !   S = (2Re[g])^{-n-1/2} Gamma[n+1/2]
    ! So, normalization term N becomes
    !   N = (2Re[g])^{n/2+1/4} Sqrt(1/Gamma[n+1/2])
    use Mod_const, only : PI
    use Mod_math, only : gamma_half_int
    integer, intent(in) :: ndiff, n
    double precision, intent(in) :: gr
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    double precision  ga(0:n)
    double precision nn, c
    
    ierr = 0
    call gamma_half_int(n, ga, ierr); CHK_ERR(ierr)
    nn = n*0.5d0+0.25d0
    c  = 2**nn * sqrt(1/ga(n))
    select case(ndiff)
    case(0)
       res = gr**nn * c
    case(1)
       res = nn*gr**(nn-1) * c
    case default
       MSG_ERR("invalid value: ndiff")
       ierr = 1; return
    end select
  end subroutine calc_nterm
  subroutine PWGTO_overlap(this,   ibra, iket, res, ierr)
    ! compute overlap.
    !
    ! Inputs
    ! ------
    !  ibra : integer : index of operator for bra
    !  iket : integer : index of operator for ket
    !
    ! Returns
    ! -------
    !  res : [[complex]] : resultant matrix
    type(Obj_PWGTO), intent(in) :: this
    integer, intent(in) :: ibra, iket
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    integer :: A, B, i, j, nAi, nBj
    complex(kind(0d0)) :: acc
    ierr = 0
    call check_matrix(this, res, ierr); CHK_ERR(ierr)
    do A = 1, this%num
       do B = 1, this%num
          acc = 0
          do i = 1, this%ops_num(ibra,A)
             do j = 1, this%ops_num(iket,B)
                nAi = this%ops_ns(ibra, A, i)
                nBj = this%ops_ns(iket, B, j)
                acc = acc + conjg(this%ops_cs(ibra,A,i))*this%ops_cs(iket,B,j) &
                     * this%d(A,B,nAi,nBj,0) 
             end do
          end do
          res(A,B) = acc * this%eAB(A,B) * this%hAB(A,B)
       end do       
    end do
  end subroutine PWGTO_overlap
  subroutine PWGTO_at(this, iop, cs, xs, res, ierr)
    use Mod_const, only : II
    type(Obj_PWGTO) :: this
    integer, intent(in) :: iop
    complex(kind(0d0)), intent(in) :: cs(:)
    double precision, intent(in) :: xs(:)
    complex(kind(0d0)), intent(out) :: res(:)
    integer, intent(out) :: ierr
    integer ix, i, A
    double precision d
    complex(kind(0d0)) tmp
    
    ierr = 0
    if(iop<1.or.this%numops<iop) then
       MSG_ERR("invalid iop")
       ierr = 1; return
    end if
    if(size(cs)<this%num) then
       MSG_ERR("size(cs)<num")
       ierr = 1; return
    end if
    if(size(xs).ne.size(res)) then
       MSG_ERR("size(xs)<size(res)")
       ierr = 1; return
    end if

    res(:) = 0
    do ix = 1, size(xs)
       tmp = 0
       do A = 1, this%num
          d = xs(ix)-this%Rs(A)
          do i = 1, this%ops_num(iop, A)
             tmp = tmp + cs(A)*this%ops_cs(iop,A,i) * d**this%ops_ns(iop,A,i) * &
                  exp(-this%gs(A)*d**2 + II*this%Ps(A)*d)
          end do
       end do
       if(abs(tmp) > 1.0d-14) then
          res(ix) = tmp
       end if
    end do
    
  end subroutine PWGTO_at
  ! -- utils --
  subroutine check_matrix(this, M, ierr)
    ! check the matrix size and raise error if it does not fit basis size.
    type(Obj_PWGTO), intent(in) :: this
    complex(kind(0d0)), intent(in) :: M(:,:)
    integer, intent(out) :: ierr

    ierr = 0
    if(size(M,1) .ne. this%num .or. size(M,2) .ne. this%num) then
       MSG_ERR("invalid size")
       ierr = 1; return
    end if
    
  end subroutine check_matrix
  function PWGTO_nterm(this,A) result(res)
    ! gives normalization term of A th basis.
    type(Obj_PWGTO) :: this
    integer, intent(in) :: A
    double precision :: res
    res = real(this%ops_cs(1,A,1))
  end function PWGTO_nterm
  ! -- calc matrix element --
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
    ! batch calculation of hermitian coefficient (see Tachikawa(2008))
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
end module Mod_PWGTO2
