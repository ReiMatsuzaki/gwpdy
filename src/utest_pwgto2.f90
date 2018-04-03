#include "macros.fpp"
#include "macros_utest.fpp"

module Mod_UTestPWGTO
  use Mod_UTest
  use Mod_PWGTO2
  implicit none
contains
  subroutine UTestPWGTO_run
    use Mod_Timer
    type(Obj_Timer) :: timer
    integer :: ierr
    
    call Timer_new(timer, "UtestPWGTO", .true., ierr); CHK_ERR(ierr)
    write(*,*) 
    write(*,*) "UTestGWP begin"
    write(*,*)

    call Timer_begin(timer, "run", ierr)
    call test_run
    call Timer_end(timer, "run", ierr)

    call Timer_begin(timer, "coef_d", ierr)
    call test_coef_d
    call Timer_end(timer, "coef_d", ierr)

    call Timer_begin(timer, "prod", ierr)
    call test_prod
    call Timer_end(timer, "prod", ierr)    

    call Timer_begin(timer, "overlap", ierr)
    call test_overlap(ierr)
    call Timer_end(timer, "overlap", ierr)

    call Timer_begin(timer, "multipole0", ierr)
    call test_multipole_0
    call Timer_end(timer, "multipole0", ierr)
    
    call Timer_begin(timer, "kinetic", ierr)
    call test_kinetic
    call Timer_end(timer, "kinetic", ierr)

    !call Timer_begin(timer, "compare", ierr)
    !    call test_compare_gwp
    !    call Timer_end(timer, "compare", ierr)            
    
  end subroutine UTestPWGTO_run
  subroutine test_run
    type(Obj_PWGTO) :: gwp
    integer, parameter :: num = 3
    integer, parameter :: nf = 1
    integer, parameter :: maxnd = 3
    integer, parameter :: numops = 1
    integer ierr
    call PWGTO_new(gwp, num, nf, maxnd, numops, ierr); CHK_ERR(ierr)
    call PWGTO_setup(gwp, ierr);         CHK_ERR(ierr)
    call PWGTO_delete(gwp, ierr);        CHK_ERR(ierr)
  end subroutine test_run
  subroutine test_coef_d
    integer, parameter :: maxnA = 1
    integer, parameter :: maxnB = 0
    integer, parameter :: maxm = 3
    complex(kind(0d0)) :: zP = (1.0d0, 0.0d0)
    complex(kind(0d0)) :: wP = (0.0d0, 0.0d0)
    double precision   :: RA = 0.0d0
    double precision   :: RB = 0.0d0
    complex(kind(0d0)) :: d(0:maxnA, 0:maxnB, 0:maxnA+maxnB)
    complex(kind(0d0)) :: ref
    integer nA, nB, Nk, ierr

    call hermite_coef_d(zP, wP, RA, RB, maxnA, maxnB, d)
    do nA = 0, maxnA
       do nB = 0, maxnB
          do Nk = 0, maxnA+maxnB
             call hermite_coef_d_0(zP, wP, RA, RB, nA, nB, Nk, ref)
             EXPECT_EQ_C(ref, d(nA, nB, Nk), ierr)
          end do
       end do
    end do
    
  end subroutine test_coef_d
  subroutine test_prod
    use Mod_const, only : II
    complex(kind(0d0)) :: zA, zB, zP
    double precision :: RA, RB, PA, PB, gA, gB
    complex(kind(0d0)) :: RP, eP
    complex(kind(0d0)) ref, calc
    double precision :: R
    integer ierr
    zA = (1.1d0, 0.3d0)
    zB = (1.2d0, -0.5d0)
    RA = 0.1d0
    gA = 0.3d0
    RB = 0.2d0
    PA = 0.2d0
    PB = 0.4d0
    gB = 0.4d0
    call prod_gauss(zA, RA, PA, gA,     zB, RB, PB, gB,    zP, RP, eP)
    R = 0.4d0
    ref = -conjg(zA)*(R-RA)**2 - zB*(R-RB)**2 &
         -ii*PA*(R-RA) +ii*PB*(R-RB) &
         +ii*(-gA+gB)
    calc = -zP*(R-RP)**2 + log(eP)
    EXPECT_EQ_C(ref, calc, ierr)
  end subroutine test_prod
  subroutine test_overlap(ierr)
    use Mod_math, only : gtoint
    type(Obj_PWGTO) :: g
    integer, parameter :: n = 4
    integer, parameter :: numops = 1
    integer, parameter :: maxnd = 5
    complex(kind(0d0)) :: calc, ref
    complex(kind(0d0)) :: S(n,n)
    complex(kind(0d0)) :: gg(0:20)
    integer A, B
    integer ierr
    complex(kind(0d0)) gP, RP, eP
    
    call PWGTO_new(g, n, 1, maxnd, numops, ierr); CHK_ERR(ierr)
    g%ns(1)=2; g%gs(1)=1.2d0;  g%Ps(1) = 10.0d0; g%Ps(1) = 1.0d0
    g%ns(2)=3; g%gs(2)=(0.9d0, -0.8d0); 
    g%ns(3)=4; g%gs(3)=1.0d0; g%Rs(3) = 0.1d0;
    g%ns(4)=2; g%gs(4)=1.1d0; g%Rs(4) = 0.1d0;
    call PWGTO_setup(g, ierr); CHK_ERR(ierr)
    
    call PWGTO_overlap(g, 1, 1, S, ierr); CHK_ERR(ierr)    

    A = 1
    B = 1
    calc = S(A, A) / (abs(PWGTO_nterm(g,A))**2)
    call gtoint(10, conjg(g%gs(A))+g%gs(B), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 2
    B = 2
    calc = S(A, A) / (abs(PWGTO_nterm(g,A))**2)
    call gtoint(10, conjg(g%gs(A))+g%gs(B), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 3
    B = 3
    calc = S(A, B)  / (abs(PWGTO_nterm(g,A))**2)
    call prod_gauss(&
         g%gs(A), g%Rs(A), g%Ps(A), g%thetas(A), &
         g%gs(B), g%Rs(B), g%Ps(B), g%thetas(B), &
         gP, RP, eP)
    call gtoint(10, conjg(g%gs(A))+g%gs(B), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 3
    B = 4
    calc = S(A, B) / (conjg(PWGTO_nterm(g,A))*PWGTO_nterm(g,B))
    call gtoint(10, conjg(g%gs(A))+g%gs(B), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)+g%ns(B))
    EXPECT_EQ_C(ref, calc, ierr)
    
  end subroutine test_overlap
  subroutine test_multipole_0
    use Mod_Math, only : gtoint
    type(Obj_PWGTO) :: g
    integer, parameter :: n = 5
    integer, parameter :: maxnd = 4
    integer, parameter :: numops = 3
    complex(kind(0d0)) :: M1(n, n)
    complex(kind(0d0)) :: ref
    integer A, B, ierr
    complex(kind(0d0)) :: gg(0:20)
    integer, parameter :: op0=1, op1=2, op2=3
    
    call PWGTO_new(g, n, 1, maxnd, numops, ierr)
    g%ns(1)=0; g%gs(1)=1.1d0
    g%ns(2)=0; g%gs(2)=1.1d0; g%Rs(2) = 0.1d0
    g%ns(3)=2; g%gs(3)=1.1d0
    g%ns(4)=2; g%gs(4)=1.1d0
    g%ops_typ(op0) = "0"
    g%ops_typ(op1) = "1"
    g%ops_typ(op2) = "2"
    call PWGTO_setup(g, ierr)
       
    A = 1; B = A
    call gtoint(10, g%gs(A)+g%gs(B), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)+g%ns(B)) *abs(PWGTO_nterm(g,A))**2
    call PWGTO_overlap(g, op0, op0, M1, ierr)    
    EXPECT_EQ_C(ref, M1(A,B), ierr)

    ref = gg(g%ns(A)+g%ns(B)+1)*abs(PWGTO_nterm(g,A))**2
    call PWGTO_overlap(g, op0, op1, M1, ierr)
    EXPECT_EQ_C(ref, M1(A,B), ierr)
    call PWGTO_overlap(g, op1, op0, M1, ierr)    
    EXPECT_EQ_C(ref, M1(A,B), ierr)

    ref = gg(g%ns(A)+g%ns(B)+2)*abs(PWGTO_nterm(g,A))**2
    call PWGTO_overlap(g, op1, op1, M1, ierr)    
    EXPECT_EQ_C(ref, M1(A,A), ierr)
    call PWGTO_overlap(g, op0, op2, M1, ierr)    
    EXPECT_EQ_C(ref, M1(A,B), ierr)
    call PWGTO_overlap(g, op2, op0, M1, ierr)    
    EXPECT_EQ_C(ref, M1(A,B), ierr)
 
    A = 2
    call gtoint(20, g%gs(A)+g%gs(B)+1, gg(:), ierr); CHK_ERR(ierr)
    call PWGTO_overlap(g, op0, op1, M1, ierr)
    ref = g%Rs(A) * gg(g%ns(A)*2) * (abs(PWGTO_nterm(g,A))**2)
    EXPECT_EQ_C(ref, M1(A,B), ierr)

  end subroutine test_multipole_0
  subroutine test_kinetic
    use Mod_math, only : gtoint
    ! d/dx Exp[-zx^2]  = -2zx Exp[]
    ! d/dx xExp[-zx^2] = Exp[]   -2zx^2Exp[]
    ! <d/dx(Exp[])|d/dx(Exp[])> =  4zz<2>
    ! <d/dx(Exp[])|d/dx(xExp[])> = 0

    ! d/dx Exp[-zx^2 + ipx]  = -2zx Exp[] +ip Exp[]
    ! <d/dx(Exp[]) | d/dx(Exp[])> = 4zz<2> + p^2<0>

    ! d/dx xExp[-zx^2 + ipx]  = Exp[] -2z(x2) Exp[] +ip xExp[]
    ! <d/dx(xExp[]) | d/dx(xExp[])> = <0> +4zz(4) -4z<2> + p2<2>
    ! <d/dx(Exp[])  | d/dx(xExp[])> = -ip<0> +2ipz<2> -2izp<2>
    
    type(Obj_PWGTO) :: g    
    integer, parameter :: n=4, maxnd=3, numops=2
    complex(kind(0d0)) :: calc, ref, z, zz, cz
    complex(kind(0d0)) :: T(n, n)
    double precision :: pA, pB
    integer ierr
    complex(kind(0d0)) :: gg(0:20)
    integer, parameter :: op0=1, opP2=2

    z = (1.0d0, 0.2d0)
    cz = conjg(z)
    zz = z + cz
    pA = 0.1d0
    pB = 0.2d0

    call PWGTO_new(g, n, 1, maxnd, numops, ierr)
    g % ns(1) = 0; g % gs(1) = z; g % Rs(1) = 0.1d0
    g % ns(2) = 1; g % gs(2) = z; g % Rs(2) = 0.1d0
    g % ns(3) = 0; g % gs(3) = z; g % Ps(3) = pA;    g % Rs(3) = 0.2d0
    g % ns(4) = 1; g % gs(4) = z; g % Ps(4) = pB;    g % Rs(4) = 0.2d0
    g % ops_typ(op0) = "0"
    g % ops_typ(opP2) = "P2"
    call PWGTO_setup(g, ierr); CHK_ERR(ierr)

    call PWGTO_overlap(g, op0, opP2, T, ierr); CHK_ERR(ierr)
    T = T/2

    call gtoint(4, zz, gg, ierr); CHK_ERR(ierr)
    ref = -0.5d0 * (-2*z*gg(0) + 4*z*z*gg(2))
    calc = T(1, 1) / (abs(PWGTO_nterm(g, 1))**2)
    EXPECT_EQ_C(ref, calc, ierr)

    ref = 0.5d0 * (gg(0) +4*cz*z*gg(4) - 2*zz*gg(2))
    calc = T(2, 2) / (abs(PWGTO_nterm(g,2))**2)
    EXPECT_EQ_C(ref, calc, ierr)

    EXPECT_EQ_C(T(3,1), conjg(T(1,3)), ierr)
    EXPECT_EQ_C((0.0d0,0.0d0), T(2,1), ierr)
    EXPECT_EQ_C((0.0d0,0.0d0), T(1,2), ierr)
    
    
    ref = 0.5d0 * (4*z*cz*gg(2) + pA*pA*gg(0))
    calc = T(3, 3) / (abs(PWGTO_nterm(g,3))**2)
    EXPECT_EQ_C(ref, calc, ierr)

    ref = 0.5d0 * (gg(0) + 4*z*cz*gg(4) + (pB*pB-2*z-2*cz)*gg(2))
    calc = T(4, 4) / (abs(PWGTO_nterm(g,4))**2)
    EXPECT_EQ_C(ref, calc, ierr)
    
    call PWGTO_delete(g, ierr); CHK_ERR(ierr)
    
  end subroutine test_kinetic
end module Mod_UTestPWGTO

program main
  use Mod_UTestPWGTO
  call UTestPWGTO_run
end program main
