#include "macros.fpp"
#include "macros_utest.fpp"

module UTestPWGTO
  use Mod_UTest
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

    call Timer_begin(timer, "multipole1", ierr)
    call test_multipole_1
    call Timer_end(timer, "multipole1", ierr)
    
    call Timer_begin(timer, "kinetic", ierr)
    call test_kinetic
    call Timer_end(timer, "kinetic", ierr)            
    
  end subroutine UTestPWGTO_run
  subroutine test_run
    use Mod_PWGTO
    type(Obj_PWGTO) :: gwp
    integer, parameter :: num = 3
    integer, parameter :: numNCs = 1
    integer, parameter :: maxnd = 3
    integer ierr
    call PWGTO_new(gwp, num, numNCs, maxnd, ierr); CHK_ERR(ierr)
    call PWGTO_setup(gwp, ierr);         CHK_ERR(ierr)
    call PWGTO_delete(gwp, ierr);        CHK_ERR(ierr)
  end subroutine test_run
  subroutine test_coef_d
    use Mod_PWGTO
    integer, parameter :: maxnA = 1
    integer, parameter :: maxnB = 0
    integer, parameter :: maxm = 3
    complex(kind(0d0)) :: zP = (1.0d0, 0.0d0)
    complex(kind(0d0)) :: wP = (0.0d0, 0.0d0)
    double precision   :: RA = 0.0d0
    double precision   :: RB = 0.0d0
    complex(kind(0d0)) :: d(0:maxnA, 0:maxnB, 0:maxnA+maxnB)
    complex(kind(0d0)) :: h(0:maxm, 0:maxm)
    complex(kind(0d0)) :: ref
    integer nA, nB, Nk, m, ierr

    call hermite_coef_d(zP, wP, RA, RB, maxnA, maxnB, d)
    do nA = 0, maxnA
       do nB = 0, maxnB
          do Nk = 0, maxnA+maxnB
             call hermite_coef_d_0(zP, wP, RA, RB, nA, nB, Nk, ref)
             EXPECT_EQ_C(ref, d(nA, nB, Nk), ierr)
          end do
       end do
    end do

    call hermite_1drm(zP, wP, maxm, h)
    do m = 0, maxm
       do Nk = 0, maxnA+maxnB
          call hermite_1drm_0(zP, wP, m, Nk, ref, ierr)
          EXPECT_EQ_C(ref, h(m, Nk), ierr)
       end do
    end do
    
  end subroutine test_coef_d
  subroutine test_prod
    use Mod_const, only : II
    use Mod_PWGTO, only : prod_gauss
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
    use Mod_PWGTO
    type(Obj_PWGTO) :: g
    integer, parameter :: n = 4
    complex(kind(0d0)) :: calc, ref
    complex(kind(0d0)) :: S(n,n)
    complex(kind(0d0)) :: gg(0:20)
    integer A, B
    integer ierr
    
    call PWGTO_new(g, n, 1, 4, ierr); CHK_ERR(ierr)
    g%ns(1)=2; g%zs(1)=1.2d0;  g%Ps(1) = 10.0d0
    g%ns(2)=3; g%zs(2)=(0.9d0, -0.8d0); 
    g%ns(3)=4; g%zs(3)=1.0d0; g%Rs(3) = 0.1d0;
    g%ns(4)=2; g%zs(4)=1.1d0; g%Rs(4) = 0.1d0;
    g%ncs(1)%typ = "0"
    call PWGTO_setup(g, ierr); CHK_ERR(ierr)
    
    call PWGTO_overlap(g, 1, 1, S, ierr); CHK_ERR(ierr)    

    A = 1
    
    calc = S(A, A) / (abs(g%ncs(1)%cs(A,1))**2)
    !    ref = gtoint(g%ns(A)*2, g%zP(A,A))
    call gtoint(10, g%zp(A,A), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 2
    calc = S(A, A) / (abs(g%ncs(1)%cs(A,1))**2)
    !ref = gtoint(g%ns(A)*2, g%zP(A,A))
    call gtoint(10, g%zp(A,A), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 3
    calc = S(A, A) / (abs(g%ncs(1)%cs(A,1))**2)
    !    ref = gtoint(g%ns(A)*2, conjg(g%zs(A))+g%zs(A))
    call gtoint(10, conjg(g%zs(A))+g%zs(A), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 3
    B = 4
    !    calc = S(A, B) / (conjg(g%ncs%cs(A,1)) * g%ncs%cs(B,1))
    calc = S(A, B) / (conjg(PWGTO_nterm(g,A))*PWGTO_nterm(g,B))
    ! ref = gtoint(g%ns(A)+g%ns(B), conjg(g%zs(A))+g%zs(B))
    call gtoint(10, conjg(g%zs(A))+g%zs(B), gg(:), ierr); CHK_ERR(ierr)
    ref = gg(g%ns(A)+g%ns(B))
    EXPECT_EQ_C(ref, calc, ierr)
    
  end subroutine test_overlap
  subroutine test_multipole_0
    use Mod_PWGTO
    use Mod_Math, only : gtoint
    type(Obj_PWGTO) :: g
    integer, parameter :: n = 5
    complex(kind(0d0)) :: Mm(n, n)
    complex(kind(0d0)) :: calc, ref
    integer A, m, ierr
    complex(kind(0d0)) :: gg(0:20)
    
    call PWGTO_new(g, n, 1, 4, ierr)
    g%ns(1)=0; g%zs(1)=1.1d0
    g%ns(2)=0; g%zs(2)=1.1d0; g%Rs(2) = 0.1d0
    g%ns(3)=2; g%zs(3)=1.1d0
    g%ns(4)=2; g%zs(4)=1.1d0
    call PWGTO_setup(g, ierr)
       
    A = 1
    call gtoint(10, g%zs(A)*2, gg(:), ierr); CHK_ERR(ierr)
    do m = 0, 2
       call PWGTO_multipole(g, 1,  m, 1, MM, ierr)
       calc = MM(A, A) / (abs(PWGTO_nterm(g,A))**2)
       !ref = gtoint(g%ns(A)*2 + m, g%zs(A)*2)
       ref = gg(g%ns(A)*2+m)
       EXPECT_EQ_C(ref, calc, ierr)
    end do
    
    A = 2
    m = 1
    call gtoint(20, g%zs(A)*2, gg(:), ierr); CHK_ERR(ierr)
    call PWGTO_multipole(g, 1, m, 1, MM, ierr)
    !calc = MM(A,A) / (abs(g%nc0%cs(A,1))**2)
    calc = MM(A,A) / (abs(PWGTO_nterm(g,A))**2)
    ref = g%Rs(A) * gg(g%ns(A)*2)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 4
    call gtoint(20, g%zs(A)*2, gg(:), ierr); CHK_ERR(ierr)
    call PWGTO_multipole(g, 1, 10, 1, MM, ierr)
    ! calc = MM(A,A) / (abs(g%nc0%cs(A,1))**2)
    calc = MM(A,A) / (abs(PWGTO_nterm(g,A))**2)
    ! ref = gtoint(14, g%zs(A)*2)
    ref = gg(14)
    EXPECT_NEAR_C(ref, calc, 1.0d-12, ierr)
    
  end subroutine test_multipole_0
  subroutine test_multipole_1
    use Mod_PWGTO
    type(Obj_PWGTO) :: g
    integer, parameter :: n = 2
    complex(kind(0d0)) :: calc, ref
    complex(kind(0d0)) :: S(n, n), R2(n, n), R0(n, n)
    integer A, B, ierr

    call PWGTO_new(g, n, 1, 2, ierr)
    g%ns(1)=1; g%zs(1)=0.4d0; g%Rs(1) = 0.1d0
    g%ns(2)=2; g%zs(2)=0.5d0; g%Rs(2) = 0.1d0
    call PWGTO_setup(g, ierr)

    call PWGTO_overlap(g,   1,    1, S, ierr)
    call PWGTO_multipole(g, 1, 0, 1, R0, ierr)
    call PWGTO_multipole(g, 1, 2, 1, R2, ierr)

    A = 1
    B = 1
    ref = S(A, A)
    calc= R0(A, A)
    EXPECT_EQ_C(ref, calc, ierr)

    A = 1
    B = 2
    calc = R2(A, B) / (conjg(PWGTO_nterm(g,A))*PWGTO_nterm(g,B))
    ref = 0.345987112123329
    EXPECT_NEAR_C(ref, calc, 1d-9, ierr)
  end subroutine test_multipole_1
  subroutine test_kinetic
    use Mod_PWGTO
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
    integer, parameter :: n = 4
    complex(kind(0d0)) :: calc, ref, z, zz, cz
    complex(kind(0d0)) :: T(n, n)
    double precision :: pA, pB
    integer ierr
    complex(kind(0d0)) :: gg(0:20)

    z = (1.0d0, 0.2d0)
    cz = conjg(z)
    zz = z + cz
    pA = 0.1d0
    pB = 0.2d0

    call PWGTO_new(g, n, 1, 2, ierr); CHK_ERR(ierr)
    g % ns(1) = 0; g % zs(1) = z; g % Rs(1) = 0.1d0
    g % ns(2) = 1; g % zs(2) = z; g % Rs(2) = 0.1d0
    g % ns(3) = 0; g % zs(3) = z; g % Ps(3) = pA;    g % Rs(3) = 0.2d0
    g % ns(4) = 1; g % zs(4) = z; g % Ps(4) = pB;    g % Rs(4) = 0.2d0
    call PWGTO_setup(g, ierr); CHK_ERR(ierr)

    call PWGTO_kineticP2(g, 1, 1, T, ierr); CHK_ERR(ierr)
    T = T/2

    call gtoint(4, zz, gg, ierr); CHK_ERR(ierr)
    !    ref = -0.5d0 * (-2*z*gtoint(0,zz) + 4*z*z*gtoint(2,zz))
    ref = -0.5d0 * (-2*z*gg(0) + 4*z*z*gg(2))
    calc = T(1, 1) / (abs(PWGTO_nterm(g, 1))**2)
    EXPECT_EQ_C(ref, calc, ierr)

    ref = 0.5d0 * (gg(0) +4*cz*z*gg(4) - 2*zz*gg(2))
    !    calc = T(2, 2) / (abs(g%nc0%cs(2,1))**2)
    calc = T(2, 2) / (abs(PWGTO_nterm(g,2))**2)
    EXPECT_EQ_C(ref, calc, ierr)

    EXPECT_EQ_C((0.0d0,0.0d0), T(1,2), ierr)
    EXPECT_EQ_C((0.0d0,0.0d0), T(2,1), ierr)
    
    ref = 0.5d0 * (4*z*cz*gg(2) + pA*pA*gg(0))
    !    calc = T(3, 3) / (abs(g%nc0%cs(3,1))**2)
    calc = T(3, 3) / (abs(PWGTO_nterm(g,3))**2)
    EXPECT_EQ_C(ref, calc, ierr)

    ref = 0.5d0 * (gg(0) + 4*z*cz*gg(4) + (pB*pB-2*z-2*cz)*gg(2))
    !    calc = T(4, 4) / (abs(g%nc0%cs(4,1))**2)
    calc = T(4, 4) / (abs(PWGTO_nterm(g,4))**2)
    EXPECT_EQ_C(ref, calc, ierr)
    
    call PWGTO_delete(g, ierr); CHK_ERR(ierr)
    
  end subroutine test_kinetic
end module UTestPWGTO

program main
  use UTestPWGTO
  call UTestPWGTO_run
end program main
