module Mod_const
  ! == mathematics ==
  double precision, parameter :: PI = 3.14159265358979323846d0
  complex(kind(0d0)), parameter :: II = (0.0d0, 1.0d0)
  ! == chemical ==
  double precision, parameter :: NABOGADRO = 6.022140857d23
  ! == atomic unit ==
  double precision, parameter :: MELE_KG = 9.1093826d-31
  double precision, parameter :: ABOHR_M = 5.291772108d-11
  double precision, parameter :: HBAR_JS = 1.05457168d-34
  double precision, parameter :: EHARTREE_J = 4.35974417d-18
  double precision, parameter :: EELECTRON_C = 1.60217653d-19
  ! == unit conversion ==
  double precision, parameter :: AU_FS = HBAR_JS/EHARTREE_J * 1.0d15
  double precision, parameter :: AU_ANG = ABOHR_M * 1.0d10
  double precision, parameter :: AU_ARW = MELE_KG*1000*NABOGADRO
  double precision, parameter :: AU_EV  = EHARTREE_J/EELECTRON_C

! 2.4188843172419448E-002
! 0.52917721080000002     
! 5.4857985137504890E-004
! 27.211384565719481     
end module Mod_const
