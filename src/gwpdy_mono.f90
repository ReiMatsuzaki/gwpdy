#include "macros.fpp"

module Mod_GWPDyMono
  ! molecular wave function \Psi are expanded as 
  !     \Psi(r,Q,t) = G(Q,t) \sum_I c_I(t) \Phi_I(r;Q) .
  ! G is auss wave packet
  use Mod_GWP
  implicit none
  integer :: nf_        ! number of freedom
  integer :: nebasis_   ! number of electron Hilbert space
  type(Obj_GWP) :: gwp_
  complex(kind(0d0)), allocatable :: c_(:)
contains
  subroutine GWPDyMono_run
    
  end subroutine GWPDyMono_run
  subroutine GWPDyMono_new(nf, nebasis, ierr)
    integer, intent(in) :: nf, nebasis
    integer, intent(out) :: ierr
    ierr = 0
    nf_  = nf
    nebasis_ = nebasis
    call GWP_new(gwp_, nf, 1, 'c', ierr) ; check_err(ierr)
    allocate(c_(nebasis))
    c_ = 0
    
    
  end subroutine GWPDyMono_new
  subroutine GWPDyMono_setup(ierr)
    integer, intent(out) :: ierr
    double precision :: norm2
    integer i
    double precision, parameter :: tol = 1.0d-15
    ierr = 0
    call GWP_setup(gwp_, ierr); check_err(ierr)
    norm2 = sum(abs(c_(i))**2)
    if(norm2<tol) then
       begin_err("norm is too small")
       ierr = 1
       write(1,*) "norm2:", norm2
       return
    end if
    c_(:) = c_(:) / sqrt(norm2)
  end subroutine GWPDyMono_setup
  subroutine GWPDyMono_delete(ierr)
    integer, intent(out) :: ierr
    ierr = 0
    call GWP_delete(gwp_)
    deallocate(c_)
  end subroutine GWPDyMono_delete
end module Mod_GWPDyMono
