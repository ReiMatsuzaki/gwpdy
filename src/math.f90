#include "macros.fpp"

module Mod_math
  implicit none
contains
  function vmv(a, S, b) result(res)
    complex(kind(0d0)), intent(in) :: a(:), S(:,:), b(:)
    complex(kind(0d0)) res
    res = dot_product(a, matmul(S, b))
  end function vmv
  subroutine lapack_zheev(n, H, w, U, ierr)
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: H(n,n)
    double precision, intent(out) :: w(n)
    complex(kind(0d0)), intent(out) :: U(n,n)
    integer, intent(out) :: ierr
    double precision, parameter :: eps = 1.0d-12

    complex(kind(0d0)) :: work(2*n)
    double precision :: rwork(3*n)
    integer :: info

    info = 0
    U(:,:) = H(:,:)
    call ZHEEV('V', 'U', n, U, n, w, work, 2*n, rwork, info)

    if(info .ne. 0) then
       MSG_ERR("Error on ZHEEV")
       ierr = 1
       write(0,*) "info:", info
       write(0,*) "n:", n
       write(*,*) "H:"
       if(size(H,1) > 3) then
          write(0,*) H(1:3,1:3)
          write(0,*) H(1:3,1:3)
          write(0,*) H(1:3,1:3)
       else
          write(*,*) H(:,:)
       end if
       return
    end if
    
  end subroutine lapack_zheev
  subroutine lapack_zggev(n, H, S, w, UL, UR, ierr)
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: H(n, n), S(n, n)
    complex(kind(0d0)), intent(out) :: w(n)
    complex(kind(0d0)), intent(out) :: UL(n,n), UR(n,n)
    integer, intent(out) :: ierr
    integer info
    complex(kind(0d0)) :: work(n*2)
    double precision   :: rwork(8*n)
    complex(kind(0d0)) :: a(n), b(n)
    complex(kind(0d0)) :: HH(n, n), SS(n, n)
    integer i, n0
    complex(kind(0d0)) norm2
    double precision, parameter :: eps = 1.0d-12
    
    ierr = 0
    HH = H    
    SS = S
    info = 0

    call ZGGEV('V', 'V', n, HH, n, SS, n, a, b,&
         UL, n, UR, n, work, 2*n, rwork, info)
    if(info .ne. 0) then
       MSG_ERR("Error on ZGGEV")
       ierr = 1
       write(0,*) "info:", info
       write(0,*) "Inputs:"
       write(0,*) "n:", n
       write(0,*) "w:", a/b
       
       if(n>3) then
          n0 = 3
       else
          n0 = n
       end if
       write(0,*) "H:"
       do i = 1, n0
          write(0,*) H(i,:n0)
       end do
       write(0,*) "S:"
       do i = 1, n0
          write(0,*) S(i,:n0)
       end do
       if(info.eq.n+1) then
          write(0,*) "unexpected error from ZHGEQZ"
       else if(info.eq.n+2) then
          write(0,*) "unexpected error from ZTGEVC"
       else if(info < 0) then
          write(0,*) -info, "argument had an illegal value"
       end if
       return
    end if

    w = a/b
    do i = 1, n
       norm2 = vmv(UL(:,i), S, UR(:,i))
       UL(:,i) = UL(:,i)/conjg(sqrt(norm2))
       UR(:,i) = UR(:,i)/sqrt(norm2)
    end do

  end subroutine lapack_zggev
  subroutine lapack_zggev_shift(n, H, S, h0, w, UL, UR, ierr)
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: H(n, n), S(n, n)
    complex(kind(0d0)), intent(in) :: h0
    complex(kind(0d0)), intent(out) :: w(n)
    complex(kind(0d0)), intent(out) :: UL(n,n), UR(n,n)
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: HH(n,n)
    integer I

    HH=H
    do I = 1, n
       HH(I,I) = H(I,I) - h0
    end do
    call lapack_zggev(n, HH, S, w, UL, UR, ierr); CHK_ERR(ierr)
    w(:) = w(:) + h0
    
  end subroutine lapack_zggev_shift
end module Mod_math
