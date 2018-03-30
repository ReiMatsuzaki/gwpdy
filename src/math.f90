#include "macros.fpp"

module Mod_math
  implicit none
contains
  ! ==== Linear Algebra ====
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
  subroutine lapack_zgesv(n, A, X, res, ierr)
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: A(:,:)    
    complex(kind(0d0)), intent(in) :: X(:,:)
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    integer info, ipiv(n)
    complex(kind(0d0)) :: AA(n,n)
    
    ierr = 0
    AA(:,:) = A(:,:)
    res(:,:) = X(:,:)

    call ZGESV(n, n, AA, n, ipiv, res, n, info)
    
  end subroutine lapack_zgesv
  subroutine lapack_zgesv_1(n, A, x, res, ierr)
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: A(:,:)    
    complex(kind(0d0)), intent(in) :: X(:)
    complex(kind(0d0)), intent(out) :: res(:)
    integer, intent(out) :: ierr
    integer info, ipiv(n)
    complex(kind(0d0)) :: AA(n,n)
    complex(kind(0d0)) :: tmp(n,1)

    ierr = 0
    AA(:,:) = A(:,:)
    tmp(:,1) = x(:)

    call ZGESV(n, 1, AA, n, ipiv, tmp, n, info)
    res(:) = tmp(:,1)
    
  end subroutine lapack_zgesv_1
  subroutine intet_diag(n, H, dt, c, ierr)
    use Mod_const, only : II
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: H(n,n)
    double precision, intent(in) :: dt
    complex(kind(0d0)), intent(inout) :: c(n)
    integer, intent(out) :: ierr
    double precision :: w(n)
    complex(kind(0d0)) :: U(n,n), UH(n,n)

    call lapack_zheev(n, H, w, U, ierr); CHK_ERR(ierr)
    UH(:,:) = conjg(transpose(U(:,:)))
    c(:)   = matmul(UH(:,:),    c(:))
    c(:)   = exp(-II*w(:)*dt) * c(:)
    c(:)   = matmul(U(:,:),     c(:))
    
  end subroutine intet_diag
  ! ==== mathematical function ====
  subroutine gtoint(maxn, z, res, ierr)
    ! gives the integrations :  { Int_0^oo x^{n}exp[-zx^2] dx | n=0,...,maxn}
    use Mod_const, only : pi
    integer, intent(in) :: maxn
    complex(kind(0d0)), intent(in) :: z
    complex(kind(0d0)), intent(out) :: res(0:)
    integer, intent(out) ::  ierr
    integer n

    ierr = 0
    if(maxn < 0) then
       MSG_ERR("maxn must be 0 or positive")
       ierr = 1
       return 
    end if

    if(size(res)<maxn+1) then
       MSG_ERR("invalid size")
       ierr = 1
       write(0,*) "size(res):", size(res)
       write(0,*) "maxn:", maxn
       return
    end if

    res(0) = sqrt(pi/z)
    if(maxn .eq. 0) then
       return
    end if

    res(1) = 0.0d0
    if(maxn .eq. 1) then
       return
    end if

    do n = 2, maxn
       res(n) = (n-1)/(2*z) * res(n-2)
    end do
    
  end subroutine gtoint
  ! ==== IO ====
  subroutine ivec2csv(x, fn, ierr)
    use Mod_sys, only : open_w
    integer, intent(in) :: x(:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("val")')
    do i = 1, size(x, 1)
       write(ifile,'(I0)') x(i)
    end do
    close(ifile)
    
  end subroutine ivec2csv
  subroutine dten2csv(x, fn, ierr)
    use Mod_sys, only : open_w
    double precision, intent(in) :: x(:,:,:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i, j, k
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,j,k,val")')
    do i = 1, size(x, 1)
       do j = 1, size(x, 2)
          do k = 1, size(x, 3)
             write(ifile,'(I0,",",I0,",",I0,",",F20.10)') i,j,k,x(i,j,k)
          end do
       end do
    end do
    close(ifile)
    
  end subroutine dten2csv
  subroutine cten2csv(x, fn, ierr)
    use Mod_sys, only : open_w
    complex(kind(0d0)), intent(in) :: x(:,:,:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i, j, k
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,j,k,re,im")')
    do i = 1, size(x, 1)
       do j = 1, size(x, 2)
          do k = 1, size(x, 3)
             write(ifile,'(I0,",",I0,",",I0,",",F20.10,F20.10)') i,j,k,real(x(i,j,k)),aimag(x(i,j,k))
          end do
       end do
    end do
    close(ifile)
    
  end subroutine cten2csv
  subroutine dvec2csv(x, fn, ierr)
    use Mod_sys, only : open_w
    double precision, intent(in) :: x(:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,val")')
    do i = 1, size(x, 1)
       write(ifile,'(I0,",",F20.10)') i,x(i)
    end do
    close(ifile)
    
  end subroutine dvec2csv
  subroutine cvec2csv(x, fn, ierr)
    use Mod_sys, only : open_w
    complex(kind(0d0)), intent(in) :: x(:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    integer i
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,re,im")')
    do i = 1, size(x, 1)
       write(ifile,'(I0,",",F20.10,F20.10)') i,real(x(i)),aimag(x(i))
    end do
    close(ifile)
    
  end subroutine cvec2csv
  subroutine cmat2csv(x, fn, ierr)
    use Mod_sys, only : open_w
    complex(kind(0d0)), intent(in) :: x(:,:)
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 16231
    complex(kind(0d0)) y
    integer i ,j
    
    call open_w(ifile, fn, ierr); CHK_ERR(ierr)    
    write(ifile, '("i,j,re,im")')
    do i = 1, size(x, 1)
       do j = 1, size(x, 2)
          y = x(i,j)
          write(ifile,'(I0,",",I0,",",F20.10,F20.10)') i,j,real(y),aimag(y)
       end do
    end do
    close(ifile)
    
  end subroutine cmat2csv
end module Mod_math
