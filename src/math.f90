! gwpdy/src/math.f90
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
  subroutine lapack_zgesv_1(n, A, res, ierr)
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: A(:,:)    
    complex(kind(0d0)), intent(inout) :: res(:)
    integer, intent(out) :: ierr
    integer info, ipiv(n)
    complex(kind(0d0)) :: AA(n,n)
    complex(kind(0d0)) :: tmp(n,1)

    ierr = 0
    AA(:,:) = A(:,:)
    tmp(:,1) = res(:)
    call ZGESV(n, 1, AA, n, ipiv, tmp, n, info)
    res(:) = tmp(:,1)
    
  end subroutine lapack_zgesv_1
  subroutine lapack_dgesv_1(n, A, x_res, ierr)
    integer, intent(in) :: n
    double precision, intent(in) :: A(:,:)    
    double precision, intent(inout) :: x_res(:)
    integer, intent(out) :: ierr
    integer info, ipiv(n)
    double precision :: AA(n,n)
    double precision :: tmp(n,1)

    if(size(A,1).ne.n) then
       MSG_ERR("invalid size")
       write(0,*) "A:", size(A,1), size(A,2)
       write(0,*) "n:", n
       ierr = 1; return
    end if

    ierr = 0
    AA(:,:) = A(:,:)
    tmp(:,1) = x_res(:)
    call DGESV(n, 1, AA, n, ipiv, tmp, n, info)
    if(info.ne.0) then
       MSG_ERR("error on DFESV")
       write(0,*) "info:", info
       ierr=1; return
    end if
    x_res(:) = tmp(:,1)
        
  end subroutine lapack_dgesv_1
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
  subroutine intet_gdiag(n, S, H, dt, c, ierr)
    use Mod_const, only : II
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: S(n,n), H(n,n)
    double precision, intent(in) :: dt
    complex(kind(0d0)), intent(inout) :: c(n)
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: w(n)
    complex(kind(0d0)) :: UL(n,n), UR(n,n)

    ! call lapack_zggev_shift(n, H, S, H(1,1), w, UL, UR, ierr); CHK_ERR(ierr)
    call lapack_zggev(n, H, S, w, UL, UR, ierr); CHK_ERR(ierr)
    c(:) = matmul(transpose(conjg(UL)), matmul(S,c))
    c(:) = exp(-II*w(:)*dt) * c(:)
    c(:) = matmul(UR(:,:), c(:))
    
  end subroutine intet_gdiag
  subroutine sort_zggev(n, w, UL, UR, ierr)
    ! sort results of lapack_zggev
    integer, intent(in) :: n
    complex(kind(0d0)), intent(inout) :: w(:), UL(:,:), UR(:,:)
    integer, intent(out) :: ierr
    complex(kind(0d0)) :: tmp_uR(n), tmp_uL(n), tmp_w
    integer A, B, num

    num = size(w)
    ierr = 0
    
    do A = 1, num
       do B = num, A+1, -1
          if(real(w(B-1)) < real(w(B))) then
             
             tmp_w  = w(B-1)
             tmp_uR = UR(:,B-1)
             tmp_uL = UL(:,B-1)

             w(B-1)    = w(B)
             UR(:,B-1) = UR(:,B)
             UL(:,B-1) = UL(:,B)

             w(B)    = tmp_w
             UR(:,B) = tmp_uR(:)
             UL(:,B) = tmp_uL(:)
             
          end if
       end do
    end do
    
  end subroutine sort_zggev
  ! ==== mathematical function ====
  subroutine dfact(maxn, res, ierr)
    ! gives double factorial
    integer, intent(in) :: maxn
    integer, intent(out) :: res(0:maxn)
    integer, intent(out) :: ierr
    integer n
    ierr = 0
    res(0) = 1
    res(1) = 1
    do n = 2, maxn
       res(n) = n*res(n-2)
    end do
  end subroutine dfact
  subroutine gamma_half_int(maxn, res, ierr)
    ! gives {Gamma(n+1/2) | n = 0,...,maxn}
    use Mod_const, only : PI
    integer, intent(in) :: maxn
    double precision, intent(out) :: res(0:maxn)
    integer, intent(out) :: ierr
    integer n
    
    ierr = 0
    res(0) = sqrt(PI)
    do n = 0, maxn-1
       ! Gamma(n+1+1/2) = (n+1/2)Gamma(n+1/2)
       res(n+1) = (n+0.5d0) * res(n)
    end do
    
  end subroutine gamma_half_int
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
  subroutine gtoint_shift(maxn, a, b, q0, res, ierr)
    ! compute integration
    !      J_n = Int_{-oo}^{+oo} dq (q-q0)^n Exp[-a (q-q0)^2 - b q^2]
    ! by the recursive formula.
    !      (1+b/a) J_n = (n-1)/(2a) J_{n-2} - b.q_0/a. J_{n-1}
    ! See 2018/4/16/qpbranch for detail
    
    use Mod_const, only : PI
    integer, intent(in) :: maxn
    complex(kind(0d0)), intent(in) :: a, b, q0
    complex(kind(0d0)), intent(out) :: res(0:)
    integer, intent(out) :: ierr
    integer n

    ierr = 0
    if(maxn<0) then
       MSG_ERR("maxn<0");
       write(0,*) "maxn:", maxn
       ierr = 1; return       
    end if
    
    res(0) = sqrt(PI) * exp(-(a*b*q0**2)/(a+b)) / sqrt(a+b)
    if(maxn.eq.0) return
    
    res(1) = sqrt(PI) * b*(-q0)*exp(a*q0**2*(a/(a+b)-1))/((a+b)**1.5d0)
    if(maxn.eq.1) return
    
    do n = 2, maxn
       res(n) = ((n-1)/(2*a)*res(n-2) - b*q0/a*res(n-1)) / (1+b/a)
    end do

  end subroutine gtoint_shift
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
