#include "macros.fpp"

module Mod_math
contains
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
end module Mod_math
