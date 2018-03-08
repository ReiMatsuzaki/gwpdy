#include "macros.fpp"
module Mod_StrUtil
  implicit none
contains
  subroutine str2i(str, res, ierr)
    implicit none
    character(*), intent(in) ::  str
    integer, intent(out) :: res
    integer, intent(out) :: ierr
    double precision  :: xd
    
    integer ierr_d, ierr_i, idx
 
    double precision eps

    ierr = 0
    eps = 1.0D-10
    read(str, *, iostat=ierr_i) res
    read(str, *, iostat=ierr_d) xd    
    
    if(ierr_i .ne. 0 .or. ierr_d.ne.0) then
       MSG_ERR("failed to convert")
       ierr = 1
       write(0,*) "str: ", str
       return
    end if

    idx = index(str, ".")
    if(abs(xd-res) > eps .or. idx.ne.0) then
       MSG_ERR("str is float not integer")
       ierr = 2
       write(0,*) "str: ", str
       return
    end if
      
  end subroutine str2i
  subroutine str2d(str, res, ierr)
    implicit none
    character(*), intent(in) :: str
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    
    read(str, *, iostat=ierr) res
    if(ierr.ne.0) then
       MSG_ERR("failed to convert")
       ierr = 1
       write(0,*) "str:", str
       return
    end if
    
  end subroutine str2d
  function is_i(str) result(res)
    character(*), intent(in) :: str
    logical :: res
    integer :: i, ierr
    call str2i(str, i, ierr)
    res = ierr.eq.0
  end function is_i
  function is_d(str) result(res)
    character(*), intent(in) :: str
    logical :: res
    double precision d
    integer :: ierr
    call str2d(str, d, ierr)
    res = ierr.eq.0    
  end function is_d
  subroutine str_split(line, sep, n, res, ierr)
    character(*), intent(in) :: line
    character(*), intent(in) :: sep
    integer, intent(out) :: n
    character(*), intent(out) :: res(:)
    integer, intent(out) :: ierr
    integer :: i, i0
    ierr = 0
    
    n = 1
    i0 = 1
    do i = 1, len_trim(line)-1
       if(line(i:i).eq.trim(sep)) then
          res(n) = line(i0:i-1)
          i0 = i+1
          n = n + 1
       end if
    end do
    res(n) = line(i0:len_trim(line))
!    if(n==1) then
!       res(1) = line
!    end if
    
  end subroutine str_split
  subroutine str2vec(str, nx, xs, ierr)
    character(100), intent(in) :: str
    integer, intent(out) :: nx
    double precision, intent(out) :: xs(:)
    integer, intent(out) :: ierr
    character(50) :: lines(10)
    double precision :: x0, x1
    integer :: n, i

    call str_split(str, ",", n, lines, ierr); CHK_ERR(ierr)
    if(n < 1) then
       MSG_ERR("invalid argument")
       ierr = 1
       return
    end if
    select case(lines(1))
    case("linspace")
       if(n .ne. 4) then
          MSG_ERR("invalid argument")
          ierr = 1
          return
       end if       
       call str2i(lines(4), nx, ierr)
       if(ierr.ne.0) then
          MSG_ERR("invalid argument")
          ierr = 1
          write(0,*) "str:", str
          write(0,*) "line[4]:", lines(4)
          return
       end if
       call str2d(lines(2), x0, ierr); CHK_ERR(ierr)
       call str2d(lines(3), x1, ierr); CHK_ERR(ierr)
       do i = 1, nx
          xs(i) = x0 + (i-1)*(x1-x0)/(nx-1)
       end do
    case("scalar")       
       nx = 1
       call str2d(lines(2), xs(1), ierr); CHK_ERR(ierr)
    case default
       MSG_ERR("not supported")
       ierr = 1
       write(0,*) "str: ", str
       return
    end select
  end subroutine str2vec
end module Mod_StrUtil
