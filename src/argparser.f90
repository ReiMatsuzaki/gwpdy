#include "macros.fpp"

module Mod_ArgParser
  implicit none
contains
  subroutine throw_err_failed_find_name(name, ierr)
    character(*), intent(in) :: name
    integer, intent(out) :: ierr
    begin_err("failed to find option")
    ierr = 1
    write(0,*) "name:", name
    return
  end subroutine throw_err_failed_find_name
  subroutine arg_parse_i(name, x, ierr)
    use Mod_StrUtil, only : str2i
    character(*), intent(in) :: name
    integer, intent(out) :: x
    integer, intent(out) :: ierr
    integer i, num
    character(100) :: ele
    ierr = 0
    num = iargc()
    do i = 1, num
       call getarg(i, ele)       
       if(name.eq.ele) then
          call getarg(i+1, ele)
          call str2i(ele, x, ierr); check_err(ierr)
          return
       end if
    end do
    
    call throw_err_failed_find_name(name, ierr) ; check_err(ierr)
  end subroutine arg_parse_i
  subroutine arg_parse_d(name, x, ierr)
    use Mod_StrUtil, only : str2d
    character(*), intent(in) :: name
    double precision, intent(out) :: x
    integer, intent(out) :: ierr
    integer i, num
    character(100) :: ele
    ierr = 0
    num = iargc()
    do i = 1, num
       call getarg(i, ele)
       if(name.eq.ele) then
          call getarg(i+1, ele)
          call str2d(ele, x, ierr); check_err(ierr)
          return
       end if
    end do
    call throw_err_failed_find_name(name, ierr); check_err(ierr)
  end subroutine arg_parse_d
  subroutine arg_parse_s(name, x, ierr)
    character(*), intent(in) :: name
    character(*), intent(out) :: x
    integer, intent(out) :: ierr
    integer i, num
    character(100) :: ele
    ierr = 0
    num = iargc()
    do i = 1, num
       call getarg(i, ele)
       if(name.eq.ele) then
          call getarg(i+1, ele)
          x = ele
          return
       end if
    end do
    call throw_err_failed_find_name(name,ierr); check_err(ierr)
  end subroutine arg_parse_s
  subroutine arg_parse_dvec(name, xs, ierr)
    use Mod_StrUtil, only : str2d
    character(*), intent(in) :: name
    double precision, intent(out) :: xs(:)
    integer, intent(out) :: ierr
    integer i, num, j
    character(100) :: ele

    num = iargc()
    do i = 1, num
       call getarg(i, ele)
       if(name.eq.ele) then
          do j = 1, num-i
             call getarg(i+j, ele)
             !if(ele(1:1).eq."-") then
             !   return
             !end if
             call str2d(ele, xs(j), ierr); check_err(ierr)
          end do
          return
       end if
    end do
    call throw_err_failed_find_name(name, ierr); check_err(ierr)    
  end subroutine arg_parse_dvec
  function arg_parse_exist(name) result(res)
    character(*), intent(in) :: name
    logical :: res
    integer i, num
    character(100) :: ele
    
    num = iargc()    
    do i = 1, num
       call getarg(i, ele)
       if(trim(name).eq.trim(ele)) then
          res = .true.
          return
       end if
    end do

    res = .false.
    
  end function arg_parse_exist
end module Mod_ArgParser
