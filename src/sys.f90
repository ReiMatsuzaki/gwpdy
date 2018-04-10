! gwpdy/src/sys.f90
#include "macros.fpp"
module mod_sys
  implicit none
contains
  subroutine mkdir(path)
    character(*), intent(in) :: path
    character(100) cmd
    write(cmd, "( 'mkdir ', A )") trim(path)
    call system(cmd)
  end subroutine mkdir
  subroutine mkdir_if_not(path)
    character(*), intent(in) :: path
    character(1000) cmd
    write(cmd, "('ls ', A, ' > /dev/null 2>&1 || mkdir ', A )") trim(path), trim(path)
    call system(cmd)
  end subroutine mkdir_if_not
  subroutine mkdirp_if_not(path)
    character(*), intent(in) :: path
    character(1000) cmd
    write(cmd, "('ls ', A, ' > /dev/null 2>&1 || mkdir -p ', A )") trim(path), trim(path)
    call system(cmd)
  end subroutine mkdirp_if_not
  subroutine open_w(ifile, fn, ierr)
    integer, intent(in) :: ifile
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    ierr = 0
    open(ifile, file=fn, status='replace', err=999)
    return
999 continue
    MSG_ERR("failed open file")
    ierr = 1
    write(0,*) "filename:", trim(fn)
    return
  end subroutine open_w
  subroutine open_r(ifile, fn, ierr)
    integer, intent(in) :: ifile
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    ierr = 0
    open(ifile, file=fn, status='old', err=999)
    return
999 continue
    MSG_ERR("failed open file")
    ierr = 1
    write(0,*) "filename:", trim(fn)
    return
  end subroutine open_r
end module mod_sys

