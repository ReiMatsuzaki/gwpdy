#include "macros.fpp"

module Mod_Main
  use Mod_DyMono
  implicit none
contains
  subroutine Main_run
    call new
  end subroutine Main_run
  subroutine new
    use Mod_ArgParser
    integer :: ierr
    write(*,*) 
    write(*,*) "Gauss Wave Packet DYnamics (GWPDY)"
    write(*,*)
    
    write(*,*) "new begin"
    call arg_parse_i("-nf", nf_, ierr); check_err(ierr)
    if(nf_.ne.1) then
       begin_err("only nf=1 is supported")
       ierr = 1; return
    end if
    call arg_parse_i("-ne", ne_, ierr); check_err(ierr)

    call arg_parse_s("-hel", , ierr); check_err()
    
    call DyMono_new(nf, ne, ierr); check_err(ierr)
    write(*,*) "nf:", nf_
    write(*,*) "ne:", ne_

    
    
    write(*,*) "new end"
  end subroutine new
end module Mod_Main

program main
  use Mod_Main
  call Main_run
end program main
