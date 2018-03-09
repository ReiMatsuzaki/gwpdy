#include "macros.fpp"

module Mod_MoleFit
  ! 1 Dimensional system represented by Spline function
  use Mod_Spline
  implicit none
  private
  integer ne_
  type(Obj_Spline), allocatable :: spl_HeIJ_(:,:)
  type(Obj_Spline), allocatable :: spl_XkIJ_(:,:)
  public :: MoleFit_new, MoleFit_new_file, MoleFit_delete, MoleFit_H_X
contains
  subroutine MoleFit_new(ne, xs, HeIJ, XkIJ, ierr)
    integer, intent(in) :: ne
    double precision, intent(in) :: xs(:)
    double precision, intent(in) :: HeIJ(:,:,:)
    double precision, intent(in) :: XkIJ(:,:,:)
    integer, intent(out) :: ierr
    integer nx, i,j

    ne_ = ne
    nx = size(xs)
    
    if(  nx.ne.size(HeIJ,1) .or. nx.ne.size(XkIJ,1) .or. &
         ne.ne.size(HeIJ,2) .or. ne.ne.size(HeIJ,3) .or. &
         ne.ne.size(XkIJ,2) .or. ne.ne.size(XkIJ,3)) then
       MSG_ERR("size mismatch")
       write(0,*) "shape(xs):", shape(xs)
       write(0,*) "shape(HeIJ):", shape(HeIJ)
       write(0,*) "shape(XkIJ):", shape(XkIJ)
       ierr = 1; return
    end if
    
    allocate(spl_HeIJ_(ne,ne))
    allocate(spl_XkIJ_(ne,ne))

    do i = 1, ne_
       do j = 1, ne_
          call Spline_new(spl_HeIJ_(i,j), nx, xs, HeIJ(:,i,j), "error", ierr)
          CHK_ERR(ierr)
          call Spline_new(spl_XkIJ_(i,j), nx, xs, XkIJ(:,i,j), "error", ierr)
          CHK_ERR(ierr)
       end do
    end do    
  end subroutine MoleFit_new
  subroutine MoleFit_new_file(fn_xs, fn_HeIJ, fn_XkIJ, ierr)
    use Mod_sys, only : open_r
    character(*), intent(in) :: fn_xs, fn_HeIJ, fn_XkIJ
    integer, intent(out) :: ierr
    integer, parameter :: ifile_xs = 2341, ifile_h = 2342
    integer, parameter :: ifile_xk = 2343
    integer nx, ne, ix, i, j, k
    character(256) :: line
    double precision, allocatable :: HeIJ(:,:,:), XkIJ(:,:,:), xs(:)
    double precision :: v
    ierr = 0

    call open_r(ifile_xs, fn_xs,   ierr); CHK_ERR(ierr)
    call open_r(ifile_h,  fn_HeIJ, ierr); CHK_ERR(ierr)
    call open_r(ifile_xk, fn_XkIJ, ierr); CHK_ERR(ierr)

    ! -- determine data size --    
    read(ifile_h, *) line
    nx = 0
    ne = 0
    do
       read(ifile_h, *, end=100) ix, i, j, v
       if(nx<ix) then
          nx = ix
       end if
       if(ne<max(i,j)) then
          ne = max(i,j)
       end if
    end do
100 continue
    rewind(ifile_h)

    ! -- allocation --
    allocate(xs(nx), HeIJ(nx,ne,ne), XkIJ(nx,ne,ne))

    ! -- set data --
    read(ifile_xs,*) line
    do ix = 1, nx
       read(ifile_xs, *) v
       xs(ix) = v
    end do
    close(ifile_xs)

    HeIJ(:,:,:) =0
    read(ifile_h, *) line
    do
       read(ifile_h, *, end=101) ix, i, j, v
       HeIJ(ix,i,j) = v
       HeIJ(ix,j,i) = v
    end do
101 continue    
    close(ifile_h)

    XkIJ(:,:,:) = 0
    read(ifile_xk, *) line
    do
       read(ifile_xk, *, end=102) ix, k, i, j, v
       XkIJ(ix,i,j) = v
       if(i.ne.j) then
          XkIJ(ix,j,i) = -v
       end if

       if(k.ne.1) then
          MSG_ERR("more than one dim is not supported")
          ierr = 1; return          
       end if       
    end do
102 continue
    close(ifile_xk)
    call MoleFit_new(ne, xs, HeIJ, XkIJ, ierr); CHK_ERR(ierr)
    
  end subroutine MoleFit_new_file
  subroutine MoleFit_delete
    integer i, j
    do i = 1, ne_
       do j = 1, ne_
          call Spline_delete(spl_HeIJ_(i,j))
          call Spline_delete(spl_XkIJ_(i,j))
       end do
    end do
    deallocate(spl_HeIJ_, spl_XkIJ_)
  end subroutine MoleFit_delete
  subroutine MoleFit_H_X(Q, HeIJ, XkIJ, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: HeIJ(:,:)
    complex(kind(0d0)), intent(out) :: XkIJ(:,:,:)
    integer, intent(out) :: ierr
    double precision h, x
    integer i, j
    ierr = 0
    do i = 1, ne_
       do j = 1, ne_
          call Spline_at(spl_HeIJ_(i,j), Q(1), h, ierr); CHK_ERR(ierr)
          HeIJ(i,j) = h
          call Spline_at(spl_XkIJ_(i,j), Q(1), x, ierr); CHK_ERR(ierr)
          XkIJ(1,i,j) = x
       end do
    end do
  end subroutine MoleFit_H_X
end module Mod_MoleFit
