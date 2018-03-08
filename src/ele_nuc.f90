#include "macros.fpp"

module Mod_EleNuc1D
  ! 1 Dimensional system represented by Spline function
  use Mod_Spline
  implicit none
  private
  integer ne_
  type(Obj_Spline), allocatable :: spl_HeIJ_(:,:)
  type(Obj_Spline), allocatable :: spl_XkIJ_(:,:)
  public :: EleNuc1D_new, EleNuc1D_new_file, EleNuc1D_delete, EleNuc1D_H_X
contains
  subroutine EleNuc1D_new(ne, xs, HeIJ, XkIJ, ierr)
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
  end subroutine EleNuc1D_new
  subroutine EleNuc1D_new_file(fn, ierr)
    use Mod_sys, only : open_r
    use Mod_StrUtil, only : str2d
    character(*), intent(in) :: fn
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 2341
    integer nx, ix, ne
    character(256) :: line
    character(30) :: lines(5)
    integer, allocatable :: idxHe(:,:), idxXk(:,:)
    double precision, allocatable :: HeIJ(:,:,:), XkIJ(:,:,:), xs(:)
    ierr = 0
    ne = 2

    allocate(idxHe(ne,ne), idxXk(ne,ne))
    idxHe(:,:) = 0
    idxXk(:,:) = 0
    idxHe(1,1) = 2
    idxHe(1,2) = 3
    idxHe(2,2) = 4
    idxXk(1,2) = 5
    
    call open_r(ifile, fn, ierr); CHK_ERR(ierr)
    read(ifile, *) line
    ix = 0
    do
       read(ifile, *, end=100)
       ix = ix+1
    end do
100 continue
    nx = ix    

    allocate(xs(nx), HeIJ(nx,ne,ne), XkIJ(nx,ne,ne))
    HeIJ(:,:,:) = 0
    XkIJ(:,:,:) = 0
    
    rewind(ifile)
    read(ifile, *) line

    do ix = 1, nx
       read(ifile, *) lines
       call str2d(lines(1), xs(ix),       ierr); CHK_ERR(ierr)
       call str2d(lines(2), HeIJ(ix,1,1), ierr); CHK_ERR(ierr)
       call str2d(lines(3), HeIJ(ix,1,2), ierr); CHK_ERR(ierr)
       call str2d(lines(4), HeIJ(ix,2,2), ierr); CHK_ERR(ierr)
       call str2d(lines(5), XkIJ(ix,1,2), ierr); CHK_ERR(ierr)

       HeIJ(:,2,1) = HeIJ(:,1,2)
       XkIJ(:,2,1) = -XkIJ(:,1,2)
    end do    
    close(ifile)

    call EleNuc1d_new(ne, xs, HeIJ, XkIJ, ierr); CHK_ERR(ierr)
    
  end subroutine EleNuc1D_new_file
  subroutine EleNuc1D_delete
    integer i, j
    do i = 1, ne_
       do j = 1, ne_
          call Spline_delete(spl_HeIJ_(i,j))
          call Spline_delete(spl_XkIJ_(i,j))
       end do
    end do
    deallocate(spl_HeIJ_, spl_XkIJ_)
  end subroutine EleNuc1D_delete
  subroutine EleNuc1D_H_X(Q, HeIJ, XkIJ, ierr)
    double precision, intent(in) :: Q(:)
    complex(kind(0d0)), intent(out) :: HeIJ(:,:)
    complex(kind(0d0)), intent(out) :: XkIJ(:,:,:)
    integer, intent(out) :: ierr
    double precision h, x
    integer i, j
    ierr = 0
    do i = 1, ne_
       do j = 1, ne_
          call Spline_at(spl_HeIJ_(i,j), Q(1), h, ierr)
          HeIJ(i,j) = h
          call Spline_at(spl_XkIJ_(i,j), Q(1), x, ierr)
          XkIJ(1,i,j) = x
       end do
    end do
  end subroutine EleNuc1D_H_X
end module Mod_EleNuc1D
