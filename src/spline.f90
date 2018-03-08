#include "macros.fpp"
module Mod_spline
  use mod_const
  implicit none
  type Obj_Spline
     integer :: n
     double precision, allocatable :: x(:), y(:), b(:), c(:), d(:)
     integer :: outvalue
  end type Obj_Spline
contains
  subroutine make_spline (x, y, b, c, d, n)
  !======================================================================
  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !  for  x(i) <= x <= x(i+1)
  !  Alex G: January 2010
  !----------------------------------------------------------------------
  !  input..
  !  x = the arrays of data abscissas (in strictly increasing order)
  !  y = the arrays of data ordinates
  !  n = size of the arrays xi() and yi() (n>=2)
  !  output..
  !  b, c, d  = arrays of spline coefficients
  !  comments ...
  !  spline.f90 program is based on fortran version of program spline.f
  !  the accompanying function fspline can be used for interpolation
  !======================================================================
  implicit none
  integer n
  double precision x(n), y(n), b(n), c(n), d(n)
  integer i, j, gap
  double precision h

  gap = n-1
  ! check input
  if ( n < 2 ) return
  if ( n < 3 ) then
     b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
     c(1) = 0.
     d(1) = 0.
     b(2) = b(1)
     c(2) = 0.
     d(2) = 0.
     return
  end if
  !
  ! step 1: preparation
  !
  d(1) = x(2) - x(1)
  c(2) = (y(2) - y(1))/d(1)
  do i = 2, gap
     d(i) = x(i+1) - x(i)
     b(i) = 2.0*(d(i-1) + d(i))
     c(i+1) = (y(i+1) - y(i))/d(i)
     c(i) = c(i+1) - c(i)
  end do
  !
  ! step 2: end conditions
  !
  b(1) = -d(1)
  b(n) = -d(n-1)
  c(1) = 0.0
  c(n) = 0.0
  if(n /= 3) then
     c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
     c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
     c(1) = c(1)*d(1)**2/(x(4)-x(1))
     c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
  end if
  !
  ! step 3: forward elimination
  !
  do i = 2, n
     h = d(i-1)/b(i-1)
     b(i) = b(i) - h*d(i-1)
     c(i) = c(i) - h*c(i-1)
  end do
  !
  ! step 4: back substitution
  !
  c(n) = c(n)/b(n)
  do j = 1, gap
     i = n-j
     c(i) = (c(i) - d(i)*c(i+1))/b(i)
  end do
  !
  ! step 5: compute spline coefficients
  !
  b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
  do i = 1, gap
     b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
     d(i) = (c(i+1) - c(i))/d(i)
     c(i) = 3.*c(i)
  end do
  c(n) = 3.0*c(n)
  d(n) = d(n-1)
  end subroutine make_spline
  function ispline(u, x, y, b, c, d, n)
  !======================================================================
  ! function ispline evaluates the cubic spline interpolation at point z
  ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  !----------------------------------------------------------------------
  ! input..
  ! u       = the abscissa at which the spline is to be evaluated
  ! x, y    = the arrays of given data points
  ! b, c, d = arrays of spline coefficients computed by spline
  ! n       = the number of data points
  ! output:
  ! ispline = interpolated value at point u
  !=======================================================================
  implicit none
  double precision ispline
  integer n
  double precision  u, x(n), y(n), b(n), c(n), d(n)
  integer i, j, k
  double precision dx

  ! if u is ouside the x() interval take a boundary value (left or right)
  if(u <= x(1)) then
     ispline = y(1)
     return
  end if
  if(u >= x(n)) then
     ispline = y(n)
     return
  end if

  !*
  !  binary search for for i, such that x(i) <= u <= x(i+1)
  !*
  i = 1
  j = n+1
  do while (j > i+1)
     k = (i+j)/2
     if(u < x(k)) then
        j=k
     else
        i=k
     end if
  end do
  !*
  !  evaluate spline interpolation
  !*
  dx = u - x(i)
  ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
  end function ispline
  subroutine Spline_new(this, n, x, y, out_value, ierr)
    type(Obj_Spline) :: this
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n)
    character(*), intent(in) :: out_value
    integer, intent(out) :: ierr
    
    this % n = n
    allocate(this%x(n), this%y(n), this%b(n), this%c(n), this%d(n))
    this % x = x
    this % y = y
    call make_spline(x, y, this%b, this%c, this%d, n)
    select case(out_value)
    case("edge")
       this%outvalue = 0
    case("const0")
       this%outvalue = 1
    case("error")
       this%outvalue = 2
    case default
       MSG_ERR("invalid option")
       ierr = 1
       write(0,*) "out_value:", out_value
       return
    end select
  end subroutine Spline_new
  subroutine Spline_new_file(this, fn, out_value, ierr)
    use Mod_sys, only : open_r
    type(Obj_Spline) :: this
    character(*), intent(in) :: fn
    character(*), intent(in) :: out_value
    integer, intent(out) :: ierr
    integer, parameter :: ifile = 292911
    integer num, i
    character(100) line
    double precision, allocatable :: xs(:), ys(:)

    ! -- count data num --
    call open_r(ifile, fn, ierr); CHK_ERR(ierr)
    read(ifile, *)
    num = 0
    do while(.true.)
       read(ifile, *, end=999) line
       if(len_trim(line).eq.0) goto 999
       num = num + 1
    end do
999 continue

    rewind(ifile)

 !   write(0,*) "spline.fn", filename

    ! -- read data --
    allocate(xs(num), ys(num))
    read(ifile, *)
    do i = 1, num
       read(ifile, *) xs(i), ys(i)
       !       write(0,*) xs(i),ys(i)
    end do
    close(ifile)

    ! -- build spline --
    call Spline_new(this, num, xs, ys, out_value, ierr); CHK_ERR(ierr)
    
  end subroutine Spline_new_file
  subroutine Spline_delete(this)
    type(Obj_Spline) :: this
    this % n = 0
    deallocate(this%x, this%y, this%b, this%c, this%d)
  end subroutine Spline_delete
  subroutine Spline_dump(this)
    type(Obj_Spline) :: this
    integer i
    write(*,*) "x, y, b, c, d"
    do i = 1, this % n
       write(*,*) this%x(i), this%y(i), this%b(i), this%c(i), this%d(i)
    end do
  end subroutine Spline_dump
  function search_index(this, u) result (res)
    type(Obj_Spline) :: this
    double precision, intent(in) :: u
    integer res
    integer i, j, k
    
    i = 1
    j = this%n+1
    do while (j > i+1)
       k = (i+j)/2
       if(u < this%x(k)) then
          j=k
       else
          i=k
       end if
    end do

    res = i
    
  end function search_index
  subroutine check_out_of_range(this, u, ierr)
    type(Obj_Spline) :: this
    double precision, intent(in) :: u
    integer, intent(out) :: ierr
    if(u < this%x(1) .or. u > this%x(this%n)) then
       MSG_ERR("out of range")
       ierr = 1
       write(0,*) "x0:", this%x(1)
       write(0,*) "xN:", this%x(this%n)
       write(0,*) "u:", u
       return 
    end if
  end subroutine check_out_of_range
  subroutine Spline_at(this, u, res, ierr)
    type(Obj_Spline) :: this
    double precision, intent(in) :: u
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    double precision dx
    integer i

    res = 0.0d0    
    if(u<this%x(1)) then
       select case(this%outvalue)
       case(0)
          res = this%y(1)
       case(1)
          res = 0
       case(2)
          call check_out_of_range(this, u, ierr); CHK_ERR(ierr)
       end select
    else if(u>this%x(this%n)) then
       select case(this%outvalue)
       case(0)
          res = this%y(this%n)
       case(1)
          res = 0
       case(2)
          call check_out_of_range(this, u, ierr); CHK_ERR(ierr)
       end select
    else
       i = search_index(this, u)
       dx = u - this%x(i)
       res = this%y(i) + dx*(this%b(i) + dx*(this%c(i) + dx*this%d(i)))
    end if
    
    !res = 0.0d0    
    !    call check_out_of_range(this, u); CHK_ERR()
    !i = search_index(this, u); CHK_ERR()
    !dx = u - this%x(i)
    !res = this%y(i) + dx*(this%b(i) + dx*(this%c(i) + dx*this%d(i)))
    
  end subroutine Spline_at
  subroutine Spline_deriv_at(this, n, u, res, ierr)
    type(Obj_Spline) :: this
    integer, intent(in) :: n
    double precision, intent(in) :: u
    double precision, intent(out) :: res
    integer, intent(out) :: ierr
    double precision dx
    integer i

    res = 0.0d0
    call check_out_of_range(this, u, ierr); CHK_ERR(ierr)

    i = search_index(this, u)
    dx = u - this%x(i)
    if(n.eq.0) then
       call Spline_at(this, u, res, ierr); CHK_ERR(ierr)
    else if(n.eq.1) then
       res = this%b(i) + dx*(2*this%c(i) + 3*dx*this%d(i))
    else if(n.eq.2) then
       res = 2*this%c(i) + 6*dx*this%d(i)
    else
       MSG_ERR("Spline_deriv_at supports only for n=0,1,2")
       ierr = 1; return       
    end if
    
  end subroutine Spline_deriv_at
end module Mod_spline
