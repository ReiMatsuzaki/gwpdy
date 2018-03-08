! store calculating timer and print as table

#include "macros.fpp"
module Mod_Timer
  implicit none
  private
  type Obj_Record
     character(100) :: key
     double precision :: tsum
     integer :: t0
     integer :: num
     logical :: act
  end type Obj_Record
  type Obj_Timer
     character(10) :: name
     logical :: print_key
     integer :: size, capacity
     type(Obj_Record), allocatable :: rec(:)
     integer t0     
  end type Obj_Timer  
  public Timer_new, Timer_delete, Timer_time, Timer_begin, Timer_end, Timer_result, Obj_Timer
contains
  subroutine init(this, num)
    type(Obj_Timer) :: this
    integer, intent(in) :: num
    integer i
    allocate(this%rec(num))
    do i = 1, num
       this%rec(i)%key = ""
       this%rec(i)%tsum = 0
       this%rec(i)%t0 = 0
       this%rec(i)%num = 0
       this%rec(i)%act = .false.
    end do
  end subroutine init
  subroutine extend(this)
    type(Obj_Timer), intent(inout) :: this
    integer num0, num1
    type(Obj_Record), allocatable :: tmp(:)
    
    num0 = this%capacity
    num1 = 2*num0

    allocate(tmp(num0))
    tmp(:) = this%rec(:)

    deallocate(this%rec)
    call init(this, num1)
    this%rec(:num0) = tmp(:)
    
    deallocate(tmp)
    this%capacity = num1
    
  end subroutine extend
  function find(this, key) result(res)
    type(Obj_Timer) :: this
    character(*), intent(in) :: key
    integer :: res
    integer i

    res = 0
    do i = 1, this%size
       if(trim(this%rec(i)%key)==trim(key)) then
          res = i
          return
       end if
    end do
    
  end function find
  function Timer_time(this, key) result(res)
    type(Obj_Timer) :: this
    character(*), intent(in) :: key
    double precision res
    integer idx

    idx = find(this, trim(key))

    res = this%rec(idx)%tsum
    
  end function Timer_time
  subroutine Timer_new(this, name, print_key, ierr)
    type(Obj_Timer) :: this
    character(*), intent(in) :: name
    logical, intent(in) :: print_key
    integer, intent(out) :: ierr
    integer, parameter :: NUM0=1024
    ierr = 0
    this%name = name
    this%size = 0
    this%capacity = NUM0
    this%print_key = print_key
    call system_clock(this%t0)
    call init(this, NUM0)
  end subroutine Timer_new
  subroutine Timer_delete(this, ierr)
    type(Obj_Timer) :: this
    integer, intent(out) :: ierr
    ierr = 0
    deallocate(this%rec)
  end subroutine Timer_delete
  subroutine Timer_begin(this, key, ierr)
    type(Obj_Timer) :: this
    character(*), intent(in) :: key
    integer, intent(out) :: ierr
    integer idx
    
    if(this%print_key) then       
       call timestamp(key // " begin")
    end if
    idx = find(this, key)
    if(idx==0) then
       if(this%size+1 > this%capacity) then
          call extend(this)
       end if       
       this%size = this%size+1
       this%rec(this%size)%key = key
       this%rec(this%size)%tsum = 0
       call system_clock(this%rec(this%size)%t0)
       this%rec(this%size)%num = 0
       this%rec(this%size)%act = .true.
    else
       if(this%rec(idx)%act) then
          MSG_ERR("Timer_start is called for active timer key.")
          ierr = 1
          write(0,*) "key:", key
          return
       end if
       call system_clock(this%rec(idx)%t0)
       this%rec(idx)%act = .true.
    end if
  end subroutine Timer_begin
  subroutine Timer_end(this, key, ierr)
    type(Obj_Timer) :: this
    character(*), intent(in) :: key
    integer, intent(out) :: ierr
    integer idx
    integer dt, t0, t1, t_rate, t_max

    idx = find(this, key)
    if(idx==0) then
       MSG_ERR("Timer_end is called for not started record")
       ierr = 1
       write(0,*) "key:", key
       return
    else
       if(.not.this%rec(idx)%act) then
          MSG_ERR("Timer_end is called for non active started record")
          ierr = 1
          write(0,*) "key:", key
          return
       end if
       t0 = this%rec(idx)%t0
       call system_clock(t1, t_rate, t_max)
       if(t1<t0) then
          dt = (t_max-t0) + t1 + 1
       else
          dt = t1-t0
       end if
       this%rec(idx)%tsum = this%rec(idx)%tsum + (1.0d0*dt)/t_rate
       this%rec(idx)%num  = this%rec(idx)%num+1
       this%rec(idx)%act  = .false.
    end if

    if(this%print_key) then
       call timestamp(key // " end")
    end if
    
  end subroutine Timer_end
  subroutine Timer_result(this, ifile, ierr)
    type(Obj_Timer) :: this    
    integer, intent(in) :: ifile
    integer, intent(out) :: ierr
    integer i
    integer t1, t_rate, t_max, dt
    
    write(ifile,'(A,":Time_result begin")') trim(this%name)
    call system_clock(t1, t_rate, t_max)

    if(t1<this%t0) then
       dt = (t_max-this%t0) + t1 + 1
    else
       dt = t1-this%t0
    end if
    write(ifile,*) "total time;", (1.0d0*dt)/t_rate
       
    do i = 1, this%size
       if(this%rec(i)%act) then
          MSG_ERR("active record is found.")
          ierr = 1
          write(0,*) "key:", this%rec(i)%key
          return
       end if
       write(ifile,'(A30, I10, f20.5)') trim(this%rec(i)%key), this%rec(i)%num, this%rec(i)%tsum
    end do
    write(ifile,'(A,":Time_result end")') trim(this%name)
  end subroutine Timer_result
  subroutine pad0(x, res)
    integer, intent(in) :: x
    character(2), intent(out) :: res
    if(x>10) then
       write(res, '(i1, i1)') x/10, mod(x, 10)
    else
       write(res, '("0", i1)') mod(x, 10)
    end if
  end subroutine pad0
  subroutine timestamp(message)
    character(*), intent(in) :: message
    integer :: date_time(8), h, m, s
    character(10) :: sys_time(3)
    character(8) :: tlabel
    character(2) :: hh, mm, ss
    call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    h = date_time(5)
    m = date_time(6)
    s = date_time(7)
    tlabel = "00:00:00"
    call pad0(s, ss)
    call pad0(m, mm)
    call pad0(h, hh)
    tlabel = hh // ":" // mm // ":" // ss
    write(*,'(A, 2x, A, 2x, A)') trim(sys_time(1)), trim(tlabel), trim(message)
  end subroutine timestamp
end module Mod_timer
