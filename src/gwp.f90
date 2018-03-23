#include "macros.fpp"
module Mod_GWP
  ! Gauss wave packet sets for arbitary dimension.
  !    G_i(Q) = N exp[ -(Q-R)^Tg(Q-R) + iP(Q-R) + ic ]
  ! nf : positive inter
  !      number of freedom
  ! exp_type : character
  !     m=> g is matrix, d=> g is diagonal, c=> g is (constant)x(id)
  ! N,c : real
  ! g   : [[complex]]
  ! R,P : [real]
  implicit none
  type Obj_GWP
     integer :: nf
     integer :: num
     character :: exp_type  
     complex(kind(0d0)), allocatable :: g(:,:,:)
     double precision ,  allocatable :: N(:), c(:), R(:,:), P(:,:)
  end type OBJ_GWP
contains
  ! ==== Constructors ====
  subroutine GWP_new(this, nf, num, exp_type, ierr)
    type(Obj_GWP), intent(out) :: this
    integer, intent(in) :: nf, num
    character, intent(in) :: exp_type
    integer, intent(out) :: ierr
    integer k
    ierr = 0

    if(exp_type.ne.'m'.and.exp_type.eq.'d'.and.exp_type.eq.'c') then       
       MSG_ERR("Invalid argument")
       ierr = 1
       write(1,*) "exp_type must be m,d,u"
       write(*,*) "exp_type:", exp_type
       return
    end if
    
    allocate(this%g(num,nf,nf))
    allocate(this%N(num))
    allocate(this%c(num))
    allocate(this%R(num,nf))
    allocate(this%P(num,nf))
    
    this%nf  = nf
    this%num = num
    this%exp_type = exp_type
    this%g = 0
    this%N = 0
    this%c = 0
    this%R = 0
    this%P = 0

    do k = 1, nf
       this%g(:,k,k) = 1
    end do
    
  end subroutine GWP_new
  subroutine GWP_setup(this, ierr)
    use Mod_const, only : PI
    type(Obj_GWP), intent(inout) :: this
    integer, intent(out) :: ierr
    
    double precision, parameter :: tol = 1.0d-10    
    integer i, j, a
    double precision :: tmp
    
    ierr = 0
    if(this%exp_type.eq."g") then
       ierr = 1
       write(1,*) "exp_type=g is not supported now."
       return
    end if
    if(this%exp_type.eq."d" .or. this%exp_type.eq."c") then
       do i = 1, this%nf
          do j = 1, this%nf
             if(i.ne.j) then
                if(maxval(abs(this%g(:,i,j))) > tol) then
                   MSG_ERR("invalid g")
                   ierr = 1; return                   
                end if
             else
                if( minval(real(this%g(:,i,i))) < tol ) then
                   MSG_ERR("invalid g")
                   ierr = 1; return                   
                end if
             end if
          end do
       end do
    end if
    if(this%exp_type.eq."c") then
       do i = 1, this%nf
          do a = 1, this%num
             if(abs(this%g(a,i,i)-this%g(a,1,1)) > tol) then
                ierr = 1
                write(*,*) "invalid g for exp_type=c."
                return
             end if
          end do
       end do
    end if

    do a = 1, this%num
       tmp = 1
       do i = 1, this%nf
          tmp = tmp * (2*real(this%g(a,i,i))/PI)**0.25d0
       end do
       this%N(a) = tmp
    end do
    
  end subroutine GWP_setup
  subroutine GWP_delete(this, ierr)
    type(Obj_GWP), intent(inout) :: this
    integer, intent(out) :: ierr
    ierr = 1
    deallocate(this%g)
    deallocate(this%N)
    deallocate(this%c)
    deallocate(this%R)
    deallocate(this%P)
    ierr = 0
  end subroutine GWP_delete
  ! ==== Calculation ====
  subroutine GWP_overlap(this, res, ierr)
    ! grive the overlap matrix of same width parameter GWP sets
    ! {g(i), R(i,:), P(i,:), theta(i) } : GWP
    use Mod_const, only : ii
    type(Obj_GWP), intent(in) :: this
    complex(kind(0d0)), intent(out) :: res(:,:)
    integer, intent(out) :: ierr
    double precision:: dR, dP, PP
    complex(kind(0d0)) :: g
    integer i, j, k
    complex(kind(0d0)) tmp

    ierr = 0

    if(this%exp_type.ne."d" .and. this%exp_type.ne.'c') then
       MSG_ERR("only exp_type==c,d is supported")
       write(0,*) "exp_type:", this%exp_type
       ierr = 1; return       
    end if

    if(size(res,1)<this%num .or. size(res,2)<this%num) then
       MSG_ERR("invalid size")
       write(0,*) "shape(res):", shape(res)
       write(0,*) "this%num:  ", this%num
       ierr = 1; return
    end if

    do i = 1, this%num
       do j = 1, this%num
          tmp = 0
          do k = 1, this%nf
             g = this%g(i,k,k)
             dR = this%R(i,k) - this%R(j,k)
             dP = this%P(i,k) - this%P(j,k)
             PP = this%P(i,k) + this%P(j,k)             
             tmp = tmp -g/2*dR**2 - 1/(4*g)*dP**2 + ii/2*dR*PP
          end do
          res(i,j) = exp(tmp -ii*this%c(i) +ii*this%c(j))
       end do
    end do
    
  end subroutine GWP_overlap
  subroutine GWP_p2(this, res, ierr)
    ! grive square of momentum operator of same width GWP sets
    ! {g(i), R(i,:), P(i,:), theta(i) } : GWP
    use Mod_const, only : ii
    type(Obj_GWP), intent(in) :: this
    complex(kind(0d0)), intent(out) :: res(:,:,:)
    integer, intent(out) :: ierr
    integer i, j, k
    complex(kind(0d0)) tmp
    double precision dR, PP
    complex(kind(0d0)) g
    complex(kind(0d0)) :: S(this%num, this%num)
    
    ierr = 0

    if(this%exp_type.ne."d" .and. this%exp_type.ne."c") then
       MSG_ERR("only exp_type==c,d is supported")
       return
    end if

    call GWP_overlap(this, S, ierr); CHK_ERR(ierr)

    do i = 1, this%num
       do j = 1, this%num
          tmp = 1
          do k = 1, this%nf
             g = this%g(i,k,k)
             dR = this%R(i,k) - this%R(j,k)
             PP = this%P(i,k) + this%P(j,k)
             tmp = tmp * ( (PP/2)**2 + g - g*g*dR**2 - II*g*dR*PP)
             res(k,i,j) = tmp * S(i,j) 
          end do
       end do
    end do
    
  end subroutine GWP_p2
  subroutine GWP_r1(this, res, ierr)
    use Mod_const, only : ii
    type(Obj_GWP), intent(in) :: this
    complex(kind(0d0)), intent(out) :: res(:,:,:)
    integer, intent(out) :: ierr
    integer i, j, k
    ierr = 0

    if(this%exp_type.ne."d" .and. this%exp_type.ne."c") then
       MSG_ERR("only exp_type==c,d is supported")
       return
    end if

    do i = 1, this%num
       do j = 1, this%num
          do k = 1, this%nf
             res(k,i,j) = 1
          end do
       end do
    end do
  end subroutine GWP_r1
end module Mod_GWP
