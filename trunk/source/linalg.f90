module cFromFortran
use,intrinsic  :: iso_c_binding
  implicit none
  private

  public :: ComplexDotProdv1
  public :: DoubleGEMMv1
  public :: ComplexGEMMv1
  public :: DiagonaliseSymmetricalv1
  public :: DiagonaliseHermitianv1
  public :: SolveSystemEqs

  type, public, bind(c) :: DoubleMatrix
    integer(c_int) :: iRows,iCols
    type(c_ptr) :: a
  end type DoubleMatrix
  type, public, bind(c) :: ComplexMatrix
    integer(c_int) :: iRows,iCols
    type(c_ptr) :: a
  end type ComplexMatrix

  integer, parameter, public :: kpr=kind(1.0d0)
  real(kpr), parameter, public :: ktol=1.0d-10

contains


  subroutine ComplexDotProdv1(c,n,a,inca,b,incb) bind(c, name="ComplexDotProdv1")
    integer(c_int), intent(in) :: n,inca,incb
    type(c_ptr), intent(in) :: a,b
    complex(c_double_complex), pointer :: ca(:)=>null(),cb(:)=>null()
    complex(c_double_complex),intent(inout) :: c
    complex(kpr), external :: zdotc


    call c_f_pointer(a,ca,(/n/))
    call c_f_pointer(b,cb,(/n/))
    c = zdotc(n, ca, inca, cb, incb )
    nullify(ca,cb)
  end  subroutine ComplexDotProdv1


  subroutine DoubleGEMMv1(transa,transb,alpha,a,b,beta,c) bind(c, name="DoubleGEMMv1")
    character(len=*), parameter :: sMyName="DoubleGEMMv1"
    type(c_ptr),intent(in) :: a,b
    type(c_ptr), intent(inout) :: c
    real(c_double), intent(in) :: alpha,beta
    character(kind=c_char), intent(in) :: transa,transb

    type(DoubleMatrix), pointer :: sa=>null(),sb=>null(),sc=>null()
    type(c_ptr),pointer :: pa=>null(),pb=>null(),pc=>null()
    real(c_double), pointer :: ca(:,:)=>null(),cb(:,:)=>null(),cc(:,:)=>null()


!     call c_f_pointer(a,sa,(/1/))
!     call c_f_pointer(sa%a,pa,(/1/))
! this is an  workaround as specifying the size for a scalar would generate seg fault
     call c_f_pointer(a,sa)
     call c_f_pointer(sa%a,pa)

    call c_f_pointer(pa,ca,(/sa%iCols,sa%iRows/))
!     write(*,'(a,3Z)')"sa,pa,ca ",transfer(c_loc(sa),0_C_INTPTR_T),transfer(c_loc(pa),0_C_INTPTR_T), transfer(c_loc(ca),0_C_INTPTR_T)

    call c_f_pointer(b,sb)
    call c_f_pointer(sb%a,pb)
    call c_f_pointer(pb,cb,(/sb%iCols,sb%iRows/))
!     write(*,'(a,3Z)')"sb,pb,cb ",transfer(c_loc(sb),0_C_INTPTR_T),transfer(c_loc(pb),0_C_INTPTR_T), transfer(c_loc(cb),0_C_INTPTR_T)

    call c_f_pointer(c,sc)
    call c_f_pointer(sc%a,pc)
    call c_f_pointer(pc,cc,(/sc%iCols,sc%iRows/))
!     write(*,'(a,3Z)')"sc,pc,cc ",transfer(c_loc(sc),0_C_INTPTR_T),transfer(c_loc(pc),0_C_INTPTR_T), transfer(c_loc(cc),0_C_INTPTR_T)

    if ((sa%iCols == sb%iRows) .and. (sc%iRows == sa%iRows) .and. (sc%iCols == sb%iCols)) then
       call dgemm(transa, transb, sb%iCols, sa%iRows, sb%iRows, alpha, cb, sb%iCols, ca, sa%iCols, beta, cc, sc%iCols)
    else
      write(*,'(a)')trim(sMyName)//"The dimensions of matrices prevent multiplication. Please check them."
      stop
    endif
     nullify(sa,pa,ca,sb,pb,cb,sc,pc,cc)
  end subroutine DoubleGEMMv1

  subroutine ComplexGEMMv1(transa,transb,alpha,a,b,beta,c) bind(c, name="ComplexGEMMv1")
    character(len=*), parameter :: sMyName="ComplexGEMMv1"
    type(c_ptr),intent(in) :: a,b
    type(c_ptr), intent(inout) :: c
    complex(c_double_complex), intent(in) :: alpha,beta
    character(kind=c_char), intent(in) :: transa,transb

    type(ComplexMatrix), pointer :: sa=>null(),sb=>null(),sc=>null()
    type(c_ptr),pointer :: pa=>null(),pb=>null(),pc=>null()
    complex(c_double_complex), pointer :: ca(:,:)=>null(),cb(:,:)=>null(),cc(:,:)=>null()
    character(len=1):: ta,tb


    call c_f_pointer(a,sa)
    call c_f_pointer(sa%a,pa)
    call c_f_pointer(pa,ca,(/sa%iCols,sa%iRows/))
!     write(*,'(a,3Z)')"sa,pa,ca ",transfer(c_loc(sa),0_C_INTPTR_T),transfer(c_loc(pa),0_C_INTPTR_T), transfer(c_loc(ca),0_C_INTPTR_T)

    call c_f_pointer(b,sb)
    call c_f_pointer(sb%a,pb)
    call c_f_pointer(pb,cb,(/sb%iCols,sb%iRows/))
!     write(*,'(a,3Z)')"sb,pb,cb ",transfer(c_loc(sb),0_C_INTPTR_T),transfer(c_loc(pb),0_C_INTPTR_T), transfer(c_loc(cb),0_C_INTPTR_T)

    call c_f_pointer(c,sc)
    call c_f_pointer(sc%a,pc)
    call c_f_pointer(pc,cc,(/sc%iCols,sc%iRows/))
!     write(*,'(a,3Z)')"sc,pc,cc ",transfer(c_loc(sc),0_C_INTPTR_T),transfer(c_loc(pc),0_C_INTPTR_T), transfer(c_loc(cc),0_C_INTPTR_T)

    if ((sa%iCols == sb%iRows) .and. (sc%iRows == sa%iRows) .and. (sc%iCols == sb%iCols)) then
      call zgemm(transa, transb, sb%iCols, sa%iRows, sb%iRows, alpha, cb, sb%iCols, ca, sa%iCols, beta, cc, sc%iCols)
    else
      write(*,'(a)')trim(sMyName)//"The dimensions of matrices prevent multiplication. Please check them."
      stop
    endif
    nullify(sa,pa,ca,sb,pb,cb,sc,pc,cc)
  end subroutine ComplexGEMMv1



  ! info=-400 means not enough memory for work array
! info=-401 means error deallocating array work
! ! info=-402 means not enough memory for issupz array
! info=-403 means error deallocating array isspux
 ! info=-404 means not enough memory for iwork array
! info=-405 means error deallocating array iwork
  subroutine DiagonaliseSymmetricalv1(a,w,eig,job,uplo,info) bind(c,name="DiagonaliseSymmetricalv1")
    character(len=*),parameter :: sMyName="DiagonaliseSymmetricalv1"
    type(c_ptr), intent(inout) :: a,eig
    type(c_ptr), intent(inout) :: w
    character(kind=c_char,len=1), intent(in) :: uplo,job
    integer(c_int) :: info

    type(DoubleMatrix), pointer :: ss=>null(),se=>null()
    type(c_ptr),pointer :: ps=>null(),pe=>null()
    real(c_double), pointer :: cs(:,:)=>null(),ce(:,:)=>null()
    real(c_double), pointer :: cw(:)=>null()
    real(kpr) :: tmp(1)
    real(kpr), allocatable :: work(:)
    integer :: lwork,liwork,ierror
    character(len=1) :: ta,rnge
    real(kpr) :: vl,vu,abstol
    integer :: il,iu,itmp(1),m
    integer, allocatable :: iwork(:),issupz(:)

    abstol=ktol
    if ((uplo=='U') .or. (uplo=='u')) then
      ta="L"
    else
      ta="U"
    endif
    rnge="A"

    call c_f_pointer(a,ss)
    call c_f_pointer(ss%a,ps)
    call c_f_pointer(ps,cs,(/ss%iCols,ss%iRows/))

    call c_f_pointer(eig,se)
    call c_f_pointer(se%a,pe)
    call c_f_pointer(pe,ce,(/se%iCols,se%iRows/))

!     write(*,'(a,3Z)')'ss,ps,cs ',transfer(c_loc(ss),0_C_INTPTR_T),transfer(c_loc(ps),0_C_INTPTR_T), transfer(c_loc(cs),0_C_INTPTR_T)
!
    call c_f_pointer(w,cw,(/ss%iRows/))
!     write(*,'(a,Z)')"cw", transfer(c_loc(cw),0_C_INTPTR_T)
    allocate(issupz(2*ss%iCols),stat=ierror)
    if (ierror /= 0) then
        info = -402
        return
     endif
    lwork=-1
    liwork=-1
    call dsyevr(job,rnge, ta, ss%iCols, cs, ss%iCols,vl,vu,il,iu,abstol,m, cw,&
       ce,ss%iCols,issupz,tmp, lwork,itmp,liwork, info)
    lwork=int(tmp(1))
    liwork=itmp(1)
     allocate(work(lwork),stat=ierror)
     if (ierror /= 0) then
        info = -400
        return
     endif
     allocate(iwork(liwork),stat=ierror)
     if (ierror /= 0) then
        info = -404
        return
     endif
     call dsyevr(job,rnge, ta, ss%iCols, cs, ss%iCols,vl,vu,il,iu,abstol,m, cw, ce,se%iCols,issupz,work, lwork,iwork,liwork, info)
     deallocate(work,stat=ierror)
     if (ierror /= 0) then
        info = -401
        return
     endif
     deallocate(issupz,stat=ierror)
     if (ierror /= 0) then
        info = -403
        return
     endif
     deallocate(iwork,stat=ierror)
     if (ierror /= 0) then
        info = -405
        return
     endif
     nullify(ss,ps,cs,se,pe,ce,cw)
  end subroutine DiagonaliseSymmetricalv1


    ! info=-400 means not enough memory for work array
! info=-401 means error deallocating array work
! ! info=-402 means not enough memory for issupz array
! info=-403 means error deallocating array isspux
 ! info=-404 means not enough memory for iwork array
! info=-405 means error deallocating array iwork
  subroutine DiagonaliseHermitianv1(a,w,eig,job,uplo,info) bind(c,name="DiagonaliseHermitianv1")
    character(len=*),parameter :: sMyName="DiagonaliseHermitianv1"
    type(c_ptr), intent(inout) :: a,eig
    type(c_ptr), intent(inout) :: w
    character(kind=c_char,len=1), intent(in) :: uplo,job
    integer(c_int) :: info

    type(DoubleMatrix), pointer :: ss=>null(),se=>null()
    type(c_ptr),pointer :: ps=>null(),pe=>null()
    real(c_double_complex), pointer :: cs(:,:)=>null(),ce(:,:)=>null()
    real(c_double), pointer :: cw(:)=>null()
    complex(kpr) :: wtmp(1)
    complex(kpr), allocatable :: work(:)
    integer :: lwork,liwork,lrwork,ierror
    character(len=1) :: ta,rnge
    real(kpr) :: vl,vu,abstol
    integer :: il,iu,itmp(1),m
    integer, allocatable :: iwork(:),issupz(:)
    real(kpr) :: rtmp(1)
    real(kpr), allocatable :: rwork(:)


    abstol=ktol
    if ((uplo=='U') .or. (uplo=='u')) then
      ta="L"
    else
      ta="U"
    endif
    rnge="A"

    call c_f_pointer(a,ss)
    call c_f_pointer(ss%a,ps)
    call c_f_pointer(ps,cs,(/ss%iCols,ss%iRows/))

    call c_f_pointer(eig,se)
    call c_f_pointer(se%a,pe)
    call c_f_pointer(pe,ce,(/se%iCols,se%iRows/))

!     write(*,'(a,3Z)')'ss,ps,cs ',transfer(c_loc(ss),0_C_INTPTR_T),transfer(c_loc(ps),0_C_INTPTR_T), transfer(c_loc(cs),0_C_INTPTR_T)
!
    call c_f_pointer(w,cw,(/ss%iRows/))
!     write(*,'(a,Z)')"cw", transfer(c_loc(cw),0_C_INTPTR_T)
    allocate(issupz(2*ss%iCols),stat=ierror)
    if (ierror /= 0) then
        info = -402
        return
     endif
    lwork=-1
    liwork=-1
    lrwork=-1
    call zheevr(job,rnge,ta, ss%iCols, cs, ss%iCols,vl,vu,il,iu,abstol,m, cw,&
       ce,ss%iCols,issupz,wtmp, lwork,rtmp,lrwork,itmp,liwork, info)
    lwork=int(wtmp(1))
    liwork=itmp(1)
    lrwork=int(rtmp(1))
     allocate(work(lwork),stat=ierror)
     if (ierror /= 0) then
        info = -400
        return
     endif
     allocate(iwork(liwork),stat=ierror)
     if (ierror /= 0) then
        info = -404
        return
     endif
     allocate(rwork(lrwork),stat=ierror)
     if (ierror /= 0) then
        info = -406
        return
     endif
    call zheevr(job,rnge,ta, ss%iCols, cs, ss%iCols,vl,vu,il,iu,abstol,m, cw,&
       ce,ss%iCols,issupz,work, lwork,rwork,lrwork,iwork,liwork, info)
     deallocate(work,stat=ierror)
     if (ierror /= 0) then
        info = -401
        return
     endif
     deallocate(issupz,stat=ierror)
     if (ierror /= 0) then
        info = -403
        return
     endif
     deallocate(iwork,stat=ierror)
     if (ierror /= 0) then
        info = -405
        return
     endif
     deallocate(rwork,stat=ierror)
     if (ierror /= 0) then
        info = -407
        return
     endif
     nullify(ss,ps,cs,se,pe,ce,cw)
  end subroutine DiagonaliseHermitianv1


!> \brief The routine solves for X the system of linear equations A*X = B,
!> \details  where A is an n-by-n matrix, the columns of matrix B are individual right-hand sides, and the columns of X are the corresponding solutions.
!> \remarks info=-400  error allocating temporary buffers
!> info=-401  error deallocating temporary buffers
  subroutine SolveSystemEqs(a,b,info) bind(c,name="SolveSystemEqs")
    character(len=*),parameter :: sMyName="SolveSystemEqs"
    type(c_ptr), intent(inout) :: a,b
    integer(c_int) :: info

    type(DoubleMatrix), pointer :: sa=>null()
    type(c_ptr),pointer :: pa=>null()
    real(c_double), pointer :: ca(:,:)=>null()
    real(c_double), pointer :: cb(:)=>null()
    real(kpr), allocatable :: work(:),x(:)
    real, allocatable :: swork(:)
    integer, allocatable :: ipiv(:)
    integer :: nrhs,iter,ierror
    integer :: i,j

    nrhs=1
    call c_f_pointer(a,sa,(/1/))
    call c_f_pointer(sa%a,pa,(/1/))
    call c_f_pointer(pa,ca,(/sa%iCols,sa%iRows/))

    call c_f_pointer(b,cb,(/sa%iCols/))

    allocate(work(sa%iCols),swork(sa%iCols),ipiv(sa%iCols),x(sa%iCols),stat=ierror)
    if (ierror /= 0) then
        info = -400
        return
    endif
!    call dsgesv(sa%iCols,nrhs,ca,sa%iCols,ipiv,cb,sa%iCols,x,sa%iCols,work,swork,iter,info)
    cb=x

    deallocate(work,swork,ipiv,x,stat=ierror)
    if (ierror /= 0) then
       info = -401
       return
    endif
!     print *,"hello"
    nullify(pa,sa,cb)
  end subroutine SolveSystemEqs

end module cFromFortran
