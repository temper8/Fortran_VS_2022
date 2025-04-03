module MKL_helper
    use, intrinsic :: iso_fortran_env,    only: real32,real64,real128
    
    implicit none
   
    integer,parameter :: I4B = selected_int_kind(9)
    !integer,parameter,public :: wp = real32     !! single precision reals
    integer,parameter,public :: DPT = real64  
    real(DPT),parameter ::  ZERO = 0.0_DPT
    real(DPT),parameter ::  HALF = 0.5_DPT
    real(DPT),parameter ::   ONE = 1.0_DPT
    real(DPT),parameter ::   TWO = 2.0_DPT
    real(DPT),parameter :: THREE = 3.0_DPT
    real(DPT),parameter ::  FOUR = 4.0_DPT
    complex(DPT),parameter :: ZEROC = (ZERO,ZERO)
    complex(DPT),parameter ::  ONEC = (ONE,ZERO)
    
    contains
    
    function sys_time()
       implicit none
       real(8) sys_time
       integer count, count_rate, count_max
       call system_clock(count, count_rate, count_max)
       sys_time = count*1.0/count_rate
       return
    end   
    
    function matmul_MKL(a,b) result(c)
        integer n
	    complex(DPT),  intent(IN)  :: a(:,:), b(:,:)
        complex(DPT)  :: c(size(a,1),size(a,1))
        complex(DPT) alpha
        !print *, size(a,1), size(a,2)
        !print *, size(b,2), size(b,2)
        alpha = 1.0
        n = size(a,1)
	    c=0.d0
        call zgemm('N', 'N', n, n, n, alpha, a, n, b, n, 0.0d0, c, n)
    end   

    subroutine inverse_matrix_mkl(A)  
        implicit none
	    complex(DPT)  :: A(:,:)
        complex(DPT), allocatable :: work(:)
        integer :: ipiv(size(A,1))           ! Массив перестановок
        integer :: info, lwork, i
        integer n
    
        n = size(A,1)

        ! Шаг 1: LU-факторизация
        call zgetrf(n, n, A, n, ipiv, info)
        if (info /= 0) then
            print *, 'Ошибка LU-разложения. Код:', info
            stop
        end if

        ! Шаг 2: Определение оптимального размера рабочего массива
        lwork = -1
        allocate(work(1))
        call zgetri(n, A, n, ipiv, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! Шаг 3: Вычисление обратной матрицы
        call zgetri(n, A, n, ipiv, work, lwork, info)
        if (info /= 0) then
            print *, 'Ошибка обращения матрицы. Код:', info
            stop
        end if

        deallocate(work)
    end subroutine 
    
    subroutine MatInv(Amat)
        !use Constants, only: ZEROC,ONEC
        !use ElapsedTime, only: InitElapsedTime,PrintElapsedTime
        implicit none
        complex(DPT),dimension(:,:),intent(inout)::Amat
        integer(I4B)::j,Nmax
        integer(I4B),dimension(size(Amat,1))::Ipp
        real(DPT)::d
        complex(DPT),dimension(size(Amat,1),size(Amat,2))::Cmat
        ! ----------------------------------------------------------------------
        Nmax=size(Amat,1)
        Cmat=ZEROC ; forall (j=1:Nmax) Cmat(j,j)=ONEC
        !call InitElapsedTime (4)
        call ludcmp7_z (Amat,Ipp,d)
        !call ludcmp9_z (Amat,Ipp,d)
        !call PrintElapsedTime (4,'             LU decomposition:          ')
        do j=1,Nmax
         call lubksb_z (Amat,Ipp,Cmat(:,j))
        end do
        !call PrintElapsedTime (4,'             LU backsubstitution:       ')
        Amat=Cmat
    end subroutine MatInv
    
    
    function assert_eq4(n1,n2,n3,n4,string)
        implicit none
        character(len=*),intent(in) :: string
        integer,intent(in) :: n1,n2,n3,n4
        integer :: assert_eq4
        character(len=*),parameter::prname='assert_eq4'
        character(len=255)::errmessage
        ! ----------------------------------------------------------------------
        if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
         assert_eq4=n1
        else
         write(errmessage,'(a)') 'assertion failed in: '//string
         print *,  prname, trim(errmessage)
        end if
    end function assert_eq4
    
    subroutine lubksb_z(a,indx,b)
        !use Constants, only: ZERO
        !use LowLevLib, only: assert_eq
        implicit none
        character(len=*),parameter::prname='lubksb_z'
        complex(DPT),dimension(:,:),intent(in)::a
        integer(I4B),dimension(:),intent(in)::indx
        complex(DPT),dimension(:),intent(inout)::b
        integer(I4B)::i,n,ii,ll
        complex(DPT)::summ
        ! ----------------------------------------------------------------------
        n=assert_eq4(size(a,1),size(a,2),size(b),size(indx),prname)
        ii=0
        do i=1,n
         ll=indx(i)
         summ=b(ll)
         b(ll)=b(i)
         if (ii /= 0) then
          summ=summ-sum(a(i,ii:i-1)*b(ii:i-1)) ! dot_product would require 2 excessive conjg ops
         else if (abs(summ) /= ZERO) then
          ii=i
         end if
         b(i)=summ
        end do
        do i=n,1,-1
         b(i) = (b(i)-sum(a(i,i+1:n)*b(i+1:n)))/a(i,i)
        end do
    end subroutine lubksb_z
    
    function assert_eq3(n1,n2,n3,string)
        implicit none
        character(len=*),intent(in) :: string
        integer,intent(in) :: n1,n2,n3
        integer :: assert_eq3
        character(len=*),parameter::prname='assert_eq3'
        character(len=255)::errmessage
        ! ----------------------------------------------------------------------
        if (n1 == n2 .and. n2 == n3) then
         assert_eq3=n1
        else
         write(errmessage,'(a)') 'assertion failed in: '//string
         print *,  prname, trim(errmessage)
        end if
    end function assert_eq3
    
    subroutine swap_zv(a,b)
        implicit none
        complex(DPT),dimension(:),intent(inout) :: a,b
        complex(DPT),dimension(size(a)) :: dum
        ! ----------------------------------------------------------------------
        dum=a
        a=b
        b=dum
    end subroutine swap_zv
    
    subroutine ludcmp7_z (a,indx,d)
        !use Constants, only: ZERO,ONE
        !use LowLevLib, only: assert_eq,swap
        implicit none
        character(len=*),parameter::prname='ludcmp7_z'
        complex(DPT),dimension(:,:),intent(inout)::a
        integer(I4B),dimension(:),intent(out)::indx
        real(DPT),intent(out)::d
        integer(I4B)::n,i,j,imax
        real(DPT),dimension(size(a,1))::vv
        ! ----------------------------------------------------------------------
        n=assert_eq3(size(a,1),size(a,2),size(indx),prname)
        d=ONE
        vv=maxval(abs(a),dim=2)
        if (any(vv == ZERO)) print *, prname,'singular matrix'
        vv=ONE/vv
        do j=1,n
         do i=2,j-1
          a(i,j)=a(i,j)-sum(a(i,1:i-1)*a(1:i-1,j))
         end do
         do i=j,n
          a(i,j)=a(i,j)-sum(a(i,1:j-1)*a(1:j-1,j))
         end do
         imax=(j-1)+sum(maxloc(vv(j:n)*abs(a(j:n,j)))) ! sum used to get a scalar from maxloc(1)
         if (j /= imax) then
          call swap_zv(a(imax,:),a(j,:))
          d=-d
          vv(imax)=vv(j)
         end if
         indx(j)=imax
         if (abs(a(j,j)) == ZERO) print *, prname, 'singular matrix'
         a(j+1:n,j)=a(j+1:n,j)/a(j,j)
        end do
    end subroutine ludcmp7_z
    
end module MKL_helper