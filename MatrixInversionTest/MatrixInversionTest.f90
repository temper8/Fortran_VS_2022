!  MatrixInversionTest.f90 
!
!  FUNCTIONS:
!  MatrixInversionTest - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: MatrixInversionTest
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
    module test
        use MKL_helper
    contains
        subroutine random_complex(matrix)

            implicit none
            complex(DPT), intent(out) :: matrix(:,:)  ! Выходная комплексная матрица
            real(DPT) :: real_part(size(matrix,1), size(matrix,2))  ! Временный массив для действительной части
            real(DPT) :: imag_part(size(matrix,1), size(matrix,2))  ! Временный массив для мнимой части

            ! Генерируем случайные значения для действительной и мнимой частей
            call random_number(real_part)
            call random_number(imag_part)

            ! Объединяем части в комплексную матрицу
            matrix = cmplx(real_part, imag_part)
        end subroutine random_complex
    
        function calc_error(zmat) result(err)
            implicit none
            integer n
            real(DPT) err
            complex(DPT), dimension(:,:) :: zmat
            real,    dimension(size(zmat,1),size(zmat,1)) :: m
            
            m = abs(zmat)
            err = MAXVAL(m)
            
        end function calc_error
        
        function mkl_mat_inv_test(n) result(err)
            implicit none
            integer n
            real(DPT) err
            complex(DPT), dimension(n,n) :: a, b, c
            
            call random_complex(a)
            b = a
            call inverse_matrix_mkl(b)
            call inverse_matrix_mkl(b)
            c = a - b
            err = calc_error(c)
        end function mkl_mat_inv_test
    
        function mat_inv_test(n) result(err)
            implicit none
            integer n
            real(DPT) err
            complex(DPT),allocatable,  dimension(:,:) :: a, b, c
            allocate(a(n,n), b(n,n), c(n,n))
            call random_complex(a)
            b = a
            call MatInv(b)
            call MatInv(b)
            c = a - b
            err = calc_error(c)
        end function mat_inv_test
        
    end module test
    
    program MatrixInversionTest
    
    use test
    implicit none
    integer n
    real(DPT) mkl_err, MatInt_err

    ! Body of MatrixInversionTest
    print *, 'Hello World'
    print *, "         n  ", "        mkl error", "        MatInv error"
    do n = 500, 5000, 500
        mkl_err    = mkl_mat_inv_test(n)
        MatInt_err = mat_inv_test(n)
        print *, n, mkl_err, MatInt_err
    end do
    
    pause
    end program MatrixInversionTest
    
