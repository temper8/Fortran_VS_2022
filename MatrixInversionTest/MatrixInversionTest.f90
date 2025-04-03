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
        
        function mkl_mat_inv_test(mat) result(err)
            implicit none
            real(DPT) err
            complex(DPT),dimension(:,:) :: mat
            complex(DPT),dimension(size(mat,1), size(mat,1)) :: b, c
            b = mat
            call inverse_matrix_mkl(b)
            call inverse_matrix_mkl(b)
            c = mat - b
            err = calc_error(c)
        end function mkl_mat_inv_test
    
        function mat_inv_test(mat) result(err)
            implicit none
            real(DPT) err
            complex(DPT),dimension(:,:) :: mat
            complex(DPT),dimension(size(mat,1),size(mat,1)) :: b, c
            b = mat
            call MatInv(b)
            call MatInv(b)
            c = mat - b
            err = calc_error(c)
        end function mat_inv_test
        
    end module test
    
    program MatrixInversionTest
    
    use test
    implicit none
    integer n
    real(DPT) mkl_err, MatInt_err
    real MatInv_time, mkl_time, time0
    complex(DPT),allocatable,  dimension(:,:) :: mat

    ! Body of MatrixInversionTest
    print *, 'Matrix Inversion Test'
    print *, "         n  ", "        mkl error", "        MatInv error"
    do n = 500, 5000, 500
        allocate(mat(n,n))
        call random_complex(mat)
        time0 = sys_time()
        mkl_err    = mkl_mat_inv_test(mat)
        mkl_time = sys_time() - time0
        time0 = sys_time()
        MatInt_err = mat_inv_test(mat)
        MatInv_time = sys_time() - time0
        print *, n, mkl_err, MatInt_err
        print *, '  time=   ', mkl_time,  MatInv_time, 'Ratio= ', MatInv_time/mkl_time
        deallocate(mat)
    end do
    
    pause
    end program MatrixInversionTest
    
