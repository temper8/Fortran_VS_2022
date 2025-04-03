module MKL_helper

    contains
    
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
        integer :: ipiv(size(A,1))           ! ������ ������������
        integer :: info, lwork, i
        integer n
    
        n = size(A,1)

        ! ��� 1: LU-������������
        call zgetrf(n, n, A, n, ipiv, info)
        if (info /= 0) then
            print *, '������ LU-����������. ���:', info
            stop
        end if

        ! ��� 2: ����������� ������������ ������� �������� �������
        lwork = -1
        allocate(work(1))
        call zgetri(n, A, n, ipiv, work, lwork, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))

        ! ��� 3: ���������� �������� �������
        call zgetri(n, A, n, ipiv, work, lwork, info)
        if (info /= 0) then
            print *, '������ ��������� �������. ���:', info
            stop
        end if

        deallocate(work)
    end subroutine 
end module MKL_helper