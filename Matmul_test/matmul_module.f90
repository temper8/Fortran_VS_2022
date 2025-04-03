module matmul_module
    contains 
    
    subroutine matmul_no_parallel(a,b,c,n)
        implicit none
        real(8) a(n,n),b(n,n),c(n,n)
        integer i,j,k,n
        c=0.d0
        !DIR$ NOPARALLEL
        do i=1,n         ! Outer loop is parallelized.
            do j=1,n      ! inner loops are interchanged
                do k=1,n   ! new inner loop is vectorized 
                    c(j,i)=c(j,i)+a(k,i)*b(j,k)
                enddo
            enddo
        enddo
    end    
    
    subroutine matmul_parallel(a,b,c,n)
	    real(8) a(n,n),b(n,n),c(n,n)
	    c=0.d0
	    !DIR$ PARALLEL  ALWAYS
	    do i=1,n         ! Outer loop is parallelized.
	        do j=1,n      ! inner loops are interchanged
		        do k=1,n   ! new inner loop is vectorized 
			        c(j,i)=c(j,i)+a(k,i)*b(j,k)
		        enddo
	        enddo
	    enddo
    end
    !subroutine matmul_MKL(a,b,c, n)
    subroutine matmul_MKL(a,b,c)
        integer n
	    real(8) a(:,:),b(:,:),c(:,:)
	    c=0.d0
        n = size(a,1)
        call dgemm('N', 'N', n, n, n, 1.0d0, a, n, b, n, 0.0d0, c, n)
    end    
    
    function matmul_MKL_f(a,b) result(c)
        integer n
	    real(8), allocatable :: a(:,:),b(:,:), c(:,:)
        n = size(a,1)
        allocate(c(n,n))
	    c=0.d0
        call dgemm('N', 'N', n, n, n, 1.0d0, a, n, b, n, 0.0d0, c, n)
    end    
end module matmul_module