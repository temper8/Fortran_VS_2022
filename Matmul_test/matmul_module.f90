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
	
	    do i=1,n         ! Outer loop is parallelized.
	        do j=1,n      ! inner loops are interchanged
		        do k=1,n   ! new inner loop is vectorized 
			        c(j,i)=c(j,i)+a(k,i)*b(j,k)
		        enddo
	        enddo
	    enddo
    end
end module matmul_module