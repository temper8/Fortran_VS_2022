!  Matmul_test.f90 
!
!  FUNCTIONS:
!  Matmul_test - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Matmul_test
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Matmul_test
    use time_module
    use matmul_module
    implicit none
     integer nn, i, istat
     character(len=70) errmes_string
     real(8), allocatable, dimension(:,:) :: a, b, c
     real   :: all_time1, all_time2
     real   :: s_time
     real   :: time21, time22
     real   ::  cpu_t1,cpu_t2 
     real   ::  cpu_dt(5)
     real   ::  sys_dt(5)
     
     all_time1 = sys_time()
     print *, 'Hello, matmul'
     
     call print_compiler_info
     
     print *, "         n", "         no-parallel", "     parallel  ",  "       MKL     ",  "     MatMul"
     
     do nn = 500, 1500, 100
         allocate(a(nn,nn), b(nn,nn),c(nn,nn), stat=istat, errmsg=errmes_string )
         if ( istat /= 0 ) then
            write ( *,* )  errmes_string,' : istat =',istat
         endif
         s_time = sys_time()
         call cpu_time(cpu_t1) 
         CALL RANDOM_NUMBER(a)
         CALL RANDOM_NUMBER(b)    
         call cpu_time(cpu_t2)
         sys_dt(1) = sys_time() - s_time
         cpu_dt(1) = cpu_t2-cpu_t1
         !print *, sum(c) 

         s_time = sys_time()
         call cpu_time(cpu_t1) 
         do i = 1,10
            !c = MatMul(a,b)
            call matmul_no_parallel(a,b,c,nn)
            !call matmul_MKL(a,b,c,nn)
         enddo

         call cpu_time(cpu_t2)
         sys_dt(2) = sys_time() - s_time
         cpu_dt(2) = cpu_t2-cpu_t1
         !print *, sum(c)     

         s_time = sys_time()
         call cpu_time(cpu_t1) 
         do i = 1,10
            call matmul_parallel(a,b,c,nn)
         enddo
         call cpu_time(cpu_t2)
         sys_dt(3) = sys_time() - s_time
         cpu_dt(3) = cpu_t2 - cpu_t1
         !print *, sum(c)
         s_time = sys_time()
         call cpu_time(cpu_t1) 
         do i = 1,10
            !c = matmul_MKL_f(a,b)
            call matmul_MKL(a,b,c)
            !call matmul_MKL(a,b,c, nn)
         enddo
         call cpu_time(cpu_t2)
         sys_dt(4) = sys_time() - s_time
         cpu_dt(4) = cpu_t2 - cpu_t1

         s_time = sys_time()
         call cpu_time(cpu_t1) 
         do i = 1,10
            c = MatMul(a,b)
         enddo
         call cpu_time(cpu_t2)
         sys_dt(5) = sys_time() - s_time
         cpu_dt(5) = cpu_t2 - cpu_t1         
         
         print *, nn, "---------------"
         print *, ' cpu time   ', cpu_dt(2:5)
         print *, ' sys time   ', sys_dt(2:5)
         deallocate(a, b, c)
     end do
     all_time2 = sys_time()
     print *,"all time =", all_time2 - all_time1

    pause
    
    contains
        subroutine print_compiler_info
        use, intrinsic :: iso_fortran_env, only : compiler_version
        use, intrinsic :: iso_fortran_env, only : compiler_options
        implicit none
           print '(4(a,/))', &
              'This file was compiled by: ', &
              compiler_version(),           &
              'using the options: ',        &
              compiler_options()
        end subroutine print_compiler_info
    end program Matmul_test

