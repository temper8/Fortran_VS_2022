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
     integer nn
     real(8), dimension(5000,5000) :: a, b, c
     real   :: all_time1, all_time2
     real   :: s_time
     real   :: time21, time22
     real   ::  cpu_t1,cpu_t2 
     real   ::  cpu_dt1, cpu_dt2, cpu_dt3 
     real   ::  sys_dt1, sys_dt2, sys_dt3 
     
     all_time1 = sys_time()
     print *, 'Hello, matmul'
     print *, "         n", "      init","         parallel", "     no-parallel"
     
     do nn = 500, 2500, 100
         s_time = sys_time()
         call cpu_time(cpu_t1) 
         CALL RANDOM_NUMBER(a)
         CALL RANDOM_NUMBER(b)    
         call cpu_time(cpu_t2)
         sys_dt1 = sys_time() - s_time
         cpu_dt1 = cpu_t2-cpu_t1
         !print *, sum(c) 

         s_time = sys_time()
         call cpu_time(cpu_t1) 
         call matmul_no_parallel(a,b,c,nn)
         call cpu_time(cpu_t2)
         sys_dt2 = sys_time() - s_time
         cpu_dt2 = cpu_t2-cpu_t1
         !print *, sum(c)     

         s_time = sys_time()
         call cpu_time(cpu_t1) 
         call matmul_parallel(a,b,c,nn)
         call cpu_time(cpu_t2)
         sys_dt3 = sys_time() - s_time
         cpu_dt3 = cpu_t2-cpu_t1
         !print *, sum(c)

         print *, nn, cpu_dt1, cpu_dt2, cpu_dt3
         print *,  0, sys_dt1, sys_dt2, sys_dt3
     end do
     all_time2 = sys_time()
     print *,"all time =", all_time2 - all_time1

    pause
    end program Matmul_test

