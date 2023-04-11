module time_module
    Contains   
    function sys_time()
       implicit none
       real(8) sys_time
       integer count, count_rate, count_max
       call system_clock(count, count_rate, count_max)
       sys_time = count*1.0/count_rate
       return
    end   
 
 end module time_module
 