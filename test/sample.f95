 
 DO k = 1,20 
   
    mass = 0.2d0*(-10+k)
  
    call Measure(stor_Z20,stor_Z31,mass)
    
    
   
    
    sm20 = 0
    sm31 = 0
    
    do i = 1,hmc_steps
       sm20 = sm20 + stor_Z20(i)
       sm31 = sm31 + stor_Z31(i)
    end do
    
    
   !print*, "mass=",mass,"Area =", Area
   
   
   mm20 = (cdlog(sm20/float(hmc_steps))/bb)/float(2)
   mm31 = (cdlog(sm31/float(hmc_steps))/bb)/float(2)
   
   
   
   print*,"final = ",k,mm20,mm31
   
   
   
   write(20,*) mass,mm20
   write(31,*) mass,mm31
   
   
   
   
   
END DO   
   
   


       