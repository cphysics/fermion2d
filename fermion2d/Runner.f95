!====================================================================================
!This module runs HMC
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 MODULE Runner
 
 !-----------------------------------------------
 use Data_set
 use Algebra
 use Mstart
 use Updator
 use Calculator
 use Partition
 use Printer
 !-----------------------------------------------
 
 
IMPLICIT NONE
CONTAINS       
      
 
!========================================================================================
!                                           HMC                                          !
!========================================================================================


 Subroutine  HMC(mass,writeid)
            implicit none
            integer::      							   t,count,istat,k,ilf,writeid
            real*8::       							   a,b,S0,S1,S2,eds,Sa,Sb,Sc,SUG
            complex*16,dimension(2,L,L,N,N)::          UU,PI,PII,UF,PIF,PI3,U3
            real*8,intent(in)::                        mass
            
            
          
         
                 

                    
                    !call fix_Ustart(UU)
                    call cold_Ustart(UU)
                    call fix_Hstart(PI)
                  
                   
                    call Action(S0,Sa,Sb,Sc,UU,mass,PI)
                    
                    
                    if (writeid .eq.1) then
                    print*,"================= first ===============>",0, S0
                    write(555,*)"==== first ======>",0, S0
                    Write(1001,*),0,S0
                    end if
                    
                   
count = 0                    
t = 0


               
26                      UF = UU
                        call Leap_Frog(UU,mass,PI,ihmc = t,printid = 0)
                        ilf = ilf + 1
                      
                       
                        call Action(S1,Sa,Sb,Sc,UU,mass,PI)
                        eds = dexp(-(S1-S0))
                        a = min(1.d0,eds)
        			    b = rand()
        			    
        			    
        		 If (b .lt. a) then 
                        count = count+1 
                        
                        
                        
                         if (writeid .eq.1) then
                         print*,"=======  success =========>",t,S1,eds,a,b
                         write(555,*)"=======  success ======>",t,S1,S1-S0,a,b
                         end if
                         
                         
                        call Plaquette(SUG,UU,mass)
                        
                        print*,t,SUG
                        Write(1001,*) t,SUG
                        
                       
                       
                      
                        call fix_Hstart(PI3)
                        PI = PI3
                        call Action(S0,Sa,Sb,Sc,UU,mass,PI)

                        
                Else 
                       
                        call fix_Hstart(PI3)
                        PI = PI3
                        UU = UF
                        
                        
                        
                         if (writeid .eq.1) then
                         print*,"=======  failed  ========>",t,S0,eds,a,b
                         write(555,*)"=======  failed ======>", t,S0,S1-S0,a,b
                         end if
                         
                         
                         call Plaquette(SUG,UU,mass)
                         
                         
                        print*,t,SUG
                        Write(1001,*) t,SUG
                        
                        
                        call Action(S0,Sa,Sb,Sc,UU,mass,PI)
               End if
                
                
             
                
  t = t+1       
 
            if (t .lt. HMC_steps) then
             go to 26
             end if
             
             
                        
  print*,"Acceptance % =", (count/float(HMC_steps))*100 ,"%"                 
                    
            return
            end
            	
    
 !========================================================================================
!                                           Measure                                         !
!========================================================================================


 Subroutine  Measure(stor_Z20,stor_Z31,mass,massid,printid)
            implicit none
            integer::      							   tt,count,p,ilf,ihmc
            real*8::       							   a,b,S0,S1,S2,eds,Sa,Sb,Sc
            complex*16,dimension(2,L,L,N,N)::          UU,PI,PII,UF,PIF,PI3,U3
            complex*16::                               Z20,Z31,Z31A,Z31B
            complex*16,intent(out)::                   stor_Z20(hmc_steps), stor_Z31(hmc_steps)
            real*8,intent(in)::                        mass
            integer,intent(in)::                       massid,printid
            
            

                    
                    call fix_Ustart(UU)
                    !call cold_Ustart(UU)
                    call fix_Hstart(PI)
                    
                    call Action(S0,Sa,Sb,Sc,UU,mass,PI) 
                   
                    do p = 1,hmc_steps
                       stor_Z20(p) = dcmplx(0,0)
                       stor_Z31(p) = dcmplx(0,0)
                    end do
                   
                   
ilf = 0                
count = 0   
tt = 1 !total trial


               
26                      UF = UU

                        if (printid .eq.1) then
                        print*, "================================================================================" 
                        print*, "L=",L,":LF_steps=",LF_steps,":hmc_steps=",hmc_steps,":mass_N=",mass_N
                        print*, "--------------------------------------------------------------------------------"
                        print*, 'Job-tracker:MassId = ',massid,'HMC-step =',tt
                        print*, "================================================================================" 
                        end if
                        
                        call Leap_Frog(UU,mass,PI,ihmc = tt,printid = 0)
                       
                      
                       
                        call Action(S1,Sa,Sb,Sc,UU,mass,PI)
                        eds = dexp(-(S1-S0))
                        a = min(1.d0,eds)
        			    b = rand()
        			    
        			    
        		 If (b .lt. a) then 
                        count = count+1 
                        
                        
                        print*,"****************Accepted********************"
                       
                        ! Do measurement after thermalization
                       
                        call PZ20( UU,mass,Z20)
                        
                        
                        
                        !call PZ31A( UU,mass,Z31A)
                        call PZ31B( UU,mass,Z31B)   !-------------------------------
                        
                        
                       ! print*, Z31A,Z31B   !----------------------------------------
                        
                        
                        stor_Z20(count) = Z20
                        
                        
                        !stor_Z31(count) = Z31A  !-------------------
                        stor_Z31(count) = Z31B   !-------------------
                         
                         
                        call fix_Hstart(PI3)
                        PI = PI3
                        call Action(S0,Sa,Sb,Sc,UU,mass,PI)

                        
                Else 
                
                     print*," xxxxxxxxxxxx Not-Accepted  xxxxxxxxxxxxxx"
                     
                     
                        call fix_Hstart(PI3)
                        PI = PI3
                        UU = UF
                        
                      ! Do nothing
                        
                       
                        call Action(S0,Sa,Sb,Sc,UU,mass,PI)
               End if
                
   tt = tt+1               
             
                      if (count .lt. hmc_steps) then
                      go to 26
                      end if
      
 
            !if (t .lt. HMC_steps+1) then
             !go to 26
             !end if
             
             
             
             
			if (printid .eq. 1) then             
			print*, "================================================================================" 
			print*,"count =",count,"total trial = ",tt,"Acceptance % =", (count/float(tt))*100 ,"%"                 
			print*, "================================================================================"    
			end if


                 
            return
            end
            	
    
!============================================================================================================

subroutine mass_finder(mass_N)
  implicit none
  integer::                 		    i,j,k,r,p,q,s,t,printid,massid
  integer,intent(in)::                  mass_N
  complex*16::                          Tau(ln2,ln2),Z20,stor_Z20(hmc_steps),stor_Z31(hmc_steps)
  real*8,dimension(hmc_steps)::         m20,m31
  real*8::                              sm20,sm31,var_m20,var_m31,std_dev_m20,std_dev_m31
  real*8::                              mass,mm20,mm31
  
 
  
 
 DO k = 1,mass_N
   
    mass = 0.1d0*(-(mass_N)/float(2) + k)
  
    call Measure(stor_Z20,stor_Z31,mass,massid = k,printid = 1)
    
    
    
        do i = 1,hmc_steps
    		m20(i) = real(cdlog(stor_Z20(i))/bb)*(1/float(2))
  		    m31(i) = real(cdlog(stor_Z31(i))/bb)*(1/float(3))
    	end do 
    	
    !print*,m31
    
        sm20 = 0.d0
        sm31 = 0.d0
    	do i = 1,hmc_steps
       		sm20 = sm20 + (m20(i))
       		sm31 = sm31 + (m31(i))
    	end do
    
    
      !print*, "mass=",mass,"Area =", Area
   
   
  	 	mm20 = sm20/float(hmc_steps)
  	 	mm31 = sm31/float(hmc_steps)
   
  	 	var_m20 = 0.d0
   	 	var_m31 = 0.d0
  	 	do i = 1,hmc_steps
    		var_m20 = var_m20  + (m20(i) - mm20)**2
  			var_m31 = var_m31  + (m31(i) - mm31)**2
    	end do 
   
   
   		do i = 1,hmc_steps
    		std_dev_m20 = dsqrt(var_m20/float(hmc_steps))
  			std_dev_m31 = dsqrt(var_m31/float(hmc_steps)) 
  		end do 
  		
    
   
   		print*,"final = ",k,mm20,mm31,std_dev_m20,std_dev_m31
   
   
   		write(20,*) mass,mm20,std_dev_m20
   		write(31,*) mass,mm31,std_dev_m31
   
   
   
   
   
END DO   
   
   


return
end
       
                     	           	
           			     
!=====================================================================================
!                        Reversibility - Checker                                      !
!=====================================================================================           			       
           			       
                           

 Subroutine  Reversibility(mass)
            implicit none
            integer::      							   t,k,p,q,r,s,ilf,drx
            real*8::       							   a,b,S0,S1,S2
            complex*16,dimension(2,L,L,N,N)::          UU,PI,UA,UB,UC,PIA,PIB,PIC
            real*8,intent(in)::                                       mass
            
            
   
                    
                    call fix_Ustart(UU)
                    call fix_Hstart(PI)
                    UA = UU
                    PIA = PI
              
                    call Leap_Frog_test(UU,mass,PI,ilf = 200,drx = 1)
                     
                    UB = UU
                    PIB = PI

                   call Leap_Frog_test(UU,mass,PI,ilf = 200,drx = -1)
                
                    UC = UU
                    PIC = PI
                    
                    
                    
         do r = 1,2
    		do s = 1,L
    			do t = 1,L
    				do p = 1,N
   					 	do q = 1,N
    						print*, r,s,t,p,q,UA(r,s,t,p,q),UC(r,s,t,p,q)
    				 	end do
                	end do
             	end do
          	end do
    	end do
    	
    	
    	
    	
 	return
 	end    


!=====================================================================================



end module Runner









