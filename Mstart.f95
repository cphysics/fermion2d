 
 

!====================================================================================
!This module generates Link matrices for configuration UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 MODULE Mstart
 
 
 
 !-----------------
 use Data_set
 use Algebra
 use Printer
 !--------------------
 
 
        IMPLICIT NONE
        CONTAINS
        
 !========================================================================   
        !GENERAL!       
 !=======================================================================    
     subroutine cold_Ustart(UU)
     implicit none
     integer::					p,q,r,s,t
     complex*16,intent(out)::   UU(2,L,L,N,N)


      do r = 1,2
        do s = 1,L
            do t = 1,L
                do p = 1,N
                    do q = 1,N
                        if (p .eq. q) then
                            UU(r,s,t,p,q) = dcmplx(1,0)
                        else if(p.ne.q) then
                            UU(r,s,t,p,q) = dcmplx(0,0)
                        end if
                     end do
                  end do
              end do
           end do
        end do

      return
      end

 !=======================================================================    
     subroutine cold_Hstart(HH)
     implicit none
     integer::					p,q,r,s,t
     complex*16,intent(out)::   HH(2,L,L,N,N)


      do r = 1,2
        do s = 1,L
            do t = 1,L
                do p = 1,N
                    do q = 1,N
                        if (p .eq. q) then
                            HH(r,s,t,p,q) = dcmplx(0.d0,0.d0)
                        else if(p.ne.q) then
                            HH(r,s,t,p,q) = dcmplx(0.d0,0.d0)
                        end if
                     end do
                  end do
              end do
           end do
        end do

      return
      end
 !===========================================================================   
     subroutine hot_Ustart(UU)
     implicit none
     integer::                                 p,q,r,s,t
     complex*16,allocatable,dimension(:,:) ::  NU
     complex*16,intent(out):: 				   UU(2,L,L,N,N)


      Allocate(NU(N,N))
      
      do r = 1,2
        do s = 1,L
            do t = 1,L
            
                    call sun_generator(NU,N)!self
             
                    do p = 1,N
                    do q = 1,N
                       UU(r,s,t,p,q) = NU(p,q)
                    end do
                    end do
                    
              end do
           end do
        end do
        
      Deallocate(NU)
      return
      end       
       
       
!=======================================================================   
     subroutine hot_Hstart(HH)
     implicit none
     integer::									p,q,r,s,t
     complex*16,allocatable,dimension(:,:) ::   NH
     complex*16,intent(out)::                   HH(2,L,L,N,N)

      Allocate(NH(N,N))
      
      do r = 1,2
        do s = 1,L
            do t = 1,L
            
             call H_generator(NH,N)
            
                    do p = 1,N
                    do q = 1,N
                       HH(r,s,t,p,q) = NH(p,q)
                    end do
                    end do
                    
              end do
           end do
        end do
        
      Deallocate(NH)
      return
      end
!================================================================================
                !IMPORTANT!               
!================================================================================   
     subroutine fix_Ustart(UU)
     implicit none
     integer::									p,q,r,s,t
     complex*16,allocatable,dimension(:,:) ::   NU
     complex*16,intent(out):: 					UU(2,L,L,N,N)



        Allocate(NU(N,N))
      ! make all links identity------
       call cold_Ustart(UU)
    !links along space direction (r = 1)----------    
     
        do s = 1,L
            do t = 1,L
            
                call sun_generator(NU,N)
                
                
                do p = 1,N
                do q = 1,N
                       UU(1,s,t,p,q) = NU(p,q)
                end do
                end do
                  
                  
              end do
           end do
     
    !links along time direction(r = 2) at all position but time t = L
    
    
        do s = 1,L
        
            call sun_generator(NU,N)
        
                    do p = 1,N
                    do q = 1,N
                       UU(2,s,L,p,q) = NU(p,q)
                    end do
                    end do
                  
                  
              end do
     
    
      Deallocate(NU)
      return
      end
      

!=================================================================================
!================================================================================   
     subroutine fix_Hstart(HH)
     implicit none
     integer::									p,q,r,s,t
     complex*16,allocatable,dimension(:,:) ::   NH
     complex*16,intent(out):: 					HH(2,L,L,N,N)



      
      Allocate(NH(N,N))
    ! make all links null------
      call Null(HH)
        
    !links along space direction (r = 1)----------    
     
        do s = 1,L
            do t = 1,L
            
                call H_generator(NH,N)
                
                
                do p = 1,N
                do q = 1,N
                       HH(1,s,t,p,q) = NH(p,q)
                end do
                end do
                  
                  
              end do
           end do
     
    !links along time direction(r = 2) at all position but time t = L
    
    
        do s = 1,L
            call H_generator(NH,N)
                    do p = 1,N
                    do q = 1,N
                       HH(2,s,L,p,q) = NH(p,q)
                    end do
                    end do
        end do
     
    
      Deallocate(NH)
      return
      end
      

!=================================================================================
 
       subroutine sun_generator_orginal(U,dim)
        implicit none
        integer::	                  			      t,s 
        integer,intent(in)::           				  dim
   	    real*8::                                     T1,T2,T3,phi,xi,theta,pi
        complex*16,allocatable,dimension(:,:)::      SU2,SUN
        complex*16::                                 ai,a_a,b_b
        complex*16,intent(out)::                     U(dim,dim)
       
        
        

       
        pi = dacos(-1.d0)
        ai = dcmplx(0.d0,1.d0)
        
        
        
                            Allocate( SU2(2,2),SUN(dim,dim))
        
                            call IIN(U,dim,dcmplx(1.d0,0.d0))!Algebra
                            
                            


                                      s = 1
            26                        t = s+1
                              
                                            
                                       
            25          T1 = rand()
                        T2 = rand()
                        T3 = rand()
                        !get angles
                        xi = pi*(2*T1-1.d0)
				        theta = 0.5d0*dacos(2*T2-1.d0)
				        phi = pi*(2*T3-1.d0)
				        !get elements
                        a_a = dcos(theta)*(cdexp(ai*phi))
				       	b_b = dsin(theta)*(cdexp(ai*xi))
				       	!generate SU2
                        SU2(1,1) = a_a
				        SU2(1,2) = b_b
                        SU2(2,2) = dconjg(SU2(1,1))
				        SU2(2,1) = -dconjg(SU2(1,2))
                            
                            
                       

       
                         call IIN(SUN,dim,dcmplx(1.d0,0.d0))!Algebra   
                                           

                          !plug SU2 to SUN
                          SUN(s,s) =  SU2(1,1)
                          SUN(s,t) =  SU2(1,2)
                          SUN(t,s) =  SU2(2,1)
                          SUN(t,t) =  SU2(2,2)
                          
                              
                          !print*,s,t
                          
                          
                          U = MATMUL(U,SUN)
                          

        t = t+1
                         
                          if (t .lt. dim+1) then
                              go to 25
                          else if (t .eq. dim+1) then
                                s = s+1
                                if (s .lt. dim) THEN
                                    go to 26
                               end if
                         end if
                                
                                
                        Deallocate( SU2,SUN)                           
                        return
                        end
  !=================================================================================
 
       subroutine sun_generator(U,dim)
        implicit none
        integer::	                  			      t,s 
        integer,intent(in)::           				  dim
   	    real*8::                                     T1,T2,T3,phi,xi,theta,pi
        complex*16,allocatable,dimension(:,:)::      SU2,SUN
        complex*16::                                 ai
        real*8::                                     a,a0,a1,a2,a3
        complex*16,intent(out)::                     U(dim,dim)
       
        
        

       
        pi = dacos(-1.d0)
        ai = dcmplx(0.d0,1.d0)
        
        
        
                            Allocate( SU2(2,2),SUN(dim,dim))
        
                            call IIN(U,dim,dcmplx(1.d0,0.d0))!Algebra
                            
                            


                                      s = 1
            26                        t = s+1
                              
                                            
                                       
            25          T1 = rand()
                        T2 = rand()
                        T3 = rand()
                       
                       
				        theta = dacos(2*T2-1.d0)
				        phi = pi*(T3)
				        a0 = (2.d0*T1 - 1.d0 )
				        
				        a=dsqrt(1.d0-a0*a0)
				        
            			a1=a*dsin(theta)*dcos(phi)
            			a2=a*dsin(theta)*dsin(phi)
            			a3=a*dcos(theta)
            			 
            			 
            			 SU2(1,1)=a0+ai*a3
                         SU2(1,2)=a2+ai*a1
                         SU2(2,1)=-dconjg(SU2(1,2))
                         SU2(2,2)=dconjg(SU2(1,1))
				       
				       
				  
       
                         call IIN(SUN,dim,dcmplx(1.d0,0.d0))!Algebra   
                                           

                          !plug SU2 to SUN
                          SUN(s,s) =  SU2(1,1)
                          SUN(s,t) =  SU2(1,2)
                          SUN(t,s) =  SU2(2,1)
                          SUN(t,t) =  SU2(2,2)
                          
                              
                          !print*,s,t
                          
                          
                          U = MATMUL(U,SUN)
                          

        t = t+1
                         
                          if (t .lt. dim+1) then
                              go to 25
                          else if (t .eq. dim+1) then
                                s = s+1
                                if (s .lt. dim) THEN
                                    go to 26
                               end if
                         end if
                                
                                
                        Deallocate( SU2,SUN)                           
                        return
                        end
                                              
 !=======================================================================================                       
   subroutine H_generator(H,dim)
        implicit none
        
        integer,intent(in)::           dim
        integer::	                   i,j,k,l,p,q,u,s
   	    real*8::                       x,phi,a,b,pi,trc,t
        complex*16,intent(out)::       H(dim,dim)
        
        
        
        pi = dacos(-1.d0)
        
        
      
        
        			!Generating upper diagonal of the matrix
        			
        			do p = 1,dim
        			do q = p,dim
        			
        			if (p .ne. q) then
        			
        			
        			    t = rand() 
                        x = dsqrt(-dlog(t)) ! 2.d0 cancels out
                        
        				!phi = 2.d0*pi*rand()
        				phi = 2.d0*pi*(rand() - 0.5d0)
        			
        				a = x*dcos(phi)
        				b = x*dsin(phi)
        				
        				H(p,q) = dcmplx(a,b)
        				H(q,p) = dcmplx(a,-b)
        				
        			else
        			
        			    t = rand()
        				x = dsqrt(-2.d0*dlog(t))
        				phi = 2.d0*pi*(rand() - 0.5d0)
        				H(p,q) =  x*dcos(phi)
        			
        			end if
        			
        			end do
        			end do
        			
        			
        		
        			!imposing tracelessness condition
        			trc = 0.d0
        			do p = 1,dim
        			trc = trc+ H(p,p)
        			end do
        
        			do p = 1,dim
        			H(p,p) = H(p,p) - (trc/float(dim))
        			end do
        			
        			
                    !call print_matrix(H,dim)
                   
        			return 
        			end
        			

 !================================================================================== 
   
     
     
     
        
         
!==================================================================================    
   end module Mstart      
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         