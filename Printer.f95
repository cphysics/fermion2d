


!====================================================================================
!This module prints matrices from configuration UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 
 MODULE Printer
 
 
        !----use------
        use Data_set
        use Algebra
        !---------------
        
        
       IMPLICIT NONE
       CONTAINS
           
 
!==============================================

	subroutine print_links(UU)
	implicit none
	integer:: s,t,r,p,q
	complex*16,dimension(:,:,:,:,:):: UU

    	do r = 1,2
    		do s = 1,L
    			do t = 1,L
    				do p = 1,N
   					 	do q = 1,N
    						print*, r,s,t,p,q,UU(r,s,t,p,q)
    				 	end do
                	end do
             	end do
          	end do
    	end do
 	return
 	end    
        
!==============================================
 subroutine print_matrix(A,dim)
	implicit none
	integer:: p,q,dim
	complex*16:: A(dim,dim)
	
    	do p = 1,dim
   		do q = 1,dim
    			print*,p,q, A(p,q)
    	end do
        end do
             	
 	return
 	end    
  
     
      
!================================================================================== 
 subroutine Unitary_checker_links(UU)
     implicit none
     integer::         							 p,q,r,s,t
     complex*16::       						 UU(2,L,L,N,N),Det
     complex*16,allocatable,dimension(:,:)::     HA,II,U,A
     
     
     Allocate(HA(N,N),II(N,N),U(N,N),A(N,N))
    
    do r = 1,2
            
        if (r .eq. 1) then
             
            do s = 1,L
            do t = 1,L
            
            			do p = 1,N
            			do q = 1,N
               				U(p,q)  = UU(r,s,t,p,q)
               				A(p,q)  = UU(r,s,t,p,q) 
            			end do
            			end do
            
            
                      call Detm(U,Det,N)			
                      print*,"r =",r,"s =",s,"t=",t,"Det=", Det
           		
                     call getH(A,HA,N)      
    		           II = matmul(A,HA)
    		         call print_matrix(II,N)  
               
            end do
            end do
                    
                    
                    
          else if (r .eq.2) then 
          
                   
                print*,"====================================== Horizontal ========================================="		
                print*,"=========================================<><><> ==========================================="	     
                     
              do s = 1,L
              
            			do p = 1,N
            			do q = 1,N
               				U(p,q)  = UU(r,s,L,p,q)
               				A(p,q)  = UU(r,s,L,p,q) 
            			end do
            			end do
                       
            
                      call Detm(U,Det,N)
                      print*,"r =",r,"s =",s,"t=",L,"Det=", Det
           		
                     call getH(A,HA,N)       !Algebra
    		           II = matmul(A,HA)
    		         call print_matrix(II,N)  !printer
               
              end do
         
     end if            
                    
          
            end do
                
     Deallocate( HA,II,U,A)
      return
      end                        
      

                        
                        
                      
!================================================================================== 
 subroutine Hermitian_checker_links(PI)
     implicit none
     integer::         							 p,q,r,s,t
     complex*16::       						 PI(2,L,L,N,N),Trc
     complex*16,allocatable,dimension(:,:)::     HA,A
     
     
     Allocate(HA(N,N),A(N,N))
    
    do r = 1,2
            
        if (r .eq. 1) then
             
            do s = 1,L
            do t = 1,L
            
            			do p = 1,N
            			do q = 1,N
               				A(p,q)  = PI(r,s,t,p,q)
            			end do
            			end do
            
            
                      call trace(A,N,trc)			
                      print*,"r =",r,"s =",s,"t=",t,"Trc=", trc
     
                	  call getH(A,HA,N) !self
                	  print*, "==============================="
                      do p = 1,N
                      do q = 1,N
                          print*, p,q, A(p,q) , HA(p,q)
                      end do
                      end do
               
            end do
            end do
                    
                    
                    
          else if (r .eq.2) then 
          
                   
                print*,"====================================== Horizontal ========================================="		
                print*,"=========================================<><><> ==========================================="	     
                     
              do s = 1,L
              
            			
            			do p = 1,N
            			do q = 1,N
               				A(p,q)  = PI(r,s,L,p,q)
            			end do
            			end do
            
            
                      call trace(A,N,trc)			
                      print*,"r =",r,"s =",s,"t=",t,"Trc=", trc
     
                	  call getH(A,HA,N) !self
                	  print*, "==============================="
                      do p = 1,N
                      do q = 1,N
                          print*, p,q, A(p,q) , HA(p,q)
                      end do
                      end do
               
              end do
         
     end if            
                    
          
            end do
                
     Deallocate( HA,A)
      return
      end                        
      

                        
                        
!=====================================================                         
                        
     
        
                        
                        
!=====================================================            
         end module Printer
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         