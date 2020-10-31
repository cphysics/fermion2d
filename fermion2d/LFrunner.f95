!====================================================================================
!This module runs HMC
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 MODULE LFrunner
 
 !-----------------------------------------------
    use Data_set
    use Algebra
    use Mstart
    use Updator
   
 !-----------------------------------------------
 
 
IMPLICIT NONE
CONTAINS       
      
 
!========================================================================================
!                                           HMC                                          !
!========================================================================================


 Subroutine  LFchecker(mass)
            implicit none
            integer::      							   t,count,istat,k,ilf,s,p,q,r,drx
            complex*16,dimension(N,N)::                 U,H
            complex*16::                                Det,HD(lln2,lln2),IHD(lln2,lln2)
            complex*16,dimension(2,L,L,N,N)::           UU,PI,PII,UF,PIF,PI3,U3
            real*8,intent(in)::                         mass
            
            
         
   

                    
                    call fix_Ustart(UU)
                    !call cold_Ustart(UU)
                    call fix_Hstart(PI)
                    
                    !call DataReader(UU,PI)
           
                    
     				       
                 
  
 

               
                     
                   call Leap_Frog_test(UU,mass,PI,ilf=1,drx = 1)
                      
                     
                        
                        
                        
           !CHECKER             
                        
            do r = 1,2
            do s = 1,L
            do t = 1,L
            
            			do p = 1,N
            			do q = 1,N
               				U(p,q)  = UU(r,s,t,p,q)
               				H(p,q)  = PI(r,s,t,p,q) 
               				
            			end do
            			end do
            
                        call Unitary_checker(U,N,'error on unitarity after LF')
                        call Hermitian_checker(H,N,'error on Hermiticity after LF')
               
            end do
            end do
            end do
            
            
            
            
            
 
              
                    
            return
            end
            	
    
    !=====================================================================================

Subroutine  DataReader(UU,PI)
            implicit none
            integer::      							   			   t,s,i,j,r
            complex*16,dimension(L,L,N,N)::            			   U,HU
            complex*16,dimension(L,N,N)::             			   V,HV
            complex*16,dimension(2,L,L,N,N),intent(out)::          UU,PI
            
            
           
          
           
           
           
           
           OPEN (unit = 101, file = 'fort.101')
           READ(101,*) ((((U(s,t,i,j),j = 1,N),i = 1,N),t = 1,L),s = 1,L)
           CLOSE(unit = 101)
           
           
           OPEN (unit = 102, file = 'fort.102')
           READ(102,*) (((V(s,i,j),j = 1,N),i = 1,N),s = 1,L)
           CLOSE(unit = 102)
        
        
           
           OPEN (unit = 401, file = 'fort.401')
           READ(401,*) ((((HU(s,t,i,j),j = 1,N),i = 1,N),t = 1,L),s = 1,L)
           CLOSE(unit = 401)
           
           
           OPEN (unit = 402, file = 'fort.402')
           READ(402,*) (((HV(s,i,j),j = 1,N),i = 1,N),s = 1,L)
           CLOSE(unit = 402)
             
       
           
      !----------- links ------------------------------------------     
           
           
           do s = 1,L
           do t = 1,L
           do i = 1,N
           do j = 1,N
           
           UU(1,s,t,i,j) = U(s,t,i,j)
          
          end do
          end do
          end do
          end do
          
          
          
           do s = 1,L
           do t = 1,L
           
           if (t .ne. L) then
           
           do i = 1,N
           do j = 1,N
           
           if (i .eq. j) then
              UU(2,s,t,i,j) = dcmplx(1.d0,0.d0)
           else
              UU(2,s,t,i,j) = dcmplx(0.d0,0.d0)
           end if
           
           end do
           end do
           
           
           else if (t .eq. L) then
           
           do i = 1,N
           do j = 1,N
           
            UU(2,s,t,i,j) = V(s,i,j)
            
           end do
           end do
            
           end if
          
           end do
           end do
          
      !----------------- momentum ------------------------------    
           
           do s = 1,L
           do t = 1,L
           do i = 1,N
           do j = 1,N
           
           PI(1,s,t,i,j) = HU(s,t,i,j)
          
          end do
          end do
          end do
          end do
          
          
          
           do s = 1,L
           do t = 1,L
           
           if (t .ne. L) then
           
           do i = 1,N
           do j = 1,N
           
           if (i .eq. j) then
              PI(2,s,t,i,j) = dcmplx(1.d0,0.d0)
           else
              PI(2,s,t,i,j) = dcmplx(0.d0,0.d0)
           end if
           
           end do
           end do
            
           else if (t .eq. L) then
           
           do i = 1,N
           do j = 1,N
            PI(2,s,t,i,j) = HV(s,i,j)
           end do
           end do
          
           end if
       
           end do
           end do
                
           
            return
            end       
 !=====================================================================================           
           
           
   subroutine printMTX(A,dim)
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
           
  
  !=====================================================================================           
         
    end module LFrunner









