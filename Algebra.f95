 

!====================================================================================
!This module does many Algebraic calculations from configuration UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 MODULE Algebra
 
        !----------------------------
        !use Printer
        use Data_set
        !---------------------------
        
        
        IMPLICIT NONE
        
       
        
        CONTAINS
                        
                        
                        
!==================================MINOR SECTION========================================
 
 subroutine Trace(A,N,trc)
     implicit none
     integer::                 p
     integer,intent(in)::      N
     complex*16,intent(in)::  A(N,N)
     complex*16,intent(out):: trc
     
                trc = dcmplx(0.d0,0.d0)
                
                do p = 1,N
                trc = trc + A(p,p)
                end do
              

      return
      end                      
   
   !=======================================================================================                         
   subroutine IIN(A,N,nbr)
     implicit none
     integer::                 p,q
     integer,intent(in)::      N
     complex*16,intent(in)::   nbr
     complex*16,intent(out)::  A(N,N)
     
                do p = 1,N
                    do q = 1,N
                        if (p .eq. q) then
                            A(p,q) = nbr
                        else if(p.ne.q) then
                            A(p,q) = dcmplx(0.d0,0.d0)
                        end if
                     end do
                end do
              

      return
      end                      
   
    !======================================================================================= 
   subroutine Null(OO)
   
            	   implicit none
            	   integer:: 				r,s,t,p,q
            	   complex*16,intent(out):: OO(2,L,L,N,N)
            	
            	    
            	    do r = 1,2
                    do s = 1,L
                	do t = 1,L
                	do p = 1,N
                    do q = 1,N
                        OO(r,s,t,p,q) = Dcmplx(0.d0,0.d0)
                    end do
                  	end do
              		end do
              		end do
              		end do
                       
                    
                     return 
                     end
 
 
 !=======================================================================================                       
   subroutine getH(A,H,N)
     implicit none
     integer::					p,q
     integer,intent(in)::       N
     complex*16,intent(in)::   A(N,N)
     complex*16,intent(out)::  H(N,N)
     
                do p = 1,N
                    do q = 1,N
                       H(p,q) = dconjg(A(q,p))
                    end do
                end do
              

      return
      end                     
       
    
    
!==================================MAJOR SECTION====================================== 
!======================================gnl=============================================== 

    subroutine Inverse_finder(A,IA,n)
     implicit none
     integer::                                      info,lwork
     integer,intent(in)::                           n
     integer,parameter::                            lwmax = 1000000
     complex*16,intent(in),dimension(n,n)::         A
     complex*16,allocatable,dimension(:,:)::        B
     integer,allocatable,dimension(:)::             ipiv
     complex*16,allocatable,dimension(:)::          work
     complex*16,intent(out)::                       IA(n,n)
                     
                   
     
                    !Lapack inversion of A
                    
                    allocate(B(n,n),work(n),ipiv(n)) !Allocate here
                    
                    B = A
                    
                    call ZGETRF(n,n,B,n,ipiv,info) !this provide info
                    
                    if (info .eq.0) then
                    
                    !print*, "yes"
                    
    	            !lwork = -1
    	            
     		        !call ZGETRI( n, B, n, ipiv, work, lwork, info)    !lwork = n
     		        
      		        !lwork = min( lwmax, int( work( 1 ) ) )
      		        
     		        call ZGETRI( n, B, n, ipiv, work, n, info)     !lwork = n
     		        
     		        end if 
     		        
     		        if (info .ne. 0) then
     		        print*, 'Matrix inversion failed!'
     		        end if
     		        
               
                    IA = B
                    
                    Deallocate(B,work,ipiv)!deallocate here
                    
                    return
                    end


!=======================================   INVERSE  ============================================== 
 
subroutine Inverser(A,IA,n)
  
     implicit none
     integer::                                      info,p,q
     integer,intent(in)::                           n
     complex*16,intent(in),dimension(n,n)::         A
     integer,allocatable,dimension(:)::             ipiv
     complex*16,allocatable,dimension(:)::          work2,work1
     complex*16,intent(out)::                       IA(n,n)
                     
                    
                    allocate(work1(n),work2(2*n),ipiv(n))
                    
                    
                   call hermitian_checker(A,n,'error on pre-hermiticity at inverse finder')
                    
                    
                    IA = A
                    
                    call ZHETRF('U', n, IA, n, ipiv, work2, 2*n, info) 
                    
                    if (info .eq.0) then
     		        call ZHETRI('U', n, IA, n, ipiv, work1, info)     
     		       
     		       
                    else
     		        print*, 'Matrix inversion failed!'
     		        end if
                    
                    
                    
                    do p = 2,n
                    do q = 1,p-1
                       IA(p,q)  = dconjg(IA(q,p))
                    end do
                    end do
                    
                    
                    
                    
                    call hermitian_checker(IA,n,'error on post-hermiticity at inverse finder')
                    
                    
                    
                    Deallocate(work1,work2,ipiv)
                    
                    return
                    end
   
 !============================================   Eigner -general  ===============================
    
    
   subroutine Eigen_finder(A,EV,W,n)
     implicit none
     integer::                      					p,q,info,lwork
     integer,intent(in)::                               n
     integer,parameter::            					lwmax = 100000000
     complex*16,intent(in),dimension(n,n)::    			A
     complex*16,allocatable,dimension(:,:):: 			B,VL,VR
     real*8,allocatable,dimension(:)::                  rwork
     complex*16,allocatable,dimension(:)::              work 
     complex*16,intent(out)::       					EV(n,n),W(n)
     CHARACTER*1::                  					UPLO,V
     
     Allocate(B(n,n),VL(n,n),VR(n,n),rwork(3*n -2),work(lwmax))
     
  

                      B = A
                      lwork = -1
                      call ZGEEV( "V", "V", n, B, n, W, VL, n, VR, n, work, lwork, rwork, info )
                      lwork = min( lwmax, int(work( 1 ) ) )
                      call ZGEEV( "V", "V", n, B, n, W, VL, n, VR, n, work, lwork, rwork, info )
                      
                      
                      
                      if( info.gt.0 ) then
                      write(*,*)'The algorithm failed to compute eigenvalues.'
                      end if
  
                      
                     do p = 1,n
                     do q = 1,n
                        EV(p,q) = VR(p,q)
                     end do
                     end do
                     
                     
                     
                     
                      Deallocate(B,VL,VR,rwork,work)
                      
                      return
                      end 
                      
!============================================   Eigner - hermitian ===============================
    
   subroutine Eigner(A,EV,Dig,n)
     implicit none
     integer::                      					p,q,info,lwork
     integer,intent(in)::                               n
     complex*16,intent(in),dimension(n,n)::    			A
     complex*16,allocatable,dimension(:,:):: 			B,II,HEV
     real*8,allocatable,dimension(:)::                  W,rwork
     complex*16,allocatable,dimension(:)::              work 
     complex*16,intent(out)::       					EV(n,n),Dig(n,n)
     
     
                     Allocate(B(n,n),rwork(3*n -2),work(2*n),W(n),II(n,n),HEV(n,n))
  	
     
     
               
                     call hermitian_checker(A,n,'Error on hermiticity at eigner')
                     
                     
                      B = A
                       
                   
                    call ZHEEV( 'V', 'U', n, B, n, W, work, 2*n, rwork, info )
                    
                      
                    if( info.ne.0 ) then
                    write(*,*)'The algorithm failed to compute eigenvalues.'
                    end if
                        
                    
                        
                    do p = 1,n
                    do q = 1,n
                       if (p .eq. q) then
                         Dig(p,q) = W(p)
                       else
                         Dig(p,q) = dcmplx(0.d0,0.d0)
                       end if
                    end do
                    end do
                      
                      
                    do p = 1,n
                    do q = 1,n
                        EV(p,q) = B(p,q)
                    end do
                    end do
                    
                    
            !----  checker --------------------- 
               
                  
                  call Unitary_checker(EV,N,'error on eigen vector matrix at Eigner')
                  
                  
                   
                     
                Deallocate(B,rwork,work,W,II,HEV)
                
                
                return
                end                      
                      
                      
  !=========================================== gnl ==========================================
	subroutine lnDetm_gnl(A,lnDet,N)
			implicit none
			integer,intent(in)::                                N
			integer::                							i,info
			complex*16,intent(in)::             				A(N,N)
			complex*16,allocatable,dimension(:,:)::             B
			integer,allocatable,dimension(:) ::                	ipiv
			complex*16::                                         Det
			real*8,intent(out)::                            lnDet



		!------------------------------------------------------------
		!https://software.intel.com/en-us/forums/topic/309460
		!------------------------------------------------------------
		
		  Allocate(B(N,N),ipiv(N))
		
		
		        B = A

        		call ZGETRF(N,N,B,N,ipiv,info)
        		
        		if (info .eq.0) then
				!Det = dcmplx(1.d0,0.d0)
				lnDet = dcmplx(0.d0,0.d0)
				
			do i= 1,N
				
					if (ipiv(i).ne.i) then
					!print*, B(i,i)  !--------------------------------<--------
					!Det = -Det*B(i,i)
					lnDet = lnDet + cdlog(-B(i,i))
					else
					!Det = Det*B(i,i)
					lnDet = lnDet + cdlog(B(i,i))
					end if
					
					!print*, lnDet
					
					
			 end do
				
				
				end if
			


        Deallocate(B,ipiv)
		return
		end
		
		                    
 
!======================================== lnDetm - hermitian =============================================

! To calculate determinant of Hermitian matrices!
	subroutine lnDetm(A,lnDet,N)
			implicit none
			integer,intent(in)::                                N
			integer::                							i,info,flag
			complex*16,intent(in)::             				A(N,N)
			complex*16::                                        B(N,N)
			integer ::                                      	ipiv(lln2)
			complex*16::                                        Det,work(2*lln2)
			real*8,intent(out)::                                lndet

  
            call hermitian_checker(A,n,'Error at lnDetm')
		
		        B = A

        	call ZHETRF('U',N,B,N,ipiv,work,2*N,info) ! For hermitian matrix
        		
        		
     		lndet=0.d0
      		flag=0
      		do i=1,N
         		if(flag.eq.0) then
            		if(ipiv(i).lt.0) then
              		   flag=1
              		   lndet = lndet+dreal(cdlog(B(i,i)*B(i+1,i+1) &
    		                    & -B(i,i+1)*dconjg(B(i,i+1))))
           		    else
               		   lndet = lndet+dreal(cdlog(B(i,i)))
            		endif
        	   else
           		   flag=0
         	  endif
         	 
     	   enddo

        		
        	
			return
			end
!======================================= detm- general ==============================================
	subroutine Detm(A,Det,N)
			implicit none
			integer,intent(in)::                                N
			integer::                							i,info
			complex*16,intent(in)::             				A(N,N)
			complex*16,allocatable,dimension(:,:)::             B
			integer,allocatable,dimension(:) ::                	ipiv
			complex*16,intent(out):: Det

		    Allocate(B(N,N),ipiv(N))
		
		
		        B = A

        		call ZGETRF(N,N,B,N,ipiv,info) !For general complex matrix
        		
        		if (info .eq.0) then
				Det = 1.d0
				
				
				do i= 1,N
					if (ipiv(i).ne.i) then
					Det = -Det*B(i,i)
					else
					Det = Det*B(i,i)
					end if
				end do
				
				
				end if
			


        Deallocate(B,ipiv)
		return
		end
		
!==========================================gnl===========================================
	subroutine Detm_gnl(A,Det,N)
			implicit none
			integer,intent(in)::                                N
			integer::                							i,info
			complex*16,intent(in)::             				A(N,N)
			complex*16,allocatable,dimension(:,:)::             B
			integer,allocatable,dimension(:) ::                	ipiv
			complex*16,intent(out):: Det



		!------------------------------------------------------------
		!https://software.intel.com/en-us/forums/topic/309460
		!------------------------------------------------------------
		
		  Allocate(B(N,N),ipiv(N))
		
		
		        B = A

        		call ZGETRF(N,N,B,N,ipiv,info)
        		
        		if (info .eq.0) then
				Det = 1.d0
				
				
				do i= 1,N
				
					if (ipiv(i).ne.i) then
					Det = -Det*B(i,i)
					else
					Det = Det*B(i,i)
					end if
					
				
				end do
				
				
				end if
			


        Deallocate(B,ipiv)
		return
		end
!=========================================  Checkers ========================================= 

 subroutine hermitian_checker(A,N,key)
     implicit none
     integer::p,q,N
     complex*16:: A(N,N)
     character*1::   key
    
            
                
           
                do p = 1,N
                    do q = 1,N
                    
                      if(cdabs(A(p,q) - dconjg(A(q,p))) .gt. absT) then 
                      		print*, key,p,q, A(p,q) , dconjg(A(q,p))
                      end if
                      
                    end do
                end do
              

      return
      end       
      
      
!================================================================================== 
 subroutine Unitary_checker(A,N,key)
     implicit none
     integer::          p,q,N
     complex*16::       A(N,N),HA(N,N),II(N,N)
     character*1::   key
     
                
           call getH(A,HA,N)
           II = matmul(HA,A)
                  
                  
                  do p = 1,N
                  do q = 1,N
                    if (p .eq. q) then
                    	if (cdabs(II(p,q) -1.d0) .gt. absT) then
                    		print*,key,p,q,II(p,q)
                    	end if
                   else 
                    	if (cdabs(II(p,q)) .gt. absT) then
                            print*,key,p,q,II(p,q)
                        end if
                  end if  
                  
                  end do
                  end do
                     
                

      return
      end                   



!==================================================================



 
   
   end module Algebra      
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         