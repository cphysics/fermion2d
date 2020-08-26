


!====================================================================================
!This module prints matrices from configuration UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 
 MODULE Test
 
 
        !----use------
        use Algebra
        use Mstart
        !---------------
        
        
       IMPLICIT NONE
       CONTAINS
           
!=====================================================                        
Subroutine Matrix_writer(A,dim) 
     implicit none                      
     integer:: dim,p,q
     complex*16:: A(dim,dim)                   
     real*8:: rl,im                   
     
    		do p = 1,dim
    		do q = 1,dim
    			rl = real(A(p,q))
    			im = aimag(A(p,q))
       			write(111,*) p,q,rl,im
    		end do
    		end do
                         
    return
    end                     


!=====================================================                         
 subroutine eigen_writer(A,dim)!file 888
	implicit none
	integer::    p,q,dim
	real*8::      rl,im
	complex*16:: A(dim,dim)
    complex*16:: EV(dim,dim),W(dim)
    
  
    	
     	call Eigen_finder(A,EV,W,dim)!Algebra
     
        !write eigenvalues in file "fort.888"
        do p = 1,dim
        rl = real(W(p))
        im = aimag(W(p))
         write(888,*) p, rl, im
        end do
        
        return
         end       
      
 !=====================================================                         
 subroutine Hermitian_eigen_tester(dim,tmax)!file 888
	implicit none
	integer,parameter::     beans = 200
	integer::    			p,q,dim,t,tmax,count(beans)
	real*8::     			rl,im,dy,len
	complex*16:: 			A(dim,dim),trc
    complex*16:: 			EV(dim,dim),W(dim)
    
    
    len = 10.0
   
    dy = len/float(beans)
    do p = 1,beans
    count(p) = 0
    end do
  
t = 1

11    	call H_generator(A,dim)
        !call trace(A,dim,trc)
        !print*,trc
     	call Eigen_finder(A,EV,W,dim)!Algebra
     
        !write eigenvalues in file "fort.808"
        
        do p = 1,dim
        rl = real(W(p))
        
        !im = aimag(W(p))
        write(404,*) p, rl, im
        
         do q = 1,beans
            if (rl .gt. (-(len/2.0) + (dy*(q-1)))) then
                if(rl .lt. (-(len/2.0) + (dy*(q)))) then
                  count(q) = count(q) + 1
                  end if
            end if
        end do 
        end do
        
       
t = t+1
        
        
        if (t .lt. tmax) then
        go to 11
        else
        do q = 1,beans
        write(808,*) q,count(q)
        end do
        end if
        
        
        return
         end       
             
    !==================================================================    
  
    
   

	subroutine Dirac_eig_writer(DD)!file 444
	implicit none
	integer::                                p
	real*8::                                 rl,im
	complex*16,intent(in)::                  DD(lln2,lln2)
	complex*16,allocatable,dimension(:,:)::  EV
	complex*16,allocatable,dimension(:)::  W
	
	
     Allocate(EV(lln2,lln2),W(lln2))
 
     
    
        call Eigen_finder(DD,EV,W,lln2)
        
 
        do p = 1,lln2
             rl = real(W(p))
             im = aimag(W(p))
             write(444,*) p, rl, im
             !print*,p, rl, im
        end do
        
        
       Deallocate(EV,W)
     
        return 
        end


!==================================================================



subroutine Dirac_Matrix(UU,HD)
     implicit none
     integer::										p,q,j,k,r,kk,rj,cj,rd,cd
     complex*16,intent(in):: 						UU(2,L,L,N,N)
     complex*16,intent(out):: 						HD(lln2,lln2)
     complex*16::                                   ai,oo,oi,ut,utd,ux,uxd
     complex*16::                                   x1,x2,x3,x4,xd1,xd2,xd3,xd4
     complex*16::                                   t1,t2,t3,t4,td1,td2,td3,td4
     character*2::                                  key
     
     
     
                        ai = dcmplx(0.d0,1.d0)
                        oi = dcmplx(1.d0,0.d0)
                        oo = dcmplx(0.d0,0.d0)
                       
    					
     					    x1 = oi
     					    x2 = -oi
     					    x3 = -oi
     					    x4 = oi
     					    
     					    xd1 = oi
     					    xd2 = oi
     					    xd3 = oi
     					    xd4 = oi
     					    
     					    t1 = oo
     					    t2 = oo
     					    t3 = oo
     					    t4 = 2.d0*oi
     					    
     					    td1 = 2.d0*oi
     					    td2 = oo
     					    td3 = oo
     					    td4 = oo
     					    
                 
    				do p = 1,lln2
    				do q = 1,lln2
     				   HD(p,q) = dcmplx(0.d0,0.d0)
     				end do   
     				end do
                      
   
     
     do j  = 1,L
     do k  = 1,L  
                   
     				
     				kk = (k-1)*ln2 ! k -diagonal-bloc 
     				
     				!j-mini-block
     				if (j .lt. L) then
    				rj = (j-1)*N
    				rd = (j)*N
    				cj = (j)*N
    				cd = (j-1)*N
    				else if (j .eq. L) then
    				rj = (L-1)*N
    				rd = 0
    				cj = 0
    				cd = (L-1)*N
    				end if
    				
    				
    				do p = 1,N
     				do q = 1,N
     				
     				ux = UU(1,j,k,p,q)
     				
         			HD(kk + rj + p,kk + cj + q) =  x1*(-0.5d0)*ux
         			HD(kk + rj + p, kk + cj + ln + q)  =  x2*(-0.5d0)*ux
         		    HD(kk + rj + ln + p, kk + cj  + q)  =  x3*(-0.5d0)*ux
         			HD(kk + rj + ln + p, kk + cj + ln + q)  =  x4*(-0.5d0)*ux
         			
         			uxd = dconjg(UU(1,j,k,q,p))
         			
         			HD(kk + rd + p, kk + cd + q)  =   xd1*(-0.5d0)*uxd
         			HD(kk + rd + p, kk + cd + ln + q)  =  xd2*(-0.5d0)*uxd
         		    HD(kk + rd + ln + p, kk + cd  + q) =   xd3*(-0.5d0)*uxd
         			HD(kk + rd + ln + p, kk + cd + ln + q) =   xd4*(-0.5d0)*uxd
         			
    				end do
    				end do
     
    end do
    end do 
     
    				
    do j = 1,L				
     				kk = (L-1)*ln2 ! L -corner-bloc 
     				
     				!j-mini-block
     				
    				rj = (j-1)*N
    				cj = (j-1)*N
    				
    				do p = 1,N
     				do q = 1,N
     				
     				ut = UU(2,j,L,p,q)
     				HD(  kk + ln + rj + p, 0 + ln + cj + q)  =  t4*(-0.5d0)*ut
         			utd = dconjg(UU(2,j,L,q,p))
     				HD(  0 + rj + p,  kk + cj + q)  =  td1*(-0.5d0)*utd
     				
    				end do
    				end do
    		
   
     
     
    end do
   
     
     
             !Mass term in Dirac Matrix
     
            
                     
    					do p = 1,lln2
     				   		HD(p,p) = HD(p,p) + (2.d0 + mass)
     					end do
     					
     		
     				if (key .eq. 'hd') then
     				 call hermitian_checker(HD,lln2,'error on hermiticity at HD-X matrix')
     				end if
     				
     				
     				
                      
  return
  end      
                     
                        


  !==================================================================    
   
   subroutine DtoC_eigner(Tau,C,ICD)

     implicit none
     integer::					p,q,K
     complex*16:: 				EV(ln2,ln2),W(ln2),Dig(ln2,ln2)
     complex*16,intent(in):: 	Tau(ln2,ln2)
     complex*16,intent(out):: 	C(ln,ln),ICD(ln,ln)
     
    
     
     			!Find eigen value of Tau matrix
     			call Eigen_finder(Tau,EV,W,ln2)
     
    			call IIN(Dig,2*ln,dcmplx(0.d0,0.d0))
    			call IIN(C,ln,dcmplx(0.d0,0.d0))
    			call IIN(ICD,ln,dcmplx(0.d0,0.d0))
    
     					do p = 1,2*ln
        					Dig(p,p) = W(p)
     					end do    

     					do p = 1,ln
        					C(p,p) = W(p)
     					end do 
        
     					do p = ln+1,2*ln
        					ICD(p-ln,p-ln) = W(p)
     					end do    
     

    			return 
    			end
                         
  !==================================================================

 subroutine C_ICD_checker(Tau)
	 implicit none
     integer::				p,q,K
     complex*16::			EV(ln2,ln2),W(ln2),Dig(ln2,ln2),C(ln,ln),ICD(ln,ln)
     complex*16,intent(in)::Tau(ln2,ln2)
     complex*16:: 			HC(ln,ln),IHC(ln,ln)
     
   
     
     			!Find eigen value of Tau matrix
     		
     			call Eigen_finder(Tau,EV,W,ln2)
     
    			call IIN(Dig,2*ln,dcmplx(0.d0,0.d0))
    			call IIN(C,ln,dcmplx(0.d0,0.d0))
    			call IIN(ICD,ln,dcmplx(0.d0,0.d0))
    
     					do p = 1,2*ln
        					Dig(p,p) = W(p)
     					end do    

     					do p = 1,ln
        					C(p,p) = W(p)
     					end do 
        
     					do p = ln+1,2*ln
        					ICD(p-ln,p-ln) = W(p)
     					end do    
                        
                        !Checking section----------------
                        
                        call geth(C,HC,ln)
                        call Inverse_finder(HC,IHC ,ln)
                        
                        do p = 1,ln
                        print*,p,IHC(p,p),ICD(p,p)
                        end do
                        

    			return 
    			end
                         

 !==================================================================

 subroutine Tau_eigen_writer(Tau)
	implicit none
	integer::    p,q,dim
	real*8::      rl,im
	complex*16,intent(in):: Tau(ln2,ln2)
    complex*16:: EV(ln2,ln2),W(ln2),Dig(ln2,ln2)
    
    
 
    	dim = ln2
    	
    	!Find eigen values of Tau
     	call Eigen_finder(Tau,EV,W,dim)!Algebra
     
        !write eigenvalues in file "fort.888"
        do p = 1,dim
        rl = real(W(p))
        im = aimag(W(p))
        write(888,*) p, rl, im
        end do
        
        return
         end


!=============================APENDIX===================================================

 !==================================================================
 !                    Determinant of Dirac Operator                 !
 !==================================================================    
 subroutine NewDetm(UU,NetDet)

     implicit none
     integer::						  p,q,k
     complex*16,intent(in):: 		  UU(2,L,L,N,N)
     complex*16,dimension(ln,ln)::    Bk,Ck
     complex*16::                     Det,Dt
     complex*16::                     Tau(ln2,ln2),II(ln2,ln2),ITau(ln2,ln2)
     complex*16,intent(out)::        NetDet  
     
    
              
          k = 1
                      
                       NetDet = dcmplx(1.d0,0.d0)
                       
                       !Bring Bk and Ck matrix from Gamma module
    101                call BkCk(K,UU,Bk,Ck)!Gamma
                       call Detm(Bk,Det,ln)
                       NetDet = NetDet*Det*(-1.d0)**(ln)
                       
       k = k+1
       
       
                     if (k .lt. L+1) then
                        go to 101
                     end if
                     
                  !Tau - det - section   
                     
                  call Tau_finder(UU,Tau)               
                  call IIN(II,2*ln,dcmplx(1.d0,0.d0))!Algebra  
                  
                  ITau =  II - Tau
                  call Detm(ITau,Dt,ln2)
                  NetDet = NetDet*Dt
                  
                     
                                       
       return
       end     
 !=============================================================================== 
        
         end module Test
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         