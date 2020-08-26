!====================================================================================
!This module calculates Hermitian Dirac Matrix from configuration UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 
 
 MODULE hdm
 
 
 !---------------------
 use Data_set
 use Algebra
 !----------------------
 
 
IMPLICIT NONE
CONTAINS

              

!========================   model- X  ====================================================

subroutine HDX_Matrix(UU,mass,HD)
     implicit none
     integer::										p,q,j,k,r,kk,rj,cj,rd,cd,kc,kr,krd,kcd
     complex*16,intent(in):: 						UU(2,L,L,N,N)
     complex*16,intent(out):: 						HD(lln2,lln2)
     complex*16::                                   ai,oo,oi,ut,utd,ux,uxd
     complex*16::                                   x1,x2,x3,x4,xd1,xd2,xd3,xd4
     complex*16::                                   t1,t2,t3,t4,td1,td2,td3,td4
     real*8,intent(in)::                            mass
     
     
     
                        ai = dcmplx(0.d0,1.d0)
                        oi = dcmplx(1.d0,0.d0)
                        oo = dcmplx(0.d0,0.d0)
                       
                       
     					
           					x1 = ai
     					    x2 = ai
     					    x3 = -ai
     					    x4 = -ai
     					    
     					    xd1 = -ai
     					    xd2 = ai
     					    xd3 = -ai
     					    xd4 = ai
     					    
     					    t1 = oo
     					    t2 = oo
     					    t3 = -2.d0*ai
     					    t4 = oo
     					    
     					    td1 = oo
     					    td2 = 2.d0*ai
     					    td3 = oo
     					    td4 = oo
     					    
     		
                    
    				do p = 1,lln2
    				do q = 1,lln2
     				   HD(p,q) = dcmplx(0.d0,0.d0)
     				end do   
     				end do
                      
                    
                    
                    
                    
 ! k -Diagonal-Bloc ------------------------------------------------------------------- 
     
     do j  = 1,L
     do k  = 1,L  
                   
     				
     				kk = (k-1)*ln2 
     				
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
     
  !Non - Diagonal Block ---------------------------------------------------------------
     
    do k = 1,L
    
             if ( k .lt. L) then 
             
                 ! general term
             		kr = (k-1)*ln2
             		kc = k*ln2
             
             		krd = k*ln2
             		kcd = (k-1)*ln2
             		
             else if (k .eq. L) then
             
                ! L -corner-bloc 
                    kk = (L-1)*ln2 
                  
            end if
            
    
    do j = 1,L
    
             ! general term
             
             rj = (j-1)*N
             cj = (j-1)*N
             
             do p = 1,N
             do q = 1,N
             
             if (p .eq. q) then
             	HD(kr + rj + ln + p, kc + cj + q) = t3*(-0.5d0)
             	HD(krd + rj + p, kcd + cj + ln + q) = td2*(-0.5d0)
             end if
             
           
             
   			!L - corner block	
    		
     				
     			ut = UU(2,j,L,p,q)
         		HD(  kk + ln + rj + p, 0 + cj + q)  =   t3*(-0.5d0)*ut
         		utd = dconjg(UU(2,j,L,q,p))
         		HD(  0 + rj + p,  kk + ln +  cj + q)  =   td2*(-0.5d0)*utd
         				             
    		end do
    		end do 
    		
    		
    end do
    end do
    
   
!Mass term in HD  ------------------------------------------------------------------
     				    
     	
        do k  = 1,L  
     	do j  = 1,L
     				
     				kk = (k-1)*ln2 
     				
    				rj = (j-1)*N
    				cj = (j-1)*N
    				
    				do p = 1,N
     				do q = 1,N
     				
     				
     				if (p .eq. q) then
         			HD(kk + rj + p, kk + cj + ln + q)  =  HD(kk + rj + p, kk + cj + ln + q) + ai*(2.d0 + mass)
         		    HD(kk + rj + ln + p, kk + cj + q) =  HD(kk + rj + ln + p, kk + cj  + q) - ai*(2.d0 + mass)
         		    end if
         		    
         		    
    				end do
    				end do
     
    	end do
    	end do 
     				   
     				
     				
     				 call hermitian_checker(HD,lln2,'error on hermiticity at HD-X matrix')
     				
     				
                      
  return
  end                    

!======================================= model - Z =====================================================

subroutine HDZ_Matrix(UU,mass,HD)
     implicit none
     integer::										p,q,i,j,k,it,ix,ixp,ixn,is,isp,isn,itp,itn
     complex*16,intent(in):: 						UU(2,L,L,N,N)
     complex*16,intent(out):: 						HD(lln2,lln2)
     complex*16::                                   ai,ut,utd,ux,uxd
     real*8,intent(in)::                            mass
     
    




       ai = dcmplx(0.d0,1.d0)
       
       

      do i=1,lln2
         do j=1,lln2
            HD(i,j)=0.d0
         enddo
      enddo



      do i=1,lln
         HD(i,i+lln)=ai*(2.d0+mass)
         HD(i+lln,i)=-ai*(2.d0+mass)
      enddo




      do it=1,L
      do ix=1,L
      
            ixp=ix+1
            if(ix.eq.L) ixp=1
            ixn=ix-1
            if(ix.eq.1) ixn=L
            is=(((it-1)*L+ix)-1)*N
            isp=(((it-1)*L+ixp)-1)*N
            isn=(((it-1)*L+ixn)-1)*N
            
            
            do i=1,N
            do j=1,N
            
                  ux=-0.5d0*ai*UU(1,ix,it,i,j)
                  
                  HD(is+i,isp+j) = HD(is+i,isp+j)+ux
                  HD(is+i+lln,isp+j) = HD(is+i+lln,isp+j)-ux
                  HD(is+i,isp+j+lln) = HD(is+i,isp+j+lln)+ux
                  HD(is+i+lln,isp+j+lln) = HD(is+i+lln,isp+j+lln)-ux
                  
                  uxd=-0.5d0*ai*dconjg(UU(1,ixn,it,j,i))
                  
                  HD(is+i,isn+j) = HD(is+i,isn+j)-uxd
                  HD(is+i+lln,isn+j) = HD(is+i+lln,isn+j)-uxd
                  HD(is+i,isn+j+lln) = HD(is+i,isn+j+lln)+uxd
                  HD(is+i+lln,isn+j+lln) = HD(is+i+lln,isn+j+lln)+uxd
                  
            end do
            end do
            
            
      end do
      end do
      

      do it=1,L
      do ix=1,L
      
      
      
            itp=it+1
            if(it.eq.L) itp=1
            itn=it-1
            if(it.eq.1) itn=L
            is=(((it-1)*L+ix)-1)*N
            isp=(((itp-1)*L+ix)-1)*N
            isn=(((itn-1)*L+ix)-1)*N
            
            if(it.ne.L) then
            
               do i=1,N
                    HD(is+i+lln,isp+i) = HD(is+i+lln,isp+i)+ai
               enddo
               
            else
            
               do i=1,n
                  do j=1,n
                     ut = UU(2,ix,L,i,j)
                     HD(is+i+lln,isp+j) = HD(is+i+lln,isp+j)+ai*ut
                  enddo
               enddo
               
            endif
            
            if(it.ne.1) then
            
               do i=1,n
                   HD(is+i,isn+i+lln) = HD(is+i,isn+i+lln)-ai
               enddo
               
            else
            
               do i=1,n
                  do j=1,n
                     utd = dconjg(UU(2,ix,L,j,i))
                     HD(is+i,isn+j+lln) = HD(is+i,isn+j+lln)- &
                                          &ai*utd
                  enddo
               enddo
               
            endif
            
            
            
            
      end  do
      end  do

  	
     				
     call hermitian_checker(HD,lln2,'error on hermiticity at HD-Z matrix')
     				
     				
!======================================================================================    								
                      
  return
  end                    

         
         end module hdm
         
         
         
         
         
         