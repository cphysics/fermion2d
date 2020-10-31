


!=========================================================================================
!This module calculates different things: like Action,Forces ,Dynamical quantities from 
!configuration UU
!=========================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!=========================================================================================
 
 
 
 MODULE Calculator
 
 !-----------------------------------------------
 use Data_set
 use Algebra
 use hdm
 use Functionary
 use Printer
 !-----------------------------------------------
 
 
IMPLICIT NONE
CONTAINS       
      
        
        
 
!=============================CONTAINS===================================



                     
!==================================<><><>=================================================                              
!                            SECTION - Action
!==================================<><><>=================================================  



      subroutine Action(S,S1,S2,S3,UU,mass,PI)
          implicit none
          integer::						 j,k,p,q
       
          complex*16,dimension(N,N)::   U1,U2,U3,U4,HU1,HU2,HU12,HU34,U12,U34,A,AA,Mtx,AH
          complex*16::                  trc_sge,trc
          complex*16:: 					HD(lln2,lln2),DD(lln2,lln2)
          real*8:: 						SK,SGE(L,L),sx,SUG,lnDet
     
          complex*16,intent(in)::		UU(2,L,L,N,N),PI(2,L,L,N,N)
          real*8,intent(out):: 			S,S1,S2,S3
          real*8,intent(in)::                                       mass 
      
        
          SK = 0.0
          S = 0.0
        
          
          
        ! Kinetic energy part using all PI-------
          
          do j = 1,L
          do k = 1,L
          
                       do p = 1,N
                       do q = 1,N
                       A(p,q) = PI(1,j,k,p,q) 
                       end do
                       end do
                       call getH(A,AH,N)
        			   AA = Matmul(A,AH) 
         			   call trace(AA,N,trc) 
         			   SK = SK + real(trc)
         
         end do
         end do
       
         
          do j = 1,L
        
          
                       do p = 1,N
                       do q = 1,N
                       A(p,q) = PI(2,j,L,p,q) 
                       end do
                       end do
                       
                       call getH(A,AH,N)
        			   AA = Matmul(A,AH) 
        			  
         			   call trace(AA,N,trc) 
         			   SK = SK + real(trc)
         
         end do
     
       
         
        ! potential energy (I)  Gauge field part only---------
           
          
          do j = 1,L
          do k = 1,L-1
          
            		do p = 1,N
            		do q = 1,N
            		U1(p,q)=UU(1,j,k,p,q)
            		U2(p,q)=UU(1,j,modi(k),p,q)
            		end do
            		end do
            
            		call getH(U1,HU1,N)
            		call getH(U2,HU2,N)
            
            		Mtx  = matmul(U1,HU2) + matmul(U2,HU1)
            
            		call trace(Mtx,N,trc_sge)
            		SGE(j,k) = real(trc_sge)
          
          end do
          end do
          
          
          do j = 1,L
           
                   do p = 1,N
                   do q = 1,N
                   U1(p,q) = UU(1,j,L,p,q) 
                   U2(p,q) = UU(2,modi(j),L,p,q)
                   U3(p,q) = UU(2,j,L,p,q)
                   U4(p,q) = UU(1,j,1,p,q)
                   end do
                   end do
                   U12 = matmul(U1,U2)
                   U34 = matmul(U3,U4)
                   call getH(U12,HU12,N)
                   call getH(U34,HU34,N)
                   Mtx = matmul(U12,HU34) + matmul(U34,HU12)
                   call trace(Mtx,N,trc_sge)
                   SGE(j,L) =  real(trc_sge)
                   
        end do
        
        
          SUG = 0.0
          
          do j = 1,L
          do k = 1,L
                    
                       SUG  =  SUG + SGE(j,k)

           end do
           end do
        
         


          S1 = (-0.5d0)*SK
          S2 = (bb*dfloat(N)*SUG)
          
          
          if (fem .eq. 0) then
                S3 = 0.0
          else if (fem .eq. 1) then
          
          
          
          	    if (model .eq. 'X') then 
          			call HDX_matrix(UU,mass,HD)
          		else if (model .eq. 'Z') then
          			call HDZ_matrix(UU,mass,HD)
          		end if
          		
          		
          		
          		call lnDetm(HD,lnDet,lln2)
          		
                S3 = Nf*lnDet
                
           end if
           
           S = S1 + S2 + S3
           
           
          return  
          end
 
                     
!==================================<><><>=================================================                              
!                            SECTION - Force 
!==================================<><><>=================================================                               
                              


!=========================================================================================
!                  Gauge Force             
!=========================================================================================
         
Subroutine   Force_gauge(UU,Force,j,k,key)
         implicit none
         integer::                                 p,q 
         integer,intent(in)::                      j,k 
         complex*16,intent(in)::                   UU(2,L,L,N,N)
         complex*16::                              trc
         complex*16,dimension(N,N)::               U1,U2,U3,U4,U5,U6,HU1,HU2,HU3,HU4,HU5,U,F,FG
         Complex*16, intent(out)::                 Force(N,N)
         character*1::                             key
     
        
       
        if (key .eq. 'x') then    
        
                                ! Part - I (Get Fige)
                                
                   if (k .eq. 1) then
                             
                                 do p = 1,N
                                 do q = 1,N
                                 U1(p,q) = UU(1,j,2,p,q) 
                                 U2(p,q) = UU(2,modi(j),L,p,q)
                                 U3(p,q) = UU(1,j,L,p,q)
                                 U4(p,q) = UU(2,j,L,p,q)
                                 end do
                                 end do
                                 
                                 call getH(U1,HU1,N)
                                 call getH(U2,HU2,N)
                                 call getH(U3,HU3,N)
                                 
                                 F = (HU1 + matmul(HU2, matmul(HU3,U4)))
                                
                                 
                  else if (k .eq. L) then
                             
                                 do p = 1,N
                                 do q = 1,N
                                 U1(p,q) = UU(2,modi(j),L,p,q) 
                                 U2(p,q) = UU(1,j,1,p,q)
                                 U3(p,q) = UU(2,j,L,p,q) 
                                 U4(p,q) = UU(1,j,L-1,p,q)
                                 end do
                                 end do
                                 
                                 call getH(U2,HU2,N)
                                 call getH(U3,HU3,N)
                                 call getH(U4,HU4,N)
                                 
                                 F = matmul(U1,matmul(HU2,HU3)) + HU4
                               
                  else
                             
                                do p = 1,N
                                do q = 1,N
                                U1(p,q) = UU(1,j,modi(k),p,q)
                                U2(p,q) = UU(1,j,fun(k),p,q) 
                                end do
                                end do
                                
                                call getH(U1,HU1,N)
                                call getH(U2,HU2,N)
                                
                                F = HU1+HU2
                             
                 end if    
                 
                  
                            
                                 
                             	do p = 1,N
                             	do q = 1,N    
                                U(p,q) = UU(1,j,k,p,q)! this is hero - U_i(x,k)  
                             	end do
                             	end do
                         
                            
   else if (key .eq. 't') then  
   
   
   
                              ! Part - I (Get Fdge)
                             
                            	do p = 1,N
                            	do q = 1,N
                            	U1(p,q) = UU(1,j,1,p,q)
                            	U2(p,q) = UU(2,modi(j),L,p,q) 
                            	U3(p,q) = UU(1,j,L,p,q) 
                            	U4(p,q) = UU(1,fun(j),1,p,q)
                            	U5(p,q) = UU(2,fun(j),L,p,q) 
                            	U6(p,q) = UU(1,fun(j),L,p,q) 
                            	end do
                            	end do
                                 
                            	call getH(U2,HU2,N)
                            	call getH(U3,HU3,N)
                            	call getH(U4,HU4,N)
                            	call getH(U5,HU5,N)
                                 
                                 
                            	F = matmul(U1,matmul(HU2,HU3)) + matmul(HU4,matmul(HU5,U6))
                             
                             
                            
                                
                             	do p = 1,N
                             	do q = 1,N    
                                U(p,q) = UU(2,j,L,p,q)  ! this is hero - U_d(x)  
                             	end do
                             	end do
                             
                           
   
   
   
   end if
   
                                ! Part - II (Get Force)
                              
   
                              
                                !multiply by i
                                do p = 1,N
                                do q = 1,N
                                F(p,q) = dcmplx(0.d0,1.d0)*F(p,q)
                                end do
                                end do
                               
                             	F = matmul(U,F) 
                             	
                             	!normalized trace
                             	trc = dcmplx(0.d0,0.d0)
                             	
                             	do p = 1,N
                             	trc = trc + F(p,p)          
                                end do
                                
                             	trc = trc/dfloat(N)    
                             	
                             	
                             	!impose tracelessness
                             	do p = 1,N
                             	      F(p,p) = F(p,p)- trc
                                end do
                                
                                !add h.c.
                                do p = 1,N
                                do q = 1,N
                                FG(p,q) = F(p,q) + dconjg(F(q,p))
                             	end do
                             	end do
                             	
                             	Force = FG
                             
                         	
                return
                end
               






!=========================================================================================
!        Fermion Force - Model -'X'                   
!=========================================================================================
     
  


subroutine Force_fermionx(U,IHD,Force,j,k,key)

     implicit none
     integer::                                   rj,cj,kk,rr,cc,p,q
     integer,intent(in)::                        j,k
     character*1,intent(in)::                    key  
     complex*16::                  				 sm
     complex*16,intent(in)::                     IHD(lln2,lln2),U(N,N)
     complex*16,dimension(N,N) ::                F,FF,mIHD
     Complex*16, intent(out)::                   Force(N,N)
  
                   
                    
     if (key .eq. 'x') then               
     				
     				! k -diagonal-bloc 
     				
     				kk = (k-1)*ln2 
     				
     				!j-mini-block
     				if (j .lt. L) then
    				rj = (j)*N
    				cj = (j-1)*N
    				else if (j .eq. L) then
    				rj =  0
    				cj = (L-1)*N
    				end if
    				
    				
    				do p = 1,N
     				do q = 1,N
     				
         				mIHD(p,q) =   0.5d0*(IHD(kk + rj + p, 		kk + cj + q) &
         				                     & + IHD(kk + rj + ln + p, kk + cj  + q) &
         				              		 & - IHD(kk + rj + p, 		kk + cj + ln + q) &
         				             		 & - IHD(kk + rj + ln + p, kk + cj + ln + q)) 
         				             		 
         				
         				       
         				
    				end do
    				end do
    				
    else if (key .eq. 't') then
    				
    				! L -corner-bloc 
    				
     				kk = (L-1)*ln2 
     				
     				!j-mini-block
     				
    				rj = (j-1)*N
    				cj = (j-1)*N
    				
    				do p = 1,N
     				do q = 1,N
     				
         				mIHD(p,q) =  -IHD(  0 + rj + p, kk + cj + ln + q) 
         			            
    				end do
    				end do
    		
     end if   
     
     
     
     
                      F = matmul(U,mIHD)
                       
                       
                        
                        sm = 0.d0
                        do p = 1,N
                        sm = sm+F(p,p)
                        end do
                        
                        sm = sm/float(N)
                        
                        do p = 1,N
                        F(p,p) = F(p,p) - sm
                        end do
                      
                        
                       
                        do p = 1,N
                        do q = 1,N
                        FF(p,q) = F(p,q) + dconjg(F(q,p))
                        end do
                        end do
                        
                        Force = FF
                   
                  
     				return 
     				end


!=========================================================================================
!        Fermion Force - Model 'Z'                     
!=========================================================================================

     
subroutine Force_fermionz(U,IHD,Force,j,k,key)

     implicit none
     integer::                                   it,ix,ixp,is,isp,p,q
     integer,intent(in)::                        j,k
     character*1,intent(in)::                    key  
     complex*16::                  				 sm
     complex*16,intent(in)::                     IHD(lln2,lln2),U(N,N)
     complex*16,dimension(N,N) ::                F,FF,mIHD
     Complex*16, intent(out)::                   Force(N,N)
  
    
    
    
       if (key .eq. 'x') then 
            
      		it = k
      		ix = j
     
            ixp = ix+1
            
            if(ix.eq.L) ixp = 1
            is = (((it-1)*L+ix)-1)*N
            isp = (((it-1)*L+ixp)-1)*N
            
            
            
            do p=1,N
            do q=1,N
            
                  mIHD(p,q) = 0.5d0*(IHD(isp+p,is+q)&
                      			 &+IHD(isp+p+lln,is+q)&
                      			 &-IHD(isp+p,is+q+lln)&
                      			 &-IHD(isp+p+lln,is+q+lln))
                      			 
              
              
            end do
            end do
            
          
            
     else if (key .eq. 't') then
      
      
      

         ix = j
         is = (((L-1)*L+ix)-1)*N
         isp = (ix-1)*N
         
         do p = 1,n
         do q = 1,n
               mIHD(p,q) = -IHD(isp+p,is+q+lln)
         end do
         end do
         
         
   end if        
               
                    
     
     
                      F = matmul(U,mIHD)
                       
                       
                        
                        sm = dcmplx(0.d0,0.d0)
                        do p = 1,N
                        sm = sm+F(p,p)
                        end do
                        
                        sm = sm/float(N)
                        
                        do p = 1,N
                        F(p,p) = F(p,p) - sm
                        end do
                      
                        
                       
                        do p = 1,N
                        do q = 1,N
                        FF(p,q) = F(p,q) + dconjg(F(q,p))
                        end do
                        end do
                        
                        Force = FF
     
                
     				return 
     				end


 !====================================================================    

                     



      subroutine Plaquette(SUG,UU,mass)
      
          implicit none
          integer::						j,k,p,q
          complex*16,dimension(N,N)::   U1,U2,U3,U4,HU1,HU2,HU12,HU34,U12,U34,A,AA,Mtx
          complex*16::                  trc_sge,trc
          real*8:: 						SGE(L,L),sx
          complex*16,intent(in)::		UU(2,L,L,N,N)
          real*8,intent(in)::           mass 
          real*8,intent(out)::		    SUG
      
        
     
       
         
        ! potential energy (I)  Gauge field part only---------
           
          
          do j = 1,L
          do k = 1,L-1
          
            		do p = 1,N
            		do q = 1,N
            		U1(p,q)=UU(1,j,k,p,q)
            		U2(p,q)=UU(1,j,modi(k),p,q)
            		end do
            		end do
            
            		call getH(U1,HU1,N)
            		call getH(U2,HU2,N)
            
            		Mtx  = matmul(U1,HU2) + matmul(U2,HU1)
            
            		call trace(Mtx,N,trc_sge)
            		SGE(j,k) = real(trc_sge)
          
          end do
          end do
          
          
          do j = 1,L
           
                   do p = 1,N
                   do q = 1,N
                   U1(p,q) = UU(1,j,L,p,q) 
                   U2(p,q) = UU(2,modi(j),L,p,q)
                   U3(p,q) = UU(2,j,L,p,q)
                   U4(p,q) = UU(1,j,1,p,q)
                   end do
                   end do
                   U12 = matmul(U1,U2)
                   U34 = matmul(U3,U4)
                   call getH(U12,HU12,N)
                   call getH(U34,HU34,N)
                   Mtx = matmul(U12,HU34) + matmul(U34,HU12)
                   call trace(Mtx,N,trc_sge)
                   SGE(j,L) =  real(trc_sge)
                   
        end do
        
        
          SUG = 0.0
          
          do j = 1,L
          do k = 1,L
                    
                       SUG  =  SUG + SGE(j,k)

           end do
           end do
        
         

         
           
          return  
          end
 

 
 !==================================================================        
         end module Calculator
         
         

         
         