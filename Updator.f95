
!====================================================================================
!This module calculates use leapfrog approximation and HMC to update links in UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 


 MODULE Updator
 
 
 !--------------------
 use data_set
 use Algebra
 use Mstart
 use hdm
 use Calculator
 use Printer
 !---------------------
 
 
IMPLICIT NONE
CONTAINS


!=========================CONTAINS=======================================


                     
!==================================<><><>=================================================                              
!                            SECTION  - Leap - Frog  
!==================================<><><>=================================================                               
                              
             Subroutine Leap_Frog(UU,mass,PI,ihmc,printid)  
             implicit none
             real*8::                        S1,S2,S3,S,dtau,Si,Sf
             real*8::                        SGi,SGf,SKi,SKf,SFi,SFf
             integer::                   	 t,p,q
             complex*16,intent(in out):: 	 UU(2,L,L,N,N),PI(2,L,L,N,N)
             real*8,intent(in) ::            mass
             integer,intent(in) ::           ihmc,printid
          
          
          
           
             dtau = dt_lf
          
            
            
                     
t = 0

                      if (printid .eq.1) then
                      print*,'**********job-tracker: HMC-step = ', ihmc,'LF-run=' , t
                      end if
                      
                      
                      call Action(Si,SKi,SGi,SFi,UU,mass,PI)
                    
                        
t = 1
                    
                      if (printid .eq.1) then
                      print*,'**********job-tracker: HMC-step = ',ihmc,'LF-run=' , t
                      end if
                      
                      call Gforce(UU,PI,dtau*0.5d0)
                      
                      if (fem .eq. 1) then
                      call Fforce(UU,mass,PI,dtau*0.5d0)
                      end if
                      
                      
                 
                    
                     
 t = 2 
                    
                     
 
11                    call UNEW(UU,PI,dtau)

                      call Gforce(UU,PI,dtau)
                      
                      if(fem .eq. 1) then
                      call Fforce(UU,mass,PI,dtau) 
                      end if
                      
                      
                    
                      if (printid .eq.1) then
                      print*,'**********job-tracker: HMC-step = ',ihmc,'LF-run=' , t
                      end if
                 
                    
t = t+1  
               
  					 if (t .lt. LF_steps+2) then
    					go to 11
    				 end if
                  
t = LF_steps + 2  
  
                     
                      if (printid .eq.1) then
                      print*,'**********job-tracker: HMC-step = ',ihmc,'LF-run=' , t
                      end if
                      
  
                    call UNEW(UU,PI,dtau) 
                                 
                    call Gforce(UU,PI,dtau*0.5d0)
                     
                    if (fem .eq. 1) then
                    call Fforce(UU,mass,PI,dtau*0.5d0)
                    end if
                
                    
                   
                
  
  return
  end
  
    !============================  U-Updator ============================================= 
    
          
 Subroutine UNEW(UU,PI,dtau)
             implicit none
             integer::       									p,q,j,k
             real*8,intent(in) ::                               dtau
             complex*16,dimension(N,N):: 			            ExpH,tPI,iU,tU
             complex*16,intent(in):: 					        PI(2,L,L,N,N)
             complex*16,intent(in out):: 					    UU(2,L,L,N,N)
    
            
           
    

            !--------------- Vertical --- r = 1 ----------------
             do j = 1,L
             do k = 1,L
                        do p = 1,N
                        do q = 1,N
                        iU(p,q) = UU(1,j,k,p,q)  
                        tPI(p,q) = PI(1,j,k,p,q) 
                        end do
                        end do
                        
                        call U_Exponenter(tPI,ExpH,dtau)
                        tU = matmul(ExpH,iU) 
                        
                        do p = 1,N
                        do q = 1,N
                        UU(1,j,k,p,q) = tU(p,q)
                        end do
                        end do 
          end do
          end do
      
      !--------------- Horizontal -------- r = 2 ----------- 
           do j = 1,L
                        do p = 1,N
                        do q = 1,N
                        iU(p,q) = UU(2,j,L,p,q)  
                        tPI(p,q) = PI(2,j,L,p,q) 
                        end do
                        end do
                        
                        call U_Exponenter(tPI,ExpH,dtau)
                        tU = matmul(ExpH,iU)  
                        
                        do p = 1,N
                        do q = 1,N
                        UU(2,j,L,p,q) = tU(p,q)
                        end do
                        end do 
        end do
        
       return 
       end   
       
 !============================   U-Exp    ====================================== 
                            
            Subroutine U_Exponenter(tPI,ExpH,dtau)
            implicit none
            integer::                    			p,q
            real*8,intent(in) ::                    dtau
            complex*16,dimension(N,N)::             A,ExD,EV,HEV,D
            complex*16,intent(in)::                 tPI(N,N)
            complex*16,intent(out)::    			ExpH(N,N) 
                  
                  
                  
                     A = tPi
                     call Eigner(A,EV,D,N)
                      
                     do p = 1,N
                     do q = 1,N
                     
                     	if (p.eq.q) then
                        		ExD(p,q) = cdexp(dcmplx(0.d0,1.d0)*D(p,p)*dtau)
                     	else
                        		ExD(p,q) = dcmplx(0.d0,0.d0)
                     	end if
                     
                     end do   
                     end do
                    
                     call getH(EV,HEV,N) 
                  
                     ExpH = Matmul(EV,matmul(ExD,HEV))
                  
                    return
                    end               
                           
        
 !============================  U-Equation =============================================   
   
 Subroutine Gforce(UU,PI,dtau)
 
             implicit none
             integer::       									p,q,r,j,k
             real*8,intent(in) ::                               dtau
             complex*16,intent(in out):: 					    PI(2,L,L,N,N)
             complex*16,intent(in):: 					        UU(2,L,L,N,N)
             complex*16,dimension(N,N)::                		Force_g
             character*1::                              		x,t
             complex*16::                                       tst
            
           
            		
            		 do j = 1,L
            		 do k = 1,L
                               
                     				call Force_gauge(UU,Force_g,j,k,'x') 
                              
                     				do p = 1,N
                     				do q = 1,N
                     				 tst = PI(1,j,k,p,q)
                                	 PI(1,j,k,p,q) = PI(1,j,k,p,q) + (dtau*Force_g(p,q)*bb*dfloat(N))
                                	end do
                                	end do
                               
                    end do
                    end do
                    
            		do j = 1,L
            		
                     				call Force_gauge(UU,Force_g,j,L,'t') 
                     		
                     				do p = 1,N
                     				do q = 1,N
                     				 tst = PI(2,j,L,p,q)
                                	 PI(2,j,L,p,q) = PI(2,j,L,p,q) + (dtau*Force_g(p,q)*bb*dfloat(N)) 
                                	end do
                                	end do
                    end do
                    
                    
                 
       
        return 
        end    
        
        
                   
                       
!============================= PI - Equation =====================================                             
                                     
           
            Subroutine Fforce(UU,mass,PI,dtau)
            implicit none
            integer::                                  p,q,r,j,k
            complex*16,dimension(N,N)::                Force_g,Force_f,U
            complex*16,intent(in)::                    UU(2,L,L,N,N)
            complex*16::                               IHD(lln2,lln2) 
            complex*16::                               PIF(2,L,L,N,N) 
            real*8,intent(in) ::                       dtau
            complex*16,intent(in out):: 			   PI(2,L,L,N,N)
            character*1::                              x,t
            real*8,intent(in)::                        mass
         
          
                   
                     call  Inverse_IHD(UU,mass,IHD)
                   
                
            		
            		 do j = 1,L
            		 do k = 1,L
            		
            		           
            		           		do p = 1,N
                             		do q = 1,N    
                                	U(p,q) = UU(1,j,k,p,q)  
                             		end do
                             		end do
                              
                                
                     			    
                     			    if (model .eq. 'X') then 
                     					call Force_fermionx(U,IHD,Force_f,j,k,'x') 
                     				else if (model .eq. 'Z') then
                     					call Force_fermionz(U,IHD,Force_f,j,k,'x') 
                     				end if
                     				
                               
                                    do p = 1,N
                     				do q = 1,N
                                	 PI(1,j,k,p,q) = PI(1,j,k,p,q) + (dtau*Force_f(p,q)*Nf)
                                	end do
                                	end do
                    
                            
                               
                    end do
                    end do
                    
         
            		
            		do j = 1,L
            		
            		           
                                do p = 1,N
                             	do q = 1,N    
                                U(p,q) = UU(2,j,L,p,q) 
                            	end do
                             	end do
            		             
                     			     
                     			     if(model .eq. 'X') then
                     					call Force_fermionx(U,IHD,Force_f,j,L,'t')
                     				 else if(model .eq. 'Z') then
                     				    call Force_fermionz(U,IHD,Force_f,j,L,'t')
                     				 end if    
                     				 	
                     				 
                     			
                                    do p = 1,N
                     				do q = 1,N
                                	 PI(2,j,L,p,q) = PI(2,j,L,p,q) + (dtau*Force_f(p,q)*Nf)
                                	end do
                                	end do
                               
                    
                    
                     end do
                     
                   ! print*, '============================================================='
            	   
                     return
                     end
                 


                     

!==================================================================  
    
             Subroutine Inverse_IHD(UU,mass,IHD)
             implicit none
             integer::                            p,q,k
             complex*16,intent(in) ::             UU(2,L,L,N,N) 
             complex*16::						  DD(lln2,lln2), HD(lln2,lln2),cdet,ac,bc
             complex*16,intent(out):: 		      IHD(lln2,lln2)
             real*8,intent(in)::                  mass
         
                        
  				          
  				           
  				           
     				       
     				       if (model .eq. 'X') then
     				        	call HDX_matrix(UU,mass,HD)
     				       else if(model .eq. 'Z') then
     				        	call HDZ_matrix(UU,mass,HD) 
     				       end if
  				       
     				       
     				       !Check hermiticity 
     				       call hermitian_checker(HD,lln2,'error on HD before Inverse_IHD')
     				       
     				       call Inverser(HD,IHD,lln2)  
     				       
     				       call hermitian_checker(IHD,lln2,'error on IHD at Inverse_IHD')
     				       
     				      !check inverse 
     				      do p = 1,lln2
     				      do q = 1,lln2
     				         cdet = 0.d0
     				         do k = 1,lln2
     				           cdet = cdet + (HD(p,k)*IHD(k,q))
     				         end do
     				         if(p .eq. q) then
     				         cdet = cdet -1.d0
     				         if (cdabs(cdet) .gt. absT) then
     				         print*,"Error on IHD", p,q,cdet
     				         end if
     				         end if
     				      end do
     				      end do
     				  
     				    
     				    
     		return 
     		end		    
!==================================================================  
      
             Subroutine Leap_Frog_test(UU,mass,PI,ilf,drx)  
             implicit none
             real*8::                        S1,S2,S3,S,dtau,Si,Sf
             real*8::                        SGi,SGf,SKi,SKf,SFi,SFf
             integer::                   	 t,p,q,ilf,drx
             complex*16,intent(in out):: 	 UU(2,L,L,N,N),PI(2,L,L,N,N)
             real*8,intent(in) ::                mass
          
          
          
            if (drx .eq. 1) then
              dtau = dt_lf
            else if (drx .eq. -1) then
              dtau = -dt_lf
            end if
            
            
                     
t = 0

                      call Action(Si,SKi,SGi,SFi,UU,mass,PI)
                     
                      !---------------------------------------------
                      write(200,*) t,Si,SKi,SGi,SFi
                      print*,"==========",t,Si,SKi,SGi,SFi
                      !--------------------------------------------
           
           
                        
t = 1
                      call Gforce(UU,PI,dtau*0.5d0)
                      
                      if (fem .eq. 1) then
                      call Fforce(UU,mass,PI,dtau*0.5d0)
                      end if
                      
                      
                      
                     !---------------------------------------------
                      call Action(S,S1,S2,S3,UU,mass,PI)
                      write(200,*) t,S,S1,S2,S3
                      print*,"==========",t,S,S1,S2,S3
                     !--------------------------------------------
                    
                     
 t = 2 
 
11                    call UNEW(UU,PI,dtau)

                      call Gforce(UU,PI,dtau)
                      
                      if(fem .eq. 1) then
                      call Fforce(UU,mass,PI,dtau) 
                      end if
                      
                      
                      
                    !---------------------------------------------
                      call Action(S,S1,S2,S3,UU,mass,PI)
                      write(200,*) t,S,S1,S2,S3
                      print*,"==========",t,S,S1,S2,S3
                    !-----------------------------------------------
                    
t = t+1  
               
  					 if (t .lt. LF_steps+2) then
    					go to 11
    				 end if
                  
t = LF_steps + 2  
  
                    call UNEW(UU,PI,dtau) 
                                 
                    call Gforce(UU,PI,dtau*0.5d0)
                     
                    if (fem .eq. 1) then
                    call Fforce(UU,mass,PI,dtau*0.5d0)
                    end if
                
                    
                    
                     !---------------------------------------------
                      call Action(Sf,SKf,SGf,SFf,UU,mass,PI)
                      write(200,*) t,Sf,SKf,SGf,SFf
                      print*,"==========",t,Sf,SKf,SGf,SFf
                    !--------------------------------------------
                   
                   
                   
                   
                       
                  write(*,*) "Model=",model
                  write(*,*) LF_steps,"initial gauge action =", SGi,"final = ", SGf,"change = ", SGf -SGi
                  write(*,*) LF_steps,"initial kinetic action =", SKi,"final = ", SKf,"change = " , SKf -SKi
                  write(*,*) LF_steps,"initial fermion action =", SFi,"final = ", SFf,"change = ", SFf -SFi
                  write(*,*) '==============================================================================================='
                  write(*,*) "Initial Action =",Si
                  Write(*,*) "Final Action = ",Sf
                  Write(*,*) "-----------------------------------------------"
                  Write (*,*)"Change in Action = ",Sf -Si
                  
                  
                   
          
  
  return
  end
 !==================================================================     
  
  
  
                  
!==================================================================         
         end module Updator
         
         
         
         
         
         