



!====================================================================================
!This module calculates Tau matrix and its eigen values (C- matrix) from conf. UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 



 MODULE Partition
 
 
 !-----------------------
 use Data_set
 use Algebra
 use Gamma
 use Printer
 !------------------------
 
 
IMPLICIT NONE

 CONTAINS
 !================================================================


!-------------------------------------------------------------------------------------------
!This subroutine construct Tau matrix by making product of all Tk matrices and a TD matrix.
! These Tk and TD matrices have double the size of Bk,Ck Td matrices
!-------------------------------------------------------------------------------------------



subroutine Tau_finder(UU,mass,Tau)

     implicit none
     integer::						  p,q,K
     complex*16,intent(in):: 		  UU(2,L,L,N,N)
     complex*16,dimension(ln,ln)::    Tk11,Tk12,Tk21,Tk22,MTd,MTdD,II,Bk,Ck,IBk
     complex*16::                     Tk(L,ln2,ln2),TkT(ln2,ln2),TD(ln2,ln2)
     complex*16,intent(out)::         Tau(ln2,ln2)
     real*8,intent(in)::              mass
    
              
                        k = 1
                      
                      
        !STEP - I: Construct Tk matrix for for different K.
                       
                       
                       
                       !Bring Bk and Ck matrix from Gamma module
    101                call BkCk(K,UU,mass,Bk,Ck)!Gamma
                  
                 
                       !We need inverse of Bk 
                       call Inverse_finder(Bk,IBk,ln)!Algebra
                       !II = matmul(IBk,Bk)
                       !call print_matrix(II,ln)
                       
                       
                       !construction of four blocks 
                       TK11 = IBk
                       TK12 = matmul(IBk,Ck)
                       TK21 = matmul(Ck,IBk)
                       TK22 = Bk - matmul(Ck,matmul(IBk,Ck))
                
                
                      do p = 1,2*ln
                      do q = 1,2*ln
                            
                                    !Block - I
                                   if (p .lt. ln+1)then
                                       if (q .lt. ln+1)then
                                           TK(k,p,q) = TK11(p,q)
                                           end if
                                           end if
                                           
                                   !Block - II        
                                   if (p .lt. ln+1)then
                                       if (q .gt. ln)then
                                           TK(k,p,q) = -TK12(p,q-ln)
                                           end if
                                           end if
                                           
                                   !Block - III
                                  if (p .gt. ln)then
                                       if (q .lt. ln+1)then
                                           TK(k,p,q) = TK21(p-ln,q)
                                            end if
                                            end if
    
                                   !Block - IV
                                   if (p .gt. ln)then
                                       if (q .gt. ln)then
                                           TK(k,p,q) = TK22(p-ln , q - ln)
                                           end if
                                           end if
                     
                     end do
                     end do  
                     
                   
                     k = k+1
                     
                     if (k .lt. L+1) then
                        go to 101
                        
                     end if
                     
                   !Bring MTD and MTdD here from Gamma module
                   call  TdD(UU,MTd,MTdD)!Gamma
                   
                    
                   !expand MTd to TD
                   call IIN(TD,2*ln,dcmplx(0.d0,0.d0))!Algebra
                   
                   do p = 1,2*ln
                   do q = 1,2*ln
                                      
                                 !Diagonal - block - I
                                 if (p .lt. ln+1) then
                                       if (q .lt. ln+1) then
                                           TD(p,q) = MTd(p,q)
                                           end if
                                           end if
                                 !Diagonal block -II
                                 if (p .gt. ln) then
                                       if (q .gt. ln) then
                                           TD(p,q) = MTd(p - ln, q - ln)
                                            end if
                                            end if
                                          
                     
                     end do
                     end do
                     
                     
                     !----TD - det - check--------------------
                     !call Detm(TD,Det,2*ln)
                     !print*,Det
                     !------------------------------------
                     
    !STEP - II :  Multiply Tk for all k and from 1 to L and with TD  
                   
                   !Initiate Tau matrix for multiplication
                   call IIN(Tau,2*ln,dcmplx(1.d0,0.d0))!Algebra 
                   
            k = 1
                    
102                 do p = 1,2*ln
                    do q = 1,2*ln
                         !pick TK for specific k
                         TkT(p,q) = Tk(k,p,q)
                         
                    end do
                    end do  
                    
                     !---- TKT - det - check--------------------
                     !call Detm(TKT,Dt,2*ln)
                     !print*,Dt
                     !------------------------------------
                     
                     
                     !---- TKT - Hermiticity - check--------------------
                     !call Hermitian_checker(TKT,2*ln)
                     !------------------------------------
                     
                       
                    Tau = matmul(Tau,TkT)   
          k = k+1
                          
                    if (k .lt. L+1) then
                    go to 102   
                    else if(k .eq. L+1) then
                    !Finally multiply with TD
                    Tau = matmul(Tau,TD)
                    end if   
                       
                       
                       
                  
                       
                       
       return
       end                
                       
 
  !==================================================================    
   
   subroutine Tau_eigner(Tau,D)

     implicit none
     integer::					p,q,K
     complex*16:: 				EV(ln2,ln2),W(ln2)!,Dig(ln2,ln2)
     complex*16,intent(in):: 	Tau(ln2,ln2)
     complex*16,intent(out):: 	D(ln2)
     
    	
     			call Eigen_finder(Tau,EV,W,ln2)
     
    			
     			do p = 1,2*ln
        			D(p) = W(p)
     			end do    

     				
     

    			return 
    			end
                         
                         
                         
!==================================================================
subroutine  PZ20(UU,mass,Z20)
     
      implicit none
      integer::                        i,j,k,p,q
      complex*16,intent(in):: 		   UU(2,L,L,N,N)
      complex*16:: 	                   Tau(ln2,ln2),D(ln2)
      complex*16,intent(out):: 	       Z20
      real*8,intent(in)::              mass



       call Tau_finder(UU,mass,Tau)
       call Tau_eigner(Tau,D)

      
              Z20 = dcmplx(0,0)
       
             do p = 1,ln2
             do q = 1,ln2
                  Z20 = Z20 + 1/(D(p)*dconjg(D(q)))
             end do
             end do    
      
            Z20  = 2.d0*Z20

     return
     end
!======================================================================

subroutine  PZ31A(UU,mass,Z31ALL)
     
      implicit none
      integer::                        i,j,k,p,q,r1,r2,r3
      complex*16,intent(in):: 		   UU(2,L,L,N,N)
      complex*16:: 	                   Tau(ln2,ln2),D(ln2)
      complex*16:: 	                   Z31I,Z31II  ! ,Z31I_1,Z31I_2,Z31I_3,Z31I_net,Z31II_1,Z31II_cub
      complex*16,intent(out):: 	       Z31ALL
      real*8,intent(in)::              mass




       call Tau_finder(UU,mass,Tau)
       call Tau_eigner(Tau,D)


      
              Z31I = dcmplx(0,0)
              Z31II = dcmplx(0,0)
              !Z31I_1 = dcmplx(0,0)
              !Z31I_2 = dcmplx(0,0)
              !Z31I_3 = dcmplx(0,0)
              !Z31I_net = dcmplx(0,0)
              !Z31II_1 = dcmplx(0,0)
              !Z31II_cub = dcmplx(0,0)
              
              
 !---------------------------------Z_31_[111]-------------------------------------------------
       
        do r1 = 1,ln2
        
                    !Z31I_3 = Z31I_3 + 1/((dconjg(D(r1))**3))
        
        do r2 = 1,ln2
        
                    !Z31I_2 = Z31I_2 + 1/((dconjg(D(r1))**2)*dconjg(D(r2)))
        
        do r3 = 1,ln2
               
               
                    !Z31I_1 = Z31I_1 + 1/(dconjg(D(r1))*dconjg(D(r2))*dconjg(D(r3)))       
                  
               
               
                  if (r1 .ne. r2 .and. r1 .ne. r3 .and. r2 .ne. r3) then
                       Z31I = Z31I + 1/(dconjg(D(r1))*dconjg(D(r2))*dconjg(D(r3)))
                  end if
                  
                   
                  
                  
        end do
        end do
        end do   
        
        
         !Z31I_net =  Z31I_1 - (3.d0*Z31I_2) + (2.d0*Z31I_3)     
        
         !print*, Z31I_net,Z31I
        
         !Print*, "Ds_sq-from A:", Z31I_2
        
        
!-------------------------------Z_31_[1,2]----------------------------------
 
 
        do r1 = 1,ln2
            
               !Z31II_cub = Z31II_cub  +  1/(dconjg(D(r1))**3)
                
        do r2 = 1,ln2
         
               !Z31II_1 = Z31II_1 + 1/(dconjg(D(r1))*(dconjg(D(r2))**2))
          
           
              if (r1 .ne. r2) then
               Z31II = Z31II + 1/(dconjg(D(r1))*(dconjg(D(r2))**2))
               !Z31II_1 = Z31II_1 + 1/(dconjg(D(r2))*(dconjg(D(r1))**2))
              end if         
             
             
                  
        end do
        end do    
      
      
      
       !print*, Z31II,Z31II_1 - Z31II_cub
 !-------------------------------------------------------------------- 
  
        Z31ALL =  -8.d0/6.d0*Z31I -  2.d0*Z31II
  
        !print*, " Z31all-from A:", Z31ALL

     return
     end

!======================================================================

subroutine  PZ31B(UU,mass,Z31ALL)
     
      implicit none
      integer::                        i,j,k,p,q,r1,r2,r3
      complex*16,intent(in):: 		   UU(2,L,L,N,N)
      complex*16:: 	                   Tau(ln2,ln2),D(ln2)
      complex*16:: 	                   Ds,Ds_sq,Ds_cub
      complex*16,intent(out):: 	       Z31ALL
      real*8,intent(in)::              mass




       call Tau_finder(UU,mass,Tau)
       call Tau_eigner(Tau,D)
       

       Ds = dcmplx(0,0)
       Ds_sq = dcmplx(0,0)
       Ds_cub = dcmplx(0,0)
            
       
        do r1 = 1,ln2
        Ds = Ds +  1/(dconjg(D(r1)))
        Ds_sq = Ds_sq + 1/(dconjg(D(r1))**2)
        Ds_cub = Ds_cub + 1/(dconjg(D(r1))**3)
        end do   
        
         
         !print*,"Ds_sq-from B:", Ds_sq*Ds
         
          
       
        Z31ALL = (-(4.0/3.0)*Ds**3) + (2.0*Ds_sq*Ds) - (2.0/3.0)*(Ds_cub)
  
       ! print*, " Z31all-from B:", Z31ALL






     return
     end


 
!==================================================================
         
         end module Partition
         
         
         
         
     
     
     
     
     
     
     
     
     
     
     
         
         