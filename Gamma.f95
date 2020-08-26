
!====================================================================================
!This module calculates Bk and Ck matrices from configuration UU
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 
 MODULE Gamma
 
 
 !---------------------
 use data_set
 use Algebra
 use Mstart
 !---------------------
 
IMPLICIT NONE
CONTAINS


 !================================================================


    subroutine BkCk(K,UU,mass,Bk,Ck)
     implicit none
     integer::                                   p,q
     integer,intent(in)::                        K
     real*8,intent(in)::                         mass
     complex*16,intent(in)::                     UU(2,L,L,N,N)
     complex*16,allocatable,dimension(:,:) ::    IMD, MTi_k, MTi_kD
     complex*16,intent(out)::                    Bk(ln,ln),Ck(ln,ln)
     
     
     allocate( IMD(ln,ln), MTi_k(ln,ln),MTi_kD(ln,ln))
     
     !Construct II matrix with mass m and dimension d
     call IIN(IMD,ln,dcmplx((2.d0 + mass),0.d0))
     !construct MTk and MTkD matrix for given time -  K
     call TTi(K,UU,MTi_k,MTi_kD)
     
    
                do p = 1,ln
                    do q = 1,ln
                       Bk(p,q) = IMD(p,q) - (0.5d0)*(MTi_k(p,q) + MTi_kD(p,q))
                       Ck(p,q) = (0.5d0)*(MTi_k(p,q) - MTi_kD(p,q))
                    end do
               end do
               
    Deallocate( IMD,MTi_k,MTi_kD)    
    return
    end 
     
     
           
 !====================================================================
 !TTi calculates  T_i(k) and T_i(k)^D matrices for specific time K. Only links 
 !from space direction contributes
 !=====================================================================  
 
     subroutine TTi(K,UU,MTi_k,MTi_kD)
     implicit none
     integer::                 					  p,q,lr,lc
     integer,intent(in)::                         K
     complex*16,intent(in)::             		  UU(2,L,L,N,N)
     complex*16,allocatable,dimension(:,:)::      O
     complex*16,allocatable,dimension(:,:,:,:)::  Ti_k
     complex*16,intent(out)::                     MTi_k(ln,ln),MTi_kD(ln,ln)



     Allocate(Ti_k(L,L,N,N),O(N,N))
     
     call IIN(O,N,dcmplx(0.d0,0.d0))
    
    
      !Initiate multi-dim-array Ti_k, each (lr,lc) block includes a NxN matrix
      !There are K = L number of Ti
      
        do lr = 1,L
            do lc = 1,L
                do p = 1,N
                    do q = 1,N
                       
                            Ti_k(lr,lc,p,q) = O(p,q)
                        
                     end do
                  end do
              end do
        end do
        
        
      !Construct multi-dim-array Tk
      !Upper diagonal exist like (12,23,34,45......)
      
       do lr = 1,L
            do lc = 1,L
            
                 if(lc .eq. lr+1) then
                    do p = 1,N
                    do q = 1,N
                    Ti_k(lr,lc,p,q) = UU(1,lr,K,p,q)!K is time slice
                    end do
                    end do
                end if
                  
              end do
      end do 
           
              
       !use boundary conditions
       
        do p = 1,N
        do q = 1,N
        Ti_k(L,1,p,q) = UU(1,L,K,p,q)!left - lower corner(lr= L,lc = 1)
        end do
        end do
        
        
      ! changing array to matrix     
      call array_to_matrix(Ti_k,MTi_k)!self
      !getting hermitian conjugate
      call getH(MTi_k,MTi_kD,ln)!Algebra

      !MTi_k and MTi_kD have  size lnxln
      Deallocate(Ti_k,O)
      return
      end

!===================================================================================
 !TdD calculates  MTd and MTdD matrices for specific space x. Only links 
 !from time direction contributes.Due to gauge fixing All Td at other (x<L)
 !are Identity except Td at (x = L).
 !==================================================================================
     
 subroutine TdD(UU,MTd,MTdD)
     implicit none
     integer::                                    p,q,lr,lc
     complex*16,intent(in)::                      UU(2,L,L,N,N)
     complex*16,allocatable,dimension(:,:)::      O
     complex*16,allocatable,dimension(:,:,:,:)::  Td
     complex*16,intent(out)::                     MTd(ln,ln),MTdD(ln,ln)  
     
     
    
         Allocate( Td(L,L,N,N),O(N,N) )     
                      
                    
           call IIN(O,N,dcmplx(0.d0,0.d0))!Algebra
           
           !Initiate multi-dim-array Td, each (lr,lc) block includes a NxN matrix
           !There is single Td
           
            do lr = 1,L
              do lc = 1,L
                do p = 1,N
                    do q = 1,N
                       
                            Td(lr,lc,p,q) = O(p,q)
                        
                     end do
                  end do
              end do
          end do
          
            
                    
        !Construct array Td for (x = L)
        !Elements exis only in diagonal position (11,22,33...) in LxL structure
        
         do lr = 1,L
                   
                    do p = 1,N
                    do q = 1,N
                  
                    Td(lr,lr,p,q) = UU(2,lr,L,p,q)! lc = L because  Td is non identity at lc = L
                   
                    end do
                    end do
                 
            
        end do 
                    
      ! changing array to matrix     
      call array_to_matrix(Td,MTd)!self
      !getting hermitian conjugate
      call getH(MTd,MTdD,ln)!Algebra

      ! All MTD and MTdD have size lnxln
      
       Deallocate( Td,O )       
      return
      end


!=============================================================================

      subroutine array_to_matrix(array,matrix)
      implicit none
      integer::                  p,q,lr,lc,u,v
      complex*16,intent(in)::    array(L,L,N,N)
      complex*16,intent(out)::   matrix(ln,ln)
      
     
      
                  !initiate the Matrix
                  call IIN(Matrix,ln,dcmplx(0.d0,0.d0))!Algebra
                  
                     p = 0
   14                q = 0 
                 
   15                      lr = int(p/float(N)) ! lr and lc finds which  (lr,lc)-block to go
                           lc = int(q/float(N))
                           u = mod(p,N)        ! u and v find which element of matrix inside (r,c) - block
                           v = mod(q,N)
                           
                           !print*,p+1,q+1,lr+1,lc+1,u+1,v+1,array(lr+1,lc+1,u+1,v+1)
                           
                    matrix(p+1,q+1) = array(lr+1,lc+1,u+1,v+1)
                    
                          !print*, p+1,q+1,Matrix(p+1,q+1)
                    
                    
                    q = q+1
                    
                    if (q .lt. (ln)) then
                         go to 15
                    else if (q .eq. (ln)) then
                         if (p .lt. (ln-1)) then 
                                p = p+1
                                go to 14
                         end if
                    end if
                    
                    
                           
                  return
                  end 
 

!=======================================================



!==================================================================
         
         end module Gamma
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         