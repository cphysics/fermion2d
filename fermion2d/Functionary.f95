 
!====================================================================================
!This module contains several functions
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================
 
 
 
 MODULE Functionary


!----use--------
use Data_set
!---------------- 

 
IMPLICIT NONE
CONTAINS       
      
        
        
 
!=============================CONTAINS===================================

!======================================================================================
real*8 function min(a,b)
        implicit none
		real*8::    a,b  
	  		if (a .lt. b) then
	    		min = a
	  		else if( b .lt. a) then
	    		min = b
	  		end if
			return
			end 
						
!==================================================================

integer function fun(s)
        implicit none
		integer:: s
		     
	  		if (s .eq. 1) then
	    		fun = L ! why N ?, big mistake
	  		else if( s .ne. 1) then
	    		fun = s-1
	  		end if
			return
			end 
!==================================================================
			
integer function modi(s)
        implicit none
		integer:: s,md
		
		md = mod(s+1,L)
	  		if (md .eq. 0) then
	    		modi = L
	  		else if( s .ne. 0) then
	    		modi = md
	  		end if
			return
			end 			




!==================================================================



!==================================================================
         
         
         end module Functionary
         
         
         
         
         
         