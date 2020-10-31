

!====================================================================================
!This is a Master Program which controls all modules
!====================================================================================
!Date: May22,2015        by: Dibakar sigdel        sub: Lattic_gauge _theory/fermion   
!====================================================================================




!========================================= Data - Set  ===================================
Module Data_set

       implicit none
       Save
       Integer,parameter ::    L = 4,   N = 3 ,ns = 2, LF_steps = 200, hmc_steps = 100, mass_N = 40
       real*8,parameter::      dt_lf = 1.0d0/dfloat(LF_steps),absT = 1.d-10, Area = 5/dfloat(2)
       real*8,parameter::      Nf = 2.d0 , bb = 1.0d0*(dfloat(L*L)/(4.d0*Area))! , mass = 0.0d0/dfloat(L)
       integer,parameter::     ln = L*N, lln = L*L*N, lln2 = L*L*N*ns, ln2 =L*N*ns,Nns = N*ns,ll = L*L
       
       
       
       character*1,parameter::    model ='X'
       
       
       
       integer,parameter::        fem = 1
      
End module Data_set
!=========================================================================================
  

  
  
  PROGRAM GOOD_START
  
  ! -------------------------  use modules  ---------------------------------------------
  use Data_set
  use Mstart
  use Algebra 
  use hdm
  use Updator
  use Calculator
  use Functionary
  use LFrunner
  use Runner
  use Printer
  use Gamma
  use Partition
  !use Test
 
  !-------------------------------------------------------------------------------------

  IMPLICIT NONE
  
  !-----------------------------Permanent declaration-----------------------------------
  integer::                 		    i,j,k,r,p,q,s,t,seed,lseed
  complex*16,dimension(2,L,L,N,N)::     UU,PI
  complex*16,dimension(lln2,lln2) ::    HD,DD,IHD,III
  complex*16::                          Tau(ln2,ln2),Z20,stor_Z20(hmc_steps),stor_Z31(hmc_steps),Z31ALL
  
  real*8,dimension(hmc_steps)::         m20,m31
  real*8::                              sm20,sm31,var_m20,var_m31,std_dev_m20,std_dev_m31
  real*8::                              mass,mm20,mm31,SUG
  
 
   mass = 0.1d0*(-(2.d0)/float(2) + 2.d0)
  
  !-----------------------------call rand()-----------------------------------
       
        seed = 7686545
        call srand(seed)
        !t = rand()
  !-------------------------------------------------------------------
  
  
!============================= Test =======================================     
  
    !call fix_Ustart(UU)   
    !call Tau_finder(UU,Tau)
    
    !call C_ICD_checker(Tau)
    !call  PZ20(UU,Z20)
    !print*, Z20
   
    !call HMC()
    
    !call Measure()
    
 !mass = 0.5

 !call LFchecker(mass)
 
 !call PZ31A(UU,mass,Z31ALL)
 !call PZ31B(UU,mass,Z31ALL)
 
 
  !call  Plaquette(SUG,UU,mass)
  !print*,SUG

!============================= MAIN =======================================    

  !call HMC(mass)

 
  call  mass_finder(mass_N)
  
  
!====================================================================    

 END PROGRAM GOOD_START
 	
 	
 	
 	
!  plot "fort.31" using 1:2 with lines'\
!        "fort.20" using 1:2 with lines

      
      
      
      
      
      
      
      
      
      
      
      
      
      
       









       
