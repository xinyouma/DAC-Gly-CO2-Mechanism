      implicit none
      integer:: nline
      integer:: i,j
      integer:: max_loc
      real*8::  prefac,h,kBT,sumI,num
      real*8::  Zn,temp,rate,tau
      real*8::  beta
      real*8::  maximum,cndag,Wdag
      real*8::  cn_min_i,cn_min_f,kappa
      real*8,allocatable,dimension(:):: cn
      real*8,allocatable,dimension(:):: W
      
      read(*,*) nline ! number of lines
      read(*,*) Zn,temp,cn_min_i,cn_min_f ! Unit of Zn=1/(nm^2 * g/mol) 
      read(*,*) kappa ! transmission coefficient
      beta=1.0/(0.0019872041*temp) !1/kcal/mol
      prefac=sqrt(Zn/(2.0*(22.0/7.0)*beta*(1E-03/4184.0))) 
      
      allocate(cn(nline))
      allocate(W(nline))
      write(*,*) "Reading PMF.dat: cn dimensionless, W in kcal/mol"
         write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
      open(unit=20,FILE="PMF.dat",status="old") 
      do i=1,nline
           read(20,*) cn(i),W(i)
      end do
      close(20)

      maximum=0.0
      do i=1,nline
         if((cn(i)>cn_min_i).and.(cn(i)<cn_min_f)) then
           if(W(i)>maximum) then
              maximum=W(i)
              max_loc=i
           endif
         endif  
      end do
      Wdag=W(max_loc)
      cndag=cn(max_loc)
      write(*,*) " Barrier | location and value  "
      write(*,"(2F12.6,1X)") cn(max_loc),W(max_loc)
       
      !!!!!!! calculate TST Forward rate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      num=exp(-beta*Wdag)
      sumI=0.5*(exp(-beta*W(1))+exp(-beta*W(max_loc)))
      do i=2,max_loc-1
         sumI=sumI+exp(-beta*W(i))
      end do
      sumI=sumI*(cn(2)-cn(1))  
      rate=prefac*(num/sumI)*(1.0/1E-09) ! per s
      rate=rate*1E-12 ! per ps
      tau=1.0/rate ! ps

      write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
      write(*,*) "Forward TST rate and time (1/ps and ps)"
      write(*,"(2(F12.6,1X))") rate,tau 
      write(*,*) "Forward: kappa-corrected rate and time (1/ps and ps)"
      write(*,"(2(F12.6,1X))") kappa*rate,tau/kappa
      write(*,*) "Forward: kappa-corrected rate (per s)"
      write(*,*) kappa*rate*(1E12) 
      write(*,*) "Forward: Mass weighted Reactant volume"
      write(*,"(2F12.6,1X)") sumI/sqrt(Zn)
      write(*,*) "Forward: Reactant volume"
      write(*,"(2F12.6,1X)") sumI
      write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"   

      !!!!!!! calculate TST Backward rate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sumI=0.5*(exp(-beta*W(nline))+exp(-beta*W(max_loc)))
      do i=max_loc+1,nline-1
         sumI=sumI+exp(-beta*W(i))
      end do
      sumI=sumI*(cn(2)-cn(1))  
      rate=prefac*(num/sumI)*(1.0/1E-09) ! per s
      rate=rate*1E-12 ! per ps
      tau=1.0/rate ! ps
      write(*,*) "Backward TST rate and time (1/ps and ps)"
      write(*,"(2(F12.6,1X))") rate,tau 
      write(*,*) "Backward: kappa-corrected rate and time (1/ps and ps)"
      write(*,"(2(F12.6,1X))") kappa*rate,tau/kappa
      write(*,*) "Backward: kappa-corrected rate (per s)"
      write(*,*) kappa*rate*(1E12)
      write(*,*) "Backward: Mass-weighted reactant vol(nm(g/mol)^(1/2))"
      write(*,"(2F12.6,1X)") sumI/sqrt(Zn)
      write(*,*) "Backward: Reactant volume"
      write(*,"(2F12.6,1X)") sumI

      deallocate(W,cn)

      end
