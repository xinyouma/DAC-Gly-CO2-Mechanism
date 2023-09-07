        implicit none
        integer:: i,min_loc,nline,max_loc
         real*8:: beta,num,Wdag,rdag
         real*8:: mu,prefac,minimum,maximum
         real*8:: rate,tau,sumI,temp,diff
         real*8:: r_min_i,r_min_f,kappa,sumF
         real*8,allocatable,dimension(:)::r,W

         read(*,*) nline ! number of lines
         read(*,*) mu,temp,r_min_i,r_min_f !reduced mass,temp,min1,min2
         read(*,*) kappa                   ! transmission coefficient
         beta=1.0/(0.0019872041*temp) !1/kcal/mol
         prefac=1.0/sqrt(2.0*(22.0/7.0)*beta*mu*(1E-03/4184.0))

         allocate(r(nline))   ! in nm
         allocate(W(nline))   ! in kcal/mol
         write(*,*) "Reading PMF.dat: r in nm, W in kcal/mol"
         write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" 
         open(unit=20,FILE="PMF.dat",status="old") 
         do i=1,nline
            read(20,*) r(i),W(i)
         end do
         close(20)

         maximum=0.0
         do i=1,nline
           if((r(i)>r_min_i).and.(r(i)<r_min_f)) then
             if(W(i)>maximum) then
                maximum=W(i)
                max_loc=i
             endif
           endif
         end do
         Wdag=W(max_loc)
         rdag=r(max_loc)
         write(*,*) " Barrier | location and value  "
         write(*,"(2F12.6,1X)") r(max_loc),W(max_loc) 

         num=(rdag**2)*exp(-beta*Wdag)
         ! FW reaction
         sumI=0.5*(((r(1)**2)*exp(-beta*W(1)))+((r(max_loc)**2)*exp(-beta*W(max_loc))))
         do i=2,max_loc-1
            sumI=sumI+((r(i)**2)*exp(-beta*W(i)))
         end do
         sumI=sumI*(r(2)-r(1))
         rate=prefac*(num/sumI)*(1.0/1E-09) ! per s
         
         write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
         rate=rate*1E-12 !per ps
         tau=1.0/rate !ps
         write(*,*) "Forward: TST rate and time (ps)"
         write(*,"(2F12.6,1X)") rate,tau 
         write(*,*) "Forward: kappa corrected TST rate and time (ps)"
         write(*,"(2F12.6,1X)") kappa*rate,tau/kappa
         write(*,*) "Forward: kappa-corrected TST rate (per s)"
         write(*,*) kappa*rate*(1E12)
         write(*,*) "Reduced mass, reactant V, mass weighted reactant V"
         write(*,*) "(g/mol), (nm^3), sqrt(g/mol)*(nm^3)"
         write(*,"(3F12.6,1X)") mu,sumI,sqrt(mu)*sumI
        
         write(*,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
         ! BW reaction
         sumI=0.5*(((r(nline)**2)*exp(-beta*W(nline)))+((r(max_loc)**2)*exp(-beta*W(max_loc))))
         do i=max_loc+1,nline-1
            sumI=sumI+((r(i)**2)*exp(-beta*W(i)))
         end do
         sumI=sumI*(r(2)-r(1))
         rate=prefac*(num/sumI)*(1.0/1E-09) ! per s
         rate=rate*1E-9 !per ns
         tau=1.0/rate !ns
         write(*,*) "Backward: TST rate and time (ns)"
         write(*,"(2F12.6,1X)") rate,tau
         write(*,*) "Backward: kappa corrected TST rate and time (ns)"
         write(*,"(2F12.6,1X)") kappa*rate,tau/kappa
         write(*,*) "Forward: kappa-corrected TST rate (per s)"
         write(*,*) kappa*rate*(1E9)
         write(*,*) "Reduced mass, reactant V, mass weighted reactant V"
         write(*,*) "(g/mol), (nm^3), sqrt(g/mol)*(nm^3)"
         write(*,"(3F12.6,1X)") mu,sumI,sqrt(mu)*sumI
         deallocate(r)
         deallocate(W)
         end
