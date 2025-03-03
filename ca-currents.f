	
      
	subroutine binevol(nt,nx,alpha,beta,dt,iseed,ndeltap,ndeltam)                  
	implicit double precision (a-h,o-z)

	! note code unstable when alpha drops too low.  
	! so just don't call a number if the forward rate is very small. 

	if(alpha.lt.0.00001) then
	ndeltap=0
	goto 200
	endif 

	na=nt-nx  ! this is the number of available clusters
	
	if(na.gt.0) then	
	xrate=alpha*dt ! this is the probability

	call binom(na,xrate,ndeltap,iseed)

	else

	ndeltap=0

	endif 
	
200	if(nx.gt.0) then
	xrate=beta*dt

	call binom(nx,xrate,ndeltam,iseed)

	else
	ndeltam=0
	endif 
	
	return
	end

c **************************************************************

	subroutine uptake(ci,vup,xup)
	implicit double precision (a-h,o-z)

	double precision, parameter::Ki=0.30d0
	double precision, parameter::Knsr=800.0d0
	double precision, parameter::HH=3.00d0

	xup=vup*ci**HH/(Ki**HH+ci**HH)
	
	return
	end

c ***************************************************************

	subroutine total(ci,cit)
	implicit double precision (a-h,o-z)
       		
      	bcal=24.0d0
      	xkcal=7.0d0
      	srmax=47.0d0
      	srkd=0.6d0
       
      	bix=bcal*ci/(xkcal+ci)   !; calmodulin 
      	six=srmax*ci/(srkd+ci)    !; SR buffering

      	cit=ci+bix+six

	return
 	end 

c ***************************************************************

	subroutine xfree(cit,ci)
	implicit double precision (a-h,o-z)
       		
	a=2.23895
	b=52.0344
	c=0.666509

	y=cit

	xa=(b+a*c-y)**2+4.0*a*c*y
	ci=(-b-a*c+y+dsqrt(xa))/(2.0*a)

	return
 	end 

c **************************************************************

 	subroutine inaca(v,frt,xnai,xnao,cao,ci,xinacaq)
        implicit double precision (a-h,o-z)
       	
	cim=ci/1000.0d0 ! convert to milli molars
	
	zw3a=xnai**3*cao*dexp(v*0.35*frt)
	zw3b=xnao**3*cim*dexp(v*(0.35-1.)*frt)

	zw3=zw3a-zw3b
      	zw4=1.0d0+0.2d0*dexp(v*(0.35d0-1.0d0)*frt)

      	xkdna=0.3d0 ! micro M
      	aloss=1.0d0/(1.0d0+(xkdna/ci)**3)

      	xmcao=1.3d0
      	xmnao=87.5d0
      	xmnai=12.3d0
      	xmcai=0.0036d0

c       Note: all concentrations are in mM

      	yz1=xmcao*xnai**3+xmnao**3*cim
      	yz2=xmnai**3*cao*(1.0d0+cim/xmcai)
      	yz3=xmcai*xnao**3*(1.0d0+(xnai/xmnai)**3)
      	yz4=xnai**3*cao+xnao**3*cim

      	zw8=yz1+yz2+yz3+yz4

      	xinacaq=aloss*zw3/(zw4*zw8)

	return

	end 


c **************************************************************

	subroutine ica(v,frt,cao,ci,pox,rca, xicaq)
	implicit double precision (a-h,o-z)
      
        xf=96.485d0        ! Faraday's constant
      	pca=0.00054d0 ! constant from Luo-Rudy
      	za=v*2.0d0*frt
	
	factor1=4.0d0*pca*xf*frt
      	factor=v*factor1

      	cim=ci/1000.0d0 !; convert micro M to mM
             
c       compute driving force here 
      
        if(dabs(za).lt.0.001d0) then
        rca=factor1*(cim*dexp(za)-0.341d0*(cao))/(2.0d0*frt)
        else   
        rca=factor*(cim*dexp(za)-0.341d0*(cao))/(dexp(za)-1.0d0)         
        endif 

	xicaq=rca*pox
	
	return
	end 

c *****************************************************************
                   
	   subroutine markov(hode,v,ci,c1,c2,xi1,xi2,po,
     + c1s,c2s,xi1s,xi2s,pos,alpha,bts,zxr)
          implicit double precision (a-h,o-z)
      
	! Ca independent rates here

	  a23=0.3
	  a32=3.0
          a42=0.00224d0
  	 
	  vth=1.0  ! these parameters determine reactivation of ica
	  s6=7.0     
     
          poinf=1.0d0/(1.0d0+dexp(-(v-vth)/s6))
	  taupo=1.0d0   

          a12=poinf/taupo
          a21=(1.0d0-poinf)/taupo

	  vy=-40.0d0
          sy=4.0d0
          prv=1.0d0-1.0d0/(1.0d0+dexp(-(v-vy)/sy))

	  vyr=-40.0d0
	  syr=10.0d0 

          recovx=10.0d0+4954.0d0*dexp(v/15.6d0)
	  recov=1.5*recovx

          tauba=(recov-450.0d0)*prv+450.0d0
	  tauba=tauba/2.0

          poix=1.0d0/(1.0d0+dexp(-(v-vyr)/syr))

          a15=poix/tauba
          a51=(1.0d0-poix)/tauba

	  vx=-40.0d0
          sx=3.0d0
          tau3=3.0d0
          poi=1.0d0/(1.0d0+dexp(-(v-vx)/sx))
          a45=(1.0d0-poi)/tau3
	
          ! Ca dependent rates here
	
	  cat=0.8d0

          fca=1.0d0/(1.0d0+(cat/ci)**2)
          a24=0.00413d0+zxr*fca
          a34=0.00195+zxr*fca

	  a43=a34*(a23/a32)*(a42/a24)
          a54=a45*(a51/a15)*(a24/a42)*(a12/a21)

          fcax=1.0d0
          a24s=0.00413d0+zxr*fcax
          a34s=0.00195+zxr*fcax

	  a43s=a34s*(a23/a32)*(a42/a24s)
          a54s=a45*(a51/a15)*(a24s/a42)*(a12/a21)

	  ! state dynamics	

	  dpo= a23*c1+a43*xi1-(a34+a32)*po-alpha*po+bts*pos
          dc2= a21*c1+a51*xi2-(a15+a12)*c2+bts*c2s
          dc1= a12*c2+a42*xi1+a32*po-(a21+a23+a24)*c1+bts*c1s

          dxi1=a24*c1+a54*xi2+a34*po-(a45+a42+a43)*xi1+bts*xi1s

	  dpos= a23*c1s+a43s*xi1s-(a34s+a32)*pos+alpha*po-bts*pos
          dc2s= a21*c1s+a51*xi2s-(a15+a12)*c2s-bts*c2s
          dc1s= a12*c2s+a42*xi1s+a32*pos-(a21+a23+a24s)*c1s-bts*c1s

 	  dxi1s=a24s*c1s+ a54s*xi2s+a34s*pos-(a45+a42+a43s)*xi1s-bts*xi1s
          dxi2s=a45*xi1s+ a15*c2s -(a51+a54s)*xi2s-bts*xi2s


	  po=po+dpo*hode
          c1=c1+dc1*hode
          c2=c2+dc2*hode
          xi1=xi1+dxi1*hode

	  pos=pos+dpos*hode
          c1s=c1s+dc1s*hode
          c2s=c2s+dc2s*hode
          xi1s=xi1s+dxi1s*hode
          xi2s=xi2s+dxi2s*hode	

	 xi2=1.0d0-c1-c2-po-xi1-pos-c1s-c2s-xi1s-xi2s	
         
	return
	end



