
        subroutine ina(hode,v,frt,xh,xj,xm,xnai,xnao,xina) ! sodium current
       	implicit double precision (a-h,o-z)

       	gna=6.0d0
       	XKMCAM=0.3d0
       	deltax=-0.18d0
       
       	ena = (1.0d0/frt)*dlog(xnao/xnai) ! sodium reversal potential

       	am = 0.32d0*(v+47.13d0)/(1.0d0-dexp(-0.1d0*(v+47.13d0)))
       	bm = 0.08d0*dexp(-v/11.0d0)

c	camfact= 1.0d0/(1.0d0+(XKMCAM/caM)**4)
c	vshift = -3.25*camfact

	camfact=0.	
	vshift=0.
	
	vx=v-vshift		
	 	
        if(vx.lt.(-40.0d0)) then
         ah=0.135*dexp((80.0+vx)/(-6.8d0))
         bh=3.56*dexp(0.079*vx)+310000.0d0*dexp(0.35d0*vx)

	 aj1a=-127140.0*dexp(0.2444*vx)
	 aj1b=0.00003474d0*dexp(-0.04391d0*vx)
	 aj1c=(vx+37.78)/(1.0d0+dexp(0.311*(vx+79.23)))

	aj=(1.0d0+camfact*deltax)*(aj1a-aj1b)*aj1c
  	bj=(0.1212*dexp(-0.01052*vx))/(1.0+dexp(-0.1378d0*(vx+40.14d0)))

        else

         ah=0.0d0
         bh=1.0d0/(0.13d0*(1.0d0+dexp((vx+10.66)/(-11.1d0))))
         aj=0.0d0

	 bja=0.3*dexp(-0.0000002535d0*vx)
	 bjb=1.0+dexp(-0.1d0*(vx+32.0d0))

	bj=bja/bjb
        endif          
	
	tauh=1.0d0/(ah+bh)  
      	tauj=1.0d0/(aj+bj)
      	taum=1.0d0/(am+bm)

	xina=gna*xh*xj*xm*xm*xm*(v-ena)

        xh = ah/(ah+bh)-((ah/(ah+bh))-xh)*dexp(-hode/tauh)
        xj = aj/(aj+bj)-((aj/(aj+bj))-xj)*dexp(-hode/tauj)
        xm = am/(am+bm)-((am/(am+bm))-xm)*dexp(-hode/taum)

	return
	end 


c *****************************************************************************

	subroutine ikr(hode,v,frt,xko,xki,xr,xikr) ! IKR current
	implicit double precision (a-h,o-z)


      	ek = (1.0d0/frt)*dlog(xko/xki) ! K reversal potential
    
      	gss=dsqrt(xko/5.40)
      	xkrv1=0.00138d0*(v+7.0d0)/( 1.0-dexp(-0.123*(v+7.0d0))  )
      	xkrv2=0.00061d0*(v+10.0d0)/(dexp( 0.145d0*(v+10.0d0))-1.0d0)
      	taukr=1.0d0/(xkrv1+xkrv2)

      	xkrinf=1.0d0/(1.0d0+dexp(-(v+50.0d0)/7.5d0))

      	rg=1.0d0/(1.0d0+dexp((v+33.0d0)/22.4d0))
              
      	gkr=0.007836d0 ! Ikr conductance 
      	xikr=gkr*gss*xr*rg*(v-ek)

      	xr=xkrinf-(xkrinf-xr)*dexp(-hode/taukr)       
    

	return
	end 

c *****************************************************************************

	subroutine iks(hode,v,frt,ci,xnao,xnai,xko,xki,xs1,qks,xiks)  ! IKS current
	implicit double precision (a-h,o-z)

       prnak=0.018330d0
      
       qks_inf=0.6d0*ci
       tauqks=1000.0d0
      
      eks=(1.0d0/frt)*dlog((xko+prnak*xnao)/(xki+prnak*xnai))
      xs1ss=1.0/(1.0+dexp(-(v-1.50d0)/16.70d0))
      xs2ss=xs1ss

      tauxs=1.0d0/(0.0000719*(v+30.0d0)/(1.0d0-dexp(
     +-0.148d0*(v+30.0)))+0.000131d0
     + *(v+30.0d0)/(dexp(0.0687d0*(v+30.0d0))-1.0d0))
      
       gksx=0.200d0*0.7! Iks conductance 

       xiks=gksx*qks*xs1**2*(v-eks) 

       xs1=xs1ss-(xs1ss-xs1)*dexp(-hode/tauxs)
       qks=qks+hode*(qks_inf-qks)/tauqks

	return
	end 

c *****************************************************************************

	subroutine ik1(hode,v,frt,xki,xko,xik1)  ! IK1 current
	implicit double precision (a-h,o-z)

	ek = (1.0d0/frt)*dlog(xko/xki)     

      gkix=0.60d0*1.0  ! Ik1 conductance reduced in Grandi model
      gki=gkix*(dsqrt(xko/5.4))
      aki=1.02/(1.0+dexp(0.2385*(v-ek-59.215)))
      bki=(0.49124*dexp(0.08032*(v-ek+5.476))+dexp(0.061750
     +*(v-ek-594.31)))/(1.0+dexp(-0.5143*(v-ek+4.753)))
      xkin=aki/(aki+bki)
      xik1=gki*xkin*(v-ek)

	return
	end 

c *****************************************************************************

	! Ito current 
       	subroutine ito(hode,v,frt,xki,xko,xtof,ytof,xtos,ytos,
     &  xito,gtof,gtos)
     
	implicit double precision (a-h,o-z)

	ek = (1.0d0/frt)*dlog(xko/xki)

	rt1=-(v+3.0)/15.0d0
	rt2=(v+33.5)/10.0d0 
	rt3=(v+60.0d0)/10.0d0 
	xtos_inf=1.0d0/(1.0+dexp(rt1))
        ytos_inf=1.0d0/(1.0d0+dexp(rt2))

	rs_inf=1.0d0/(1.0d0+dexp(rt2))

	txs=9.0d0/(1.0d0+dexp(-rt1)) + 0.5d0
	tys=3000.0d0/(1.0+dexp(rt3)) + 30.0d0 
	
	xitos=gtos*xtos*(ytos+0.5d0*rs_inf)*(v-ek) ! ito slow

	xtos = xtos_inf-(xtos_inf-xtos)*dexp(-hode/txs)
        ytos = ytos_inf-(ytos_inf-ytos)*dexp(-hode/tys)

c --------- Ito fast following Shannon et. al. 2005 -----------

	xtof_inf=xtos_inf
	ytof_inf=ytos_inf 

	rt4=-(v/30.0d0)*(v/30.0d0)
	rt5=(v+33.5d0)/10.0d0
	txf=3.5d0*dexp(rt4)+1.5d0
	tyf=20.0/(1.0+dexp(rt5))+20.0d0	
	
	xitof=gtof*xtof*ytof*(v-ek)
	
        xtof = xtof_inf-(xtof_inf-xtof)*dexp(-hode/txf)
        ytof = ytof_inf-(ytof_inf-ytof)*dexp(-hode/tyf)

	xito=xitos+xitof ! total ito current 


	return
	end 	
	

c *****************************************************************************

	! Inak current
	subroutine inak(v,frt,xko,xnao,xnai,xinak)
	implicit double precision (a-h,o-z)

c -------  Inak (sodium-potassium exchanger) following Shannon --------------

	 xibarnak=1.50d0

         xkmko=1.5d0 ! these are Inak constants adjusted to fit
c                  ! the experimentally measured dynamic restitution curve
         xkmnai=12.0d0
         hh=1.0d0  ! Na dependence exponent

         sigma = (dexp(xnao/67.3d0)-1.0d0)/7.0d0  
      fnak = 1.0d0/(1+0.1245*dexp(-0.1*v*frt)
     +           +0.0365*sigma*dexp(-v*frt)) 
      xinak = xibarnak*fnak*(1./(1.+(xkmnai/xnai)**hh))*xko/(xko+xkmko)

	return
	end 


	


