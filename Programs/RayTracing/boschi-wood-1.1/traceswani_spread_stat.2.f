c--JWKB ray tracing for a surface wave is a shooting problem where the only
c--dependent variable whose value we can change is the wavevector angle zeta
c--at the source. (See Larson et al., Geophys. J. Int. 132 (3), 654-666, 1998, 
c--and Numerical Recipes in Fortran, 1992 edition, sections 17.0, 17.1.)
c------------------------------------------------------------------------------
c---this version integrates equations (3.21) and (3.22) of Larson et al. to 
c---calculate the geometrical spreading.
c------------------------------------------------------------------------------
c-----this version allows to calculate orbits other than 1.
c------------------------------------------------------------------------------
c------this version includes calculation of phase anomaly
	dimension urot(3,3),urott(3,3)
	PARAMETER(nmax=5,KMAXX=10000)
	real y(nmax)
	REAL xp(KMAXX),yp(NMAX,KMAXX)
	REAL xpsav(KMAXX),ypsav(KMAXX)

c--working storage
	parameter(LMX=55)
	double precision d(2*LMX+1,2*LMX+1)
	complex*8 vec1(2*(2*LMX+1)),vec2(2*(2*LMX+1))

c----the flwng line will have to change for waves other than R100
c	parameter(c0=4.08957,prd=100.)    !vel. of 100s R waves in PREM
	parameter(thremin=1.e-2) !threshold for convergence

c--storage of generalized spherical harmonic coefficients
	dimension dlc0((LMX+1)**2)
     &     ,dlc2(2*(LMX+1)**2-2*2**2)
     &     ,dlc4(2*(LMX+1)**2-2*4**2)
	dimension dlc0r((LMX+1)**2)
     &     ,dlc2r(2*(LMX+1)**2-2*2**2)
     &     ,dlc4r(2*(LMX+1)**2-2*4**2)

	dimension zeta0min(50),yp30min(50),phas_min(50)

	COMMON /path/ kmax,kount,dxsav,xp,yp
	common /expansion/dlc0r,dlc2r,dlc4r,lmax0,lmax2,lmax4
	common /rearth/rearth,omega,pi
	character chiray*4,finishlat*80,sostat*80,chnmin*1
c       character chisave*5
	character filename0*80,filename2*80,filename4*80
	external derivs

	rearth=6371.
	radian=90./asin(1.)
        pi=2.*asin(1.)
	iray=0

	print*,"reference phase velocity?"
	read*,c0
	print*,"period?"
	read*,prd

	omega=2.*pi/prd !angular frequency, needed for phase anomaly

c--input reference model
	print*,"file with isotropic harmonic coefficients?"
	read*,filename0
	print*,"file with 2psi harmonic coefficients?"
	read*,filename2
	print*,"file with 4psi harmonic coefficients?"
	read*,filename4
c-----------------------choose orbit
	print*,"what orbit?"
	read*,iorbit
c--output file
	print*,"file to store spreading factors in?"
	read*,finishlat
c--file with source station locations
	print*,"file with source station locations?"
	read*,sostat

c----------------------------------------------------read in model coefficients
	open(1,file=filename0,status="old")
	read(1,*) lmax0,(dlc0(i),i=1,(lmax0+1)**2)
	close(1)
	open(1,file=filename2,status="old")
	read(1,*) lmax2,(dlc2(i),i=1,2*(lmax2+1)**2-2*2**2)
	close(1)
	open(1,file=filename4,status="old")
	read(1,*) lmax4,(dlc4(i),i=1,2*(lmax4+1)**2-2*4**2)
	close(1)
	if(lmax0.gt.lmx) stop 'lmax of 0psi model too large'
	if(lmax2.gt.lmx) stop 'lmax of 2psi model too large'
	if(lmax4.gt.lmx) stop 'lmax of 4psi model too large'
	write(6,'(''lmax0='',i3,''  lmax2='',i3,''  lmax4='',i3)')lmax0,lmax2,lmax4
	len0=(lmax0+1)**2
	len2=2*(lmax2+1)**2-2*2**2
	len4=2*(lmax4+1)**2-2*4**2

c----from relative perturbation get absolute velocity (in what units?).
c----units should not matter as in ray eqs. only derivatives of ln of velocity
c----appear.
	dlc0(1)=dlc0(1)*c0+c0*sqrt(8.*asin(1.))
	do i=2,len0
	   dlc0(i)=dlc0(i)*c0
	enddo
	do i=1,len2
	   dlc2(i)=dlc2(i)*c0
	enddo
	do i=1,len4
	   dlc4(i)=dlc4(i)*c0
	enddo
c---------------all this assumes reference model is isotropic

	isave=0
c---------loop over a set of latitudes and longitudes, to do statistics
c	open(121,file="random_locs.txt",status="old")
	open(121,file=sostat,status="old")
	do k=1,80
	   if(finishlat(k:k).eq." ")goto2
	enddo
 2	open(21,file=finishlat(1:k-1)//".nomulti")
	open(41,file=finishlat(1:k-1)//".multi")
	open(31,file=finishlat(1:k-1)//".nomulti_sr_st")
 1	read(121,*,end=122)eplo,epla,stlo,stla

c---find Euler angles and rotation matrix
        print*,eplo,epla,stlo,stla
	call pole(epla,eplo,stla,stlo,xlatp,xlonp,azmp,delta)
	alph=xlonp/radian
	beta=(90.-xlatp)/radian
	gama=(180.-azmp-.5*delta)/radian
	call setrot(alph,beta,gama,urot)
c---the inverse of a rotation matrix equals its transpose
	do i=1,3
	do j=1,3
	urott(i,j)=urot(j,i)
	enddo
	enddo

	print*,'epicentral distance=',delta
	del=delta/radian
c--define finishing longitude based on orbit
	if(mod(iorbit,2).ne.0)then
	   finish=(pi*(float(iorbit-1)))+del
	   h1=1./10000.
	   zeta0mid=0.
	else
	   finish=(pi*(float(iorbit)))-del
	   finish=-finish
	   h1=-1./10000.
	   zeta0mid=pi
	endif
	phase0=(omega/c0)*(abs(finish))*rearth

c---find expansion of model in generalized spherical harmonics defined
c---in the rotated reference frame
	call rotvcani(
     &     dlc0,dlc2,dlc4,lmax0,lmax2,lmax4
     &     ,alph,beta,gama,d,vec1,vec2
     &     ,dlc0r,dlc2r,dlc4r)

	print*,"start loop over initial azimuth"

	zeta00=100.
	f0=100.
	fmin=100.
	nmin=0
	iray=iray+1
	icounter=0
	finlat0=100.
c-------------------------------loop over initial value of azimuth
	do zeta0=zeta0mid-0.1,zeta0mid+0.1,0.002
	   icounter=icounter+1

c--initial conditions of ray-tracing eqs.
	y(1)=zeta0
	y(2)=asin(1.)
c----initial conditions of eqs. (3.21) and (3.22) of Larson, Tromp and Ekstrom 98.
	y(3)=0.
	y(4)=1.
	y(5)=0. ! initial value for phase: should not matter
c	iray=iray+1
	write(chiray,'(i4.4)')iray
	dxsav=0.000001 !min interval between saved samples of ray (radians)
	kmax=0
c	h1=del/10000.
c	h1=1./10000.
c	hmin=0.
c	call odeint(y,4,0.,del,1.e-6,h1,0.,nok,nbad,nstprk,derivs)
	call odeint(y,5,0.,finish,1.e-6,h1,0.,nok,nbad,nstprk,derivs)
	call score(y,f,nmax)
	finlat=yp(2,kount)-asin(1.)
c	if(f.gt.fmin.and.fmin.eq.f0.and.fmin.lt.thremin)then
	if((finlat*finlat0.lt.0..or.finlat0.eq.0.).and.(icounter.gt.1))then
	   if(f.lt.f0)then
	      kountsav=kount
	      do i=1,kountsav
		 ypsav(i)=yp(2,i)
		 xpsav(i)=xp(i)
	      enddo
	      zeta00=zeta0
	      yp30=yp(3,kount)	! spreading factor
	      phas0=yp(5,kount)-phase0  ! phase anomaly
	      print*,"found a minimum at take-off angle=",zeta00*90./asin(1.)
	   else
	      print*,"found a minimum at take-off angle=",zeta00*90./asin(1.)
	   endif

c	   fmin=thremin ! rough fix lapo 9/11/05
	   nmin=nmin+1
c	   write(21,"(6(f11.6,1x))")eplo,epla,stlo,stla,zeta00*90./asin(1.),yp30
	   zeta0min(nmin)=zeta00
	   yp30min(nmin)=yp30
	   phas_min(nmin)=phas0
c	   isave=isave+1
c	   write(chisave,"(i5.5)")isave
	   write(chnmin,"(i1.1)")nmin
	   if(nmin.eq.1)then
	   open(97,file="rays/gca.txt"//chiray)
	   write(97,*)eplo,epla
	   write(97,*)stlo,stla
	   close(97)
	   endif
	   write(*,*)"save raypath as file rays/ray.rotated"//chiray//"_"//chnmin
	   open(99,file="rays/ray.rotated"//chiray//"_"//chnmin)
	   do i=1,kountsav
	      yin=90.-ypsav(i)*radian
	      xin=xpsav(i)*radian
	      call rotll(yin,xin,dla,dlo,urott)
	      if(dlo.lt.0.)dlo=dlo+360.
	      write(99,*)dlo,dla,ypsav(i)
	   enddo
	   close(99)
c-------------------
	endif
	f0=f
	finlat0=finlat
	zeta00=zeta0
	yp30=yp(3,kount) ! spreading factor
	phas0=yp(5,kount)-phase0 ! phase anomaly
	
c-------------also save ray each time
	kountsav=kount
	do i=1,kountsav
	   ypsav(i)=yp(2,i)
	   xpsav(i)=xp(i)
	enddo
c-------------ray saved
	enddo
c---------------------------------end of loop over initial azimuth
	if(nmin.eq.0)print*,"did not find a minimum ",eplo,epla,stlo,stla
	if(nmin.eq.1)then
c---implement eq. 41 of Woodhouse and Wong 1986
	   ampan=abs(sin(del)/yp30min(1))
	   write(21,"(8(f11.6,1x))")eplo,epla,stlo,stla,delta,zeta0min(1)*90./asin(1.)
     &,ampan,phas_min(1)
	   write(31,*)eplo,epla,stlo,stla
	else
	   do kmin=1,nmin
c---implement eq. 41 of Woodhouse and Wong 1986. yp30min is spreading factor
              ampan=abs(sin(del)/yp30min(kmin))
	      write(41,"(8(f11.6,1x))")
     &  eplo,epla,stlo,stla,delta,zeta0min(kmin)*90./asin(1.),ampan,
     &  phas_min(kmin)
	   enddo
	endif

	goto1
 122	close(21)
	close(31)
	close(41)
	end
c-------------------------------------------------------------------

	subroutine score(y,f,nmax)
c------how well did we do?
	real y(nmax)
	f=y(2)-asin(1.)
c	f=f*f
	f=abs(f)
	return
	end

c-------------------------------------------------------------------
	subroutine derivs(x,y,dydx)
	real x,dydx(5),y(5)
	parameter(lmx=55)
	dimension dlc0r((LMX+1)**2)
     &     ,dlc2r(2*(LMX+1)**2-2*2**2)
     &     ,dlc4r(2*(LMX+1)**2-2*4**2)
	common /expansion/dlc0r,dlc2r,dlc4r,lmax0,lmax2,lmax4
c	common /phase/phas,omega,pi
	common /rearth/rearth,omega,pi

	dimension yth2((lmx+1)**2),ythph((lmx+1)**2)
	dimension wk1((lmx+1)*4),wk2((lmx+1)*4),wk3((lmx+1)*4)

	dimension y0((lmx+1)**2),yph((lmx+1)**2),yth((lmx+1)**2)
        dimension ypp((lmx-1)*(2*lmx+6))
     &  ,ytp((lmx-1)*(2*lmx+6))
     &  ,ypppp((lmx-3)*(2*lmx+10))
     &  ,ytppp((lmx-3)*(2*lmx+10))
     &  ,ypp_ph((lmx-1)*(2*lmx+6))
     &  ,ytp_ph((lmx-1)*(2*lmx+6))
     &  ,ypp_th((lmx-1)*(2*lmx+6))
     &  ,ytp_th((lmx-1)*(2*lmx+6))
     1         ,ypp_thph((lmx-1)*(2*lmx+6))
     1         ,ytp_thph((lmx-1)*(2*lmx+6))
     1         ,ypp_th2((lmx-1)*(2*lmx+6))
     1         ,ytp_th2((lmx-1)*(2*lmx+6))
     &  ,ypppp_ph((lmx-3)*(2*lmx+10))
     &  ,ytppp_ph((lmx-3)*(2*lmx+10))
     &  ,ypppp_th((lmx-3)*(2*lmx+10))
     &  ,ytppp_th((lmx-3)*(2*lmx+10))
     1         ,ypppp_thph((lmx-3)*(2*lmx+10))
     1         ,ytppp_thph((lmx-3)*(2*lmx+10))
     1         ,ypppp_th2((lmx-3)*(2*lmx+10))
     1         ,ytppp_th2((lmx-3)*(2*lmx+10))

        double precision wk(9*(2*lmx+1))
	dimension sigma(2,2,2,2),tau(2,2),xnu(2),xnu_az(2),xnu_az2(2)
	dimension tau_th(2,2),tau_ph(2,2),sigma_th(2,2,2,2),sigma_ph(2,2,2,2)
	dimension tau_th2(2,2),tau_thph(2,2),sigma_th2(2,2,2,2),sigma_thph(2,2,2,2)
cTEST
c	common/cTEST/dlncdaz,dlncdthp,dlncdphp,dlncdth2,dlncdthph,dlncdaz2,dlncdphaz,dlncdthaz

	radian=90./asin(1.)
	len0=(lmax0+1)**2
	len2=2*(lmax2+1)**2-2*2**2
	len4=2*(lmax4+1)**2-2*4**2

c---in this version x is incremental epicentral distance, need to find longitude.

c---y(1)=zeta (azimuth from East); y(2)=theta; x=phi
        xlat=90.-(y(2)*radian)
        xlon=x*radian
c--------------it is faster to find scalar harmonics separately if lmax0 >> lmax2 and lmax4
c---this is probably the most time-consuming part of the process
	call apkyd(lmax0,x,y(2),y0,yph,yth,yth2,ythph,wk1,wk2,wk3)
	call ylmv4_der(xlat,xlon,max(lmax2,lmax4),y0,ypp,ytp,yph,yth,
     &  ypp_th,ypp_ph,ytp_ph,ytp_th,
     &  ypppp,ytppp,ypppp_ph,ypppp_th,ytppp_ph,ytppp_th,
     &	ypp_thph,ypp_th2,ytp_thph,ytp_th2,
     &	ypppp_thph,ytppp_thph,ypppp_th2,ytppp_th2,wk)
c--------------------------------------------------evaluate model at this location
c--careful at difference in definition of azimuth between T&W (from North) and LT&E (from East).
	psi=90./radian-y(1)

cTEST
c	print*,psi*radian, y(1)*radian,xlon,xlat
c	pause

c---------------contribution of 0psi term
	c=sdot(len0,dlc0r,1,y0,1)
c---------------contribution of 2psi term
	xnu(1)=-cos(psi)
	xnu(2)=+sin(psi)

	tau(1,1)=-sdot(len2,dlc2r,1,ypp,1)
	tau(1,2)=+sdot(len2,dlc2r,1,ytp,1)
	tau(2,2)=-tau(1,1)
	tau(2,1)=+tau(1,2)
	call vecmatvec(xnu,tau,xnu,2,2,delc)
	c=c+delc
c---------------contribution of 4psi term
	do 2 i=1,2
	   do 2 j=1,2
	      do 2 k=1,2
	         do 2 l=1,2
2	            sigma(i,j,k,l)=0.
c--------theta,theta,theta,theta and phi,phi,phi,phi
	sigma(1,1,1,1)=sdot(len4,dlc4r,1,ypppp,1)
	sigma(2,2,2,2)=+sigma(1,1,1,1)
c--------theta,theta,theta,phi and permutations
	sigma(1,1,1,2)=-sdot(len4,dlc4r,1,ytppp,1)
	sigma(1,1,2,1)=+sigma(1,1,1,2)
	sigma(1,2,1,1)=+sigma(1,1,1,2)
	sigma(2,1,1,1)=+sigma(1,1,1,2)
c--------theta,theta,phi,phi and permutations
	sigma(1,1,2,2)=-sigma(1,1,1,1)
	sigma(1,2,2,1)=-sigma(1,1,1,1)
	sigma(1,2,1,2)=-sigma(1,1,1,1)
	sigma(2,2,1,1)=-sigma(1,1,1,1)
	sigma(2,1,1,2)=-sigma(1,1,1,1)
	sigma(2,1,2,1)=-sigma(1,1,1,1)
c--------phi,phi,phi,theta and permutations
	sigma(2,2,2,1)=-sigma(1,1,1,2)
	sigma(2,2,1,2)=-sigma(1,1,1,2)
	sigma(2,1,2,2)=-sigma(1,1,1,2)
	sigma(1,2,2,2)=-sigma(1,1,1,2)
	call vecvecmatvecvec(xnu,xnu,sigma,xnu,xnu,2,2,2,2,delc)
	c=c+delc

c---------------------calculate derivative of model at this location wrt azimuth y(1)
	dcdaz=0.
	xnu_az(1)=+sin(psi)
	xnu_az(2)=+cos(psi)
	call vecmatvec(xnu_az,tau,xnu,2,2,deldcdaz)
	dcdaz=dcdaz+2.*deldcdaz
	call vecvecmatvecvec(xnu_az,xnu,sigma,xnu,xnu,2,2,2,2,deldcdaz)
	dcdaz=dcdaz+4.*deldcdaz
c--------I found the psi-derivative. Turn it into zeta-derivative:
        dcdaz=-dcdaz

c------------------------second derivative with respect to azimuth
c---------exploits symmetry
	xnu_az2(1)=+cos(psi)
	xnu_az2(2)=-sin(psi)
	call vecmatvec(xnu_az,tau,xnu_az,2,2,dcdaz2_2)
	call vecmatvec(xnu,tau,xnu_az2,2,2,deldcdaz)
	dcdaz2_2=dcdaz2_2+deldcdaz
	call vecvecmatvecvec(xnu_az2,xnu,sigma,xnu,xnu,2,2,2,2,dcdaz2_4)
	call vecvecmatvecvec(xnu_az,xnu_az,sigma,xnu,xnu,2,2,2,2,deldcdaz)
	dcdaz2_4=dcdaz2_4+3.*deldcdaz
	dcdaz2=2.*dcdaz2_2+4.*dcdaz2_4
	dlncdaz2=-((dcdaz*dcdaz)/(c*c))+(dcdaz2/c)

c-----------calculate derivatives of model at this location wrt colatitude th and longitude ph
c---------------contribution of 0psi term
	dcdth=sdot(len0,dlc0r,1,yth,1)
	dcdph=sdot(len0,dlc0r,1,yph,1)
	dcdth2=sdot(len0,dlc0r,1,yth2,1)
	dcdthph=sdot(len0,dlc0r,1,ythph,1)

c---------------contribution of 2psi term
	tau_th(1,1)=-sdot(len2,dlc2r,1,ypp_th,1)
	tau_th(1,2)=+sdot(len2,dlc2r,1,ytp_th,1)
	tau_ph(1,1)=-sdot(len2,dlc2r,1,ypp_ph,1)
	tau_ph(1,2)=+sdot(len2,dlc2r,1,ytp_ph,1)
	tau_th(2,2)=-tau_th(1,1)
	tau_th(2,1)=+tau_th(1,2)
	tau_ph(2,2)=-tau_ph(1,1)
	tau_ph(2,1)=+tau_ph(1,2)
	tau_th2(1,1)=-sdot(len2,dlc2r,1,ypp_th2,1)
	tau_th2(1,2)=+sdot(len2,dlc2r,1,ytp_th2,1)
	tau_th2(2,2)=-tau_th2(1,1)
	tau_th2(2,1)=+tau_th2(1,2)
	tau_thph(1,1)=-sdot(len2,dlc2r,1,ypp_thph,1)
	tau_thph(1,2)=+sdot(len2,dlc2r,1,ytp_thph,1)
	tau_thph(2,2)=-tau_thph(1,1)
	tau_thph(2,1)=+tau_thph(1,2)
	call vecmatvec(xnu,tau_th,xnu,2,2,deldcdth)
	call vecmatvec(xnu,tau_ph,xnu,2,2,deldcdph)
	call vecmatvec(xnu,tau_th2,xnu,2,2,deldcdth2)
	call vecmatvec(xnu,tau_thph,xnu,2,2,deldcdthph)
	dcdth=dcdth+deldcdth
	dcdph=dcdph+deldcdph
	dcdth2=dcdth2+deldcdth2
	dcdthph=dcdthph+deldcdthph

c---------------contribution of 4psi term
	do 1 i=1,2
	   do 1 j=1,2
	      do 1 k=1,2
	         do 1 l=1,2
	            sigma_th(i,j,k,l)=0.
1	            sigma_ph(i,j,k,l)=0.
	sigma_th(1,1,1,1)=+sdot(len4,dlc4r,1,ypppp_th,1)
	sigma_th(1,1,1,2)=-sdot(len4,dlc4r,1,ytppp_th,1)
	sigma_th(2,2,2,2)=+sigma_th(1,1,1,1)
	sigma_th(1,1,2,1)=+sigma_th(1,1,1,2)
	sigma_th(1,2,1,1)=+sigma_th(1,1,1,2)
	sigma_th(2,1,1,1)=+sigma_th(1,1,1,2)
	sigma_th(1,1,2,2)=-sigma_th(1,1,1,1)
	sigma_th(1,2,2,1)=-sigma_th(1,1,1,1)
	sigma_th(1,2,1,2)=-sigma_th(1,1,1,1)
	sigma_th(2,2,1,1)=-sigma_th(1,1,1,1)
	sigma_th(2,1,1,2)=-sigma_th(1,1,1,1)
	sigma_th(2,1,2,1)=-sigma_th(1,1,1,1)
	sigma_th(2,2,2,1)=-sigma_th(1,1,1,2)
	sigma_th(2,2,1,2)=-sigma_th(1,1,1,2)
	sigma_th(2,1,2,2)=-sigma_th(1,1,1,2)
	sigma_th(1,2,2,2)=-sigma_th(1,1,1,2)
	sigma_ph(1,1,1,1)=+sdot(len4,dlc4r,1,ypppp_ph,1)
	sigma_ph(1,1,1,2)=-sdot(len4,dlc4r,1,ytppp_ph,1)
	sigma_ph(2,2,2,2)=+sigma_ph(1,1,1,1)
	sigma_ph(1,1,2,1)=+sigma_ph(1,1,1,2)
	sigma_ph(1,2,1,1)=+sigma_ph(1,1,1,2)
	sigma_ph(2,1,1,1)=+sigma_ph(1,1,1,2)
	sigma_ph(1,1,2,2)=-sigma_ph(1,1,1,1)
	sigma_ph(1,2,2,1)=-sigma_ph(1,1,1,1)
	sigma_ph(1,2,1,2)=-sigma_ph(1,1,1,1)
	sigma_ph(2,2,1,1)=-sigma_ph(1,1,1,1)
	sigma_ph(2,1,1,2)=-sigma_ph(1,1,1,1)
	sigma_ph(2,1,2,1)=-sigma_ph(1,1,1,1)
	sigma_ph(2,2,2,1)=-sigma_ph(1,1,1,2)
	sigma_ph(2,2,1,2)=-sigma_ph(1,1,1,2)
	sigma_ph(2,1,2,2)=-sigma_ph(1,1,1,2)
	sigma_ph(1,2,2,2)=-sigma_ph(1,1,1,2)
	sigma_th2(1,1,1,1)=+sdot(len4,dlc4r,1,ypppp_th2,1)
	sigma_th2(1,1,1,2)=-sdot(len4,dlc4r,1,ytppp_th2,1)
	sigma_th2(2,2,2,2)=+sigma_th2(1,1,1,1)
	sigma_th2(1,1,2,1)=+sigma_th2(1,1,1,2)
	sigma_th2(1,2,1,1)=+sigma_th2(1,1,1,2)
	sigma_th2(2,1,1,1)=+sigma_th2(1,1,1,2)
	sigma_th2(1,1,2,2)=-sigma_th2(1,1,1,1)
	sigma_th2(1,2,2,1)=-sigma_th2(1,1,1,1)
	sigma_th2(1,2,1,2)=-sigma_th2(1,1,1,1)
	sigma_th2(2,2,1,1)=-sigma_th2(1,1,1,1)
	sigma_th2(2,1,1,2)=-sigma_th2(1,1,1,1)
	sigma_th2(2,1,2,1)=-sigma_th2(1,1,1,1)
	sigma_th2(2,2,2,1)=-sigma_th2(1,1,1,2)
	sigma_th2(2,2,1,2)=-sigma_th2(1,1,1,2)
	sigma_th2(2,1,2,2)=-sigma_th2(1,1,1,2)
	sigma_th2(1,2,2,2)=-sigma_th2(1,1,1,2)
	sigma_thph(1,1,1,1)=+sdot(len4,dlc4r,1,ypppp_thph,1)
	sigma_thph(1,1,1,2)=-sdot(len4,dlc4r,1,ytppp_thph,1)
	sigma_thph(2,2,2,2)=+sigma_thph(1,1,1,1)
	sigma_thph(1,1,2,1)=+sigma_thph(1,1,1,2)
	sigma_thph(1,2,1,1)=+sigma_thph(1,1,1,2)
	sigma_thph(2,1,1,1)=+sigma_thph(1,1,1,2)
	sigma_thph(1,1,2,2)=-sigma_thph(1,1,1,1)
	sigma_thph(1,2,2,1)=-sigma_thph(1,1,1,1)
	sigma_thph(1,2,1,2)=-sigma_thph(1,1,1,1)
	sigma_thph(2,2,1,1)=-sigma_thph(1,1,1,1)
	sigma_thph(2,1,1,2)=-sigma_thph(1,1,1,1)
	sigma_thph(2,1,2,1)=-sigma_thph(1,1,1,1)
	sigma_thph(2,2,2,1)=-sigma_thph(1,1,1,2)
	sigma_thph(2,2,1,2)=-sigma_thph(1,1,1,2)
	sigma_thph(2,1,2,2)=-sigma_thph(1,1,1,2)
	sigma_thph(1,2,2,2)=-sigma_thph(1,1,1,2)
	call vecvecmatvecvec(xnu,xnu,sigma_th,xnu,xnu,2,2,2,2,deldcdth)
	call vecvecmatvecvec(xnu,xnu,sigma_ph,xnu,xnu,2,2,2,2,deldcdph)
	call vecvecmatvecvec(xnu,xnu,sigma_th2,xnu,xnu,2,2,2,2,deldcdth2)
	call vecvecmatvecvec(xnu,xnu,sigma_thph,xnu,xnu,2,2,2,2,deldcdthph)
	dcdth=dcdth+deldcdth
	dcdph=dcdph+deldcdph
	dcdth2=dcdth2+deldcdth2
	dcdthph=dcdthph+deldcdthph

c-----------find derivatives of ln(model)
	dlncdthp=dcdth/c
	dlncdphp=dcdph/c
	dlncdaz=dcdaz/c
	dlncdth2=-((dcdth*dcdth)/(c*c))+(dcdth2/c)
	dlncdthph=-((dcdth*dcdph)/(c*c))+(dcdthph/c)

c-----------multiple derivatives of model with respect also to azimuth
	call vecmatvec(xnu_az,tau_th,xnu,2,2,deldcdthaz)
	dcdthaz=2.*deldcdthaz
	call vecvecmatvecvec(xnu_az,xnu,sigma_th,xnu,xnu,2,2,2,2,deldcdthaz)
	dcdthaz=dcdthaz+4.*deldcdthaz
	call vecmatvec(xnu_az,tau_ph,xnu,2,2,deldcdphaz)
	dcdphaz=2.*deldcdphaz
	call vecvecmatvecvec(xnu_az,xnu,sigma_ph,xnu,xnu,2,2,2,2,deldcdphaz)
	dcdphaz=dcdphaz+4.*deldcdphaz

c--------Turn psi-derivative into zeta-derivative:
        dcdthaz=-dcdthaz
	dcdphaz=-dcdphaz
	dlncdthaz=-((dcdth*dcdaz)/(c*c))+(dcdthaz/c)
	dlncdphaz=-((dcdph*dcdaz)/(c*c))+(dcdphaz/c)

c-------------ray eqs (3.16) and (3.15) from Larson et al., GJI 1998
c-------------now use their definition of azimuth y(1) (from East)
c--multiply them by cos(az) at num and den, to make them stable
	dydx(1)=-cos(y(2))*cos(y(1))+sin(y(2))*cos(y(1))*dlncdthp 
     &    +sin(y(1))*dlncdphp
	dydx(1)=dydx(1)/(cos(y(1))-sin(y(1))*dlncdaz)
	dydx(2)=-sin(y(2))*(sin(y(1))+dlncdaz*cos(y(1)))
	dydx(2)=dydx(2)/(cos(y(1))-sin(y(1))*dlncdaz)
c--eqs. (3.21) and (3.22) of Larson et al., to calculate geometrical spreading
	deno=1./(1.-tan(y(1))*dlncdaz)
	facto=tan(y(1))+dlncdaz
	coszm2=1./(cos(y(1))**2)
	facto2=sin(y(2))*dlncdthp+tan(y(1))*dlncdphp-cos(y(2))
c-----------------------------------(3.21) ! gives focusing
	coef3=cos(y(2))*facto
	coef3=coef3+ sin(y(2))*(1.+deno*tan(y(1))*facto)*dlncdthaz ! nonzero only if aniso
	coef4=coszm2+dlncdaz2 
	coef4=coef4+facto*deno*(coszm2*dlncdaz+tan(y(1))*dlncdaz2)
	coef4=coef4*sin(y(2))
	dydx(3)=-deno*(coef3*y(3)+coef4*y(4))  ! coef4 is zero in isotropic case
c-----------------------------------(3.22)
	coef3=cos(y(2))*dlncdthp+sin(y(2))*dlncdth2+tan(y(1))*dlncdthph+sin(y(2))
	coef3=coef3+ deno*tan(y(1))*dlncdthaz*facto2  ! nonzero only if aniso
	coef4=coszm2*dlncdphp
	coef4=coef4+sin(y(2))*dlncdthaz+tan(y(1))*dlncdphaz ! nonzero only if aniso
	coef4=coef4+facto2*deno*(coszm2*dlncdaz+tan(y(1))*dlncdaz2)
	dydx(4)=deno*(coef3*y(3)+coef4*y(4))
c---------(3.18) needed to find phase
	crad=c/rearth
	waven=omega/crad
	phas=sin(y(2))*deno/cos(y(1))
	phas=waven*phas
	dydx(5)=phas
	return
	end

c-------------------------------------------------------------------

      subroutine odeint
     & (ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,nstp,derivs)
c--integration of the ray-tracing equations between source and receiver
c--for a given value ystart(1) of zeta at the source 
c      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
c      real*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      real eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs
      PARAMETER (nmax=5,KMAXX=10000,maxstp=kmaxx,TINY=1.e-30)
      INTEGER i,kmax,kount,nstp
c      real*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
      real dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     *yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
c      common/integration_method/imet

	kmax=KMAXX
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0

      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue

      if (kmax.gt.0) xsav=x-2.*dxsav
c-------------------------------------------------main loop
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
c	do 12 i=1,nvar
c          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
c12      continue
c-------------------the following is convenient since values of y() are bounded
	yscal(1)=1.
	yscal(2)=3.
c--------------I am not sure what to do with y(3) and y(4) so stick to Press et al.
	do 12 i=3,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue

c-------------------------store path
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif

        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x

cTEST	
c	print*,"call rkqs with x=",x

	call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)

cTEST
c	write(*,"(9(e11.5,1x))")x,y,dydx

        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif

        if((x-x2)*(x2-x1).ge.0.)then  ! i.e. if x outside (x1,x2) 
          do 14 i=1,nvar
            ystart(i)=y(i)   ! solution at or near x2
14        continue

c-----------------------------store path (last point)
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
cTEST
c	  print*,"integrated after ",nok,nstp," runge kutta steps"
          return
        endif

        if(abs(hnext).lt.hmin) pause
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      print*,"too many steps in odeint"
	pause

      return
      END

c-------------------------------------------------------------------

      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c	implicit real*8 (a-h,o-z)
      INTEGER n,NMAX
c      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (nmax=5)
CU    USES derivs,rkck
      INTEGER i
c      REAL*8 errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
      REAL errmax,h,htemp,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     *PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)

cTEST
c	common/cTEST/dcdaz,dcdt,dcdp
c	common/cTEST/dcdaz,dcdt,dcdp,dcdt2,dcdtp,dcdaz2,dcdpaz,dcdtaz
cTEST
c	print*,"in rkqs"

      h=htry
1      call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))

cTEST
c	print*,i,yerr(i),yscal(i)

11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x)then
	   print*,'stepsize underflow in rkqs'
c	   pause
	   return	! Lapo 04 oct 03
	endif
        goto 1
      else

cTEST
c	write(99,"(e16.9,1x,e12.5,1x,e16.9,12(1x,e10.3))")x,ytemp(1),ytemp(2),dydx(1),dydx(2),
c     &yerr(1),yerr(2),dcdaz,dcdt,dcdp,dcdt2,dcdtp,dcdaz2,dcdpaz,dcdtaz

        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
12      continue
        return
      endif
      END

c-------------------------------------------------------------------

      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
c	implicit real*8 (a-h,o-z)
      INTEGER n,NMAX
c      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      REAL h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (nmax=5)
CU    USES derivs
      INTEGER i
c      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
      REAL ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)
      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      END      

c------------------------------------------------------
c set up rotation matrix corresponding to Euler angles
c alph, beta, gama

      subroutine setrot(alph,beta,gama,urot)
      dimension urot(3,3)
      sal=sin(alph)
      cal=cos(alph)
      sbe=sin(beta)
      cbe=cos(beta)
      sga=sin(gama)
      cga=cos(gama)
      urot(1,1)=cga*cbe*cal-sga*sal
      urot(1,2)=cga*cbe*sal+sga*cal
      urot(1,3)=-cga*sbe
      urot(2,1)=-sga*cbe*cal-cga*sal
      urot(2,2)=-sga*cbe*sal+cga*cal
      urot(2,3)=sga*sbe
      urot(3,1)=sbe*cal
      urot(3,2)=sbe*sal
      urot(3,3)=cbe
      return
      end


c------------------------------------------------------
c rotate lat x lony to new values using rotation matrix urot.

      subroutine rotll(x,y,x1,y1,urot)
      dimension urot(3,3),v(3),v1(3)
      data radian/57.29578/
      a=x/radian
      b=y/radian
      sa=sin(a)
      ca=cos(a)
      sb=sin(b)
      cb=cos(b)
      v(1)=ca*cb
      v(2)=ca*sb
      v(3)=sa
      do 10 i=1,3
      v1(i)=0.
      do 10 j=1,3
   10 v1(i)=v1(i)+urot(i,j)*v(j)
      x1=atan2(v1(3),sqrt(abs(1.-v1(3)*v1(3))))*radian
      y1=atan2(v1(2),v1(1))*radian
      return
      end

c----------------------------------------------------------------
c find the pole of the path xlatp, xlonp and the azimuth
c at the pole of the epicentre (I thnk).

      subroutine pole(epla,eplo,stla,stlo,xlap,xlop,azmp,delta)
      implicit double precision (a-h,o-z)
      real epla,eplo,stla,stlo,xlap,xlop,azmp,delta
      data radian/57.295779513082d0/
      th1=(90.d0-epla)/radian
      th2=(90.d0-stla)/radian
      ph1=eplo/radian
      ph2=stlo/radian
      cth1=dcos(th1)
      cth2=dcos(th2)
      sth1=dsin(th1)
      sth2=dsin(th2)
      cph1=dcos(ph1)
      cph2=dcos(ph2)
      sph1=dsin(ph1)
      sph2=dsin(ph2)
      cph21=cph1*cph2+sph1*sph2
      sph21=sph2*cph1-sph1*cph2
      cdel=sth1*sth2*cph21+cth1*cth2
      sdel=dsqrt(1.d0-cdel**2)
      delta=datan2(sdel,cdel)*radian
      ccth=sth1*sth2*sph21/sdel
      scth=dsqrt(1.d0-ccth**2)
      cth=datan2(scth,ccth)
      coco=1.d0/(sdel*scth)
      scph=coco*(cth1*sth2*cph2-cth2*sth1*cph1)
      ccph=coco*(sth1*cth2*sph1-sth2*cth1*sph2)
      cph=datan2(scph,ccph)
      xlap=90.d0-cth*radian
      xlop=cph*radian
      v1=sth1*cph1+sth2*cph2
      v2=sth1*sph1+sth2*sph2
      v3=cth1+cth2
      vn=dsqrt(v1**2+v2**2+v3**2)
      v1=v1/vn
      v2=v2/vn
      v3=v3/vn
      cthm=v3
      sthm=dsqrt(1.d0-v3**2)
      cphm=v1/sthm
      sphm=v2/sthm
      cphmp=ccth*sthm*(cphm*ccph+sphm*scph)-scth*cthm
      sphmp=sthm*(sphm*ccph-cphm*scph)
      phmp=datan2(sphmp,cphmp)
      azmp=180.d0-phmp*radian
      return
      end

c-------------------------------------------------------------------
      subroutine ylm(xlat,xlon,lmax,y,wk1,wk2,wk3)
c
      complex temp,fac,dfac
      dimension wk1(1),wk2(1),wk3(1),y(1)
c
c     wk1,wk2,wk3 should be dimensioned at least (lmax+1)*4
c
      data radian/57.2957795/
c
      theta=(90.-xlat)/radian
      phi=xlon/radian
c
      ind=0
      lm1=lmax+1
      do 10 il1=1,lm1
      l=il1-1
      call legndr(theta,l,l,wk1,wk2,wk3)
c
      fac=(1.,0.)
      dfac=cexp(cmplx(0.,phi))
c
      do 20 im=1,il1
      temp=fac*cmplx(wk1(im),0.)
      ind=ind+1
      y(ind)=real(temp)
      if(im.eq.1) goto 20
      ind=ind+1
      y(ind)=aimag(temp)
   20 fac=fac*dfac
c
   10 continue
      return
      end
c-------------------------------------------------------------------
      SUBROUTINE LEGNDR(THETA,L,M,X,XP,XCOSEC)
      DIMENSION X(*),XP(*),XCOSEC(*)
      DOUBLE PRECISION SMALL,SUM,COMPAR,CT,ST,FCT,COT,FPI,X1,X2,X3,
     1F1,F2,XM,TH,DFLOAT
      DATA FPI/12.56637062D0/
      DFLOAT(I)=FLOAT(I)
      SUM=0.D0
      LP1=L+1
      TH=THETA
      CT=DCOS(TH)
      ST=DSIN(TH)
      MP1=M+1
      FCT=DSQRT(DFLOAT(2*L+1)/FPI)
      SFL3=SQRT(FLOAT(L*(L+1)))
      COMPAR=DFLOAT(2*L+1)/FPI
      DSFL3=SFL3
      SMALL=1.D-16*COMPAR
      DO 1 I=1,MP1
      X(I)=0.
      XCOSEC(I)=0.
    1 XP(I)=0.
      IF(L.GT.1.AND.ABS(THETA).GT.1.E-5) GO TO 3
      X(1)=FCT
      IF(L.EQ.0) RETURN
      X(1)=CT*FCT
      X(2)=-ST*FCT/DSFL3
      XP(1)=-ST*FCT
      XP(2)=-.5D0*CT*FCT*DSFL3
      IF(ABS(THETA).LT.1.E-5) XCOSEC(2)=XP(2)
      IF(ABS(THETA).GE.1.E-5) XCOSEC(2)=X(2)/ST
      RETURN
    3 X1=1.D0
      X2=CT
      DO 4 I=2,L
      X3=(DFLOAT(2*I-1)*CT*X2-DFLOAT(I-1)*X1)/DFLOAT(I)
      X1=X2
    4 X2=X3
      COT=CT/ST
      COSEC=1./ST
      X3=X2*FCT
      X2=DFLOAT(L)*(X1-CT*X2)*FCT/ST
      X(1)=X3
      X(2)=X2
      SUM=X3*X3
      XP(1)=-X2
      XP(2)=DFLOAT(L*(L+1))*X3-COT*X2
      X(2)=-X(2)/SFL3
      XCOSEC(2)=X(2)*COSEC
      XP(2)=-XP(2)/SFL3
      SUM=SUM+2.D0*X(2)*X(2)
      IF(SUM-COMPAR.GT.SMALL) RETURN
      X1=X3
      X2=-X2/DSQRT(DFLOAT(L*(L+1)))
      DO 5 I=3,MP1
      K=I-1
      F1=DSQRT(DFLOAT((L+I-1)*(L-I+2)))
      F2=DSQRT(DFLOAT((L+I-2)*(L-I+3)))
      XM=K
      X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
      SUM=SUM+2.D0*X3*X3
      IF(SUM-COMPAR.GT.SMALL.AND.I.NE.LP1) RETURN
      X(I)=X3
      XCOSEC(I)=X(I)*COSEC
      X1=X2
      XP(I)=-(F1*X2+XM*COT*X3)
    5 X2=X3
      RETURN
      END
c-------------------------------------------------------------------
cprog rotmx2
cxref
      subroutine rotmx2(nmax,l,theta,d,id1,id2)
      implicit double precision (a-h,o-z)
      double precision d,theta
      dimension d(id1,id2)
c     data big,small,dlbig,dlsml/1.d35,1.d-35,35.d0,-35.d0/
      data big,small,dlbig,dlsml/1.d25,1.d-25,25.d0,-25.d0/
      data pi/3.14159265358979d0/
      dfloat(n)=n
      th=theta
      if((th.gt.pi).or.(th.lt.0.d0)) stop 'illegal arg in rotmx2'
      if(l.ne.0) goto 350
      d(1+nmax,l+1)=1.d0
      return
350   isup=1
      if(th.le.pi/2.d0) goto 310
      th=pi-th
      isup=-1
310   nm=2*l+1
      nmp1=nm+1
      lp1=l+1
      lm1=l-1
      lp2=l+2
      nrow=2*nmax+1
      nmaxp1=nmax+1
      lmn=l-nmax
      if(th.ne.0.d0) goto 320
      do 330 im1ct=1,nrow
      im1=im1ct+lmn
      do 330 im2=lp1,nm
      d(im1ct,im2)=0.d0
      if(im1.eq.im2) d(im1ct,im2)=1.d0
330   continue
      goto 400
320   continue
c
c     zero l.h.s. of matrix
c
      do 340 im1=1,nrow
      do 340 im2=1,lp1
340   d(im1,im2)=0.d0
c
c        set up parameters
c
      shth=dsin(0.5d0*th)
      chth=dcos(0.5d0*th)
      sth=2.d0*shth*chth
      cth=2.d0*chth*chth-1.d0
      dlogf=dlog10(chth/shth)
      dlogs=dlog10(shth)
c
c       iterate from last column using 1. as starting value
c
      do 10 im1ct=1,nrow
      im1=im1ct+lmn
      m1=im1-lp1
      rm1=m1
      nm2=min0(im1-1,nm-im1)
      d(im1ct,nm)=1.d0
      if(nm2.eq.0) goto 10
      do 20 nit=1,nm2
      m2=l-nit
      im2=m2+lp1
      if(m2.ne.lm1) goto 70
      t1=0.d0
      goto 30
70    t1=-dsqrt(dfloat((im2+1)*(l-m2-1)))*d(im1ct,im2+2)
30    d(im1ct,im2)=t1-(2.d0/sth)*(cth*dfloat(m2+1)-rm1)
     1    *d(im1ct,im2+1)
      d(im1ct,im2)=d(im1ct,im2)/dsqrt(dfloat(im2*(l-m2)))
      temp=d(im1ct,im2)
      rmod=dabs(temp)
      if(rmod.lt.big) goto 20
      if(nit.eq.nm2) goto 20
      d(im1ct,nit+1)=dlbig
      d(im1ct,im2)=d(im1ct,im2)/big
      d(im1ct,im2+1)=d(im1ct,im2+1)/big
20    continue
10    continue
c
c        set up normalization for rightmost column
c
      t1=dfloat(2*l)*dlogs
      if(lmn.eq.0) goto 720
      do 710 i=1,lmn
      m1=i-l
      t1=dlogf+0.5d0*dlog10(dfloat(lp1-m1)/dfloat(l+m1))+t1
710   continue
720   d(1,1)=t1
      if(nrow.eq.1) goto 730
      do 110 im1ct=2,nrow
      m1=im1ct-nmaxp1
110   d(im1ct,1)=dlogf+0.5d0*dlog10(dfloat(l-m1+1)/dfloat(l+m1))
     1     +d(im1ct-1,1)
730   sgn=-1.d0
      if((lmn/2)*2.ne.lmn) sgn=1.d0
c
c       renormalize rows
c
      do 120 im1ct=1,nrow
      im1=im1ct+lmn
      sgn=-sgn
      csum=d(im1ct,1)
      mult=1
520   if(dabs(csum).lt.dlbig) goto 510
      mult=mult*2
      csum=0.5*csum
      goto 520
510   fac=10.d0**csum
      sfac=small/fac
      nm2=min0(im1-1,nm-im1)
      nm2p1=nm2+1
      do 130 im2=1,nm2p1
      if((d(im1ct,im2+1).eq.0.d0).or.(im2.ge.nm2)) goto 250
      csum=csum*dfloat(mult)+d(im1ct,im2+1)
      mult=1
220   if(dabs(csum).lt.dlbig) goto 210
      mult=mult*2
      csum=0.5d0*csum
      goto 220
210   fac=10.d0**csum
      sfac=small/fac
250   in2=nmp1-im2
      do 270 i=1,mult
      temp=d(im1ct,in2)
      rmod=dabs(temp)
      if(rmod.gt.sfac) goto 260
      d(im1ct,in2)=0.d0
      goto 130
260   d(im1ct,in2)=d(im1ct,in2)*fac
270   continue
      d(im1ct,in2)=sgn*d(im1ct,in2)
130   continue
120   continue
c
c       fill rest of matrix
c
400   if(isup.gt.0) goto 410
      sgn=-1.d0
      if((lmn/2)*2.ne.lmn) sgn=1.d0
      do 420 im1ct=1,nrow
      sgn=-sgn
      im1=im1ct+lmn
      nm2=min0(im1,nmp1-im1)
      do 420 in2=1,nm2
      im2=nmp1-in2
420   d(im1ct,in2)=sgn*d(im1ct,im2)
      do 430 im1ct=1,nrow
      im1=im1ct+lmn
      in1=nmp1-im1
      in1ct=in1-lmn
      sgn=-1.d0
      nm2=min0(im1,in1)
      do 440 nit=1,nm2
      sgn=-sgn
      im2=1+nm2-nit
      in2=nmp1-im2
      im2ct=im2-lmn
      in2ct=in2-lmn
      d(in1ct,in2)=sgn*d(im1ct,im2)
      if(in2ct.gt.nrow) goto 440
      d(im2ct,im1)=d(in1ct,in2)
      d(in2ct,in1)=d(im1ct,im2)
440   continue
430   continue
      return
410   do 450 im1ct=1,nrow
      im1=im1ct+lmn
      in1=nmp1-im1
      in1ct=in1-lmn
      sgn=-1.d0
      nm2=min0(im1,in1)
      do 460 nit=1,nm2
      sgn=-sgn
      im2=nm-nm2+nit
      in2=nmp1-im2
      im2ct=im2-lmn
      in2ct=in2-lmn
      d(in1ct,in2)=sgn*d(im1ct,im2)
      if(im2ct.gt.nrow) goto 460
      d(im2ct,im1)=d(in1ct,in2)
      d(in2ct,in1)=d(im1ct,im2)
460   continue
450   continue
      return
      end
c-------------------------------------------------------------------
c        subroutine ylmv4_der()
c
c  evaluates at xlon, xlat, generalized spherical harmonics needed to describe
c  completely symmetric, trace-free, second (gsh N=2) and fourth (gsh N=4) rank 
c  tensors (Trampert and Woodhouse, GJI, 2003), and places them in the
c  vectors ypp,ytp, of length lenv=(lmax-1)*(2*lmax+6)
c  and vectors ypppp,ytppp, of length lenv4=(lmax-3)*(2*lmax+10).
c  Also evaluates  their first derivatives with 
c  respect to theta, phi: compare with subroutine ylmv4() which 
c  does not do derivatives.
c
c  As in ylmv4, generalized spherical harmonics are computed making use
c  of the coincidence between generalized legendre polynomial P_Nlm(cos(beta))
c  and rotation tensor D^(l)_Nm(beta) (e.g., Dahlen and Tromp, 1998, eq. C.255).
c  This is why the subroutine rotmx2, that calculates rotation matrices,
c  is used to find generalized legendre polynomials.
c
c  derivatives of ytp, ypp wrt theta, phi are ytp_th, ytp_ph, ypp_th, ypp_ph
c  derivatives of ypppp, ytppp are ypppp_th, etc.
c
c  Scalar spherical harmonics are also calculated and
c  placed in y(1) -- y( (lmax+1)**2 ).
c  Their derivatives are placed in yth, yph.
c
c  all derivatives (N=0,2,4) are found by means of Dahlen and Tromp's eq. C.119
c  (appendix C, page 899 of the 1998 edition).  
c
c  d is double precision workspace.
      subroutine ylmv4_der(xlat,xlon,lmax,y,ypp,ytp,yph,yth,
     &  ypp_th,ypp_ph,ytp_ph,ytp_th,
     &  ypppp,ytppp,ypppp_ph,ypppp_th,ytppp_ph,ytppp_th,
     &  ypp_thph,ypp_th2,ytp_thph,ytp_th2,
     &  ypppp_thph,ytppp_thph,ypppp_th2,ytppp_th2,d)
	parameter(lmx=55)
      dimension y((lmx+1)**2)
     1         ,ypp((lmx-1)*(2*lmx+6))
     1         ,ytp((lmx-1)*(2*lmx+6))
     1         ,ypppp((lmx-3)*(2*lmx+10))
     1         ,ytppp((lmx-3)*(2*lmx+10))
     1         ,ypp_ph((lmx-1)*(2*lmx+6))
     1         ,ytp_ph((lmx-1)*(2*lmx+6))
     1         ,ypp_th((lmx-1)*(2*lmx+6))
     1         ,ytp_th((lmx-1)*(2*lmx+6))
     1         ,ypp_thph((lmx-1)*(2*lmx+6))
     1         ,ytp_thph((lmx-1)*(2*lmx+6))
     1         ,ypp_th2((lmx-1)*(2*lmx+6))
     1         ,ytp_th2((lmx-1)*(2*lmx+6))
     1         ,ypppp_ph((lmx-3)*(2*lmx+10))
     1         ,ytppp_ph((lmx-3)*(2*lmx+10))
     1         ,ypppp_th((lmx-3)*(2*lmx+10))
     1         ,ytppp_th((lmx-3)*(2*lmx+10))
     1         ,ypppp_thph((lmx-3)*(2*lmx+10))
     1         ,ytppp_thph((lmx-3)*(2*lmx+10))
     1         ,ypppp_th2((lmx-3)*(2*lmx+10))
     1         ,ytppp_th2((lmx-3)*(2*lmx+10))
     &         ,yth((lmx+1)**2),yph((lmx+1)**2)
      double precision theta,d(9*(2*lmx+1))
      complex cfac,dfac,cfac0

      data radian/57.2957795/,rr4pi/0.28209479/

      theta=dble((90.-xlat)/radian)
      dfac=cexp(cmplx(0.,xlon/radian))
      k=0
      ka=0
      ka4=0
      do l=0,lmax
        if(l.lt.2) then
c-------------for l less than 2 we only calculate scalar spherical harmonics
c-------------as no generalized spherical harmonics with N=2 or 4 (the only
c-------------ones we need) are defined. First argument of rotmx2 is value of N.
          call rotmx2(0,l,theta,d,1,2*l+1)
          ind=l
          cfac=rr4pi*sqrt(2.*l+1)
          k0=k+1
          ind0=ind+1	! this is okay lapo 13 dec 03
c          ind0=ind
          cfac0=cfac
          do m=0,l
            k=k+1
            ind=ind+1
            y(k)=d(ind)*real(cfac)
            yph(k)=-float(m)*d(ind)*aimag(cfac)
            if(m.ne.0) then
              dplm0_th=-float(m)*(cos(theta)/sin(theta))*d(ind) 
     &             -sqrt(float((l+m)*(l-m+1)))*d(ind-1)
              yth(k)=dplm0_th*real(cfac)
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              yth(k)=dplm0_th*aimag(cfac)
              yph(k)=float(m)*d(ind)*real(cfac)
            endif
            cfac=cfac*dfac
          enddo
          if(l.eq.0)then
             dplm0_th=0.
          else
             if(ind0.gt.(2*l+1))stop" error in ylmv4_der"
             dplm0_th=sqrt(float(l*(l+1)))*d(ind0+1)
          endif
          yth(k0)=dplm0_th*real(cfac0)
        else if(l.lt.4) then
c-------------for l=2,3, calculate scalar spherical harmonics and generalized
c-------------spherical harmonics with N=2.
          call rotmx2(2,l,theta,d,5,2*l+1)
          ind=5*l+3
          indp=ind+2
          indm=indp
          cfac=rr4pi*sqrt(2.*l+1)
          k0=k+1
          ind0=ind
          cfac0=cfac
	  ka0=ka+1
	  indp0=indp
          do m=0,l
            k=k+1
            y(k)=d(ind)*real(cfac)
            yph(k)=-float(m)*d(ind)*aimag(cfac)
            ka=ka+1
            ypp(ka)=-d(indp)*real(cfac)
            ytp(ka)=-d(indp)*aimag(cfac)
            ypp_ph(ka)=float(m)*d(indp)*aimag(cfac)
            ytp_ph(ka)=-float(m)*d(indp)*real(cfac)
            ka=ka+1
            ypp(ka)=+d(indp)*aimag(cfac)
            ytp(ka)=-d(indp)*real(cfac)
            ypp_ph(ka)=float(m)*d(indp)*real(cfac)
            ytp_ph(ka)=float(m)*d(indp)*aimag(cfac)
            if(m.ne.0) then
              dplm0_th=-float(m)*(cos(theta)/sin(theta))*d(ind) 
     &             -sqrt(float((l+m)*(l-m+1)))*d(ind-5)
              yth(k)=dplm0_th*real(cfac)
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              yth(k)=dplm0_th*aimag(cfac)
              yph(k)=float(m)*d(ind)*real(cfac)
              dplm2p_th=-float(m)*(cos(theta)/sin(theta))*d(indp) 
     &   	   +2.*d(indp)/sin(theta)
     &             -sqrt(float((l+m)*(l-m+1)))*d(indp-5)
	      dplm2p_th2=(
     &           (( (2.-float(m)*cos(theta))**2 +float(m)-2.*cos(theta))
     &	            /(sin(theta)*sin(theta)))
     &	           *d(indp))
c     &	           +(sqrt(float((l+m)*(l-m+1)))*(4.-float(2*m-1)
     &	           -(sqrt(float((l+m)*(l-m+1)))*(4.-float(2*m-1)
     &	           *cos(theta))/sin(theta))*d(indp-5)
     &	           +sqrt(float((l-m+1)*(l-m+2)*(l+m-1)*(l+m)))*d(indp-10)

	      ypp_th(ka-1)=-dplm2p_th*real(cfac)
	      ytp_th(ka-1)=-dplm2p_th*aimag(cfac)
	      ypp_th(ka)=+dplm2p_th*aimag(cfac)
	      ytp_th(ka)=-dplm2p_th*real(cfac)
c
	ypp_th2(ka-1)=-dplm2p_th2*real(cfac)
	ytp_th2(ka-1)=-dplm2p_th2*aimag(cfac)
	ypp_th2(ka)=+dplm2p_th2*aimag(cfac)
	ytp_th2(ka)=-dplm2p_th2*real(cfac)

	ypp_thph(ka-1)=+float(m)*dplm2p_th*aimag(cfac)
	ytp_thph(ka-1)=-float(m)*dplm2p_th*real(cfac)
	ypp_thph(ka)=+float(m)*dplm2p_th*real(cfac)
	ytp_thph(ka)=+float(m)*dplm2p_th*aimag(cfac)
c
              dplm2m_th=-float(m)*(cos(theta)/sin(theta))*d(indm)
     &             -2.*d(indm)/sin(theta)
     &             +sqrt(float((l+m)*(l-m+1)))*d(indm+5)

	      dplm2m_th2=(( (-2.-float(m)*cos(theta)) **2
     &	           +float(m)+2.*cos(theta)) / (sin(theta)*sin(theta))
     &	           *d(indm))
     &	           +sqrt(float((l+m)*(l-m+1)))*(-4.+float(-2*m+1)
     &	           *cos(theta))/sin(theta)*d(indm+5)
     &	           +sqrt(float((l-m+1)*(l-m+2)*(l+m-1)*(l+m)))*d(indm+10)

              ka=ka+1
              ypp(ka)=-d(indm)*real(cfac)
              ytp(ka)=+d(indm)*aimag(cfac)
	      ypp_ph(ka)=+float(m)*d(indm)*aimag(cfac)
	      ytp_ph(ka)=+float(m)*d(indm)*real(cfac)
	      ypp_th(ka)=-dplm2m_th*real(cfac)
	      ytp_th(ka)=+dplm2m_th*aimag(cfac)
c
	ypp_th2(ka)=-dplm2m_th2*real(cfac)
	ytp_th2(ka)=+dplm2m_th2*aimag(cfac)
	ypp_thph(ka)=+float(m)*dplm2m_th*aimag(cfac)
	ytp_thph(ka)=+float(m)*dplm2m_th*real(cfac)
c
              ka=ka+1
              ypp(ka)=-d(indm)*aimag(cfac)
              ytp(ka)=-d(indm)*real(cfac)
	      ypp_ph(ka)=-float(m)*d(indm)*real(cfac)
	      ytp_ph(ka)=+float(m)*d(indm)*aimag(cfac)
	      ypp_th(ka)=-dplm2m_th*aimag(cfac)
	      ytp_th(ka)=-dplm2m_th*real(cfac)
c
	ypp_th2(ka)=-dplm2m_th2*aimag(cfac)
	ytp_th2(ka)=-dplm2m_th2*real(cfac)
	ypp_thph(ka)=-float(m)*dplm2m_th*real(cfac)
	ytp_thph(ka)=+float(m)*dplm2m_th*aimag(cfac)
c
            endif
            ind=ind+5
            indp=indp+5
            indm=indm-5
            cfac=cfac*dfac
          enddo
          dplm0_th=sqrt(float(l*(l+1)))*d(ind0+5)
          yth(k0)=dplm0_th*real(cfac0)
          dplm2p_th=-2.*d(indp0)/sin(theta)
     &          +sqrt(float(l*(l+1)))*d(indp0+5)


	  dplm2p_th2=((4.
     &	        +2.*cos(theta))/(sin(theta)*sin(theta))
     &          *d(indp0))
     &          +sqrt(float((l)*(l+1)))*(-4.+
c     &          -sqrt(float((l)*(l+1)))*(-2.+
     &          cos(theta))/sin(theta)*d(indp0+5)
     &          +sqrt(float((l+1)*(l+2)*(l-1)*l))*d(indp0+10)

	  ypp_th(ka0)=-dplm2p_th*real(cfac0)
	  ytp_th(ka0)=-dplm2p_th*aimag(cfac0)
	  ypp_th(ka0+1)=+dplm2p_th*aimag(cfac0)
	  ytp_th(ka0+1)=-dplm2p_th*real(cfac0)
c
	ypp_th2(ka0)=-dplm2p_th2*real(cfac0)
	ytp_th2(ka0)=-dplm2p_th2*aimag(cfac0)
c	ypp_thph(ka0)=+float(m)*dplm2p_th*aimag(cfac0)
c	ytp_thph(ka0)=-float(m)*dplm2p_th*real(cfac0)
	ypp_thph(ka0)=0.
	ytp_thph(ka0)=0.
	ypp_th2(ka0+1)=+dplm2p_th2*aimag(cfac0)
	ytp_th2(ka0+1)=-dplm2p_th2*real(cfac0)
c	ypp_thph(ka0+1)=+float(m)*dplm2p_th*real(cfac0)
c	ytp_thph(ka0+1)=+float(m)*dplm2p_th*aimag(cfac0)
	ypp_thph(ka0+1)=0.
	ytp_thph(ka0+1)=0.
c
        else
c-------------for l=4 or larger calculate scalar spherical harmonics and generalized
c-------------spherical harmonics with N=2 and N=4.
          call rotmx2(4,l,theta,d,9,2*l+1)
          ind=9*l+5
          indp=ind+2
          indm=indp
          indp4=ind+4
          indm4=indp4
          cfac=rr4pi*sqrt(2.*l+1)
          k0=k+1
          ind0=ind
          cfac0=cfac
	  ka0=ka+1
	  indp0=indp
	  ka40=ka4+1
	  indp40=indp4
          do m=0,l
            k=k+1
            y(k)=d(ind)*real(cfac)
            yph(k)=-float(m)*d(ind)*aimag(cfac)
            ka=ka+1
            ypp(ka)=-d(indp)*real(cfac)
            ytp(ka)=-d(indp)*aimag(cfac)
            ypp_ph(ka)=float(m)*d(indp)*aimag(cfac)
            ytp_ph(ka)=-float(m)*d(indp)*real(cfac)
            ka=ka+1
            ypp(ka)=+d(indp)*aimag(cfac)
            ytp(ka)=-d(indp)*real(cfac)
            ypp_ph(ka)=+float(m)*d(indp)*real(cfac)
            ytp_ph(ka)=+float(m)*d(indp)*aimag(cfac)
            ka4=ka4+1
            ypppp(ka4)=-d(indp4)*real(cfac)
            ytppp(ka4)=-d(indp4)*aimag(cfac)
	    ypppp_ph(ka4)=+float(m)*d(indp4)*aimag(cfac)
	    ytppp_ph(ka4)=-float(m)*d(indp4)*real(cfac)
            ka4=ka4+1
            ypppp(ka4)=+d(indp4)*aimag(cfac)
            ytppp(ka4)=-d(indp4)*real(cfac)
	    ypppp_ph(ka4)=+float(m)*d(indp4)*real(cfac)
	    ytppp_ph(ka4)=+float(m)*d(indp4)*aimag(cfac)
            if(m.ne.0) then
              dplm0_th=-float(m)*(cos(theta)/sin(theta))*d(ind) 
     &             -sqrt(float((l+m)*(l-m+1)))*d(ind-9)
              yth(k)=dplm0_th*real(cfac)
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              yth(k)=dplm0_th*aimag(cfac)
              yph(k)=float(m)*d(ind)*real(cfac)
              dplm2p_th=-float(m)*(cos(theta)/sin(theta))*d(indp) 
     &   	   +2.*d(indp)/sin(theta)
     &             -sqrt(float((l+m)*(l-m+1)))*d(indp-9)

	      dplm2p_th2=(((2.-float(m)*cos(theta))**2
     &	           +float(m)-2.*cos(theta))/(sin(theta)*sin(theta))
     &	           *d(indp))
c     &	           +sqrt(float((l+m)*(l-m+1)))*(4.-float(2*m-1)
     &	           -sqrt(float((l+m)*(l-m+1)))*(4.-float(2*m-1)
     &	           *cos(theta))/sin(theta)*d(indp-9)
     &	           +sqrt(float((l-m+1)*(l-m+2)*(l+m-1)*(l+m)))*d(indp-18)

	      ypp_th(ka-1)=-dplm2p_th*real(cfac)
	      ytp_th(ka-1)=-dplm2p_th*aimag(cfac)
	      ypp_th(ka)=+dplm2p_th*aimag(cfac)
	      ytp_th(ka)=-dplm2p_th*real(cfac)
c
	ypp_th2(ka-1)=-dplm2p_th2*real(cfac)
	ytp_th2(ka-1)=-dplm2p_th2*aimag(cfac)
	ypp_th2(ka)=+dplm2p_th2*aimag(cfac)
	ytp_th2(ka)=-dplm2p_th2*real(cfac)
	ypp_thph(ka-1)=+float(m)*dplm2p_th*aimag(cfac)
	ytp_thph(ka-1)=-float(m)*dplm2p_th*real(cfac)
	ypp_thph(ka)=+float(m)*dplm2p_th*real(cfac)
	ytp_thph(ka)=+float(m)*dplm2p_th*aimag(cfac)
c
              dplm2m_th=-float(m)*(cos(theta)/sin(theta))*d(indm) 
     &             -2.*d(indm)/sin(theta)
     &             +sqrt(float((l+m)*(l-m+1)))*d(indm+9)

	      dplm2m_th2=(( (-2.-float(m)*cos(theta)) **2
     &	           +float(m)+2.*cos(theta)) / (sin(theta)*sin(theta))
     &	           *d(indm))
     &	           +sqrt(float((l+m)*(l-m+1)))*(-4.+float(-2*m+1)
     &	           *cos(theta))/sin(theta)*d(indm+9)
     &	           +sqrt(float((l-m+1)*(l-m+2)*(l+m-1)*(l+m)))*d(indm+18)

              ka=ka+1
              ypp(ka)=-d(indm)*real(cfac)
              ytp(ka)=+d(indm)*aimag(cfac)
	      ypp_ph(ka)=+float(m)*d(indm)*aimag(cfac)
	      ytp_ph(ka)=+float(m)*d(indm)*real(cfac)
	      ypp_th(ka)=-dplm2m_th*real(cfac)
	      ytp_th(ka)=+dplm2m_th*aimag(cfac)
c
	ypp_th2(ka)=-dplm2m_th2*real(cfac)
	ytp_th2(ka)=+dplm2m_th2*aimag(cfac)
	ypp_thph(ka)=+float(m)*dplm2m_th*aimag(cfac)
	ytp_thph(ka)=+float(m)*dplm2m_th*real(cfac)
c
              ka=ka+1
              ypp(ka)=-d(indm)*aimag(cfac)
              ytp(ka)=-d(indm)*real(cfac)
	      ypp_ph(ka)=-float(m)*d(indm)*real(cfac)
	      ytp_ph(ka)=+float(m)*d(indm)*aimag(cfac)
	      ypp_th(ka)=-dplm2m_th*aimag(cfac)
	      ytp_th(ka)=-dplm2m_th*real(cfac)
c
	ypp_th2(ka)=-dplm2m_th2*aimag(cfac)
	ytp_th2(ka)=-dplm2m_th2*real(cfac)
	ypp_thph(ka)=-float(m)*dplm2m_th*real(cfac)
	ytp_thph(ka)=+float(m)*dplm2m_th*aimag(cfac)
c
              dplm4p_th=-float(m)*(cos(theta)/sin(theta))*d(indp4) 
     &   	   +4.*d(indp4)/sin(theta)
     &             -sqrt(float((l+m)*(l-m+1)))*d(indp4-9)

	      dplm4p_th2=(((4.-float(m)*cos(theta))**2
     &	           +float(m)-4.*cos(theta))/(sin(theta)*sin(theta))
     &	           *d(indp4))
c     &	           +sqrt(float((l+m)*(l-m+1)))*(8.-float(2*m-1)
     &	           -sqrt(float((l+m)*(l-m+1)))*(8.-float(2*m-1)
     &	           *cos(theta))/sin(theta)*d(indp4-9)
     &	          +sqrt(float((l-m+1)*(l-m+2)*(l+m-1)*(l+m)))*d(indp4-18)

	      ypppp_th(ka4-1)=-dplm4p_th*real(cfac)
	      ytppp_th(ka4-1)=-dplm4p_th*aimag(cfac)
	      ypppp_th(ka4)=+dplm4p_th*aimag(cfac)
	      ytppp_th(ka4)=-dplm4p_th*real(cfac)
c
	ypppp_th2(ka4-1)=-dplm4p_th2*real(cfac)
	ytppp_th2(ka4-1)=-dplm4p_th2*aimag(cfac)
	ypppp_th2(ka4)=+dplm4p_th2*aimag(cfac)
	ytppp_th2(ka4)=-dplm4p_th2*real(cfac)
	ypppp_thph(ka4-1)=+float(m)*dplm4p_th*aimag(cfac)
	ytppp_thph(ka4-1)=-float(m)*dplm4p_th*real(cfac)
	ypppp_thph(ka4)=+float(m)*dplm4p_th*real(cfac)
	ytppp_thph(ka4)=+float(m)*dplm4p_th*aimag(cfac)
c
              dplm4m_th=-float(m)*(cos(theta)/sin(theta))*d(indm4) 
     &             -4.*d(indm4)/sin(theta)
     &             +sqrt(float((l+m)*(l-m+1)))*d(indm4+9)

	      dplm4m_th2=(( (-4.-float(m)*cos(theta)) **2
     &	           +float(m)+4.*cos(theta)) / (sin(theta)*sin(theta))
     &	           *d(indm4))
     &	           +sqrt(float((l+m)*(l-m+1)))*(-8.+float(-2*m+1)
     &	           *cos(theta))/sin(theta)*d(indm4+9)
     &	           +sqrt(float((l-m+1)*(l-m+2)*(l+m-1)*(l+m)))*d(indm4+18)

              ka4=ka4+1
              ypppp(ka4)=-d(indm4)*real(cfac)
              ytppp(ka4)=+d(indm4)*aimag(cfac)
	      ypppp_ph(ka4)=+float(m)*d(indm4)*aimag(cfac)
	      ytppp_ph(ka4)=-float(m)*d(indm4)*real(cfac)
	      ypppp_th(ka4)=-dplm4m_th*real(cfac)
	      ytppp_th(ka4)=+dplm4m_th*aimag(cfac)
c
	ypppp_th2(ka4)=-dplm4m_th2*real(cfac)
	ytppp_th2(ka4)=+dplm4m_th2*aimag(cfac)
	ypppp_thph(ka4)=+float(m)*dplm4m_th*aimag(cfac)
	ytppp_thph(ka4)=+float(m)*dplm4m_th*real(cfac)
c
              ka4=ka4+1
              ypppp(ka4)=-d(indm4)*aimag(cfac)
              ytppp(ka4)=-d(indm4)*real(cfac)
	      ypppp_ph(ka4)=-float(m)*d(indm4)*real(cfac)
	      ytppp_ph(ka4)=+float(m)*d(indm4)*aimag(cfac)
	      ypppp_th(ka4)=-dplm4m_th*aimag(cfac)
	      ytppp_th(ka4)=-dplm4m_th*real(cfac)
c
	ypppp_th2(ka4)=-dplm4m_th2*aimag(cfac)
	ytppp_th2(ka4)=-dplm4m_th2*real(cfac)
	ypppp_thph(ka4)=-float(m)*dplm4m_th*real(cfac)
	ytppp_thph(ka4)=+float(m)*dplm4m_th*aimag(cfac)
c
            endif
            ind=ind+9
            indp=indp+9
            indm=indm-9
            indp4=indp4+9
            indm4=indm4-9
            cfac=cfac*dfac
          enddo
          dplm0_th=sqrt(float(l*(l+1)))*d(ind0+9)
          yth(k0)=dplm0_th*real(cfac0)
          dplm2p_th=-2.*d(indp0)/sin(theta)
     &          +sqrt(float(l*(l+1)))*d(indp0+9)

	      dplm2p_th2=((4.
     &	           +2.*cos(theta))/(sin(theta)*sin(theta))
     &	           *d(indp0))
     &	           +sqrt(float((l)*(l+1)))*(-4.+
c     &	           -sqrt(float((l)*(l+1)))*(-4.+
     &	           cos(theta))/sin(theta)*d(indp0+9)
     &	           +sqrt(float((l+1)*(l+2)*(l-1)*l))*d(indp0+18)

	  ypp_th(ka0)=-dplm2p_th*real(cfac0)
	  ytp_th(ka0)=-dplm2p_th*aimag(cfac0)
	  ypp_th(ka0+1)=+dplm2p_th*aimag(cfac0)
	  ytp_th(ka0+1)=-dplm2p_th*real(cfac0)
c
	ypp_th2(ka0)=-dplm2p_th2*real(cfac0)
	ytp_th2(ka0)=-dplm2p_th2*aimag(cfac0)
c	ypp_thph(ka0)=+float(m)*dplm2p_th*aimag(cfac0)
c	ytp_thph(ka0)=-float(m)*dplm2p_th*real(cfac0)
	ypp_thph(ka0)=0.
	ytp_thph(ka0)=0.
	ypp_th2(ka0+1)=+dplm2p_th2*aimag(cfac0)
	ytp_th2(ka0+1)=-dplm2p_th2*real(cfac0)
c	ypp_thph(ka0+1)=+float(m)*dplm2p_th*real(cfac0)
c	ytp_thph(ka0+1)=+float(m)*dplm2p_th*aimag(cfac0)
	ypp_thph(ka0+1)=0.
	ytp_thph(ka0+1)=0.
c
          dplm4p_th=-4.*d(indp40)/sin(theta)
     &          +sqrt(float(l*(l+1)))*d(indp40+9)

	      dplm4p_th2=((16.
     &	           +4.*cos(theta)) /(sin(theta)*sin(theta))
     &	           *d(indp40))
     &	           +sqrt(float((l)*(l+1)))*(-8.+
c     &	           -sqrt(float((l)*(l+1)))*(-8.+
     &	           cos(theta))/sin(theta) *d(indp40+9)
     &	          +sqrt(float((l+1)*(l+2)*(l-1)*l))*d(indp40+18)

	  ypppp_th(ka40)=-dplm4p_th*real(cfac0)
	  ytppp_th(ka40)=-dplm4p_th*aimag(cfac0)
	  ypppp_th(ka40+1)=+dplm4p_th*aimag(cfac0)
	  ytppp_th(ka40+1)=-dplm4p_th*real(cfac0)
c
      ypppp_th2(ka40)=-dplm4p_th2*real(cfac0)
      ytppp_th2(ka40)=-dplm4p_th2*aimag(cfac0)
c      ypppp_thph(ka40)=+float(m)*dplm4p_th*aimag(cfac0)
c      ytppp_thph(ka40)=-float(m)*dplm4p_th*real(cfac0)
      ypppp_thph(ka40)=0.
      ytppp_thph(ka40)=0.
      ypppp_th2(ka40+1)=+dplm4p_th2*aimag(cfac0)
      ytppp_th2(ka40+1)=-dplm4p_th2*real(cfac0)
c      ypppp_thph(ka40+1)=+float(m)*dplm4p_th*real(cfac0)
c      ytppp_thph(ka40+1)=-float(m)*dplm4p_th*aimag(cfac0)
      ypppp_thph(ka40+1)=0.
      ytppp_thph(ka40+1)=0.
c
        endif
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine vecmatvec(u,a,v,m,n,uav)
      dimension a(m,n)
      dimension u(m)
      dimension v(n)
      real*4 uav,av
      uav=0.
      do i=1,m
         av=0.
         do j=1,n
            av=av+a(i,j)*v(j)
         enddo
         uav=uav+u(i)*av
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine vecvecmatvecvec(u1,u2,a,u3,u4,n1,n2,n3,n4,uuauu)
      dimension a(n1,n2,n3,n4)
      dimension u1(n1),u2(n2),u3(n3),u4(n4)
      uuauu=0.
      do i=1,n1
         uauu=0.
         do j=1,n2
            auu=0.
            do k=1,n3
               au=0.
               do l=1,n4
                  au=au+a(i,j,k,l)*u4(l)
               enddo
               auu=auu+au*u3(k)
            enddo
            uauu=uauu+auu*u2(j)
         enddo
         uuauu=uuauu+uauu*u1(i)
      enddo
      return
      end
c-------------------------------------------------------------------
      function sdot(n,a,inca,b,incb)
      dimension a(inca,*),b(incb,*)
      sum=0.
      do i=1,n
        sum=sum+a(1,i)*b(1,i)
      enddo
      sdot=sum
      return
      end 
c-------------------------------------------------------------------
      subroutine rotvcani(
     1     c0,c2,c4,lmax0,lmax2,lmax4
     1     ,alph,beta,gama,d,vec1,vec2
     1     ,c0r,c2r,c4r)
c-------------rotates generalized spherical harmonic coefficients.
c-------------coefficients vector is real in calling program and complex here,
c-------------wich means that the first and second entries of input vector
c-------------are the real and imaginary parts of the first complex vector entry,
c-------------the third and fourth input entries form the second complex entry
c-------------and so on.

c Dimensions of input and output vectors. The
c dimension statement has no effect unless array bound
c checking is in effect, but is included as an aide memoire
      dimension 
     1  c0((lmax0+1)**2)
     1 ,c0r((lmax0+1)**2)
      complex
     1  c2((lmax2-1)*(lmax2+3))
     1 ,c4((lmax4-3)*(lmax4+5))
     1 ,c2r((lmax2-1)*(lmax2+3))
     1 ,c4r((lmax4-3)*(lmax4+5))

c workspace: lmax = max(lmax0,lmax2,lmax4)

      double precision d(*)         ! dimension at least (2*lmax+1)**2 d.p.
      complex vec1(0:*),vec2(0:*)   ! dimension at least (2*lmax+1) complex
      double precision dbeta
      complex ealph,egama

      lmax=max(lmax0,lmax2,lmax4)
      dbeta=beta
      do l=0,lmax
        call rotmx2(l,l,dbeta,d,2*l+1,2*l+1)
        ealph=cexp(cmplx(0.,alph))
        egama=cexp(cmplx(0.,gama))
     
        if(l.le.lmax0) then
        
          k0=l**2
          k=k0
          sgn=1
          do m=0,l
            if (m.eq.0) then
              k=k+1
              vec1(l+m)=cmplx(c0(k),0.)
            else
              k=k+1
              vec1(l+m)=cmplx(.5*c0(k),-.5*c0(k+1))
              k=k+1
              vec1(l-m)=sgn*conjg(vec1(l+m))
            endif
            sgn=-sgn
          enddo
          call crotl(d,l,ealph,egama,vec1,vec2)
          k=k0
          do m=0,l
            if(m.eq.0) then
              k=k+1
              c0r(k)=real(vec2(l+m))
            else
              k=k+1
              c0r(k)=2.*real(vec2(l+m))
              k=k+1
              c0r(k)=-2.*aimag(vec2(l+m))
            endif
          enddo
        endif

        if(l.ge.2.and.l.le.lmax2) then
          k0=l**2-4
          k=k0
          do m=0,l
            k=k+1
            vec1(l+m)=c2(k)
            if(m.ne.0) then
              k=k+1
              vec1(l-m)=c2(k)
            endif
          enddo
          call crotl(d,l,ealph,egama,vec1,vec2)
          k=k0
          do m=0,l
            k=k+1
            c2r(k)=vec2(l+m)
            if(m.ne.0) then
              k=k+1
              c2r(k)=vec2(l-m)
            endif
          enddo
        endif


        if(l.ge.4.and.l.le.lmax4) then
          k0=l**2-16
          k=k0
          do m=0,l
            k=k+1
            vec1(l+m)=c4(k)
            if(m.ne.0) then
              k=k+1
              vec1(l-m)=c4(k)
            endif
          enddo
          call crotl(d,l,ealph,egama,vec1,vec2)
          k=k0
          do m=0,l
            k=k+1
            c4r(k)=vec2(l+m)
            if(m.ne.0) then
              k=k+1
              c4r(k)=vec2(l-m)
            endif
          enddo
        endif

      enddo               ! end of do loop over l
      return
      end
c-------------------------------------------------------------------
c perform spherical harmonic rotation assuming
c that the matrix d(beta) and the values
c ealph=exp( i alpha) egama=exp( i gama) have
c already been calculated. This is called by
c rotvcani()

      subroutine crotl(d,l,ealph,egama,vec1,vec2)
      complex ealph,egama,vec1(-l:l),vec2(-l:l)
      double precision d(-l:l,-l:l)
      complex fct
      fct=ealph
      do m=1,l
        vec1(m)=vec1(m)*fct
        vec1(-m)=vec1(-m)*conjg(fct)
        fct=fct*ealph
      enddo
      do m=-l,l
        vec2(m)=(0.,0.)
        do mp=-l,l
          vec2(m)=vec2(m)+d(m,mp)*vec1(mp)
        enddo
      enddo
      fct=egama
      do m=1,l
        vec2(m)=vec2(m)*fct
        vec2(-m)=vec2(-m)*conjg(fct)
        fct=fct*egama
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine apkyd(lmax,phi,theta,y,yph,yth,yth2,ythph,wk1,wk2,wk3)
      dimension y(*),yph(*),yth(*),yth2(*),ythph(*),wk1(*),wk2(*),wk3(*)
      complex temp,fac,dfac,tdph,tdth,tdth2,tdthph
      cosec=1./sin(theta) 
      cot=cos(theta)*cosec
      cosec2=cosec*cosec
      ind=0 
      lm1=lmax+1
      do il1=1,lm1 
        l=il1-1 
        fl3=float(l*(l+1))
        call legndr(theta,l,l,wk1,wk2,wk3)
        fac=(1.,0.) 
        dfac=cexp(cmplx(0.,phi))
        do im=1,il1
          xm=float(im-1)
          xm2=xm*xm 
          temp=fac*cmplx(wk1(im),0.)
          tdph=temp*cmplx(0.,xm)
          tdth=fac*cmplx(wk2(im),0.)
          tdth2=fac*cmplx(-cot*wk2(im)-(fl3-xm2*cosec2)*wk1(im),0.) 
          tdthph=tdth*cmplx(0.,xm)
          ind=ind+1 
          y(ind)=real(temp) 
          yph(ind)=real(tdph) 
          yth(ind)=real(tdth) 
          yth2(ind)=real(tdth2) 
          ythph(ind)=real(tdthph) 
          if(im.ne.1) then 
            ind=ind+1 
            y(ind)=aimag(temp)
            yph(ind)=aimag(tdph)
            yth(ind)=aimag(tdth)
            yth2(ind)=aimag(tdth2)
            ythph(ind)=aimag(tdthph)
          endif
          fac=fac*dfac
        enddo
      enddo
      return
      end 
