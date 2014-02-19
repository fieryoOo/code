	parameter(nx=180,ny=90,idsize=4*nx*ny)
	parameter(lmx=55)
	character namein*80,nameout*80,chlmax*3,chialpha*1
	dimension y((lmx+1)**2)
     &         ,ypp((lmx-1)*(2*lmx+6))
     &         ,ytp((lmx-1)*(2*lmx+6))
     &         ,ypppp((lmx-3)*(2*lmx+10))
     &         ,ytppp((lmx-3)*(2*lmx+10))
	double precision wkspc(9*(2*lmx+1))
	dimension cmod((LMX+1)**2)

	print*,"what input model?"
	read*,namein
	print*,"what alpha (0, 1, 2, 3, 4)?"
	read*,ialpha
	if(ialpha.lt.0.or.ialpha.gt.4)stop "non-existing term"
	print*,"what increment for grid?"
	read*,xincr
c--------------------following line useful e.g. if model given as
c--------------------relative perturbation and we need percent
	print*,"scale model by? (1 if no scale)"
	read*,scale

	open(1,file=namein,status='old')
	read(1,*)lmax
	print*,"maximum l=",lmax
	print*,"interrupt expansion at l=?"
	read*,lmaxuser
	if(ialpha.eq.0)then
	   ncoef=(lmax+1)**2
	elseif(ialpha.le.2)then
	   ncoef=(lmax-1)*(2*lmax+6)
	else
	   ncoef=(lmax-3)*(2*lmax+10)
	endif
	read(1,*)(cmod(i),i=1,ncoef)
	close(1)
	do k=1,80
	   if(namein(k:k).eq.' ')goto2
	enddo
2	continue
	lmax=lmaxuser		! user can filter-lapo 2004-04-03
	write(chlmax,'(i3.3)')lmax
	write(chialpha,"(i1.1)")ialpha
	nameout=namein(1:k-1)//".eps"//chialpha//"."//chlmax//'.xyz'
	open(2,file=nameout)

	igrid=0
	if(lmax.gt.lmx)stop "l too big"
	do xlat=-89.,89,xincr
	do xlon=1.,359.,xincr
	igrid=igrid+1
	if(mod(igrid,1000).eq.0)print*,igrid," grid points done"
	call ylmv4(xlat,xlon,lmx,y,ypp,ytp,ypppp,ytppp,wkspc)
	z=0.
	if(ialpha.eq.0)then
c-------------------------------scalar SH
	   ncoef=(lmax+1)**2
	   do k=1,ncoef
	      z=z+y(k)*cmod(k)
	   enddo
	elseif(ialpha.eq.1)then
c-------------------------------generalized SH for second order tensor
	   ncoef=(lmax-1)*(2*lmax+6)
	   do k=1,ncoef
	      z=z-ypp(k)*cmod(k)
	   enddo
	elseif(ialpha.eq.2)then
	   ncoef=(lmax-1)*(2*lmax+6)
	   do k=1,ncoef
	      z=z-ytp(k)*cmod(k)
	   enddo
	elseif(ialpha.eq.3)then
c-------------------------------generalized SH for fourth order tensor
c	   ncoef=(lmax-3)*(2*lmax+10)
           do k=1,ncoef
	      z=z+ypppp(k)*cmod(k)
	   enddo
	elseif(ialpha.eq.4)then
c	   ncoef=(lmax-3)*(2*lmax+10)
	   do k=1,ncoef
	      z=z+ytppp(k)*cmod(k)
	   enddo
	else
	   stop "error in variable ialpha"
	endif
	write(2,*)xlon,xlat,z*scale
	enddo
	enddo

	close(2)
	end

c---------------------------------------------------------------------------
c        subroutine ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,d)
c
c  evaluates the vectors ypp,ytp, of length lenv=(lmax-1)*(2*lmax+6)
c  and vectors ypppp,ytppp, of length lenv4=(lmax-3)*(2*lmax+10).
c
c  ypp,ytp represent contributions to the symmetric trace-free
c  second rank tensor c.
c
c  ypppp,ytppp represent contributions to the completely symmetric trace-free
c  fourth rank tensor e.
c
c    In spherical coordinates ('t'=theta=colatitiude, 'p'=phi=longitude)
c
c    c  =  ( c_tt  c_tp )
c          ( c_pt  c_pp )
c
c    and  c_pp   = -c_tt   = ypp . coefs
c         c_tp    = c_pt   = ytp . coefs
c
c     Similarly for e:
c     e_pppp =  e_tttt = -e_pptt   = ypppp . coefs
c     e_tppp = -e_tttp             = ytppp . coefs
c       
c 
c     Scalar spherical harmonics are also calculated and
c     placed in y(1) -- y( (lmax+1)**2 )
c------------------------------------------------------------------
c I believe scalar harmonics here are orthogonal but NOT
c orthonormal: a factor sqrt(2) is missing from harmonics with
c nonzero m. see Dahlen and Tromp eq. B.72
c------------------------------------------------------------------
c
c     The companion routine ylmavv4() (q.v.) calculates the
c     contribution to the average value of  c_ll (= l.c.l) and e_llll
c     along the minor arc and the complete great circle,
c     where 'l' represents the tangent to the path:
c         
c           l_t = -cos(az)      l_p = sin(az)
c 
c     and az is the local azimuth of the path.
c
c     Thus:
c            c_ll   =  - c_pp * cos(2*az)  - c_tp * sin(2*az)
c     and  
c            e_llll =  e_pppp * cos(4*az) + e_tppp * sin(4*az)
c----------------------------------------------------------------
c     'd' is d.p. workspace. (notice that it's bigger than in ylmv)
c----------------------------------------------------------------
      subroutine ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,d)
      dimension y((lmax+1)**2)
     1         ,ypp((lmax-1)*(2*lmax+6))
     1         ,ytp((lmax-1)*(2*lmax+6))
     1         ,ypppp((lmax-3)*(2*lmax+10))
     1         ,ytppp((lmax-3)*(2*lmax+10))
      double precision theta,d(9*(2*lmax+1))
      complex cfac,dfac
      data radian/57.2957795/,rr4pi/0.28209479/
      theta=(90.-xlat)/radian
      dfac=cexp(cmplx(0.,xlon/radian))
      k=0
      ka=0
      ka4=0
      do l=0,lmax
        if(l.lt.2) then
          call rotmx2(0,l,theta,d,1,2*l+1)
          ind=l
          cfac=rr4pi*sqrt(2.*l+1)
          do m=0,l
            k=k+1
            ind=ind+1
            y(k)=d(ind)*real(cfac)
            if(m.ne.0) then
              k=k+1
              y(k)=d(ind)*aimag(cfac)
            endif
            cfac=cfac*dfac
          enddo
        else if(l.lt.4) then
          call rotmx2(2,l,theta,d,5,2*l+1)
          ind=5*l+3
          indp=ind+2
          indm=indp
          cfac=rr4pi*sqrt(2.*l+1)
          do m=0,l
            k=k+1
            y(k)=d(ind)*real(cfac)
            ka=ka+1
            ypp(ka)=-d(indp)*real(cfac)
            ytp(ka)=-d(indp)*aimag(cfac)
            ka=ka+1
            ypp(ka)=+d(indp)*aimag(cfac)
            ytp(ka)=-d(indp)*real(cfac)
            if(m.ne.0) then
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*real(cfac)
              ytp(ka)=+d(indm)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*aimag(cfac)
              ytp(ka)=-d(indm)*real(cfac)
            endif
            ind=ind+5
            indp=indp+5
            indm=indm-5
            cfac=cfac*dfac
          enddo
        else
          call rotmx2(4,l,theta,d,9,2*l+1)
          ind=9*l+5
          indp=ind+2
          indm=indp
          indp4=ind+4
          indm4=indp4
          cfac=rr4pi*sqrt(2.*l+1)
          do m=0,l
            k=k+1
            y(k)=d(ind)*real(cfac)
            ka=ka+1
            ypp(ka)=-d(indp)*real(cfac)
            ytp(ka)=-d(indp)*aimag(cfac)
            ka=ka+1
            ypp(ka)=+d(indp)*aimag(cfac)
            ytp(ka)=-d(indp)*real(cfac)
            ka4=ka4+1
            ypppp(ka4)=-d(indp4)*real(cfac)
            ytppp(ka4)=-d(indp4)*aimag(cfac)
            ka4=ka4+1
            ypppp(ka4)=+d(indp4)*aimag(cfac)
            ytppp(ka4)=-d(indp4)*real(cfac)
            if(m.ne.0) then
              k=k+1
              y(k)=d(ind)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*real(cfac)
              ytp(ka)=+d(indm)*aimag(cfac)
              ka=ka+1
              ypp(ka)=-d(indm)*aimag(cfac)
              ytp(ka)=-d(indm)*real(cfac)
              ka4=ka4+1
              ypppp(ka4)=-d(indm4)*real(cfac)
              ytppp(ka4)=+d(indm4)*aimag(cfac)
              ka4=ka4+1
              ypppp(ka4)=-d(indm4)*aimag(cfac)
              ytppp(ka4)=-d(indm4)*real(cfac)
            endif
            ind=ind+9
            indp=indp+9
            indm=indm-9
            indp4=indp4+9
            indm4=indm4-9
            cfac=cfac*dfac
          enddo
        endif
      enddo
      return
      end


c---------------------------------------------------------------------------
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
