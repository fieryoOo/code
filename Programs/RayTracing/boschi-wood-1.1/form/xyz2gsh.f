c---find expansion of xyz function over generalized (N=2 or N=4) or scalar spherical
c---harmonics (N=0) by solution of least squares problem. needed to express
c---maps of azimuthal anisotropy in terms of generalized s.h. coefficients.
	parameter(lmx=29)
	parameter(nunk=(lmx-1)*(2*lmx+6))
	parameter(ilata=(nunk*(nunk+1))/2)
	parameter(nx=180,ny=90,idsize=nx*ny)  !increase idsize if spacing < 1 deg
	character namein*80,nameout*80,namein2*80
	dimension y((lmx+1)**2)
     &         ,ypp((lmx-1)*(2*lmx+6))
     &         ,ytp((lmx-1)*(2*lmx+6))
     &         ,ypppp((lmx-3)*(2*lmx+10))
     &         ,ytppp((lmx-3)*(2*lmx+10))
	double precision wkspc(9*(2*lmx+1))
	dimension dat(2*idsize),rowa(nunk),ata(ilata),atd(nunk),x(nunk)
     &         ,gstore(ilata),ytemp(nunk),a(2*idsize,nunk)
	real ataout(ilata),atdout(nunk),xmasked(nunk)
	integer maskv(nunk),maskata(ilata)
	integer*4 dampcoef(nunk)

	print*,"term 0, 2, or 4?"
	read*,ialpha
	if(ialpha.eq.0)then
	   print*,"input model?"
	   read*,namein
	   open(1,file=namein,status='old')
	   ncoef=(lmx+1)**2
	elseif(ialpha.eq.2.or.ialpha.eq.4)then
	   if(ialpha.eq.2)ncoef=(lmx-1)*(2*lmx+6)
	   if(ialpha.eq.4)ncoef=(lmx-3)*(2*lmx+10)
	   print*,"a1 or a3 input model?"
	   read*,namein
	   print*,"a2 or a4 input model?"
	   read*,namein2
	   open(1,file=namein,status='old')
	   open(2,file=namein2,status='old')
	else
	   stop "non-existing term"
	endif
	print*,"output file?"
	read*,nameout
	print*,"damping factor?"
	read*,damp

	do i=1,ilata
	ata(i)=0.
	enddo

c----------------read in all the gridpoints and increment matrices
	i=1
1	continue
	if(ialpha.eq.0)then
c-------------------------------scalar SH
	   read(1,*,end=20)xlon,xlat,dat(i)
	   call ylmv4(xlat,xlon,lmx,y,ypp,ytp,ypppp,ytppp,wkspc)
	   do k=1,ncoef
	      rowa(k)=y(k)
	      a(i,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)
	elseif(ialpha.eq.2)then
c-------------------------------generalized SH for second order tensor
	   read(1,*,end=20)xlon,xlat,dat(i)
	   call ylmv4(xlat,xlon,lmx,y,ypp,ytp,ypppp,ytppp,wkspc)
	   do k=1,ncoef
	      rowa(k)=ypp(k)
	      a(i,k)=rowa(k)
	   enddo
	   dat(i)=-dat(i)
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)
	   i=i+1
	   read(2,*)xlon1,xlat1,dat(i)
	   if(xlon1.ne.xlon.or.xlat1.ne.xlat)stop "input files not compatible"
	   do k=1,ncoef
	      rowa(k)=ytp(k)
	      a(i,k)=rowa(k)
	   enddo
	   dat(i)=-dat(i)
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)
	elseif(ialpha.eq.4)then
c-------------------------------generalized SH for fourth order tensor
	   read(1,*,end=20)xlon,xlat,dat(i)
	   call ylmv4(xlat,xlon,lmx,y,ypp,ytp,ypppp,ytppp,wkspc)
	   do k=1,ncoef
	      rowa(k)=ypppp(k)
	      a(i,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)
	   i=i+1
	   read(2,*)xlon1,xlat1,dat(i)
	   if(xlon1.ne.xlon.or.xlat1.ne.xlat)stop "input files not compatible"
	   do k=1,ncoef
	      rowa(k)=ytppp(k)
	      a(i,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)
	else
	   stop "error in variable ialpha"
	endif

	if(mod(i,1000).eq.0)print*,i," gridpoints read"
c--go back to read the next datum
	i=i+1
	goto1

20	print*,i-1,' gridpoints read total'  
	ndat=i-1
	close(1)
	if(ialpha.ne.0)close(2)

c---------------------------------------regularization
	k=0
	il0=ialpha
c----------------------MUST BE CHANGED FOR GENERALIZED SH PARAMETERIZATION
	do il=il0,lmx
	   immax=2*il+1
	   do im=1,immax
	      k=k+1
	      dampcoef(k)=float(il*(il+1))
	   enddo
	enddo
	k=0
	icount=0
	do ix=1,ncoef
	   do iy=1,ix
	      k=k+1
	      if(iy.eq.ix)then
	         icount=icount+1
c	         ata(k)=ata(k)+damp*dampcoef(icount)
c------------NON-L-DEPENDENT DAMPING FOR THE TIME BEING
	         ata(k)=ata(k)+damp
	      endif
	   enddo
	enddo

cc--mask coefficients associated with basis functions that cannot be constrained
c SHOULD NOT BE NEEDED NOW
c	print*,"mask ata matrix and atd vector..."
c	call atamask(ata,atd,ataout,atdout,ilata,nunk,maskata,maskv,nonzero)

c--find least squares solution x via single-processor cholesky factorization of ATA
	print*,"cholesky factorization..."
	call choles(ata,gstore,atd,ytemp,x,ncoef,nono)
c	call choles(ataout,gstore,atdout,ytemp,xmasked,nonzero,nono)

	if(nono.ne.0)then
	print*,"cholesky factorization failed ",nono
c	do j=0,nonzero-1
c	istart=j*(j+1)/2+1
c	j1=j+1
c	istop=j1*(j1+1)/2
cc	write(*,*)j+1,(ata(i),i=istart,istop)
cc	pause
c	enddo
	stop
	endif

c--unmask solution coefficients that were eliminated by atamask
c	print*,"unmask solution..."
c	k=0
c	do i=1,nunk
c	   if(maskv(i).eq.1)then
c	      k=k+1
c	      x(i)=xmasked(k)
c	   elseif(maskv(i).eq.0)then
c	      x(i)=0.
c	   else
c	      print*,"illegal maskv value for coeff. ",i
c	      stop
c	   endif
c	enddo

c--see how well the least squares solution fits the gridpoints
	print*,"computing variance reduction..."
	truerms=0.
	realchisq=0.
	denom=0.
	do irow=1,ndat
	   error=0.
	   syndat=0.
	   do icl=1,ncoef
	      syndat=syndat+a(irow,icl)*x(icl)
	   enddo
	   error=(syndat-dat(irow))**2
	   realchisq=error+realchisq
	   truerms=truerms+abs(syndat-dat(irow))
	   denom=denom+(dat(irow)*dat(irow))
	enddo
	truerms=truerms/float(ndat)
	varred=realchisq/denom
	varred=1.-varred
	print*,'NUMBER OF DATA:',ndat
	print*,'CHI-SQUARE:',realchisq
	print*,'RMS MISFIT OBTAINED:',truerms
	print*,'VARIANCE REDUCTION:',varred

	open(2,file=nameout)
	write(2,*)lmx
	write(2,*)(x(i),i=1,ncoef)
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
c------------------------single-processor Cholesky factorization
      subroutine choles(a,g,b,y,x,n,nono)

      implicit real*4 (a-h, o-z)
      real*4 a(1),g(1)
      real*4 b(1),y(1),x(1)
c
c        a= row-wise p.d. symm. system  n*(n+1)/2
c        g= cholesky storage
c        b= r.h.s. vector               n
c        y= temp. vector
c        x= answer vector
c        n= system dimension
c        nono .gt. 0 is the level at which p.d. failed
c
c        (a,g) and (b,y,x) may be equivalenced.
c
c----------------------------------------------------------
c-----first compute cholesky decomposition

      nono=0
      
      if(a(1).le.0.) then
      nono=1
      return
      endif
      
      g(1)=sqrt(a(1))
      y(1)=b(1)/g(1)

      do 400 i=2,n
      
      kz=(i*(i-1))/2
      g(kz+1)=a(kz+1)/g(1)
      sg=g(kz+1)**2
      y(i)=b(i)-g(kz+1)*y(1)
      
      if(i.gt.2) then
      
      jmax=i-1
      
      do 200 j=2,jmax
      
      gkz=a(kz+j)
      kj=(j*(j-1))/2
      kmax=j-1
      
      do 100 k=1,kmax
      gkz=gkz-g(kz+k)*g(kj+k)
100   continue

      g(kz+j)=gkz/g(kj+j)
      y(i)=y(i)-g(kz+j)*y(j)
      sg=sg+g(kz+j)**2
      
200   continue

      endif
      
      gkz=a(kz+i)-sg
      
      if(gkz.le.0.) then
      nono=i
      return
      endif
      
      g(kz+i)=sqrt(gkz)
      y(i)=y(i)/g(kz+i)
      
400   continue

      kz=(n*(n-1))/2
      x(n)=y(n)/g(kz+n)
      if(n.le.1) return

c-----
c     compute solution for particular rhs
      
      do 600 k=2,n
      
      i=n+1-k
      x(i)=y(i)
      jmin=i+1
      
      do 500 j=jmin,n
      kj=(j*(j-1))/2
      x(i)=x(i)-g(kj+i)*x(j)
500   continue

      kz=(i*(i+1))/2
      x(i)=x(i)/g(kz)
      
600   continue

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
c---------------------------------------------------------------------------
	subroutine atamask(a,v,aout,vout,ilata,nunk,maskata,maskv,kvout)
c---removes from ata matrix (a) to be inverted columns/rows that are
c---systematically 0. corresponding atd (v) entries are also removed.
	real a(ilata),aout(ilata),v(nunk),vout(nunk)
	integer maskv(nunk),maskata(ilata)
	if(ilata.ne.nunk*(nunk+1)/2)stop "wrong dimensions"

	kvout=0
	do i=1,ilata
	   maskata(i)=1
	enddo
	do i=1,nunk
	   index=i*(i+1)/2
	   if(a(index).eq.0.)then
	      print*,i," diagonal entry is 0"
	      maskv(i)=0
	      do j=1,i
	         index=i*(i-1)/2+j

cTEST
c	print*,i,j,index,a(index)

	         if(a(index).ne.0.)stop "this does not make sense"
	         maskata(index)=0
	      enddo
	      do j=i+1,nunk
	         index=j*(j-1)/2+i

cTEST
c	print*,i,j,index,a(index)

	         if(a(index).ne.0.)stop "this does not make sense"
	         maskata(index)=0
	      enddo

c	pause

	   else
	      maskv(i)=1
	      kvout=kvout+1
	      vout(kvout)=v(i)
	   endif
	enddo
	kaout=0
	do i=1,ilata
	   if(maskata(i).ne.0)then
	      kaout=kaout+1
	      aout(kaout)=a(i)

cTEST
c	if(a(i).eq.0.)stop

	   endif
	enddo
	return
	end
c---------------------------------------------------------------------------
	subroutine contribution_ata(rowa,dat,ncoef,nunk,ilata,ata,atd)
c----given a row of A (rowa) and the corresponding datum (dat), add their
c----contribution to the matrix A^tA and the vector A^td
	dimension rowa(nunk),ata(ilata),atd(nunk)
	k1=0
	do i1=1,ncoef
	   atd(i1)=atd(i1)+rowa(i1)*dat
	   do j1=1,i1
	      k1=k1+1 
	      ata(k1)=ata(k1)+rowa(i1)*rowa(j1)
	   enddo
	enddo
        return
        end
c---------------------------------------------------------------------------
