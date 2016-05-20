c  transferp2z.f

c  This program calculates transfer function from pressure noise to vertical displacement
c  Intent is to remove gravity wave noise from vertical to enhance signal, so input will be
c  signal free noise samples on both components.  Pick a day with no earthquake interference
c  and select multiple signal-free windows.  Can also detect phase shifts between DPG and Z.

c  Pipe in data   transferp2z  < transferp2zinp
 
      parameter (maxnfreq=400, maxpts = 1000000, maxevnts = 100)

      real*4 freq(maxnfreq)
      real*4 timebeg(maxevnts,maxevnts)
      real*4 sumcrossre(maxnfreq),sumcrossim(maxnfreq)
      real*4 sumppower(maxnfreq)
      real*4 sumzpower(maxnfreq), tlength(maxevnts)
      real*4 begp(maxevnts), begz(maxevnts),deltp(maxevnts)
      real*4 deltz(maxevnts)
      real*4 winp(maxpts),winz(maxpts),gausst(maxpts),tdataz(maxpts)
      real*4 tdatap(maxpts)
      real*4 qreal(maxnfreq),qimag(maxnfreq),coher(maxnfreq)
      real*4 epsilon(maxnfreq)
      real*4 qamp(maxnfreq), qphse(maxnfreq),ampdev(maxnfreq)
      real*4 phsedev(maxnfreq)
      real*4 preal,pimag
      
      integer*4 nfreq, nwindows(maxevnts)
       integer*4 nptsp(maxevnts), nptsz(maxevnts)
      character*70 foutput, fnp(maxevnts), fnz(maxevnts)

      pi = 3.1415928
      convdeg = 3.1415928/180.
      circ = 6371.*3.1415928/180.
      twopi = 3.1415928*2.

      nfreq = 396
      do ifreq = 1, nfreq
         freq(ifreq) = (ifreq-1)*.001 + .005
	 sumcrossre(ifreq) = 0.0
	 sumcrossim(ifreq) = 0.0
	 sumppower(ifreq) = 0.0
	 sumzpower(ifreq) = 0.0
       enddo

c  foutput file for output of transfer function and coherence 
        read(*,'(a)') foutput
      open(11, file = foutput)
      
c  read list of files to be analyzed (events or time series)
c  ASSUMES VERTICAL FILES are DISPLACEMENT with RESPONSE REMOVED
c  nogauss = 0 means don't apply gaussian window - use record as is (For already windowed
c   signals)
c  nogauss = 1 means do apply gaussian window (For noise sample)
      read(*,*) nevents, nogauss
      nobs = 0
      do iev = 1, nevents
c  read in pressure file name, then associated vertical file name - both should have
c  same beginning time and same time increments
        read(*,'(a)') fnp(iev)
        read(*,'(a)') fnz(iev)
c  read in number of windows selected and timelength of window for seismogram
	read(*,*) nwindows(iev), tlength(iev)
c  timebeg is beginning of window relative to start of record
	do iwin = 1, nwindows(iev)
	  read(*,*) timebeg(iev,iwin)
	enddo
        call rsac1(fnp(iev),tdatap,nptsp(iev),begp(iev),
     1       deltp(iev),maxpts,nerr)
        call rsac1(fnz(iev),tdataz,nptsz(iev),begz(iev),
     1       deltz(iev),maxpts,nerr)
        if (begp(iev).ne.begz(iev)) then
	  write (*,*) iev,' beginning times not identical'
	endif
	if (deltp(iev).ne.deltz(iev)) write(*,*) iev,' delt not identical'
c  set up Gaussian window function - want 3 sigma at limits of time window
c  and amplitude = 1 at center
	npts = tlength(iev)/deltp(iev) +1
        if (nogauss.eq.1) then
          sigma = tlength(iev)/6.0
    	  sigma2 = 2.0*sigma*sigma
	  tmiddle = tlength(iev)/2.0
	  do ig = 1, npts
	    tg = (ig-1)*deltp(iev) - tmiddle
	    gausst(ig) = exp(-(tg*tg)/sigma2)
	  enddo
	endif
C  now apply gaussian to windows
        do iwin = 1, nwindows(iev)
	  nobs = nobs + 1
	  iwstart = timebeg(iev,iwin)/deltp(iev)
	  do iw = 1, npts
	    if (nogauss.eq.1) then
	      winp(iw) = tdatap(iw+iwstart)*gausst(iw)
	      winz(iw) = tdataz(iw+iwstart)*gausst(iw)
	    else
	      winp(iw) = tdatap(iw+iwstart)
	      winz(iw) = tdataz(iw+iwstart)
	    endif
	  enddo
c
c  get real and imaginary (phases and amplitudes) at desired frequencies
c
          do ifreq = 1, nfreq
            call frt2(winp, freq(ifreq), preal,
     1                      pimag,
     2                      npts,deltp(iev))
            call frt2(winz, freq(ifreq), zreal,
     1                      zimag,
     2                      npts,deltz(iev))
            crossre = preal*zreal + pimag*zimag
	    crossim = preal*zimag - pimag*zreal
	    ppower = preal*preal + pimag*pimag
	    zpower = zreal*zreal + zimag*zimag
	    sumcrossre(ifreq) = sumcrossre(ifreq) + crossre
	    sumcrossim(ifreq) = sumcrossim(ifreq) + crossim
	    sumppower(ifreq) = sumppower(ifreq) + ppower
	    sumzpower(ifreq) = sumzpower(ifreq) + zpower
          enddo  
        enddo
      enddo
c  end of input loop over events

c  calculate transfer function q and coherence as function of frequency - and translate into
c  amplitude and phase and associated uncertainty
      do ifreq = 1, nfreq	
        qreal(ifreq) = sumcrossre(ifreq)/sumppower(ifreq)
        qimag(ifreq) = sumcrossim(ifreq)/sumppower(ifreq)
        coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/
     1                  (sumppower(ifreq)*sumzpower(ifreq))
c  epsilon is fractional uncertainty in amplitude of transfer function (1 stddev)
c  and also approximately 1 stddev in phase (in radians) Box and Jenkins page 317
        epsilon(ifreq) = sqrt((1.0 - coher(ifreq))/
     1                   (2.0*coher(ifreq)*nobs))
	call phase(qreal(ifreq),qimag(ifreq),qphse(ifreq))
	if (qphse(ifreq).gt.0.5) qphse(ifreq)=qphse(ifreq)-1.0
	qamp(ifreq) = sqrt(qreal(ifreq)**2 + qimag(ifreq)**2)
	ampdev(ifreq) = qamp(ifreq)*epsilon(ifreq)
c  phase deviation and phase in cycles
	phsedev(ifreq) = epsilon(ifreq)/twopi
c  output 
        write(11,100) freq(ifreq),qamp(ifreq),ampdev(ifreq),
     1     qphse(ifreq),phsedev(ifreq),coher(ifreq), 
     2     qreal(ifreq), qimag(ifreq)
 100   format(f7.4,2(1x,e10.4),3(1x,f7.4),2(1x,e11.4))      
      enddo
      close(unit = 11)
      end

c---positive fourier transform
c 
      SUBROUTINE FRT(UP,FR,ARZ,PRZ,NALL,DELT)
      dimension UP(1000000)
      DIMENSION W(1000000)
      THETA=6.283185*FR*DELT
      C=COS(THETA)*2.0
      NR1=1
      NR2=NALL
      NDR1=NR2-1
      W(1)=UP(NR2)
      W(2)=C*W(1)+UP(NDR1)
      NDATA=NR2-NR1+1
      NTR1=NDATA-1
      NTR2=NDATA-2
      DO 97 I=3,NDATA
      I1=I-1
      I2=I-2
      NDRI=NR2-I+1
      W(I)=C*W(I1)-W(I2)+UP(NDRI)
97    CONTINUE
      ZRE=(W(NDATA)-W(NTR2)+UP(NR1))*DELT/2.0
      ZIM=W(NTR1)*SIN(THETA)*DELT
      CALL PHASE(ZRE,ZIM,PRZ)
      ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)
      RETURN
      END

      SUBROUTINE FRT2(UP,FR,ZRE,ZIM,NALL,DELT)
      dimension UP(1000000)
      DIMENSION W(1000000)
      real*4 ZRE,ZIM,PRZ
      THETA=6.283185*FR*DELT
      C=COS(THETA)*2.0
      NR1=1
      NR2=NALL
      NDR1=NR2-1
      W(1)=UP(NR2)
      W(2)=C*W(1)+UP(NDR1)
      NDATA=NR2-NR1+1
      NTR1=NDATA-1
      NTR2=NDATA-2
      DO 97 I=3,NDATA
      I1=I-1
      I2=I-2
      NDRI=NR2-I+1
      W(I)=C*W(I1)-W(I2)+UP(NDRI)
97    CONTINUE
      ZRE=(W(NDATA)-W(NTR2)+UP(NR1))*DELT/2.0
      ZIM=W(NTR1)*SIN(THETA)*DELT
c      CALL PHASE(ZRE,ZIM,PRZ)
c      ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)
      RETURN
      END

      SUBROUTINE PHASE(X,Y,PHI)
      real*4 x,y,phi
      IF(X) 21,20,22
20    IF(Y) 23,24,25
23    PHI=1.5*3.141592
      GO TO 28
24    PHI=0.0
      GO TO 28
25    PHI=0.5*3.141592
      GO TO 28
21    PHI=ATAN(Y/X) +3.141592
      GO TO 28
22    IF(Y) 26,27,27
27    PHI=ATAN(Y/X)
      GO TO 28
26    PHI=ATAN(Y/X)+2.0*3.141592
      GO TO 28
28    CONTINUE
      PHI=PHI/6.283184
      PHI=PHI-AINT(PHI)
      RETURN
      END

 
