c  findtiltaxis.f

c  This program rotates horizontal components to find maximum coherence 
c  with the vertical component (using signal-free noise samples as input).
c  Uses coherence with complex transfer function to account for phase shifts
c  between vertical and horizontal components.
c  Then tilts z-axis away from that direction to minimize horizontal noise
c  on the vertical.  

c  For input, pick multiple windows (~20?) 
c  files as input, avoiding times with any earthquake signals (use maximum
c  noise periods)  


c  Pipe in data   findtiltBYtransfer  < findtiltBYtransfer.inp
 
      parameter (maxnfreq=400, maxpts = 1000000, maxevnts = 100)

      real*4 bh1(maxpts), bh2(maxpts), bhz(maxpts) 
      real*4 freq(maxnfreq)
      real*4 timebeg(maxevnts,maxevnts)
      real*4 sumcrossre(maxnfreq),sumcrossim(maxnfreq)
      real*4 sumppower(maxnfreq)
      real*4 sumzpower(maxnfreq), tlength(maxevnts)
      real*4 begh1(maxevnts), begz(maxevnts),delth1(maxevnts)
      real*4 deltz(maxevnts), begh2(maxevnts),delta(maxevnts)
      real*4 winh1(maxevnts,maxpts),winz(maxevnts,maxpts)
      real*4 winh2(maxevnts,maxpts)
      real*4 gausst(maxpts),tdataz(maxpts)
      real*4 tdatah1(maxpts), tdatah2(maxpts),delth2(maxevnts)
      real*4 qreal(maxnfreq),qimag(maxnfreq),coher(maxnfreq)
      real*4 epsilon(maxnfreq)
      real*4 qamp(maxnfreq), qphse(maxnfreq),ampdev(maxnfreq)
      real*4 phsedev(maxnfreq)
      real*4 preal,pimag
      real*4 horz(maxpts),vert(maxpts)
      
      integer*4 nfreq, nwindows(maxevnts), nptsw(maxevnts)
       integer*4 nptsh1(maxevnts), nptsh2(maxevnts), nptsz(maxevnts)
      
      character*70 foutput
      
      character*70 BH1fn(maxevnts), BH2fn(maxevnts), BHZfn(maxevnts)

      pi = 3.1415928
      convdeg = 3.1415928/180.
      circ = 6371.*3.1415928/180.
      twopi = 3.1415928*2.
 
c  initially, just do a few frequencies for which we will find average coherence
       nfreq = 7
      do ifreq = 1, nfreq
         freq(ifreq) = (ifreq-1)*.005 + .01
      enddo

c  foutput file for output of transfer function and coherence 
        read(*,'(a)') foutput
      open(11, file = foutput)

c  read list of files to be analyzed (events or time series)
c  nogauss = 0 means don't apply gaussian window - use record as is (For already windowed
c   signals)
c  nogauss = 1 means do apply gaussian window (For noise sample)
      read(*,*) nevents, nogauss
      nobs = 0
      do iev = 1, nevents
c  read in horizontal BH1 file name, then BH2, then associated vertical file name - all should have
c  same beginning time and same time increments
        read(*,'(a)') BH1fn(iev)
        read(*,'(a)') BH2fn(iev)
        read(*,'(a)') BHZfn(iev)
c  read in number of windows selected and timelength of window for seismogram
	read(*,*) nwindows(iev), tlength(iev)
c  timebeg is beginning of window relative to start of record
	do iwin = 1, nwindows(iev)
	  read(*,*) timebeg(iev,iwin)
	enddo
        call rsac1(BH1fn(iev),tdatah1,nptsh1(iev),begh1(iev),
     1       delth1(iev),maxpts,nerr)
        call rsac1(BH2fn(iev),tdatah2,nptsh2(iev),begh2(iev),
     1       delth2(iev),maxpts,nerr)
        call rsac1(BHZfn(iev),tdataz,nptsz(iev),begz(iev),
     1       deltz(iev),maxpts,nerr)
        if ((begh1(iev).ne.begz(iev)).or.(begh2(iev).ne.begz(iev))) then
	  write (*,*) iev,' beginning times not identical'
	endif
	if ((delth1(iev).ne.deltz(iev)).or.(delth2(iev).ne.deltz(iev))) then
	  write(*,*) iev,' delt not identical'
	endif
c  set up Gaussian window function - want 3 sigma at limits of time window
c  and amplitude = 1 at center
	npts = tlength(iev)/deltz(iev) +1
        if (nogauss.eq.1) then
          sigma = tlength(iev)/6.0
    	  sigma2 = 2.0*sigma*sigma
	  tmiddle = tlength(iev)/2.0
	  do ig = 1, npts
	    tg = (ig-1)*deltz(iev) - tmiddle
	    gausst(ig) = exp(-(tg*tg)/sigma2)
	  enddo
	endif
C  now apply gaussian to windows
        do iwin = 1, nwindows(iev)
	  nobs = nobs + 1
	  iwstart = timebeg(iev,iwin)/deltz(iev)
	  delta(nobs) = deltz(iev)
	  nptsw(nobs) = npts
	  do iw = 1, npts
	    if (nogauss.eq.1) then
	      winh1(nobs,iw) = tdatah1(iw+iwstart)*gausst(iw)
	      winh2(nobs,iw) = tdatah2(iw+iwstart)*gausst(iw)
	      winz(nobs,iw) = tdataz(iw+iwstart)*gausst(iw)
	    else
	      winh1(nobs,iw) = tdatah1(iw+iwstart)
	      winh2(nobs,iw) = tdatah2(iw+iwstart)
	      winz(nobs,iw) = tdataz(iw+iwstart)
	    endif
	  enddo
        enddo
      enddo
c  end of input loop over events


	      
c rotate horizontal components every 5 degrees and calculate coherence
c  with vertical
      cohermax = 0.0
      do jang = 0, 175, 5
        theta = jang*pi/180.0
	costheta = cos(theta)
	sintheta = sin(theta)
	avgcoher = 0.0

        do ifreq = 1, nfreq
	  sumcrossre(ifreq) = 0.0
	  sumcrossim(ifreq) = 0.0
	  sumppower(ifreq) = 0.0
	  sumzpower(ifreq) = 0.0
        enddo
c  cycle through all windows, all events to find coherence
        do iwind = 1, nobs
	   do i = 1, nptsw(iwind)
	     vert(i) = winz(iwind,i)
	     horz(i) = winh1(iwind,i)*costheta+winh2(iwind,i)*sintheta
	   enddo
c
c  get real and imaginary (phases and amplitudes) at desired frequencies
c
          do ifreq = 1, nfreq
            call frt2(horz, freq(ifreq), preal,
     1                      pimag,
     2                      nptsw(iwind),delta(iwind))
            call frt2(vert, freq(ifreq), zreal,
     1                      zimag,
     2                      nptsw(iwind),delta(iwind))
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
c  calculate coherence as function of frequency and find average
        do ifreq = 1, nfreq	
         coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/
     1                  (sumppower(ifreq)*sumzpower(ifreq))
         avgcoher = avgcoher + coher(ifreq)
        enddo
        avgcoher = avgcoher/nfreq
	if (avgcoher.gt.cohermax) then
	  maxtheta = jang
	  cohermax = avgcoher
	endif
      enddo
      
c  now search every 1 degree to find better estimate of maximum
c  coherence direction between horizontal and vertical	
      do jang = maxtheta-3,maxtheta+3, 1
        theta = jang*pi/180.0
	costheta = cos(theta)
	sintheta = sin(theta)
	avgcoher = 0.0

        do ifreq = 1, nfreq
	  sumcrossre(ifreq) = 0.0
	  sumcrossim(ifreq) = 0.0
	  sumppower(ifreq) = 0.0
	  sumzpower(ifreq) = 0.0
        enddo
c  cycle through all windows, all events to find coherence
        do iwind = 1, nobs
	   do i = 1, nptsw(iwind)
	     vert(i) = winz(iwind,i)
	     horz(i) = winh1(iwind,i)*costheta+winh2(iwind,i)*sintheta
	   enddo
c
c  get real and imaginary (phases and amplitudes) at desired frequencies
c
          do ifreq = 1, nfreq
            call frt2(horz, freq(ifreq), preal,
     1                      pimag,
     2                      nptsw(iwind),delta(iwind))
            call frt2(vert, freq(ifreq), zreal,
     1                      zimag,
     2                      nptsw(iwind),delta(iwind))
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
c  calculate coherence as function of frequency and find average 
        do ifreq = 1, nfreq	
         qreal(ifreq) = sumcrossre(ifreq)/sumppower(ifreq)
         qimag(ifreq) = sumcrossim(ifreq)/sumppower(ifreq)
         coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/
     1                  (sumppower(ifreq)*sumzpower(ifreq))
         avgcoher = avgcoher + coher(ifreq)
	  call phase(qreal(ifreq),qimag(ifreq),qphse(ifreq))
	  if (qphse(ifreq).gt.0.5) qphse(ifreq)=qphse(ifreq)-1.0
        enddo
        avgcoher = avgcoher/nfreq
	if (avgcoher.gt.cohermax) then
	  maxtheta = jang
	  cohermax = avgcoher
	endif
      enddo
c  check to see if better to rotate by 180 degrees for optimum direction to get phase close to zero
c  instead of close to pi
       if (abs(qphse(5)).gt. 0.3) maxtheta = maxtheta + 180.
c      phideg = phi * 180./pi
      write(*,*) maxtheta, cohermax 
      write(11,*) maxtheta, cohermax
       
c  at optimum direction, calculate transfer function
      theta = maxtheta*pi/180.0
	costheta = cos(theta)
	sintheta = sin(theta)
	avgcoher = 0.0

      nfreq = 396
        do ifreq = 1, nfreq
         freq(ifreq) = (ifreq-1)*.001 + .005
	  sumcrossre(ifreq) = 0.0
	  sumcrossim(ifreq) = 0.0
	  sumppower(ifreq) = 0.0
	  sumzpower(ifreq) = 0.0
        enddo
c  cycle through all windows, all events to find coherence
        do iwind = 1, nobs
	   do i = 1, nptsw(iwind)
	     vert(i) = winz(iwind,i)
	     horz(i) = winh1(iwind,i)*costheta+winh2(iwind,i)*sintheta
	   enddo
c
c  get real and imaginary (phases and amplitudes) at desired frequencies
c
          do ifreq = 1, nfreq
            call frt2(horz, freq(ifreq), preal,
     1                      pimag,
     2                      nptsw(iwind),delta(iwind))
            call frt2(vert, freq(ifreq), zreal,
     1                      zimag,
     2                      nptsw(iwind),delta(iwind))
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
      write(11,*) 'freq   admittance std(admit)  phase  std(phase) coher
     1    real    imag'  
      do ifreq = 1, nfreq	
        qreal(ifreq) = sumcrossre(ifreq)/sumppower(ifreq)
        qimag(ifreq) = sumcrossim(ifreq)/sumppower(ifreq)
        coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/
     1                  (sumppower(ifreq)*sumzpower(ifreq))
c  epsilon is fractional uncertainty in amplitude of transfer function (1 stddev)
c  and also approximately 1 stddev in phase (in radians) Bendat and Piersol page 317
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

 
