c---------------------------------------------------------------------
c---------------------------------------------------------------------
      implicit none
      character*8 codi(2000),neti(2000)
      integer*4   nstai
      real*4      figi(2000),fici(2000),lami(2000)
C      common /stn/ nstai,codi,figi,fici,lami,neti
c ---
      real*4 cor(500,2,2000)
C      common /trk/ cor
c ---
      real*4 elat,elatc,elon,elonc
      real*4 v(3,2000),dvdz(3,2000),tm(6)
      real*4 ampr(2000),ampl(2000),ratio(2000),qR(2000),qL(2000)

      real*4 cr(2000),ur(2000),wvr(2000)
      real*4 cl(2000),ul(2000),wvl(2000),fr(2000)
      integer idate(5)
      logical*1 applyQ
      character*1 sigR,sigL,its,chnprefix
      character*3  symb, chn
      character*2  symbo,symbik
      character*2 wavetype
      character*10 symbol
      character*255 runfile,infile,outfile
      character*256 model,fnam11,fnam12,fnam13
      logical key_compr
c ---
      integer*4 idepth,im,iq,j
      integer*4 lso,lsy,ista,mode,n1,n2pow,narg,iargc,nbase,ncor
      integer*4 nd,nf,nout,lnblnk,npoints,nt
      real*4    aM,cazz,cince,cinch,cincz
      real*4    df,drad,dt,f0
      real*4    fix_vel,cincn,depth
      real*4    pi,R0,sdate,sqr_8pi
      real*4    TMIN,TMAX,vmax,I0(2000)
C      real*4 azi, azi_back

      integer*4 nsize
      parameter (nsize=8192)
      real*4 seismz(nsize),seismn(nsize),seisme(nsize)

c ---
c      data GEO/0.993277/
      data pi/3.14159265/,sqr_8pi/5.01325655/
      data idate/1999,1,1,1,1/,sdate/0.0/,im/6/
      data cincz/0.0/,cazz/0.0/,cinch/90.0/,cincn/0.0/,cince/90.0/
      data key_compr/.false./,fix_vel/2.85/,R0/6371./
C      data const1/1.E-07/
C-----------ARGUMENTS-----------------------------------------------
C-----runfile: file with parameters for running program-------------
C-----infile: output file of SURFLEV--------------------------------
C-----outfiles: output files of SURFLEV-----------------------------
C-----wavetype: R/L/RL/LR ------------------------------------------
C-----mode:0,1,2,3,...-----------------------------------------------
C-------------INITIATION--------------------------------------------S
       print*,'to make synthetic seismogram using output of SURF_DEEP'
       print*,'vertical direction is down; displacements as output;'
       print*,'step-function source-time function for earthquake;'
       print*,'delta-function as a source-function for simple force or'
       print*,'explosion; (-pi/4) and polarization factors are Included'
C------------------------------------------------------------------------
C------notations of eigenfunctions follow the Kluwer's book------
C-------------------Check input parameters---------------------------S
      narg=iargc()
      drad=pi/180.
      if (narg.lt.6.or.narg.gt.8)then
         print *,'usage : surfsyn runfile infile outfiles wavetype mode',
     +           '  depth [-c [fix_vel]]' 
         stop
      end if

C-------------------Read in params and Initialize--------------------S
      call GETARG(1,runfile)
      call GETARG(2,infile)
      call GETARG(3,outfile)
      nout=lnblnk(outfile)
      call GETARG(4,wavetype)
      call GETARG(5,symbo)
      read(symbo,*)mode
      mode=mode+1
      symbik='  '
      if(mode.lt.10)write(symbik,'(i1)')mode
      if(mode.ge.10)write(symbik,'(i2)')mode
      lso=lnblnk(symbo)
      call GETARG(6,symbol)
      read(symbol,*)depth
      idepth=depth
      symb='   '
      if(idepth.lt.10)write(symb,'(I1)')idepth
      if(idepth.ge.10.and.idepth.lt.100)write(symb,'(I2)')idepth
      if(idepth.ge.100)write(symb,'(I3)')idepth
      lsy=lnblnk(symb)
C      n1=nout+3+lsy+lso
C      current(1:n1)=outfile(1:nout)//'_'//symb(1:lsy)//'_'//symbo(1:lso)//'_'
      if(wavetype.ne.'R '.and.wavetype.ne.'L '.and.wavetype.ne.'RL'.and.
     +wavetype.ne.'LR') then
         print*, 'WRONG WAVE TYPE'
         STOP 
      end if
      sigR='-'
      sigL='-'
      if(wavetype(1:1).eq.'R'.or.wavetype(2:2).eq.'R') sigR='+'
      if(wavetype(1:1).eq.'L'.or.wavetype(2:2).eq.'L') sigL='+'
      if(narg.gt.6)then
         call GETARG(7,symbol)
         if(symbol(1:2).eq.'-c')then
            key_compr=.true.
            print*,'nondispersive seismograms will be created'
            if(narg.eq.8)then
               call GETARG(8,symbol)
               read(symbol,*)fix_vel
            endif
            print*,'nondispersive seismograms will be created; ',
     +             'fixed velocity=',fix_vel
         else
            print*,'Wrong argument, try again'
         endif 
      endif 

      PRINT*,'Wavetype=',wavetype,' mode=',mode,' depth=',idepth

      call read_params( infile, runfile, mode, model, elat, elon, aM, dt, 
     +                  tm, im, iq, its, nd, npoints, nt, TMIN, TMAX, vmax,
     +                  nstai, codi, figi, fici, lami, neti )

C---------- trace along the GCP for each station ---------------------
      applyQ = .TRUE.
C      elatc = atan(GEO*tan(drad*elat))/drad
      elatc = elat
      call atracer(model,elatc,elon,ncor,applyQ, nstai,fici,lami, cor)

      call surfread(infile,sigR,sigL,symbik,nt,nd,depth,fr,cr,ur,wvr,
     +              cl,ul,wvl,v,dvdz,ampr,ampl,ratio,qR,qL,I0)

      PRINT*,nt,' input periods'

C----------FINDING BASE SIZE ---------------------S  
      do j=1,12
         if(2**j.ge.npoints) go to 2   
      end do
2     nbase=2**j
      n2pow=j
      PRINT*,'n2pow=',n2pow,'  nbase=',nbase,' npoints=',npoints
      df=1./(dt*nbase)
      f0=0.0 
      nf=nbase/2

C--------------- main loop on stations ----------------E  
      elatc = elatc*drad
      elonc = elon*drad
      do ista =1,nstai
         write(*,*) ista,' STA= ',codi(ista),' NET= ',neti(ista),
     +              ' STLAT= ',figi(ista),' STLON= ',lami(ista)
         call cal_synsac( ista, its, sigR, sigL, cor, TMIN, TMAX, vmax, n2pow, fix_vel, iq,
     +                    npoints, fr, df, nf, dt, nt, elatc, elonc, key_compr, qR, qL,
     +                    im, aM, tm, ampl, cr, ul, ur, wvl, wvr, v, dvdz, ratio, I0, 
     +                    fici(ista), lami(ista), seismz, seismn, seisme )
C ------------ write sac files ---------
            if( dt.ge.1 ) then
               chnprefix='L'
            else 
               chnprefix='B'
            endif

            chn = chnprefix//'HZ'

            fnam11 = 'sac/'//codi(ista)(1:lnblnk(codi(ista)))//'.'//chn//'.SAC'
            call wsac(fnam11,elat,elon,codi(ista),figi(ista),
     +                lami(ista),dt,chn,neti(ista),npoints,seismz)
            chn = chnprefix//'HN'
            fnam12 = 'sac/'//codi(ista)(1:lnblnk(codi(ista)))//'.'//chn//'.SAC'
            call wsac(fnam12,elat,elon,codi(ista),figi(ista),
     +                lami(ista),dt,chn,neti(ista),npoints,seismn)
            chn = chnprefix//'HE'
            fnam13 = 'sac/'//codi(ista)(1:lnblnk(codi(ista)))//'.'//chn//'.SAC'
            call wsac(fnam13,elat,elon,codi(ista),figi(ista),
     +                lami(ista),dt,chn,neti(ista),npoints,seisme)
C            fnam11 = 'mv *.R.spec '//'SPEC/'//codi(ista)(1:lnblnk(codi(ista)))//'.R.sp'
C            call system(fnam11)
C            fnam11 = 'mv *.Z.spec '//'SPEC/'//codi(ista)(1:lnblnk(codi(ista)))//'.Z.sp'
C            call system(fnam11)
      enddo

      stop
      END
