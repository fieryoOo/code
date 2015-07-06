c---------------------------------------------------------------------
c---------------------------------------------------------------------
      implicit none
      character*8 cod(2000),codi(2000),net(2000),neti(2000)
      integer*4   nsta,nstai,nsize
      real*4      fig(2000),lam(2000),figi(2000),fici(2000),lami(2000)
      common /stn/nsta,nstai,cod,fig,lam,codi,figi,fici,lami
c ---
      real*4 cor(500,2,2000)
      common /trk/ cor
c ---
      character*8 chn
      complex*8 br(6),bl(6),sumr,suml
      complex*8 al(2000),az(2000),ah(2000)
      complex*8 step
      real*4    elat,elatc,elon,elonc,slat,slon
      real*4 v(3,2000),dvdz(3,2000),tm(6)
      real*4 ampr(2000),ampl(2000),ratio(2000),qR(2000),qL(2000)
      real*4 amp(2048),T_curr(2048)
      parameter (nsize=4096)
      real*4 cr(2000),ur(2000),wvr(2000),t(2000),wvar(2000)
      real*4 cl(2000),ul(2000),wvl(2000),fr(2000)
      real*4 wr(nsize),wl(nsize)
      real*4 sreV(nsize),simV(nsize)
      real*4 sreR(nsize),simR(nsize)
      real*4 sreT(nsize),simT(nsize)
      real*4 vu(3),du(3),wvn(2)
      real*4 seismz(nsize),seismn(nsize),seisme(nsize)
      real*4 ul_int(nsize),ur_int(nsize)
      real*4 seismt(nsize),seismr(nsize),qL_int(nsize),qR_int(nsize)
      integer idate(5),nper(20)
      character*2 wavetype
      character*10 symbol,bred
      character*3  symb
      character*2  symbo,symbik
      character*1 sigR,sigL,its
      character*255 runfile,infile,outfile,current,stafile
      character*256 model,fnam11,fnam12,fnam13
      character*80 str,outseism(50),outtit(50),outspec(50)
      logical key_compr
c ---
      integer*4 i,idepth,ierr,ii,im,iq,k,j,k_spec,kl,l,lin,ll,lq
      integer*4 lso,lsy,m,mm,mmm,mode,n1,n2pow,narg,iargc,nbase,ncor
      integer*4 nd,NDIST,nf,nk7,nout,lnblnk,npoints,nt
      real*4    GEO,aM,azi,azi_back,cazz,cince,cinch,cincz,const1
      real*4    cs,cs_b,DEL,DELS,df,dist,dper,drad,dt,f0,fi,fi_back
      real*4    fix_vel,cincn,depth,strike,dip,slip,per1,per2
      real*4    PERMAX,PERMIN,pi,R0,sc,sc_b,sdate,spread,sqr_8pi
      real*4    TMIN,TMAX,tstart,vmax,w,I0(2000)
c ---
      data GEO/0.993277/
      data pi/3.14159/,sqr_8pi/5.01325655/,step/(1.0,0.0)/
      data idate/1999,1,1,1,1/,sdate/0.0/,im/6/
      data cincz/0.0/,cazz/0.0/,cinch/90.0/,cincn/0.0/,cince/90.0/
      data key_compr/.false./,fix_vel/2.85/,const1/1.E-09/,R0/6371./
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
C-------------------Initiation-------------------------------------------S
           narg=iargc()
           drad=pi/180.
      if (narg.lt.6.or.narg.gt.8)then
      print *,'usage : surfsyn runfile infile outfiles wavetype mode',
     +         '  depth [-c [fix_vel]]' 
                    stop
            end if
           call GETARG(1,runfile)
           call GETARG(2,infile)
           lin=lnblnk(infile)
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
           n1=nout+3+lsy+lso
      current(1:n1)=outfile(1:nout)//'_'//symb(1:lsy)//'_'//symbo(1:lso)//'_'
           if(wavetype.ne.'R '.and.wavetype.ne.'L '.and.wavetype.ne.'RL'.and.
     +wavetype.ne.'LR') then
             print*, 'WRONG WAVE TYPE'
            STOP 
            end if
            sigR='-'
            sigL='-'
            if(wavetype(1:1).eq.'R'.or.wavetype(2:2).eq.'R') sigR='+'
            if(wavetype(1:1).eq.'L'.or.wavetype(2:2).eq.'L') sigL='+'
            if(narg.gt.6)THEn
            call GETARG(7,symbol)
            if(symbol(1:2).eq.'-c')then
            key_compr=.true.
            print*,'nondispersive seismograms will be created'
            if(narg.eq.8)then
            call GETARG(8,symbol)
            read(symbol,*)fix_vel
                     endif
      print*,'nondispersive seismograms will be created; ',
     +        'fixed velocity=',fix_vel
                                   else
            print*,'Wrong argument, try again'
                                  endif 
                                  ENDif 
C-------------------Initiation-------------------------------------------E
C-----------reading phase velocity info-----S
            open(1,file=infile(1:lin)//'.phv',status='OLD')
            k=1
            m=0   
            do i=1,1000
            read(1,'(a)',end=9988)bred
            if(i.eq.1)read(bred,*)per1
            if(i.eq.2)read(bred,*)per2
            m=m+1
            if(bred(6:10).eq.'     ')then
            nper(k)=m-1
            read(1,'(1X)')
            k=k+1
            m=0       
                                     endif
            enddo
9988        print*,k-1,' modes in input'
            close(1)
            nt=nper(mode)
            dper=per2-per1
            PERMAX=(nt-1)*dper+per1
            PERMIN=per1
            write(*,*) nt,' periods for the mode ',mode,
     +                 ' per_incr=',dper,' PERMAX=',PERMAX
            if(nt.gt.1998)stop 'too many periods!'
C-----------reading phase velocity info-----E
C-----------Reading of input parameters--------------S
           open(2,file=runfile,STATUS='OLD')
C-----------nt is a number of periods------------
C-----------nd is a number of depth points at SURLEV--------
C-----------npoints is a number of points of the output seismpogram
C-----------dt is an output time increment-------
C-----------elat,elon  are event coordinates
C-----------stlat,stlon  are station coordinates
C-----------depth is a source depth (one of given at SURFLEV)-----
C-----------TMIN-the lowest untapered period---------------------------
C-----------TMAX-the highest untapered period--------------------------
C-----------iq  is a power of taper          --------------------------
C-----------vmax is a maximal velocity for starting point -------------
C-----------its is type of source:F/single force;Q/quake-tensor;E/explosion
C-----------its is type of source:C/quake-double dipole----------------
      PRINT*,'Wavetype=',wavetype,' mode=',symbo
c --- read parameter file  ---
      read(2,'(a200)')model
      read(2,'(a200)')stafile
c --  read station list  ---
      open(20,file=stafile,status='old')
      nsta = 1
    1 read(20,*,end=98) net(nsta),cod(nsta),fig(nsta),lam(nsta)
      nsta = nsta+1
      goto 1
   98 nsta = nsta-1
      close(20)
      read(2,*)nd,npoints,dt,TMIN,TMAX,iq,vmax
      if(TMAX.gt.(nt-1)*dper+per1)TMAX=PERMAX-2.*dper
      if(TMIN.lt.per1)TMIN=PERMIN+2.*dper
      print*,' ndepth= ',nd
      print*,'npoints= ',npoints, ' dt = ',dt
      print*,'depth= ',depth,' Tmin=',tmin,' Tmax=',Tmax
           read(2,'(a)') its
      if(its.ne.'Q'.and.its.ne.'E'.and.its.ne.'F'.and.its.ne.'C')
     +STOP ' WRONG TYPE OF SOURCE: Q,E or F are wanted'
      print 1200 ,          iq,vmax,its 
1200   format(' iq=',i4,' vmax=',F10.3'    type=',a1)
C--if its=Q ---time function is a step function-----------
C---in the following order:xx,yy,zz,xy,xz,yz---
C-from CMT is should be taken as
C as del-del, fi-fi, rr, -(fi-del),r-del,-(fi-r)
C--if its=E or F time source function is a delta function 
C--if its=C ---time function is a step function and moment tensor is found from angles
C--strike, dip and rake(slip), aM is the best double couple moment;
C-----time source function is a step function for earthquake--
           if(its.eq.'F'.or.its.eq.'E') im=3 
           if(its.ne.'C')read(2,*)aM,(tm(m),m=1,im)
            if(its.eq.'C') then
       read(2,*)aM,strike,dip,slip
       call angles2tensor(strike,dip,slip,tm)
                           endif
c --- read requested station list ---
       read(2,*) elat,elon
        elatc = atan(GEO*tan(drad*elat))/drad
       write(*,*) ' ELAT=',elat,' ELON=',elon
       nstai = 1
    3  read(2,'(a80)',end=4) str
       if(str(1:1).eq.' ') then
         read(str,*) figi(nstai),lami(nstai)
         fici(nstai) = atan(GEO*tan(drad*figi(nstai)))/drad
         write(codi(nstai),'("Z",i04)') nstai 
         neti(nstai) = 'ZZ'
         do i =1,lnblnk(codi(nstai))
            if(codi(nstai)(i:i).eq.' ') codi(nstai)(i:i)='0'
         enddo
         write(*,*) nstai,' STA=',codi(nstai),' NET=',neti(nstai),
     +                  ' STLAT=',figi(nstai),' STLON=',lami(nstai)
       else
         read(str,*) codi(nstai)
         do ii = 1,nsta
           if(cod(ii).eq.codi(nstai)) then
             figi(nstai) = fig(ii)
             fici(nstai) = atan(GEO*tan(drad*fig(ii)))/drad
             lami(nstai) = lam(ii)
             neti(nstai) = net(ii)
             write(*,*) nstai,' STA=',codi(nstai),' NET=',neti(nstai),
     +                  ' STLAT=',figi(nstai),' STLON=',lami(nstai)
             goto 5
           endif
         enddo
         goto 3
       endif
    5  nstai = nstai+1
       goto 3
    4  nstai = nstai-1
       if(its.eq.'F')write(*,'(5x,"aM",8x,"fx ",8x,"fy ",8x,"fzz")')
       if(its.ne.'F')write(*,1000)
       write(*,'(7e11.3)') aM,(tm(i),i=1,im)
       call atracer(model,elatc,elon,ncor)
       call surfread(infile,sigR,sigL,symbik,nt,nd,depth,t,cr,ur,wvr,
     + cl,ul,wvl,v,dvdz,ampr,ampl,ratio,qR,qL,I0)
       PRINT*,nt,' input periods'
C-----------Reading of input parameters--------------E
C----------FINDING BASE SIZE ---------------------S  
           do j=1,12
           if(2**j.ge.npoints) go to 2   
           end do
2          nbase=2**j
           n2pow=j
           PRINT*,'n2pow=',n2pow,'  nbase=',nbase,' npoints=',npoints
           df=1./dt/nbase
           f0=0.0 
           nf=nbase/2
C----------FINDING BASE SIZE ---------------------E  
           elatc=elatc*drad
           elonc = elon*drad
      do mm =1,nstai
           write(*,*) "Station: ",codi(mm)
           slat=fici(mm)*drad
           slon=lami(mm)*drad
           call AZIDL(SLAT,SLON,ELATC,ELONC,DEL,DELS,AZI_BACK,AZI)
           dist=dels
           fi=AZI/drad
           fi_back=AZI_BACK/drad
         write(*,*) nstai,' STA= ',codi(mm),' NET= ',neti(mm),
     +             ' STLAT= ',figi(nstai),' STLON= ',lami(nstai),
     +             ' DIST= ',dist,' AZI= ',fi,' BACK_AZI= ',fi_back
c ---
           cs=cos(AZI)
           sc=sin(AZI)
           cs_b=cos(AZI_BACK+pi)
           sc_b=sin(AZI_BACK+pi)
C-----------Reading of spectral data from SURF_DEEP output-----S
           do j=2,nt+1
           fr(j)=1./t(j)
           wvar(j) = pi*2.0/cor(j-1,1,mm)*fr(j)
C          PRint*,t(j),ampr(j),ratio(j),qR(j),wvr(j)
           enddo
           wvar(1) = wvr(1)
           wvar(nt+2) = wvr(nt+2)
           fr(1)=10.0
           fr(nt+2)=0.0
           qR(1)=qR(2)
           qL(1)=qL(2)
           qR(nt+2)=qR(nt+1)
           qL(nt+2)=qL(nt+1)
C-----------Reading of spectral data from SURF_DEEP  output-----E
C----------CALCULATION OF SPECTRA--------------------------------------S
         write(*,*) 'NT ======',nt
           Do j=2,nt+1
c----------Source term calculations-----------------S
           vu(1)=v(1,j)
           vu(2)=v(2,j)
           vu(3)=v(3,j)
           du(1)=dvdz(1,j)
           du(2)=dvdz(2,j)
           du(3)=dvdz(3,j)
           w=pi*2.0*fr(j)
           wvn(1)=wvr(j)
           wvn(2)=wvl(j)
           if(its.eq.'F') call force(sigR,sigL,cs,sc,vu,br,bl)
         if(its.ne.'F')then
           if(its.eq.'Q'.or.its.eq.'C') step=cmplx(0.,-1./w)
           call source(sigR,sigL,cs,sc,wvn,vu,du,br,bl)
                      end if
           suml=(0.0,0.0)
           sumr=(0.0,0.0)
          do m=1,im
           sumr= sumr+tm(m)*br(m)*step
           suml= suml+tm(m)*bl(m)*step
                  enddo
           al(j)=suml*ampl(j)*aM*const1
cMB           spread=1./sqrt(R0*sin(DEL)*wvn(1))
           spread=1./sqrt(R0*sin(DEL))
c          if(j.eq.2)PRint*,' SPREAD=',spread*sqrt(wvn(1)), ita
cMB        az(j)=sumr*ampr(j)*aM*const1*conjg(cor(j-1))
           az(j)=sumr*aM*const1*cor(j-1,2,mm)/(2.0*cr(j)*sqrt(ur(j)*I0(j)))/sqrt(2.0*pi)
cMB      write(*,*)t(j),cr(j),ur(j),wvr(j),ampr(j),cor(j-1,1),cor(j-1,2)
cMB        az(j)=sumr*ampr(j)*aM*const1
cxx        WRite(33,*)t(j),cabs(az(j))*spread
           ah(j)=az(j)*ratio(j)
                       end do
c----------Source term calculations-----------------E
C----------INTERPOLATION OF WAVENUMBERS AND Q-factors--S        
           if (sigR.eq.'+')then
cMB        call intpol(fr,wvr,nt+2,f0,df,nf,wr,ierr)
           call intpol(fr,wvar,nt+2,f0,df,nf,wr,ierr)
           if(ierr.ne.0)STOP'WRONG INTERPOLATION: RAYLEIGH WAVENUMBER'
           call intpol(fr,qR,nt+2,f0,df,nf,qR_int,ierr)
           if(ierr.ne.0)STOP'WRONG INTERPOLATION: RAYLEIGH ATTENUATION'
           call intpol(fr,ur,nt+2,f0,df,nf,ur_int,ierr)
           if(ierr.ne.0)STOP'WRONG INTERPOLATION: RAYLEIGH GR.VELOCITY'
           DO mmm=1,nf
           if(f0+df*(mmm-1).gt.fr(2))qR_int(mmm)=qR_int(mmm-1)
           enddo
                       endif
           if (sigL.eq.'+')then
           call intpol(fr,wvl,nt+2,f0,df,nf,wl,ierr)
           if(ierr.ne.0)STOP'WRONG INTERPOLATION: LOVE WVENUMBER'
           call intpol(fr,qL,nt+2,f0,df,nf,qL_int,ierr)
           if(ierr.ne.0)STOP'WRONG INTERPOLATION: LOVE ATTENUATION'
           call intpol(fr,ul,nt+2,f0,df,nf,ul_int,ierr)
           if(ierr.ne.0)STOP'WRONG INTERPOLATION: LOVE GR.VELOCITY'
           DO mmm=1,nf
           if(f0+df*(mmm-1).gt.fr(2))then
           qL_int(mmm)=qL_int(mmm-1)
           uL_int(mmm)=uL_int(mmm-1)
                                     endif
           enddo
                         endif
C----------INTERPOLATION OF WAVENUMBERS AND Q-factors--E        
C----------INTERPOLATION OF SPECTRAL AMPLITUDES-------------S  
           if (sigR.eq.'+')
     +call calspecr(nt,nf,f0,df,fr,az,ah,sreV,simV,sreR,simR,TMIN,TMAX,iq)
          if (sigL.eq.'+') 
     +call calspecl(nt,nf,f0,df,fr,al,sreT,simT,TMIN,TMAX,iq)
C----------INTERPOLATION OF SPECTRAL AMPLITUDES-------------E  
C----------CALCULATION OF SPECTRA--------------------------------------E
C----------CALCULATION OF SEISMOGRAMS--------------------S
           NDIST=1
      call namer(dels,ndist,n1,current,outseism,outtit,outspec,sigR,sigL,nk7)
C----------output filenames are defined--------->>>>>>
C----------DISTANCE LOOP-----------------------S      
           DO l=1,ndist
           lq=(l-1)*3
C--------------calculating everything and writing spectra----S
           if (sigR.eq.'+') THEN
      call syn(wr,dist,dt,f0,df,sreV,simV,ur_int,qR_int,seismz,tstart,
     +vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
      open(41,file=outspec(lq+1)(1:nk7+1))
      do kl=2,k_spec
      write(41,*)T_curr(kl),amp(kl)
      enddo
      close(41)
      call syn(wr,dist,dt,f0,df,sreR,simR,ur_int,qR_int,seismr,tstart,
     +vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
      open(41,file=outspec(lq+2)(1:nk7+1))
      do kl=2,k_spec
      write(41,*)T_curr(kl),amp(kl)
      enddo
      close(41)
            fnam11='w/'//codi(mm)(1:lnblnk(codi(mm)))//'.'//outseism(lq+1)(1:nk7)
            open(11,file=fnam11)
           ELSEIF (sigL.eq.'+') THEN
      call syn(wl,dist,dt,f0,df,sreT,simT,ul_int,qL_int,seismt,tstart,
     +vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
      open(41,file=outspec(lq+3)(1:nk7+1))
      do kl=2,k_spec
      write(41,*)T_curr(kl),amp(kl)
      enddo
      close(41)
           ENDIF
C--------------calculating everything and writing spectra----E
C--------------writing seismograms----S
          fnam12='w/'//codi(mm)(1:lnblnk(codi(mm)))//'.'//outseism(lq+2)(1:nk7)
            open(12,file=fnam12)
          fnam13='w/'//codi(mm)(1:lnblnk(codi(mm)))//'.'//outseism(lq+3)(1:nk7)
            open(13,file=fnam13)
           do ll=1,npoints
           seismn(ll)=0.0
           seisme(ll)=0.0
           if (sigR.eq.'+') then
C          seismz(ll)=-seismz(ll)
C-----------Positive 'Z' direction is 'Down'-----!!!!!!!!!!!!!!!!!!!!!!
C-----------horizontal components----S
           seismn(ll)=seismr(ll)*cs_b
           seisme(ll)=seismr(ll)*sc_b
           elseif (sigL.eq.'+') then
           seismn(ll)=seismn(ll) -seismt(ll)*sc_b
           seisme(ll)=seisme(ll) +seismt(ll)*cs_b
                            endif
C-----------horizontal components----E
C-----------vertical   component-----S
           if (sigR.eq.'+')write(11,*)tstart+dt*(ll-1),seismz(ll)
           write(12,*)tstart+dt*(ll-1),seismn(ll)
           write(13,*)tstart+dt*(ll-1),seisme(ll)
C-----------vertical   component-----E
           enddo
C--------------writing seismograms----E
           close(11)
           close(12)
          close(13)
C ------------ write sac files ---------
        fnam11 = 'sac/'//codi(mm)(1:lnblnk(codi(mm)))//'.BHZ.SAC'
        fnam12 = 'sac/'//codi(mm)(1:lnblnk(codi(mm)))//'.BHN.SAC'
        fnam13 = 'sac/'//codi(mm)(1:lnblnk(codi(mm)))//'.BHE.SAC'
        chn = 'BHZ'
        call wsac(fnam11,elat,elon,codi(mm),figi(mm),
     +            lami(mm),dt,chn,neti(mm),npoints,seismz)
        chn = 'BHN'
        call wsac(fnam12,elat,elon,codi(mm),figi(mm),
     +            lami(mm),dt,chn,neti(mm),npoints,seismn)
        chn = 'BHE'
        call wsac(fnam13,elat,elon,codi(mm),figi(mm),
     +            lami(mm),dt,chn,neti(mm),npoints,seisme)
        fnam11 = 'mv *.R.spec '//'SPEC/'//codi(mm)(1:lnblnk(codi(mm)))//'.R.sp'
        call system(fnam11)
        fnam11 = 'mv *.Z.spec '//'SPEC/'//codi(mm)(1:lnblnk(codi(mm)))//'.Z.sp'
        call system(fnam11)
           end DO
C----------DISTANCE LOOP------ ----------------E
C----------CALCULATION OF SEISMOGRAMS--------------------S
        enddo
           stop
1000  format(5x,"aM",8x,"mxx",8x,"myy",8x,"mzz",8x,"mxy",8x,"mxz",8x,"myz")
           END
