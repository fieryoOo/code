      subroutine cal_synsac( ista, its, sigR, sigL, cor, TMIN, TMAX, vmax, n2pow, fix_vel, iq,
     +                       npoints, fr, df, nf, dt, nt, elatc, elonc, key_compr, qR, qL,
     +                       im, aM, tm, ampl, cr, ul, ur, wvl, wvr, v, dvdz, ratio, I0,
     +                       sfici, slami, seismz, seismn, seisme )
         implicit none
         logical key_compr
         integer*4 j, ista, nsize, nt, nf, npoints, iq
         parameter (nsize=8192)
         integer*4 im, ierr, k_spec, kl, l, ll, lq, m, mmm, nk7, ndist, n2pow
         real*4 cs, cs_b, DEL, DELS, df, dt, dist, fi, fi_back, sc, sc_b, aM
         real*4 slat, slon, sfici, slami, tstart, w, TMIN, TMAX, vmax, fix_vel
         real*4 elat,elatc,elon,elonc
         real*4 azi, azi_back, f0
         real*4 sreV(nsize),simV(nsize)
         real*4 sreR(nsize),simR(nsize)
         real*4 sreT(nsize),simT(nsize)
         real*4 vu(3),du(3),wvn(2)
         real*4 fr(2000), wvar(2000), qR(2000), qL(2000), ampl(2000)
         real*4 amp(2048), T_curr(2048)
         real*4 cr(2000),ul(2000),ur(2000),wvr(2000),wvl(2000),ratio(2000),I0(2000)
         real*4 v(3,2000), dvdz(3,2000), cor(500,2,2000), tm(6)
         real*4 wr(nsize),wl(nsize)
         real*4 qL_int(nsize), qR_int(nsize)
         real*4 ul_int(nsize),ur_int(nsize)
         real*4 seismz(nsize),seismn(nsize),seisme(nsize)
         real*4 seismt(nsize),seismr(nsize)
         real*4 const1,pi,drad
         complex*8 al(2000),az(2000),ah(2000)
         complex*8 br(6),bl(6),sumr,suml,step
         character*1 its, sigR, sigL, chnprefix
         character*3 chn
         character*80 outseism(50),outtit(50),outspec(50)
C         character*255 current
         character*256 fnam11,fnam12,fnam13
         data const1/1.E-7/, step/(1.0,0.0)/, pi/3.14159265/
C         common /trk/ cor
C         common /stn/ nstai,codi,figi,fici,lami,neti

         drad = pi/180.
C         write(*,*) "Station: ",codi(ista)
         slat=sfici*drad
         slon=slami*drad
         call AZIDL(SLAT,SLON,ELATC,ELONC,DEL,DELS,AZI_BACK,AZI)
         dist=dels
         fi=AZI/drad
         fi_back=AZI_BACK/drad
         write(*,*) ' DIST= ',dist,' AZI= ',fi,' BACK_AZI= ',fi_back
c ---
         cs=cos(AZI)
         sc=sin(AZI)
         cs_b=cos(AZI_BACK+pi)
         sc_b=sin(AZI_BACK+pi)
C-----------Reading of spectral data from SURF_DEEP output-----S
      write(*,*) "check cor: ", cor(5, 0, 120), cor(80, 1, 22)
      write(*,*) "check qR: ", qR(33), qR(129)
      write(*,*) "check wvr: ", wvr(33), wvr(129)
      STOP "debug"
         do j=2,nt+1
C            fr(j)=1./t(j)
            wvar(j) = pi*2.0/cor(j-1,1,ista)*fr(j)
C           PRint*,t(j),ampr(j),ratio(j),qR(j),wvr(j)
         enddo
         wvar(1) = wvr(1)
         wvar(nt+2) = wvr(nt+2)
         fr(1)=10.0
         fr(nt+2)=0.0
         qR(1)=qR(2)
         qL(1)=qL(2)
         qR(nt+2)=qR(nt+1)
         qL(nt+2)=qL(nt+1)

C----------CALCULATION OF SPECTRA--------------------------------------S
         write(*,*) 'NT ======',nt
         do j=2,nt+1
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
cMB         spread=1./sqrt(R0*sin(DEL)*wvn(1))
C           spread=1./sqrt(R0*sin(DEL))
c           if(j.eq.2)PRint*,' SPREAD=',spread*sqrt(wvn(1)), ita
cMB         az(j)=sumr*ampr(j)*aM*const1*conjg(cor(j-1))
            az(j)=sumr*aM*const1*cor(j-1,2,ista)/
     +            (2.0*sqrt(cr(j)*ur(j)*I0(j)))/sqrt(2.0*pi)
C     +          (2.0*sqrt(cr(j)*ur(j)*I0(j)))/sqrt(2.0*pi)
cMB         write(*,*)t(j),cr(j),ur(j),wvr(j),ampr(j),cor(j-1,1),cor(j-1,2)
cMB         az(j)=sumr*ampr(j)*aM*const1
cxx         WRite(33,*)t(j),cabs(az(j))*spread
            ah(j)=az(j)*ratio(j)
         end do

C----------INTERPOLATION OF WAVENUMBERS AND Q-factors--S        
         f0=0.0 
         if (sigR.eq.'+')then
cMB         call intpol(fr,wvr,nt+2,f0,df,nf,wr,ierr)
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
C----------INTERPOLATION OF SPECTRAL AMPLITUDES-------------S  
         if (sigR.eq.'+')
     +      call calspecr(nt,nf,f0,df,fr,az,ah,sreV,simV,sreR,simR,TMIN,TMAX,iq)
         if (sigL.eq.'+') 
     +      call calspecl(nt,nf,f0,df,fr,al,sreT,simT,TMIN,TMAX,iq)

C----------CALCULATION OF SPECTRA--------------------------------------E

C----------CALCULATION OF SEISMOGRAMS--------------------S
         ndist=1
C         call namer(dels,ndist,n1,current,outseism,outtit,outspec,sigR,sigL,nk7)
C----------output filenames are defined--------->>>>>>
C----------DISTANCE LOOP-----------------------S      
         do l=1,ndist
            lq=(l-1)*3
C--------------calculating everything and writing spectra----S
            if (sigR.eq.'+') then
               call syn(wr,dist,dt,f0,df,sreV,simV,ur_int,qR_int,seismz,tstart,
     +         vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
C               open(41,file=outspec(lq+1)(1:nk7+1))
C               do kl=2,k_spec
C                  write(41,*)T_curr(kl),amp(kl)
C               enddo
C               close(41)
               call syn(wr,dist,dt,f0,df,sreR,simR,ur_int,qR_int,seismr,tstart,
     +                  vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
C               open(41,file=outspec(lq+2)(1:nk7+1))
C               do kl=2,k_spec
C                  write(41,*)T_curr(kl),amp(kl)
C               enddo
C               close(41)
C               fnam11='w/'//codi(ista)(1:lnblnk(codi(ista)))//'.'//outseism(lq+1)(1:nk7)
C               open(11,file=fnam11)
            elseif (sigL.eq.'+') then
               call syn(wl,dist,dt,f0,df,sreT,simT,ul_int,qL_int,seismt,tstart,
     +                  vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
C               open(41,file=outspec(lq+3)(1:nk7+1))
C               do kl=2,k_spec
C                  write(41,*)T_curr(kl),amp(kl)
C               enddo
C               close(41)
            endif
C--------------calculating everything and writing spectra----E
C--------------writing seismograms----S
C            fnam12='w/'//codi(ista)(1:lnblnk(codi(ista)))//'.'//outseism(lq+2)(1:nk7)
C            open(12,file=fnam12)
C            fnam13='w/'//codi(ista)(1:lnblnk(codi(ista)))//'.'//outseism(lq+3)(1:nk7)
C            open(13,file=fnam13)
            do ll=1,npoints
               seismn(ll)=0.0
               seisme(ll)=0.0
               if (sigR.eq.'+') then
C-----------Positive 'Z' direction is 'Down'-----!!!!!!!!!!!!!!!!!!!!!!
C-----------horizontal components----S
C                 seismz(ll)=-seismz(ll)
                  seismn(ll)=seismr(ll)*cs_b
                  seisme(ll)=seismr(ll)*sc_b
               elseif (sigL.eq.'+') then
                  seismn(ll)=seismn(ll) -seismt(ll)*sc_b
                  seisme(ll)=seisme(ll) +seismt(ll)*cs_b
               endif
C-----------horizontal components----E
C-----------vertical   component-----S
C               if (sigR.eq.'+')write(11,*)tstart+dt*(ll-1),seismz(ll)
C               write(12,*)tstart+dt*(ll-1),seismn(ll)
C               write(13,*)tstart+dt*(ll-1),seisme(ll)
C-----------vertical   component-----E
            enddo
C--------------writing seismograms----E
C            close(11)
C            close(12)
C            close(13)
         end do
      end
