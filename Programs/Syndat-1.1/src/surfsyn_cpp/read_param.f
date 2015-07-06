      subroutine read_params( infile, runfile, mode,  model, elat, elon, aM, dt,
     +                        tm, im, iq, its, nd, npoints, nt, TMIN, TMAX, vmax,
     +                        nstai, codi, figi, fici, lami, neti )
      implicit none

      character*1 its
      character*10 bred
      character*80 str
      character*255 infile, runfile, stafile
      character*256 model
      integer*4 mode, i, ii, k, m, im, iq, lin, nd, npoints, nt
      integer nper(20)
      real*4 strike,dip,slip,per1,per2,dper
      real*4 elat,elon,aM,dt,tm(6)
      real*4 PERMAX,PERMIN
      real*4 pi, drad, GEO
      real*4 TMIN,TMAX,vmax
      data pi/3.14159265/, GEO/0.993277/

      character*8 codi(2000), neti(2000)
      integer*4   nstai
      real*4      figi(2000),fici(2000),lami(2000)
C      common /stn/ nstai,codi,figi,fici,lami,neti

      drad = pi/180.
      lin=lnblnk(infile)
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
9988  print*,k-1,' modes in input'
      close(1)
      nt=nper(mode)
      dper=per2-per1
      PERMAX=(nt-1)*dper+per1
      PERMIN=per1
      write(*,*) nt,' periods for the mode ',mode,
     +              ' per_incr=',dper,' PERMAX=',PERMAX
      if(nt.gt.1998)stop 'too many periods!'

C-----------Reading param file--------------S
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
      open(2,file=runfile,STATUS='OLD')
c --- read parameter file  ---
      read(2,'(a200)')model
      read(2,'(a200)')stafile
c --  read station list  ---
      open(20,file=stafile,status='old')
      nstai = 1
    1 read(20,*,end=98) neti(nstai),codi(nstai),figi(nstai),lami(nstai)
      fici(nstai) = atan(GEO*tan(drad*figi(nstai)))/drad
      do i =1,lnblnk(codi(nstai))
         if(codi(nstai)(i:i).eq.' ') codi(nstai)(i:i)='0'
      enddo
      nstai = nstai+1
      goto 1
   98 nstai = nstai-1
      close(20)
      read(2,*)nd,npoints,dt,TMIN,TMAX,iq,vmax
      if(TMAX.gt.(nt-1)*dper+per1)TMAX=PERMAX-2.*dper
      if(TMIN.lt.per1)TMIN=PERMIN+2.*dper
      print*,' ndepth= ',nd
      print*,'npoints= ',npoints, ' dt = ',dt
      print*,'Tmin=',TMIN,' Tmax=',TMAX
      read(2,'(a)') its
      if(its.ne.'Q'.and.its.ne.'E'.and.its.ne.'F'.and.its.ne.'C')
     +STOP ' WRONG TYPE OF SOURCE: Q,E or F are wanted'
      print 1200 ,          iq,vmax,its 
1200  format(' iq=',i4,' vmax=',F10.3'    type=',a1)
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
      write(*,*) ' ELAT=',elat,' ELON=',elon
C      nstai = 1
C    3 read(2,'(a80)',end=4) str
C      if(str(1:1).eq.' ') then
C         read(str,*) figi(nstai),lami(nstai)
C         fici(nstai) = atan(GEO*tan(drad*figi(nstai)))/drad
C         write(codi(nstai),'("Z",i04)') nstai 
C         neti(nstai) = 'ZZ'
C         do i =1,lnblnk(codi(nstai))
C            if(codi(nstai)(i:i).eq.' ') codi(nstai)(i:i)='0'
C         enddo
C         write(*,*) nstai,' STA=',codi(nstai),' NET=',neti(nstai),
C     +              ' STLAT=',figi(nstai),' STLON=',lami(nstai)
C      else
C         read(str,*) codi(nstai)
C         do ii = 1,nsta
C            if(cod(ii).eq.codi(nstai)) then
C               figi(nstai) = fig(ii)
C               fici(nstai) = atan(GEO*tan(drad*fig(ii)))/drad
C               lami(nstai) = lam(ii)
C               neti(nstai) = net(ii)
C               write(*,*) nstai,' STA=',codi(nstai),' NET=',neti(nstai),
C     +                    ' STLAT=',figi(nstai),' STLON=',lami(nstai)
C               goto 5
C            endif
C         enddo
C         goto 3
C      endif
C    5 nstai = nstai+1
C      goto 3
C    4 nstai = nstai-1
      if(its.eq.'F')write(*,'(5x,"aM",8x,"fx ",8x,"fy ",8x,"fzz")')
      if(its.ne.'F')write(*,1000)
      write(*,'(7e11.3)') aM,(tm(i),i=1,im)
1000  format(5x,"aM",8x,"mxx",8x,"myy",8x,"mzz",8x,"mxy",8x,"mxz",8x,"myz")
      end
