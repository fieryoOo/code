       subroutine surfread(infile,sigR,sigL,symb,nt,nd,depth,t,
     +cr,ur,wvr,cl,ul,wvl,v,dvdz,ampr,ampl,ratio,qR,qL,I0)
C------to read SURFLEV output for producing synthetic seismograms--
C----------------INPUT ARGUMENTS-----------------------------------
C------infile -file with spectral eigenfunctions, wvnumbers ,etc---
C------sign is R or L----------------------------------------------
C------nt is number of periods-------------------------------------
C------nd is number of depths for which these functions exist------
C------depth is a wanted depth of a point source-------------------
C---------------OUTPUT ARGUMENTS-----------------------------------
C------c is a phase velocity;u -group velocity;wvn -wavenumber;
C------v(3,2000)---eigenfunction; dvdz(3,2000) - depth derivative; 
C------amp -ampitude factor; ratio is ellipticity------
                 character*80 infile,bred,vzdor,mura
                 character*1 sigR,sigL
                 character*2 symb
                 real*4 v(3,2000),dvdz(3,2000),ratio(2000)
                 real*4 t(2000),ur(2000),cr(2000),wvr(2000),ampr(2000),qR(2000)
                 real*4        ul(2000),cl(2000),wvl(2000),ampl(2000),qL(2000)
                 real*4 I0(2000)
                 data pi2/6.28318/,tlim/10000.0/,eps/1.0/
C-----------------------------------------------------------------
           lsy=lnblnk(symb)
           wvl(1)=0.0
           wvr(1)=0.0
           t(1)=0.0
           t(nt+2)=10000.
           ires=0
         nd_real=nd
C----------Reading Love stuff----------------------S
          if (sigL.eq.'+') THEN
                 open(1,file=infile,STATUS='OLD')
                 do m=1,2000000
                 read(1,'(a)',end =11),bred
                 do j=1,80
          if(bred(j:j+3).eq.'Love')then
                 if(lsy.eq.1.and.bred(j+11:j+11).eq.symb(1:1))go to 1
                 if(lsy.eq.2.and.bred(j+11:j+12).eq.symb(1:2))go to 1
                                          endif
                 end do
                 end do
                 STOP 'NO LOVE'
C----------Reading Love stuff----------------------S
1        Do k=2,nt+1
         if(nd_real.eq.nd)  read(1,'(a)',end=11)vzdor
           read(1, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
           read(1,'(a)') ,vzdor
C---------Love component---S
                  do l=1,nd
           read(1,'(a)')mura
           if(mura(2:4).eq.'@@@')then
           nd_real=l-1 
C          PRint*,'k=',k-1,' nd_real=',nd_real
           goto 5454
                                 endif
           read(mura,'(3(E14.7,2X))',end=11),dep,va,vd
            if(abs(dep-depth).lt.eps) then
           v(3,k)=va
           dvdz(3,k)=vd
                                     end if
                           end do
          nd_real=nd
C---------Love component---E
5454      continue
                                   enddo
           wvl(1)=pi2/cl(2)/tlim
           wvl(nt+2)=0.0
           ampl(nt+2)=0.0
           ires=1
                      END IF
C----------Reading Love stuff----------------------E
11                if(ires.eq.0) print *,'NO LOVE'
                   close(1)
C----------Reading Rayleigh stuff----------------------S
                 if (sigR.eq.'+') Then
                  open(1,file=infile,STATUS='OLD')
C----------------------Search for Rayleigh--------S
                 do m=1,2000000
                 read(1,'(a)',end =99),bred
                 do j=1,80
                 if(bred(j:j+3).eq.'Rayl') then
                 if(lsy.eq.1.and.bred(j+15:j+15).eq.symb(1:1))go to 2
                 if(lsy.eq.2.and.bred(j+15:j+16).eq.symb(1:2))go to 2
                                           endif
                 end do
                 end do
                 STOP 'NO RAYLEIGH'
C----------------------------------PERIOD LOOP------S
2          DO k=2,nt+1
           nd_real=nd
           read(1,'(a)',end=9797)vzdor
           read(1, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
c          read(1,'(a)') ,vzdor
           read(1,'(E14.7)') I0(k)
C----------Rayl. Horizontal component------S
                       do l=1,nd
           read(1,'(a)')mura
           if(mura(2:4).eq.'$$$')then 
           nd_real=l-1
C          PRint*,'k=',k-1,' nd_real=',nd_real,' nd=',nd
           go to 5554
                                 endif
           read(mura,'(3(E14.7,2X))'),dep,va,vd
           if(abs(dep-depth).lt.eps) then
           v(1,k)=va
           dvdz(1,k)=vd
                                     end if
                               end do
C----------Rayl. Horizontal component------E
5554       if(nd_real.eq.nd)read(1,'(a)')vzdor
C----------Rayl. Vertical   component------S
                       do l=1,nd_real
           read(1,*),dep,va,vd
           if(abs(dep-depth).lt.eps) then
           v(2,k)=va
           dvdz(2,k)=vd
                                     end if
                               end do
C----------Rayl. Vertical   component------E
           wvr(1)=pi2/cr(2)/tlim
C          PRint*,k,t(k),v(1,k),v(2,k),dvdz(1,k),dvdz(2,k),wvr(1)
           wvr(nt+2)=0.0
           ampr(nt+2)=0.0
                   end DO
C----------------------------------PERIOD LOOP------E
C----------Reading Rayleigh stuff----------------------E
9797       close(1)
C          nt=k-1
           ires=ires+1
                      End if
99         if(ires.gt.0)  return
           STOP 'NO RAYLEIGH OR LOVE'
            END
