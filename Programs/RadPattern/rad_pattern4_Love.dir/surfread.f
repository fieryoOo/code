       subroutine surfread(infile,sigR,sigL,symb,nt,nd,depth,t,
     +cr,ur,wvr,cl,ul,wvl,v,dvdz,ampl)
C------to read SURFLEV output for producing synthetic seismograms--
C----------------INPUT ARGUMENTS-----------------------------------
C------infile -file with spectral eigenfunctions, wvnumbers ,etc---
C------sign is R or L----------------------------------------------
C------nt is number of periods-------------------------------------
C------nd is number of depths for which these functions exist------
C------depth is a wanted depth of a point source-------------------
C---------------OUTPUT ARGUMENTS-----------------------------------
C------c is a phase velocity;u -group velocity;wvn -wavenumber;
C------v(3,200)---eigenfunction; dvdz(3,200) - depth derivative; 
C------amp -ampitude factor; ratio is ellipticity------
		 character*80 infile,bred,vzdor,mura
		 character*1 sigR,sigL
		 character*2 symb
		 real*4 v(3,200),dvdz(3,200),ratio(200)
		 real*4 t(200),ur(200),cr(200),wvr(200),ampr(200),qR(200)
		 real*4        ul(200),cl(200),wvl(200),ampl(200),qL(200)
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
		 do m=1,200000
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
1         ires=1 
          Do k=1,nt
         if(nd_real.eq.nd)  read(1,'(a)',end=11)vzdor
           read(1, '(6(E14.7,2X))') t(k),cl(k),ul(k),wvl(k),ampl(k),qL(k)
	   read(1,'(a)') ,vzdor
C---------Love component---S
	                  do l=1,nd
           read(1,'(a)',end=11)mura
           if(mura(2:4).eq.'@@@')then
           nd_real=l-1 
           PRint*,k,nd_real
           go to 5454
                                 endif
	   read(mura,'(3(E14.7,2X))',end=11),dep,va,vd
	   if(abs(dep-depth).lt.eps) then
	   v(3,k)=va
	   dvdz(3,k)=vd
				     end if
                           end do
6666      nd_real=nd
C---------Love component---E
5454     			   end Do
C	   ires=1
		   close(1)
		      END IF
C----------Reading Love stuff----------------------E
	          if(ires.eq.0) print *,'NO LOVE'
11         continue
C----------Reading Rayleigh stuff----------------------S
  		 if (sigR.eq.'+') Then
		  open(1,file=infile,STATUS='OLD')
C----------------------Search for Rayleigh--------S
		  do m=1,200000
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
2          DO k=1,nt
           nd_real=nd
           read(1,'(a)',end=9797)vzdor
           read(1, '(7(E14.7,2X))') t(k),cr(k),ur(k),wvr(k),ampr(k),ratio(k),qR(k)
	   read(1,'(a)') ,vzdor
C----------Rayl. Horizontal component------S
	                       do l=1,nd+2
           read(1,'(a)')mura
           if(mura(2:4).eq.'$$$')then 
           nd_real=l-1
           go to 5554
                                 endif
	   read(mura,'(3(E14.7,2X))'),dep,va,vd
	   if(abs(dep-depth).lt.eps) then
	   v(1,k)=va
	   dvdz(1,k)=vd
				     end if
                               end do
C----------Rayl. Horizontal component------E
5554       continue                                        
C----------Rayl. Vertical   component------S
      	               do l=1,nd_real
	   read(1,*),dep,va,vd
	   if(abs(dep-depth).lt.eps) then
	   v(2,k)=va
	   dvdz(2,k)=vd
				     end if
                               end do
C----------Rayl. Vertical   component------E
		   end DO
C----------------------------------PERIOD LOOP------E
C----------Reading Rayleigh stuff----------------------E
9797	   close(1)
           nt=k-1
	   ires=ires+1
		      End if
99         if(ires.gt.0)  return
	   STOP 'NO RAYLEIGH OR LOVE'
            END
