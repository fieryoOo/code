        subroutine calcul(ncount,dx,imax,idispr,idispl,number,numbel,kind,t_base,*)
C--------------------------------------------------------------------
        parameter (nsize=1000,nper=1000)
	common/d/  a(nsize),b(nsize),rho(nsize),d(nsize),qs(nsize)
        common/o/  c(nper,10),t(nper),ratio(nper,10)
        common/c/nmax,mmax,kmax,idrop,iedit,ndiv,mode,fact
        common/newly/ t1,dt,c1,dc,lstop,iq,istru
        common/rar/ depth(nsize),amprfi(nsize),ampz(nsize),strz(nsize),
     *  strrfi(nsize),mmm
        common/rar1/ dcda,dcdb,dcdr,dwx,g1(nsize),g2(nsize)
        common/rco/ ccc,cvar,ugr,wvno,rat,arle
        common/rco1/ sumi0,sumi1,sumi2,sumi3,flagr
        common/deriv/ nderiv,ndpth,dpth(nsize),dspz(nsize),dsprfi(nsize),
     *  drvz(5,nsize),drvrfi(5,nsize)
        common/log/ KEY_ATTEN,KEY_FLAT,KEY_DERIV,KEY_EIGEN,KEY_EIG_NORM,KEY_EIGEN_DER1,KEY_EIGEN_DER2
	common/ref/ a_ref(nsize),b_ref(nsize),rho_ref(nsize),d_ref(nsize),qs_ref(nsize),dpth_ref(1000)
C-------------------------------------------------------------------------
        real*8 dcda(nsize),dcdb(nsize),dcdr(nsize),dwx(nsize)
        logical KEY_DERIV,KEY_ATTEN,KEY_FLAT,KEY_EIGEN,KEY_EIG_NORM
	real*4 ugr_plus(nper),ugr_minus(nper),ugr0(nper)
	real*4 ccc_plus(nper),ccc_minus(nper),ccc0(nper)
        integer imax(10)
C-------------------------------------------------------------------------
	data pi/3.1415927/,R0/6371.0/,small_const/1.E-10/
C-----------------------FORMATS-------------------------------------------
11      format(a,/2(6x,i3,2x),6x,i2,8x,i1)
44      format(1x,/(4f10.4))
55      format(10x,3f10.4)
12      format(2(7x,i1,2x),5x,f3.0)
1000    format(5x,i3,2x,4(5x,f6.0,2x),6x,i1)
1001    format(1x,/(8f10.0))
1002    format(a)
1004    format(34x,'Love mode ',i2)
1005    format(32x,'Rayleigh mode ',i2)
1006    format(8x,i2,9x,i3)
1007    format(1x,/(7f10.0))
5000    format(1x,/1x,'t=',f6.2,2x,'c=',f6.4,2x,'cvar=',f6.4,2x,'ugr=',
     *  f6.4,2x,'wvnumb=',e11.4,2x,'al=',e11.4)
5002    format(1x,/3x,'I0 =',e11.4,5x,'I1 =',e11.4,5x,'I2 =',e11.4,
     *  5x,'L =',e11.4)
5004    format(1x,/4x,'m',5x,'depth',7x,'disp',7x,'stress',8x,
     *  'dc/db',8x,'dc/dr',10x,'g')
5005    format(2x,i3,f10.2,5e13.4)
5006    format(1x,/1x,'t=',f6.2,2x,'c=',f6.4,2x,'cvar=',f6.4,2x,'ugr=',
     *  f6.4,2x,'wvnumb=',e11.4,2x,'ar=',e11.4,/36x,'ratio =',e10.4)
5007    format(1x,'I0 =',e10.4,2x,'I1 =',e10.4,2x,'I2 =',e10.4,
     *  2x,'I3 =',e10.4,2x,'L =',e10.4)
5008    format(1x,/4x,'m',5x,'depth',7x,'dispr',8x,'dispz',
     *  7x,'stresr',7x,'stresz',8x,'dwx')
5009    format(1x,/4x,'m',5x,'depth',7x,'dcda',9x,'dcdb',9x,
     *  'dcdr',10x,'g1',11x,'g2')
5010    format(1x,/10x,
     *  'Eigenfunction V (n=0) and its derivatives with respect to z',
     *  /25x,'(n - the order of derivative)',//2x,'depth\\n',
     *  6x,'0',11x,'1',11x,'2',11x,'3',11x,'4',11x,'5')
5020    format(1x,/10x,
     *  'Eigenfunction Vr (n=0) and its derivatives with respect to z',
     *  /25x,'(n - the order of derivative)',//2x,'depth\\n',
     *  6x,'0',11x,'1',11x,'2',11x,'3',11x,'4',11x,'5')
5030    format(1x,/10x,
     *  'Eigenfunction Vz (n=0) and its derivatives with respect to z',
     *  /25x,'(n - the order of derivative)',//2x,'depth\\n',
     *  6x,'0',11x,'1',11x,'2',11x,'3',11x,'4',11x,'5')
5011    format(1x,f7.2,6e12.4)
3007    format(1x,/2x,'t / mode',3x,'1',6x,'2',6x,'3',6x,'4',6x,'5',
     *  6x,'6',6x,'7',6x,'8',6x,'9',6x,'10',/)
3008    format(1x,f6.2,3x,10f7.3)
3009    format(1x,/11x,'Mode ',i2,',  period',f6.2,
     *  'sec.  Subroutine Nevill. Too many cycles.')
      R0_invsq=1./R0**2
      PRint*, 'KEY_ATTEN=',KEY_ATTEN
      PRint*, 'KEY_FLAT =',KEY_FLAT 
      PRint*, 'KEY_EIGEN=',KEY_EIGEN
      PRint*, 'KEY_EIG_N=',KEY_EIG_NORM
      PRint*, 'KEY_DERIV=',KEY_DERIV
C--------------INITIATION-------------------------------------------
        cc=c1
                     ifunc=kind   
                                        continue
        c1=cc
c
c
        kmode=mode
                                do 2 i=1,10
2       imax(i)=0
                                do 9998 k=1,kmax
        t1=t(k)
	if (KEY_ATTEN) then
C----------------current model for a given t1---------------------
	do i=1,mmax
	b(i)=b_ref(i)*(1.+qs_ref(i)*alog(t_base/t1)/pi)
	qp=qs_ref(i)*1.33333333*b_ref(i)**2/a_ref(i)**2
	a(i)=a_ref(i)*(1.+qp*alog(t_base/t1)/pi)
	rho(i)=rho_ref(i)
	d(i)=d_ref(i)
	enddo

         if(KEY_FLAT)call flat1(d,rho,a,b,d,rho,a,b,mmax,kind)
	endif

c	write(*,*) "now itenerate!~~~~~~ ",i,mmax,t1

                                do 9997 iq=1,kmode
c        PRint*,iq,kmode
                        if(k-1) 605,605,599
599                     if(iq-2) 600,601,601
c600     c1=c(k-1,1)-0.1 
600     c1=0.90*c(k-1,1)
                                                        goto 605
601                     if(dc*(c(k-1,iq)-c(k,iq-1))) 602,603,604
602     c1=c(k,iq-1)+0.01*dc
                                                        goto 605
603     c1=c(k,iq-1)+0.01*dc
                                                        goto 605
604     c1=c(k-1,iq)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       There were c1=c(k-1,iq)-0.1
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
605                                     continue
        idrop=0
999     del1=dltar(c1,t1,ifunc)
80      c2=c1+dc
        idrop=0
        del2=dltar(c2,t1,ifunc)
                    if(sign(1.,del1).ne.sign(1.,del2))  goto 54
81	  c1=c2
        del1=del2
C       PRint*,'TUTA ',c1,b(1),b(mmax)
                        if(c1-0.8*b(1)) 250,251,251
251                     if(c1-(b(mmax)+0.3)) 252,250,250
252                                                     goto 80
54                  continue         
        call nevill (t1,c1,c2,del1,del2,ifunc,cn)
                        if(lstop.eq.0)                  go to 3006
        lstop=0
                        if(k.eq.1.and.iq.eq.1)          go to 3005
        if(ncount.eq.1)write(2,3007)
                        if(k.eq.1)                      go to 3004
                                do 3003 i=1,k-1
                                do 3000 j=1,mode
                        if(i.gt.imax(j))                go to 3001
3000                                    continue
        jmax=mode
                                                        go to 3002
3001    jmax=j-1
3002     continue
3003                                    continue
3004                    if(iq.eq.1)                     go to 3005
        jmax=iq-1
3005    continue 
	if(ncount.eq.1)write(2,*)iq,t(k)
                                                        go to 9999
3006    c1=cn
                        if(c1-b(mmax)) 121,121,250
121     c(k,iq)=c1
                                                goto (6001,6002),ifunc
6002                                    continue
        ratio(k,iq)=dltar(c1,t1,3)
                                                        goto 6003
6001                                    continue
6003                                    continue
        c1=c1+0.01*dc
        imax(iq)=k
                                                        goto 9997
250                     if(k*iq-1) 256,256,255
256     continue 
         if(ncount.eq.1)write(2,258)
258     format(1x,/22x,'Improper initial value. No zero found')
                                                        goto 9999
9997                                    continue
C--------------------search of root is finished--------------------
c	write(*,*) "search of root finished"
                                                        goto 9998
255     kmode=iq-1
                        if(kmode) 9996,9996,9998
9998                                    continue
9996                                    continue
        mmax=nmax
                                do 9995 iq=1,mode
        j=imax(iq)
                        if(j) 9995,9995,9994
9994                                            goto (7001,7002),ifunc
C------------------LOVE     WAVE PART------------------------------S
7001                                    continue
        number=-1
        if(ncount.eq.1)write(2,1004)iq
	print *,'MODE=',iq-1,';  NUMBER OF PERIODS IN OUTPUT=',j
                                do 8889 lip=1,j
	if(KEY_ATTEN)then
	do i=1,mmax
	b(i)=b_ref(i)*(1.+qs_ref(i)*alog(t_base/t(lip))/pi)
	qp=qs_ref(i)*1.33333333*b_ref(i)**2/a_ref(i)**2
	a(i)=a_ref(i)*(1.+qp*alog(t_base/t(lip))/pi)
	d(i)=d_ref(i)
	rho(i)=rho_ref(i)
	qs(i)=qs_ref(i)
	enddo
         if(KEY_FLAT)call flat1(d,rho,a,b,d,rho,a,b,mmax,kind)
	endif
        call leigen(numbel,lip,iq)
        numbel=numbel+1
C------------ATTENUATION OUTPUT---------------S
        qL_app=20000.
	if(KEY_ATTEN) then
	     skd=0.0
	     do i=1,mmm
	     skd=skd+dcdb(i+1)*b(i)*qs(i)
	     enddo
	   alphL=pi/t(lip)*skd/ccc/ccc
	   qL_app= pi/alphL/ugr/t(lip)
C       print '(2(F10.3,6X),10X,F10.3)',t(lip),alphL*1000.,qL_app
           endif
	if(ncount.eq.1)then
	write(2,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        if(.not.KEY_FLAT)write(2,'(6(E14.7,2X))')t(lip),ccc,ugr,wvno,arle,qL_app
        if(KEY_FLAT)write(2,'(6(E14.7,2X))')t(lip),ccc,ugr,sqrt(wvno**2-R0_invsq),arle,qL_app
	write(22,*)t(lip),ugr
	write(21,*)t(lip),ccc,cvar
        write(2,*)sumi0,sumi1,sumi2,flagr
        if(KEY_ATTEN)write(29,*) t(lip),qL_app
	ugr0(lip)=ugr
	ccc0(lip)=ccc
		       endif
        if(ncount.eq.2)then
	ugr_minus(lip)=ugr
	ccc_minus(lip)=ccc
		       endif
	if(ncount.eq.3)then
	   ugr_plus(lip)=ugr
	   ccc_plus(lip)=ccc
		       endif
                        if(iedit.eq.0)                  go to 8889
                        if(iedit.eq.2)                  go to 8890
C       if(ncount.eq.1)write(2,5005)(m,depth(m),amprfi(m),strrfi(m),dcdb(m),dcdr(m),g1(m),m=1,mmm)
                                                        go to 8889
8890    if(ncount.eq.1.and.KEY_EIGEN) then
        IF(KEY_EIG_NORM) then
	bmax=0.0
	do i=1,nmax
        if(abs(amprfi(i)).gt.bmax)bmax=abs(amprfi(i))
	enddo
	do i=1,nmax
	amprfi(i)=amprfi(i)/bmax
	enddo
			 endif
        write (27,*)amprfi(1), z(R0,depth(1),KEY_FLAT),' T(s)=',t(lip)
                                do  m=2,nmax 
        if(amprfi(m).eq.0.0) go to 3333
        write (27,*)amprfi(m),  z(R0,depth(m),KEY_FLAT)
                                   enddo         
        write (27,'(/)'  )
                   endif
3333			 do 5012 m=1,ndpth
        if(abs(dsprfi(m)).gt.small_const)then
        if(.not.KEY_FLAT)write(2,'(5(E14.7,2X))')dpth_ref(m),dsprfi(m),(drvrfi(lm,m),lm=1,nderiv)
	if(KEY_FLAT)write(2,'(5(E14.7,2X))')dpth_ref(m),dsprfi(m)*(1.-dpth_ref(m)/R0),drvrfi(1,m)-dsprfi(m)/R0
                                         endif
5012			 enddo
8889                                    continue
C------------------LOVE     WAVE PART------------------------------E
                                                        go to 3335
C------------------RAYLEIGH WAVE PART------------------------------S
7002                    if(ncount.eq.1)write(2,1005)iq
	print *,'MODE=',iq-1,';  NUMBER OF PERIODS IN  OUTPUT=',j
                                do 8888 lip=1,j
	if(KEY_ATTEN)then
	do i=1,mmax
	b(i)=b_ref(i)*(1.+qs_ref(i)*alog(t_base/t(lip))/pi)
	qp=qs_ref(i)*1.33333333*b_ref(i)**2/a_ref(i)**2
	a(i)=a_ref(i)*(1.+qp*alog(t_base/t(lip))/pi)
	d(i)=d_ref(i)
	rho(i)=rho_ref(i)
	qs(i)=qs_ref(i)
	enddo
         if(KEY_FLAT)call flat1(d,rho,a,b,d,rho,a,b,mmax,kind)
	endif
        call reigen(number,lip,iq)
                        if(number.eq.-1) number=0
        number=number+1
	qR_app=20000.
	if(KEY_ATTEN) then
        skd=0.0
	do i=1,mmm
	skd=skd+dwx(i+1)*qs(i)
	enddo
	alphR=pi/t(lip)*skd/ccc/ccc
	qR_app= pi/alphR/ugr/t(lip)
C	print '(2(F10.3,6X),10X,F10.3))',t(lip),alphR*1000.,qR_app
            endif
	if(ncount.eq.1)then
	write(2,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	if(.not.KEY_FLAT)write(2,'(7(E14.7,2X))')t(lip),ccc,ugr,wvno,arle,rat,qR_app
	if(KEY_FLAT)WRite(2,'(7(E14.7,2X))')t(lip),ccc,ugr,sqrt(wvno**2-R0_invsq),arle,rat,qR_app
	write(22,*)t(lip),ugr
        if(KEY_ATTEN)write(29,*)t(lip),qR_app   
	ugr0(lip)=ugr
	ccc0(lip)=ccc
		end if
	if(ncount.eq.2)then 
	ugr_minus(lip)=ugr
	ccc_minus(lip)=ccc
	         endif
        if(ncount.eq.3)then
	ugr_plus(lip)=ugr
	ccc_plus(lip)=ccc
	         endif
	if(ncount.eq.1)write(21,*)t(lip),ccc,cvar
        if(ncount.eq.1)write(2,*)sumi0,sumi1,sumi2,sumi3,flagr
                        if(iedit.eq.0)                  go to 8888
                        if(iedit.eq.2)                  go to 8891
C       if(ncount.eq.1)write(2,5005)(m,depth(m),amprfi(m),ampz(m),strrfi(m),strz(m),
C    *  dwx(m),m=1,mmm)
C       if(ncount.eq.1)write(2,5005)(m,depth(m),dcda(m),dcdb(m),dcdr(m),g1(m),
C    *  g2(m),m=1,mmm)
C                                                       go to 8888
8891    if(ncount.eq.1.and.KEY_EIGEN) then
        IF(KEY_EIG_NORM) then
	bmax=0.0
	cmax=0.0
	do i=1,nmax
        if(abs(amprfi(i)).gt.bmax)bmax=abs(amprfi(i))
        if(abs(ampz(i)).gt.cmax)cmax=abs(ampz(i))
	enddo
	do i=1,nmax
	amprfi(i)=amprfi(i)/bmax
	ampz(i)=ampz(i)/cmax
	enddo
			 endif
 	write(27,*) amprfi(1),z(R0,depth(1),KEY_FLAT), ' T(s) =',t(lip)
 	write(28,*) ampz(1),  z(R0,depth(1),KEY_FLAT), ' T(s) =',t(lip)
                                do  m=2,nmax 
	if(amprfi(m).eq.0.0.or.ampz(m).eq.0.0) go to 3334
 	write(27,*) amprfi(m),z(R0,depth(m),KEY_FLAT)
 	write(28,*) ampz(m),z(R0,depth(m),KEY_FLAT)
				enddo
3334    write(27,'(/)')
	write(28,'(/)')
                                      endif
        mdpth=ndpth
			        do 5013 m=1,ndpth
       if(abs(dsprfi(m)).gt.0.0. and.abs(dsprfi(m)).lt.small_const)then
        mdpth=m-1
        go to 6666
                                         endif
        if(.not.KEY_FLAT)write(2,'(5(E14.7,2X))')dpth_ref(m),dsprfi(m),(drvrfi(lm,m),lm=1,nderiv)
	if(KEY_FLAT)write(2,'(5(E14.7,2X))')dpth_ref(m),dsprfi(m)*(1.-dpth_ref(m)/R0),drvrfi(1,m)-dsprfi(m)/R0
5013                                    continue
6666    write(2,*)'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
                                do 5014 m=1,mdpth
C       if(abs(dspz(m)).gt.small_const)then
        if(.not.KEY_FLAT)write(2,'(5(E14.7,2X))')dpth_ref(m),dspz(m),(drvz(lm,m),lm=1,nderiv)
	if(KEY_FLAT)write(2,'(5(E14.7,2X))')dpth_ref(m),dspz(m)*(1.-dpth_ref(m)/R0),drvz(1,m)-dspz(m)/R0
C                                        endif
5014                                    continue
8888                                    continue
C------------------RAYLEIGH WAVE PART------------------------------E
3335        if(mode.ne.0)then
            write(21,'(/)')
            write(22,'(/)')
         if(KEY_ATTEN)write(29,'(/)')
                         endif
9995	continue
9999                                    continue
                 if(ncount.lt.3.and.KEY_DERIV)      return1         
777           if(KEY_DERIV) then
	      call qderiv(ugr_minus,ugr0,ugr_plus,t,kmax,dx,25)
            write(25,'(/)')
            write(26,'(/)')
	      call qderiv(ccc_minus,ccc0,ccc_plus,t,kmax,dx,23)
            write(23,'(/)')
            write(24,'(/)')
		   return 
			    endif
                                                                END
C#########################################################################
	    subroutine qderiv(u1,u2,u3,t,nper,dx,nfile)
	    real*4 u1(2),u2(2),u3(2),t(2)
	    do i=1,nper
	    du1=(u3(i)-u1(i))/2./dx
	    du2=(u1(i)+u3(i)-2.*u2(i))/dx**2
	    write(nfile,*)t(i),du1
	    write(nfile+1,*)t(i),du2
	    enddo
	    return
            end
C#########################################################################
            function z(R0,dep,KEY)
C----------to return a  real depth if flattenning was applied
	    logical KEY
	    z=dep
	    if(.not.KEY)return
	    z=R0*(1.-exp(-dep/R0))
	    return
	    end
