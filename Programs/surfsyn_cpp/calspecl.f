        subroutine calspecl(nt,nf,f0,df,fr,al,sre,sim,tmin,tmax,iq)
C-------to produce interpolated Love wave spectrum--------
        parameter(nsize=8192)
        real*4 fr(2000),sre(nsize),sim(nsize),ar(2000),ai(2000),s(nsize)
        complex*8 al(2000),rc
C------------------------------------------------------------
        n_end=int(fr(2)/df+0.5)
        n_beg=int(fr(nt+1)/df-0.5)+1
        nn=int(1./tmax/df+0.5)
        nnn=int(1./tmin/df+0.5)
        nw1=nn-n_beg
        nw2=n_end-nnn
        print *,'tapering window:',1./(n_beg-1)/df,1./(nn-1)/df,
     +          '   ',1./(nnn-1)/df,1./(n_end-1)/df
        call tapwin1(nf,nw1,nw2,n_beg,n_end,iq,s)
        do i=1,nt+2
        rc=al(i)
C-------to introduce a phase shift -pi/2 for horizontal component----
        ar(i)=aimag(rc)
        ai(i)=-real(rc)
        end do
        call intpol(fr,ar,nt+2,f0,df,nf,sre,ierr)
        call intpol(fr,ai,nt+2,f0,df,nf,sim,ierr)
        do i=1,nf
        sre(i)=sre(i)*s(i)
        sim(i)=sim(i)*s(i)
        end do
        return 
        end
