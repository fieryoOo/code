           subroutine calspecr(nt,nf,f0,df,fr,az,ah,sre1,sim1,
     + sre2,sim2,tmin,tmax,iq)
C-------to produce interpolated Rayleigh wave spectra-------
        parameter (nsize=8192)
        real*4 fr(2000),sre1(nsize),sim1(nsize),ar(2000),ai(2000) 
        real*4 sre2(nsize),sim2(nsize),s(nsize)
        complex*8 az(2000),ah(2000),rc
C-----------------------------------------------------------
        n_end=int(fr(2)/df+0.5)
        n_beg=int(fr(nt+1)/df-0.5)+1
        nn=int(1./tmax/df+0.5)
        nnn=int(1./tmin/df+0.5)
        nw1=nn-n_beg
        nw2=n_end-nnn
        PRint*,'nf=',nf
        print *,'tapering window:',1./n_beg/df,1./(nn-1)/df,'   ',1./(nnn-1)/df,1./(n_end-1)/df
        call tapwin1(nf,nw1,nw2,n_beg,n_end,iq,s)
C-------to introduce a phase shift pi/2 for horizontal component----
        do i=1,nt+2
        rc=ah(i)
        ar(i)=-aimag(rc)
        ai(i)=real(rc)
        end do
C------------------INTERPOLATION ------------------------------
        call intpol(fr,ar,nt+2,f0,df,nf,sre2,ierr)
        call intpol(fr,ai,nt+2,f0,df,nf,sim2,ierr)
        do i=1,nt+2
        rc=az(i)
        ar(i)=real(rc)
        ai(i)=aimag(rc)
        end do
        call intpol(fr,ar,nt+2,f0,df,nf,sre1,ierr)
        call intpol(fr,ai,nt+2,f0,df,nf,sim1,ierr)
        do i=1,nf
        sre1(i)=sre1(i)*s(i)
        sre2(i)=sre2(i)*s(i)
        sim1(i)=sim1(i)*s(i)
        sim2(i)=sim2(i)*s(i)
        end do
        return 
        end
