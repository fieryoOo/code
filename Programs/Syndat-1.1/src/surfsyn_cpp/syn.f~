             subroutine syn(wv,dist,dt,f0,df,sre,sim,u_int,q_int,seism,tstart,
     + vmax,n2pow,npoints,key_compr,fix_vel,amp,T_curr,tmin,tmax,k_spec)
C------------to calculate seismogram----------------------
C------------INPUT ARGUMENTS------------------------------
C---n2pow: power of 2; dist:distance; dt: time increment;vmax: max.gr.vel;
C---f0: starting frequency; df -frequency increment; npoints-n. of output points
C-----------wv: wavenumbers; u_int: gr.vel; q_int:apparent surface wave Q;
C-----------sre, sim: spectrum--------------------------------
C------------OUTPUT: seism------------------------------------
         real*4 sre(2048),sim(2048),seism(8192),wv(2048),q_int(2048),u_int(2048)
         real*4 asre(8192),asim(8192),amp(2048),T_curr(2048)
         logical key_compr
             data pi2/6.2831854/
             nbase=2**n2pow
             n=nbase/2
             tstart=dist/vmax
             do i=2,n
             f_curr=f0+df*(i-1)
             om_curr=pi2*f_curr
             arg=wv(i)*dist
             if(key_compr)arg=om_curr*dist/fix_vel
C-----dnom: geometrical spreading----
             dnom=sqrt(arg)
C---- aq: phase delay----------------
             aq=arg+pi2/8.-tstart*om_curr
cMB          aq=pi2/8.-tstart*om_curr
C-----att: attenuation factor--------
             power=dist*om_curr/u_int(i)/q_int(i)/2.
             att=0.0
             if(power.lt.20.)  att=exp(-power)
C-----full spectrum=source spectrum*propagation factor----S
             cs=cos(aq)/dnom*att
             sc=sin(-aq)/dnom*att
             sr=2.*sre(i)/dt
             si=2.*sim(i)/dt
C            sr=sre(i)/dt
C            si=sim(i)/dt
             asre(i)=sr*cs-si*sc
             asim(i)=sr*sc+si*cs
C            k=nbase-i+2
C            asre(k)=asre(i)
C            asim(k)=-asim(i)
             k=n+i
             asre(k)=0.0    
             asim(k)=0.0     
             enddo
             asre(1)=0.0
             asim(1)=0.0
             asre(n+1)=0.0
             asim(n+1)=0.0
C-----full spectrum=source spectrum*propagation factor----S
C----------seismogram outputting----------S
           k_spec=0
           do kkk=2, n
           T_c=1./(f0+df*(kkk-1))
           if(T_c.lt.tmin.or.T_c.gt.tmax)go to 2222 
           k_spec=k_spec+1
           T_curr(k_spec)=T_c
           amp(k_spec)=sqrt(asre(k_spec)**2+asim(k_spec)**2)
2222       continue
           enddo
              call FFT(n2pow,asre,asim,1)
              do i=1,npoints
              seism(i)=asre(i)
              end do
C----------seismogram outputting----------E
              return
              end
