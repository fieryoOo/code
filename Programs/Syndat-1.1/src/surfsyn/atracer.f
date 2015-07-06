c group velocity tracer for cross-correlation paths
c ---
      subroutine atracer(model,fsol,lsol,n,applyQ)
      implicit none
c ---
      character*8 cod(2000),codi(2000)
      integer*4   nsta,nstai
      real*4      fig(2000),lam(2000),figi(2000),fici(2000),lami(2000)
      common /stn/nsta,nstai,cod,fig,lam,codi,figi,fici,lami
c --- common /mdl/ ----------------------
      integer*4 ic,jc,nn,nper,nfi,nla
      real*8 fi,sfi,la,sla,pper,sper,bf,ef,bl,el,bp,ep,
     +       hfi(97),hla(201),hper(225),chfi(97)
      real*8 uw(97,201),cw(97,201)
      real*8 gw(97,201),aw(97,201)
      common /mdl/ic,jc,nn,nper,nfi,nla,fi,sfi,la,sla,pper,sper,
     +        bf,ef,bl,el,bp,ep,
     +        hfi,hla,hper,chfi,uw,cw,gw,aw
c ---
      real*4 cor(500,2,2000)
      common /trk/ cor
c ---
      character *1 applyQ
      real*4    fsol,lsol
      real*8     afi,del,per
      real *8    GEO,rad,pi2,sol(3),dst(3),trres(4)
      integer*4 ierr,k,ntr,n,m,lnblnk
c     real*8     cor(500,2,2000)
      character *256 model
      data GEO/1.0/
c     data GEO/0.993277d0/
c ---
      rad = datan(1.0d0)/45.0d0
      pi2 = datan(1.0d0)*8.0d0
      write(*,*) model(1:lnblnk(model))
c --- get event coordinates ----
      afi = datan(GEO*dtan(rad*fsol))/rad
      sol(1)=DSIN((90.0d0-afi)*rad)*DCOS(lsol*rad)
      sol(2)=DSIN((90.0d0-afi)*rad)*DSIN(lsol*rad)
      sol(3)=DCOS((90.0d0-afi)*rad)
c --- MAIN LOOP ------------
      ntr = 4
      call read_rect_model(model,0,per,ierr)
      n = 0
      do k = 1,225
        call read_rect_model(model,1,per,ierr)
        do m = 1,nstai
c     write(*,*) per,ierr
c         afi = datan(GEO*dtan(rad*figi(m)))/rad
          afi = fici(m)
          dst(1)=DSIN((90.0d0-afi)*rad)*DCOS(lami(m)*rad)
          dst(2)=DSIN((90.0d0-afi)*rad)*DSIN(lami(m)*rad)
          dst(3)=DCOS((90.0d0-afi)*rad)
          call tracer(sol,dst,.01d0,ntr,trres,del,0,ierr)
cMB    write(*,*) 'X ',per,trres,del
          if(ierr.eq.0) then  
cMB     trres(2) = trres(2)*pi2/per
cMB     cor(k) = dcmplx(dcos(trres(2)),dsin(trres(2)))*dexp(-trres(3))*trres(4)
cMB     write(*,*) 'A ', per,6371.0*del,trres(2)
            cor(k,1,m) = 6371.0*del/trres(2)
            if( applyQ.eq.'Y' ) then
               cor(k,2,m) = dexp(-trres(3))*trres(4)
            else
               cor(k,2,m) = trres(4)
            endif
CCC     write(*,*) m,per, 6371.0*del,cor(k,1,m)
c       write(*,*) per,trres,cor(k)
          else
            cor(k,1,m) = -1
            cor(k,2,m) = 0
          endif
        enddo
      enddo
      n = 225
      call read_rect_model(model,-1,per,ierr)
      end
