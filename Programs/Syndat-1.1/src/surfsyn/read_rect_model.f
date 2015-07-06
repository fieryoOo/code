c ==========================================================
c read rectangular models into memory
c ==========================================================
      subroutine read_rect_model(namea,nmod,p,ierr)
      implicit none
      integer*4 nmod,ierr
      real*8    p,geo,drad
      character*255 namea
c --- common /mdl/ ----------------------
      integer*4 ic,jc,n,nper,nfi,nla
      real*8 fi,sfi,la,sla,per,sper,bf,ef,bl,el,bp,ep,
     +       hfi(97),hla(201),hper(225),chfi(97)
      real*8 uw(97,201),cw(97,201)
      real*8 gw(97,201),aw(97,201)
      common /mdl/ic,jc,n,nper,nfi,nla,fi,sfi,la,sla,per,sper,
     +        bf,ef,bl,el,bp,ep,
     +        hfi,hla,hper,chfi,uw,cw,gw,aw
c ---
      integer*4 i
      ierr = 0
c --- read unformated file for period per--
      if(nmod.eq.0) then
        ic = 1
        jc = 1
        geo = 0.993277d0
        drad = datan(1.0d0)/45.0d0
        open(10,file=namea,form='unformatted',status='old')
        read(10) n,fi,nfi,sfi,la,nla,sla,per,nper,sper
        n =0
        do i = 1,nfi
          hfi(i) = fi+(i-1)*sfi
          chfi(i) = datan(geo*dtan(drad*hfi(i)))/drad
        enddo
        bf = chfi(1)
        ef = chfi(nfi)
        do i = 1,nla
          hla(i) = la+(i-1)*sla
        enddo
        bl = hla(1)
        el = hla(nla)
        do i = 1,nper
          hper(i) = nper+(i-1)*sper
        enddo
        bp = hper(1)
        ep = hper(nper)
        return
      else if(nmod.eq.-1) then
        close(10)
        return
      endif
      p = per+n*sper
      read(10) uw,cw,gw,aw
      n = n+1
      end
