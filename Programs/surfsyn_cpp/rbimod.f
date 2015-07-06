c****************************************************************
      subroutine rbimod(e,m,vv,dbg,ierr)
      implicit none
      integer*4 m,ierr,dbg
      real*8    e(3),vv(4)
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
      real*8 drad,dff,dfl
      real*8 hf,dhf,hl,dhl,v1,v2,v3,v4
c ---
      drad = datan(1.0d0)/45.0d0
      dfl = datan2(e(2),e(1))/drad
      ierr = 0
      if(dfl.lt.0.d0)dfl = dfl+360.d0
      dff = e(3)
      if(dabs(dff).gt.1.d0) dff = dsign(1.d0,e(3))
      dff = DASIN(dff)/drad
      if(dbg.ne.0) write(*,*) 'Point: ',dff,dfl
c --- check that e is model's internal point
      if(dff.lt.bf.or.dff.gt.ef.or.dfl.lt.bl.or.dfl.gt.el) then
        goto 9
      endif
c --- check that e is cell's internal point
    2 if(dff.lt.chfi(ic)) then
        ic = ic-1
        if(ic.eq.0) then
          ic = 1
          goto 9
        endif
        goto 2
      endif
c ---
    3 if(dff.gt.chfi(ic+1)) then
        ic = ic+1
        if(ic.gt.nfi) then
          ic = nfi
          goto 9
        endif
        goto 3
      endif
c ---
    4 if(dfl.lt.hla(jc)) then
        jc = jc-1
        if(jc.eq.0) then
          jc = 1
          goto 9
        endif
        goto 4
      endif
c ---
    5 if(dfl.gt.hla(jc+1)) then
        jc = jc+1
        if(jc.gt.nla) then
          jc = nla
          goto 9
        endif
        goto 5
      endif
c --- four maps  interpolation ---
      m = 0
      hf = chfi(ic+1)-chfi(ic)
      dhf = dff-chfi(ic)
      hl = hla(jc+1)-hla(jc)
      dhl = dfl-hla(jc)
      v1 = uw(ic,jc)
      v2 = uw(ic+1,jc)
      v3 = uw(ic,jc+1)
      v4 = uw(ic+1,jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(1) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      v1 = cw(ic,jc)
      v2 = cw(ic+1,jc)
      v3 = cw(ic,jc+1)
      v4 = cw(ic+1,jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(2) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      v1 = gw(ic,jc)
      v2 = gw(ic+1,jc)
      v3 = gw(ic,jc+1)
      v4 = gw(ic+1,jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(3) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      v1 = aw(ic,jc)
      v2 = aw(ic+1,jc)
      v3 = aw(ic,jc+1)
      v4 = aw(ic+1,jc+1)
      if(v1.lt.0.0d0.or.v2.lt.0.0d0.or.v3.lt.0.0d0.or.v4.lt.0.0d0) goto 9
      vv(4) = v1+(v2-v1)*dhf/hf+(v3-v1)*dhl/hl+
     *    (v1+v4-v2-v3)*dhf*dhl/hf/hl
      m = m+1
c ---
      return
   9  ierr = 1
      end
