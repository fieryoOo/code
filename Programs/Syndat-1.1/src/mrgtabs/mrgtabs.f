c=========================================================
c create single unformatted files for phase, group
c gamma and for receiver correction.
c Each read statment read spatial table for fixed deriod.
c=========================================================
      implicit none
      real*8 u(97,201,225),c(97,201,225)
      real*8 g(97,201,225),a(97,201,225)
      real*8 uw(97,201),cw(97,201)
      real*8 gw(97,201),aw(97,201)
      integer*4 narg,iargc
      integer*4 i,j,k,n,nn,nfi,nla,nper
      real*8    fi,sfi,la,sla,per,sper
      character*255 list,param,out_tabs,smod
c ---
      narg=iargc()
      if(narg.ne.3) then
        write(*,*) 'USAGE: mrgtab list param out_tab'
        stop
      endif
c --- read command arguments   ---
      call GETARG(1,list)
      call GETARG(2,param)
      call GETARG(3,out_tabs)
      do k = 1,97
        do j = 1,201
          do i = 1,225
            u(k,j,i) = -1.0d0
            c(k,j,i) = -1.0d0
            g(k,j,i) = -1.0d0
            a(k,j,i) = -1.0d0
          enddo
        enddo
      enddo
c ---   read parameters    ---
      open(10,file=param,status='old')
        read(10,*) fi,nfi,sfi,la,nla,sla,per,nper,sper
        write(*,'(3(f6.1,i4,f6.2))') fi,nfi,sfi,la,nla,sla,per,nper,sper
      close(10)
      n = 0
c ---      MAIN ---
      open(10,file=list,status='old')
  1     read(10,'(a20)',end=9) smod
        n = n+1
        open(11,file=smod,form='unformatted',status='old')
        read(11) nn,fi,nfi,sfi,la,nla,sla,per,nper,sper
        n = n+nn
        do k = 1,nper
          read(11) uw,cw,gw,aw
          do i = 1,nfi
            do j = 1,nla
              if(u(i,j,k).lt.0.0d0.and.uw(i,j).gt.0.0d0) u(i,j,k) = uw(i,j)
              if(c(i,j,k).lt.0.0d0.and.cw(i,j).gt.0.0d0) c(i,j,k) = cw(i,j)
              if(g(i,j,k).lt.0.0d0.and.gw(i,j).gt.0.0d0) g(i,j,k) = gw(i,j)
              if(a(i,j,k).lt.0.0d0.and.aw(i,j).gt.0.0d0) a(i,j,k) = aw(i,j)
            enddo
          enddo
        enddo
        close(11)
        goto 1
  9   close(10)
c --- write unformated file --
      open(10,file=out_tabs,form='unformatted',status='unknown')
      write(10) n,fi,nfi,sfi,la,nla,sla,per,nper,sper
      do k = 1,nper
        write(*,*) ' Period: ',per+(k-1)*sper
        do i = 1,nfi
          do j = 1,nla
            uw(i,j) = u(i,j,k)
            cw(i,j) = c(i,j,k)
            gw(i,j) = g(i,j,k)
            aw(i,j) = a(i,j,k)
          enddo
        enddo
        write(10) uw,cw,gw,aw
      enddo
      close(10)
      end
