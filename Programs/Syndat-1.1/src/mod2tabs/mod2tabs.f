c=========================================================
c create single unformatted files for phase, group
c gamma and for receiver correction.
c Each read statment read spatial table for fixed deriod.
c=========================================================
      implicit none
      real*8 u(97,201,225),c(97,201,225)
      real*8 g(97,201,225),a(97,201,225)
      real*8 uu(97,201),cc(97,201)
      real*8 gg(97,201),aa(97,201)
      common /dat/ u,c,g,a
      integer*4 narg,iargc,lnblnk
      integer*4 i,j,k,n,nfi,nla,nper,ires
      real*8    pi,fi,sfi,la,sla,per,sper,d1,d2,p,dum(2),wvr
      character*20 smod,ssmod
      character*255 list,in_data,in_tabs,param,cmd,str
c ---
      pi = datan(1.0d0)*4.0d0
      narg=iargc()
      if(narg.ne.4) then
        write(*,*) 'USAGE: mod2tab list param in_data in_tabs'
        stop
      endif
c --- read command arguments   ---
      call GETARG(1,list)
      call GETARG(2,param)
      call GETARG(3,in_data)
      call GETARG(4,in_tabs)
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
         do k = 1,20
            ssmod(k:k) = smod(k:k)
            if(smod(k:k).eq.'_') ssmod(k:k) = ' '
         enddo
         read(ssmod,*) d1,d2
         i = (d1-fi)/sfi+1.4d0
         j = (d2-la)/sla+1.4d0
         cmd = 'echo model '//smod//' >> LOG'
         ires = system(cmd)
         cmd = '../../bin/SURF_DISP '//in_data(1:lnblnk(in_data))//'/'//
     +        smod(1:lnblnk(smod))//' R R 0 0 4 60 .25 -a -f >> LOG'
         ires = system(cmd)
         n = n+1
         write(*,'(i6,1x,a15,i3)') n,smod,ires
         open(11,file='R.R.grv',status='old')
         open(12,file='R.R.phv',status='old')
         open(13,file='R.R.att',status='old')
         open(14,file='R.R',status='old')
         do k = 1,nper
           read(11,*)p,u(i,j,k)
           read(12,*)p,c(i,j,k)
           read(13,*)p,g(i,j,k)
           g(i,j,k) = pi/(g(i,j,k)*p*u(i,j,k))
         enddo
         close(11)
         close(12)
         close(13)
         k = 1
    2    read(14,'(a200)',end=99) str
         if(str(2:6).eq.'@@@@@') then
         read(14,*) p,dum,wvr
         read(14,*) a(i,j,k)
         a(i,j,k) = 1.0d-6/dsqrt(a(i,j,k)*c(i,j,k)*u(i,j,k)*wvr)
         k = k+1
         endif
         goto 2
   99    close(14)
         cmd = '/bin/rm -f R.R R.R.att R.R.grv R.R.phv'
         ires = system(cmd)
        goto 1
  9   close(10)
c --- write unformated file --
      open(10,file=in_tabs,form='unformatted',status='unknown')
      write(10) n,fi,nfi,sfi,la,nla,sla,per,nper,sper
      do k = 1,nper
        write(*,*) ' Period: ',per+(k-1)*sper
        do i = 1,nfi
          do j = 1,nla
            uu(i,j) = u(i,j,k)
            cc(i,j) = c(i,j,k)
            gg(i,j) = g(i,j,k)
            aa(i,j) = a(i,j,k)
          enddo
        enddo
        write(10) uu,cc,gg,aa
      enddo
      close(10)
      end
