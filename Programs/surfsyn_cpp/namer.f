           subroutine namer(dist,ndist,n1,current,outseism,
     *     outtit,outspec,sigR,sigL,nk7)
C----------output files are defined---------------
          real*4 dist(2)
          character*80 current, outseism(50),outtit(50),outspec(50)
          character*5  sg
          character*1  sigR,sigL,sp
C-------------------------------------------------
          sp = sigR
          sp = sigL
          do l=1,ndist
           sg='    '
           ld=int(dist(l))
           if(ld.lt.100)write(sg,'(i2)')ld
           if(ld.ge.100.and.ld.lt.1000)write(sg,'(i3)')ld
           if(ld.ge.1000.and.ld.lt.10000)write(sg,'(i4)')ld
           if(ld.ge.10000)write(sg,'(i5)')ld
           nd=lnblnk(sg)
           lq=(l-1)*3
           nqw=n1+nd-1
           nk2=nqw+2
           nk4=nqw+4
           nk7=nqw+7
           outseism(lq+1)(1:n1)=current(1:n1)
           outseism(lq+2)(1:n1)=current(1:n1)
           outseism(lq+3)(1:n1)=current(1:n1)
           outspec(lq+1)(1:n1)=current(1:n1)
           outspec(lq+2)(1:n1)=current(1:n1)
           outspec(lq+3)(1:n1)=current(1:n1)
           outseism(lq+1)(n1+1:nqw+1)= sg(1:nd)
           outseism(lq+2)(n1+1:nqw+1)= sg(1:nd)
           outseism(lq+3)(n1+1:nqw+1)= sg(1:nd)
           outspec(lq+1)(n1+1:nqw+1)= sg(1:nd)
           outspec(lq+2)(n1+1:nqw+1)= sg(1:nd)
           outspec(lq+3)(n1+1:nqw+1)= sg(1:nd)
           outseism(lq+1)(nk2:nk7)='.Z.dig'
           outseism(lq+2)(nk2:nk7)='.N.dig'
           outseism(lq+3)(nk2:nk7)='.E.dig'
           outspec(lq+1)(nk2:nk7+1)='.Z.spec'
           outspec(lq+2)(nk2:nk7+1)='.R.spec'
           outspec(lq+3)(nk2:nk7+1)='.T.spec'
           outtit(lq+1)(1:nk4)=outseism(lq+1)(1:nk4)
           outtit(lq+2)(1:nk4)=outseism(lq+2)(1:nk4)
           outtit(lq+3)(1:nk4)=outseism(lq+3)(1:nk4)
           outtit(lq+1)(1:nk4)=outseism(lq+1)(1:nk4)
           outtit(lq+2)(1:nk4)=outseism(lq+2)(1:nk4)
           outtit(lq+3)(1:nk4)=outseism(lq+3)(1:nk4)
           outtit(lq+1)(nqw+5:nk7)='tit'
           outtit(lq+2)(nqw+5:nk7)='tit'
           outtit(lq+3)(nqw+5:nk7)='tit'
           end do
           return
           end

