c--find eps1 and eps2 or eps3 and eps4 from fast azim. and ampl. of anisotropy
      character fileamp*80,fileang*80,stri1*4,stri2*4
      print*,"2psi or 4psi (2/4)"
      read*,ipsi
      if(ipsi.ne.2.and.ipsi.ne.4)stop "non existing term"
      print*,"xyz file with max amplitude"
      read*,fileamp
      print*,"xyz file with fast azimuth"
      read*,fileang
      open(1,file=fileamp,status="old")
      open(2,file=fileang,status="old")
      do k=1,80
         if(fileamp(k:k).eq." ")goto10
      enddo
 10   if(ipsi.eq.2)then 
         stri1="EPS1"
         stri2="EPS2"
      else
         stri1="EPS3"
         stri2="EPS4"
      endif
      open(91,file=fileamp(1:k-1)//"."//stri1)
      open(92,file=fileamp(1:k-1)//"."//stri2)
 1    read(1,*,end=2)x,y,amp
      read(2,*,end=3)x1,y1,ang
      if(x.ne.x1.or.y.ne.y1)stop "files are not compatible"
      if(ipsi.eq.2)then
         a2_a1=tan(2.*ang)
         a1=amp/(cos(2.*ang)+(a2_a1*sin(2.*ang)))
      else
         a2_a1=tan(4.*ang)
         a1=amp/(cos(4.*ang)+(a2_a1*sin(4.*ang)))
      endif
      a2=a1*a2_a1
      write(91,*)x,y,a1
      write(92,*)x,y,a2
      goto1
 3    print*,"end of file 2 reached"
 2    print*,"end of file 1 reached (OK)"
      close(1)
      close(2)
      close(91)
      close(92)
      end
