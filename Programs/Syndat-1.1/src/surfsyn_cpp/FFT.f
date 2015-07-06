      subroutine FFT(r,re,im,isi)
c-----------------------------------------------------------------------
      integer r,n,isi
      real*4 re(8),im(8)
      data pi2/6.2831852/
c-----------------------------------------------------------------------
      n=2**r
      lim1=n-1
      lim2=n/2
      j = 1
      do 3 i=1,lim1
      if (i.lt.j)            then
      a1=re(j)
      a2=im(j)
      re(j)=re(i)
      im(j)=im(i)
      re(i)=a1
      im(i)=a2
                            end if
      l = lim2
 2    if(l.ge.j) go to 3
      j = j-l
      l = l/2
      go to 2
3     							j = j+l
      do 4 i=1,r
      lim1 = 2**(i-1)
      lim2 = 2**(r-i)
      arg=-pi2*lim2*isi/float(n)
      cs = 1.0
      si = 0.0
      cstep = cos(arg)
      sstep = sin(arg)
      do 4 m=1,lim1
      do 5 l=1,lim2
      lim3 = (m-1)+(l-1)*2*lim1+1
      b1=re(lim3)
      b2=im(lim3)
      c1=re(lim3+lim1)
      c2=im(lim3+lim1)
      a1 = c1*cs+c2*si
      a2 = -c1*si+c2*cs
      re(lim3)=b1+a1
      im(lim3)=b2+a2
      re(lim3+lim1)=b1-a1
      im(lim3+lim1)=b2-a2
 5    continue
      cs1 = cs*cstep-si*sstep
      si1 = si*cstep+cs*sstep
      cs = cs1
      si = si1
 4    continue
c     if (isi) 11,25,20
      if (isi.lt.0) goto 11
      if (isi.eq.0) goto 25
      w=1./float(n)
      do 10 i=1,n
      re(i)=re(i)*w
      im(i)=im(i)*w
   10 continue
   11 return
   25 stop
      end
