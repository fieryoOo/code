           subroutine unwrapR(dper,per,n,tph,ph,grt,r)
c-    to unwrap phase spectrum and get group time-----
C---equal increment by T
C---d(phi)/d(omega)=d(phi)/dT**2/(-pi2)
	   real*4 ph(200),tph(200),grt(200),per(200)
           integer jump(10)
	   data pi2/6.2831854/,const/50./
C-------------saving original phase temporarily--S
	   do i=1,n
	   ph(i)=tph(i)
	   end do
C-------------saving original phase temporarily--E
cc--------------------frequency loop--------S
1	   do i=1,n-1
	   dp=ph(i+1)-ph(i)
C----------------unwrapping from PI2----S
	   if(abs(dp).gt.pi2/r) THEN
            ph(i+1)=ph(i+1)-sign(pi2,dp)
c            PRint*,dp,ph(i+1),ph(i),i
                                ENDIF
           enddo
C----------------unwrapping from PI2----E
C--------------check for phase PI jump-----S
         kj=0
		do i=2,n
        diff=ph(i)-ph(i-1)
        if(abs(diff).gt.3.00)then
c                  PRint*,' jump! ',per(i)
        kj=kj+1
        jump(kj)=i
                              endif
               enddo
C--------------check for phase PI jump-----E
C---------------source group time-----S
C-----q=d(T)/d(omega)-----------
             do i=1,n
        q=-per(i)**2/pi2
        if(i.gt.1.and.i.lt.n) grt(i)=(ph(i+1)-ph(i-1))/2./dper*q 
        if(i.eq.1)grt(1)=(ph(2)-ph(1))/dper*q                        
        if(i.eq.n)grt(n)=(ph(n)-ph(n-1))/dper*q                        
               enddo
              m=0
          DO k=1,n
          qqq=grt(k)
          if(abs(grt(k)).gt.const.and.kj.ne.0)then
             do l=1,kj
                if(jump(l).eq.k) then
                   grt(k)=sign(const,qqq)
                endif
             enddo
              m=m+1
                                  endif
          endDO
C---------------source group time-----E
                return
 	end
