           subroutine unwrap(dper,per,n,tph,ph,grt,r)
c-    to unwrap phase spectrum and get group time-----
C---equal increment by T
C---d(phi)/d(omega)=d(phi)/dT**2/(-pi2)
	   real*4 ph(n),tph(n),grt(n),per(n),dper
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
C-Ye------------source group time-----S
      do i=1, n
         if(i.eq.1) then
            dp=ph(2)-ph(1)
         else if(i.eq.n) then
            dp=ph(n)-ph(n-1)
         else 
            dp=(ph(i+1)-ph(i-1))/2.
         endif
         if(abs(dp).gt.3.00) then
            grt(i)=-123456
         else
            q=-per(i)**2/pi2
            grt(i)=dp/dper*q
         endif
      enddo
C-Ye------------source group time-----E
                return
 	end
