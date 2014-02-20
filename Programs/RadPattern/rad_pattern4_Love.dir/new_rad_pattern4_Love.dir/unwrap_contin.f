C#########################################################
           subroutine unwrap(dper,per,n,tph,unph,grt,r,jkl)
c-    to unwrap phase spectrum and get group time-----
C---equal increment by T
C---d(phi)/d(omega)=d(phi)/dT**2/(-pi2)
	   real*4 ph(200),tph(200),unph(200),grt(200),per(200)
	   data pi2/6.2831854/,const/25./
           r=1.2
	   do i=1,n
	   ph(i)=tph(i)
           write(19,*),per(i),ph(i)
	   end do
cc--------------------frequency loop--------S
1          m=0
 	   do i=1,n-1
	   dp=abs(ph(i+1)-ph(i))
C----------------unwrapping from PI2----S
	   if(dp.gt.pi2/r) THEN
           m=m+1
	        if(ph(i+1).gt.ph(i)) then
		ph(i+1)=ph(i+1)-pi2
	        do j=i+2,n
	        ph(j)=ph(j)-pi2
	        end do
C----------------unwrapping from PI2----S
	    else 
	        ph(i+1)=ph(i+1)+pi2
	        do j=i+2,n
	        ph(j)=ph(j)+pi2
	        end do
		end if
		    END IF
C----------------unwrapping from PI2----S
		end do
            if(m.ne.0)goto 1
		do i=1,n  
		unph(i)=ph(i)
		end do
cc--------------------frequency loop--------S
C--------------check for phase PI jump-----S
		do i=2,n
        diff=unph(i)-unph(i-1)
        if(abs(diff).gt.3.00) then
         unph(i)=unph(i)-pi2/2.*sign(1,diff)
C        PRint*,' jump! ',per(i)
                              endif
               enddo
C--------------check for phase PI jump-----E
C---------------source group time-----S
               do i=1,n
C-----q=d(T)/d(omega)-----------
        q=-per(i)**2/pi2
	if(i.gt.1.and.i.lt.n) grt(i)=(unph(i+1)-unph(i-1))/2./dper*q              
	if(i.eq.1)grt(1)=(unph(2)-unph(1))/dper*q                        
	if(i.eq.n)grt(n)=(unph(n)-unph(n-1))/dper*q                        
        if(abs(grt(i)).gt.const)grt(i)=const*sign(1,grt(i))
        if(i.eq.1)write(75,*)i,unph(i),grt(i),jkl
        if(i.gt.1)write(75,*)i,unph(i),grt(i)
		end do
        write(75,'(1X)')
C---------------source group time-----E
                return
		end
