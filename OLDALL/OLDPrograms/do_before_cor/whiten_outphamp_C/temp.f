       PROGRAM TEST 
       do i=1,10
	 if (i.eq.5) goto 10
	 write (*,*) i
10     continue
       enddo
       end
