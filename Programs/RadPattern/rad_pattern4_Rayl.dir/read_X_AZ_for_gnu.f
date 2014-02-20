	     character*80 infile,outname,perlist
	     character*8 symbol         
	     character*2 sym
             real*4,per(200),az(181),T(181,200),period(20)
	     data marg/5/,eps/0.0001/
	     narg=iargc()
             PRINT*,'to extract gnuplot files from T_AZ files'
	     if(narg.ne.marg)STOP'USAGE: read_X_AZ_for_gnu infile outname perlist nper naz '
             call GETARG(1,infile)
             call GETARG(2,outname)
             lin=lnblnk(outname)
             call GETARG(3,perlist)
             call GETARG(4,symbol)
             read(symbol,*)nper
             call GETARG(5,symbol)
             read(symbol,*)naz 
	     open(1,file=infile,status='OLD')
             open(3,file=perlist,status='OLD')
             do i=1,20
                 read(3,*,end=999)period(i) 
             enddo
999          np=i-1
             do j=1,nper
             do i=1,naz
             read(1,*,end=99)az(i),per(j),T(i,j)
             enddo
             enddo 
99           continue 
             do i=1,np
             do k=1,nper
             if(abs(period(i)-per(k)).lt.eps)then
             sym='  '
             iper=per(k)
             if(iper.lt.10)write(sym,'(I1)')iper
             if(iper.ge.10)write(sym,'(I2)')iper
             lys=lnblnk(sym) 
	     open(2,file=outname(1:lin)//'_'//sym(1:lys))
	     do j=1,naz    
	     write(2,*) az(j),T(j,k)
	     enddo
                                              endif
	     enddo
	     enddo
             end
