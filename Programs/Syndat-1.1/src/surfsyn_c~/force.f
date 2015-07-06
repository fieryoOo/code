            subroutine force(sigR,sigL,cs,sc,v,br,bl)
C--to calculate single force spectral function bl(3),br(3)
C----------INPUT ARGUMENTS-----------------------------------------
C----code is R or L------------------------------------------------
C----v(3) is an eigenfunction, dvdz is its depth derivative dvdz(3)
C---------OUTPUT ARGUMENTS-----------------------------------------
C--
C----------------------------------------------------------------------
            character*1 sigR,sigL
            real*4 v(3)
            complex*8 bl(3),br(3)
            do i=1,3
            bl(i)=(0.0,0.0)
            br(i)=(0.0,0.0)
            end do
            if(sigR.eq.'+') then   
C-----------Rayleigh wave--------------------------------
C-----------xx------------------------------------------- 
            br(1)=cmplx(0.0,-v(1))*cs
C-----------yy------------------------------------------ 
            br(2)=cmplx(0.0,-v(1))*sc
C-----------zz------------------------------------------ 
            br(3)=cmplx(v(2),0.0) 
                end if  
            if(sigL.eq.'+') then
C-----------Love     wave--------------------------------
C-----------xx------------------------------------------- 
            bl(1)=cmplx(0.0,v(3))*sc
C-----------yy------------------------------------------ 
            bl(2)=cmplx(0.0,-v(3))*cs
C-----------zz------------------------------------------ 
            bl(3)=(0.0,0.0)                
                            end if
                    return
                  end
                            
