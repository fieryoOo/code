c ==========================================================
c subroutine dmultifft computes cross-correlation function
c            with the double precision accuracy.
c input:
c nlen   - double the length of the input (integer*4)
c n1    - number of samples in arrays f1, (integer*4)
c f1amp - first input arrays amp (real*4)
c f1phase  - first input arrays phase (real*4)
c n2    - number of samples in arrays f2, (integer*4)
c f2amp - second input arraysamp (real*4)
c f2phase- second input arraysphase (real*4)
c lag   - lag the cross-correlation function, (integer*4)
c       - lag < max(n1,n2)/2
c output:
c cor   - cross-correlation array with lag*2+1 samples, (real*4)
c ==========================================================
      subroutine dmultifft(nlen,amp1,pha1,amp2,pha2,seis_out,ns)
      implicit none
      integer*4 nlen, ns
c      integer, intent(in) :: nlen
      real*4    amp1(nlen),pha1(nlen),amp2(nlen),pha2(nlen)
c      real*4    amp1(10000000),pha1(10000000)
      include 'fftw3.h'
      integer*4 i,k
c      real*4    cor(10001)
      real*8    plan3
      double complex temp2(nlen*2),temp3(nlen*2),czero
      real*4  seis_out(nlen*2)
c      common    /core/amp1,pha1

c read amp and phase into complex temp2
c multiply to get correlation in freq domain
 
      do i = 1,nlen
         temp3(i) = amp1(i)*amp2(i)*
     *              exp((0.0d0,-1.0d0)*(pha2(i)-pha1(i)))
      enddo

      czero = (0.0d0,0.0d0)
      ns = (nlen-1)*2
      do k = nlen+1,ns
         temp3(k) = czero
      enddo


c compute ifft 
      call dfftw_plan_dft_1d(plan3,ns,temp3,temp2,
     *                         FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan3)
      call dfftw_destroy_plan(plan3)
c extract correlatio

      do k = 1,ns
        seis_out(k) = real(dreal(temp2(k)))
      enddo

      return

      end
