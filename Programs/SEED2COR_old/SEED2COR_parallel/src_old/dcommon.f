c ==========================================================
c subroutine dcommon computes cross-correlation function
c            with the double precision accuracy.
c input:
c n1    - number of samples in arrays f1, (integer*4)
c len   - double the maximum length of either input array(integer*4)
c amp   - real part of the input of FFT (real*4)
c phase - imag part of the input of FFT (real*4)
c ==========================================================
      subroutine dcommon(nlen, samp, spha)
      implicit  none
      integer*4 i,nlen
      real*4    samp(nlen),spha(nlen)
      real*4    amp1(10000000),pha1(10000000)
      common    /core/amp1,pha1

c read amp and phase files into complex double-temp1

      do i = 1,nlen
	  amp1(i) = samp(i)
	  pha1(i) = spha(i)
c          temp1(i) = dcmplx(sreal(i),simag(i))
      enddo
      return
      end
