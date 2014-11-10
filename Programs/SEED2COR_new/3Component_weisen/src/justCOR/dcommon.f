c ==========================================================
c subroutine dcommon computes cross-correlation function
c            with the double precision accuracy.
c input:
c n1    - number of samples in arrays f1, (integer*4)
c len   - double the maximum length of either input array(integer*4)
c amp   - real part of the input of FFT (real*4)
c phase - imag part of the input of FFT (real*4)
c ==========================================================
      subroutine dcommon(len, amp, phase)
      implicit  none
      integer*4 i,len
      real*4    amp(400000),phase(400000)
      double    complex temp1(400000)
      common    /core/temp1  

c read amp and phase files into complex double-temp1

      do i = 1,len

          temp1(i) = dcmplx(amp(i)*cos(phase(i)),
     1                      amp(i)*sin(phase(i)))
      enddo
      return
      end
