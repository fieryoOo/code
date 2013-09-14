c****** CELLULAR PAAMETERS, (anisotr. part) ***********************
      integer*4 NCELLA
c     parameter (NCELLA=64000)
      parameter (NCELLA=700000)
c---------------COMMON /mata/---------------------------------------
      integer*4    nnda,nda
      real*8       dta(601),dra(602)
      common /mata/ nnda,nda,dta,dra
c---------------COMMON /wina/---------------------------------------
      integer*4    n0a,m0a,nnna,mmma
      common /wina/ n0a,m0a,nnna,mmma
c---------------COMMON /trasa/--------------------------------------
      integer*4     nalla,nma,icrpnta,ilocra,ioutra
      real*8        dpxa
      common /trasa/ nalla,nma,icrpnta(NCELLA),ilocra(NCELLA),
     +ioutra(NCELLA),dpxa(4,NCELLA)
c*******************************************************************
