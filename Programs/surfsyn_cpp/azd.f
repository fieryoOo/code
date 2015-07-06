c####################################################################
      SUBROUTINE AZIDL(SLAT,SLON,ELAT,ELON,DEL,DELS,AZIS,AZIE)
C     AZIS IS AZIMUTH FROM STATION, AZIE - FROM EPICENTER
C     DELS IS THE DISTANCE OVER THE SURFACE OF A SPHERICAL EARTH 
C     ALL IS IN RADIANS 
      DATA R/6371./,PI/3.14159265359/
      COLAT(A)=1.57079632679-A
C------------------------------------------------------------------
      SCOLAT=COLAT(SLAT)
      ECOLAT=COLAT(ELAT)
      C=ELON-SLON
      SC5=COS(C)
      SC6=SIN(C)
      SC1=SIN(SCOLAT) 
      SC2=COS(SCOLAT) 
      SC3=SIN(SLON)
      SC4=COS(SLON)
      EC3=SIN(ELON)
      EC4=COS(ELON)
      EC1=SIN(ECOLAT) 
      EC2=COS(ECOLAT) 
      AE=EC1*EC4
      BE=EC1*EC3
C________AZIMUTHS CALCULATION__________________________________
      AZI1=(AE-SC3)**2+(BE+SC4)**2+EC2*EC2-2.
      AZI2=(AE-SC2*SC4)**2+(BE-SC2*SC3)**2+(EC2+SC1)**2-2.
      IF(AZI2.EQ.0.) GO TO 10
      AZIS=ATAN2(AZI1,AZI2)
      GO TO 20
10    AZIS=3.141592653-SIGN(1.570796326,AZI1)
20    CONTINUE
      AS=SC1*SC4
      BS=SC1*SC3
      AZI1=(AS-EC3)**2+(BS+EC4)**2+SC2*SC2-2.
      AZI2=(AS-EC2*EC4)**2+(BS-EC2*EC3)**2+(SC2+EC1)**2-2.
      IF(AZI2.EQ.0.) GO TO 30
      AZIE=ATAN2(AZI1,AZI2)
      GO TO 40
30    AZIE=3.141592653-SIGN(1.570796326,AZI1)
40    CONTINUE
C__________DELTA CALCULATION_____________________________________
      COSD=SC2*EC2+SC1*EC1*SC5
      DEL= ACOS(COSD)
      DELS=R*DEL
      RETURN
      END 
