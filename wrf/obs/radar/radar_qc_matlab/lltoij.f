C                                                                       
C========================================================== K.KURIHARA  
      SUBROUTINE LLTOIJ (FI,FJ, FLAT,FLON,   IM,JM,                     
     \                   NPROJC,DELS,SLON, XI,XJ,XLAT,XLON)             
C     ---- (LAT,LON) TO (I,J) CONVERSION ----                           
      DIMENSION  FI (IM,*),   FJ (IM,*),                                
     \          FLAT(IM,*),  FLON(IM,*)                                 
      CHARACTER * 4 NPROJC                                              
                                                                        
      IF (NPROJC(1:2) .EQ. 'LL') THEN                                   
          RDELS =   1.D0 / DELS                                         
          CYCLEI= 360.D0 / DELS                                         
          DO 10 I=1,IM*JM                                               
                FI(I,1) = XI + RDELS*(FLON(I,1) - XLON)                 
                FJ(I,1) = XJ - RDELS*(FLAT(I,1) - XLAT)                 
 10       CONTINUE                                                      
          DO 20 I=1,IM*JM                                               
                IF (FI(I,1).LT.   0.0) FI(I,1) = FI(I,1) + CYCLEI       
                IF (FI(I,1).GT.CYCLEI) FI(I,1) = FI(I,1) - CYCLEI       
 20       CONTINUE                                                      
      ELSE                                                              
          DO 30 J=1,JM                                                  
          DO 30 I=1,IM                                                  
                CALL RLTLN (FI(I,J),FJ(I,J), FLAT(I,J),FLON(I,J),       
     \                      NPROJC,DELS,SLON, XI,XJ,XLAT,XLON)          
 30       CONTINUE                                                      
      END IF                                                            
                                                                        
C     IF (NPROJC(1:3) .EQ. 'MER') THEN                                  
C         SLAT = 22.5                                                   
C         RA   =6371.D3         !RADIUS OF EARTH!                       
C         PAI  =3.14159D0                                               
C         RAD  = PAI /180.D0                                            
C         CYCLEI= 2.D0*PAI*RA * COS(SLAT*RAD) / DELS                    
C         DO 50 I=1,IM*JM                                               
C               IF (FI(I,1).LT.   0.0) FI(I,1) = FI(I,1) + CYCLEI       
C               IF (FI(I,1).GT.CYCLEI) FI(I,1) = FI(I,1) - CYCLEI       
C50       CONTINUE                                                      
C     END IF                                                            
                                                                        
      RETURN                                                            
      END                                                               
C                                                                       
C============================================================ K.KURIHARA
      SUBROUTINE RLTLN(GI,GJ,RLAT,RLON,                                 
     \                 NPROJC,DELS,SLON, XI,XJ,XLAT,XLON )              
*********************************************************************** 
*    GRID POINT OF LATITUDE AND LONGITUDE                             * 
*                                                                     * 
*OUT-PUT  (GI,GJ)      GRID POINT                                     * 
*                                                                     * 
*IN-PUT   RLAT         LATITUDE   (DEGREE)   -90.<   < 90.            * 
*         RLON         LONGITUDE  (DEGREE)     0.<   <360.            * 
*         NPROJC       MAP PROJECTION                                 * 
*            'PSN ' POLAR STEREO     SLAT= 60 N                       * 
*            'PSS ' POLAR STEREO     SLAT=-60 N                       * 
*            'MER ' MERCATOR         SLAT= 22.5 N                     * 
*            'LMN ' LAMBERT          SLAT= 30 N , 60 N                * 
*            'LMS ' LAMBERT          SLAT=-30 N ,-60 N                * 
*         DELS         GRID PARAM (     M)                            * 
*         SLON         STANDARD LONGITUDE  ; Y-COORDINATE  0.<   <360.* 
*                                                                     * 
*         (XI,XJ) <----------> (XLAT,XLON)                            * 
*                      STANDARD POINT                                 * 
*********************************************************************** 
C                                                                       
      CHARACTER*4 NPROJC 
      DOUBLE PRECISION GI,GJ,RLAT,RLON,DELS,SLON,XI,XJ,XLAT,XLON                                             
C                                                                       
      RA=6371.E3             !RADIUS OF EARTH!                          
      PAI=3.14159                                                       
      RAD=PAI /180.                                                     
C-----------------------------------POLAR STEREO----------------------- 
      IF(NPROJC(1:3).EQ.'PSN') THEN                                     
         SLAT=60.                   !STANDARD LATITUDE ;MAP FACTOR=1!   
C                                                                       
            SLAT1=SLAT*RAD                                              
            SLON1=SLON*RAD                                              
            XLAT1=XLAT*RAD                                              
            XLON1=XLON*RAD                                              
         AL0=RA*(1.+SIN(SLAT1))*TAN(0.5*(0.5*PAI-XLAT1))                
         X0=AL0*SIN(XLON1-SLON1)                                        
         Y0=AL0*COS(XLON1-SLON1)                                        
C                                                                       
            RLAT1=RLAT*RAD                                              
            RLON1=RLON*RAD                                              
         ALI=RA*(1.+SIN(SLAT1))*TAN(0.5*(0.5*PAI-RLAT1))                
         DI=ALI*SIN(RLON1-SLON1)                                        
         DJ=ALI*COS(RLON1-SLON1)                                        
C                                                                       
         GI=(DI-X0)/DELS+XI                                             
         GJ=(DJ-Y0)/DELS+XJ                                             
      END IF                                                            
C                                                                       
      IF(NPROJC(1:3).EQ.'PSS') THEN                                     
         SLAT=-60.                  !STANDARD LATITUDE ;MAP FACTOR=1!   
C                                                                       
            SLAT1=-SLAT*RAD                                             
            SLON1= SLON*RAD                                             
            XLAT1=-XLAT*RAD                                             
            XLON1= XLON*RAD                                             
         AL0=RA*(1.+SIN(SLAT1))*TAN(0.5*(0.5*PAI-XLAT1))                
         X0=AL0*SIN(XLON1-SLON1)                                        
         Y0=AL0*COS(XLON1-SLON1)                                        
C                                                                       
            RLAT1=-RLAT*RAD                                             
            RLON1= RLON*RAD                                             
         ALI=RA*(1.+SIN(SLAT1))*TAN(0.5*(0.5*PAI-RLAT1))                
         DI=ALI*SIN(RLON1-SLON1)                                        
         DJ=ALI*COS(RLON1-SLON1)                                        
C                                                                       
         GI=(DI-X0)/DELS+XI                                             
         GJ=(Y0-DJ)/DELS+XJ                                             
      END IF                                                            
C-----------------------------------MERCATOR--------------------------- 
      IF(NPROJC(1:3).EQ.'MER') THEN                                     
         SLAT=22.5                  !STANDARD LATITUDE ;MAP FACTOR=1!   
C                                                                       
            SLAT1=SLAT*RAD                                              
            SLON1=SLON*RAD                                              
            XLAT1=XLAT*RAD                                              
            XLON1=XLON*RAD                                              
            AC=RA*COS(SLAT1)                                            
C                                                                       
            X0=0.                                                       
            Y0=AC*LOG((1.+SIN(XLAT1))/COS(XLAT1))                       
C                                                                       
            RLAT1=RLAT*RAD                                              
            RLON1=RLON*RAD                                              
C           ----------------------------                                
            DI = RLON - XLON                                            
            DI = MOD (DI+900., 360.) - 180.                             
            DI = AC*(DI*RAD)                                            
C           ----------------------------                                
C ----            DI=AC*(RLON1-XLON1)                                   
            DJ=AC*LOG((1.+SIN(RLAT1))/COS(RLAT1))                       
         GI=(DI-X0)/DELS+XI                                             
         GJ=(Y0-DJ)/DELS+XJ                                             
      END IF                                                            
C-----------------------------------LAMBERT---------------------------- 
      IF(NPROJC(1:3).EQ.'LMN') THEN                                     
         SLATA=30.                  !STANDARD LATITUDE ;MAP FACTOR=1!   
         SLATB=60.                                                      
C                                                                       
         PI4=PAI*0.25                                                   
            SLATA1=SLATA*RAD                                            
            SLATB1=SLATB*RAD                                            
            SLAT1=PI4-SLATA1*0.5                                        
            SLAT2=PI4-SLATB1*0.5                                        
            SLON1=SLON*RAD                                              
            XLAT1=XLAT*RAD                                              
            XLON1=XLON*RAD                                              
         CK=LOG(COS(SLATA1)/COS(SLATB1))/LOG(TAN(SLAT1)/TAN(SLAT2))     
         ACN=RA*COS(SLATA1)/CK                                          
         R0=ACN/(TAN(SLAT1))**CK                                        
                                                                        
         AL0=R0*(TAN(PI4-XLAT1*0.5))**CK                                
C        -----------------------------------                            
         XSLON = XLON-SLON                                              
         XSLON = MOD (XSLON+900., 360.) - 180.                          
         XSLON = XSLON * RAD * CK                                       
         X0=AL0*SIN(XSLON)                                              
         Y0=AL0*COS(XSLON)                                              
C        -----------------------------------                            
C ----        X0=AL0*SIN((XLON1-SLON1)*CK)                              
C ----        Y0=AL0*COS((XLON1-SLON1)*CK)                              
         RLAT1=RLAT*RAD                                                 
         RLON1=RLON*RAD                                                 
!        ALI=R0*(TAN(PI4-RLAT1*0.5))**CK                                
!
!     To prevent abnormal end on UNIX  2001.1.25 K.ONOGI
         ARG=MAX((PI4-RLAT1*0.5),0.0)                                
         ALI=R0*(TAN(ARG))**CK                                
C        -----------------------------------                            
         RXLON = RLON-XLON                                              
         RXLON = MOD (RXLON+900., 360.) - 180.                          
         RXLON = RXLON * RAD * CK                                       
         DI=ALI*SIN(RXLON + XSLON)                                      
         DJ=ALI*COS(RXLON + XSLON)                                      
C        -----------------------------------                            
C ----             DI=ALI*SIN((RLON1-SLON1)*CK)                         
C ----             DJ=ALI*COS((RLON1-SLON1)*CK)                         
         GI=(DI-X0)/DELS+XI                                             
         GJ=(DJ-Y0)/DELS+XJ                                             
      END IF                                                            
C                                                                       
      IF(NPROJC(1:3).EQ.'LMS') THEN                                     
         SLATA=-30.                  !STANDARD LATITUDE ;MAP FACTOR=1!  
         SLATB=-60.                                                     
C                                                                       
         PI4=PAI*0.25                                                   
            SLATA1=-SLATA*RAD                                           
            SLATB1=-SLATB*RAD                                           
            SLAT1=PI4-SLATA1*0.5                                        
            SLAT2=PI4-SLATB1*0.5                                        
            SLON1= SLON*RAD                                             
            XLAT1=-XLAT*RAD                                             
            XLON1= XLON*RAD                                             
         CK=LOG(COS(SLATA1)/COS(SLATB1))/LOG(TAN(SLAT1)/TAN(SLAT2))     
         ACN=RA*COS(SLATA1)/CK                                          
         R0=ACN/(TAN(SLAT1))**CK                                        
         AL0=R0*(TAN(PI4-XLAT1*0.5))**CK                                
C        -----------------------------------                            
         XSLON = XLON-SLON                                              
         XSLON = MOD (XSLON+900., 360.) - 180.                          
         XSLON = XSLON * RAD * CK                                       
         X0 = AL0*SIN(XSLON)                                            
         Y0 = AL0*COS(XSLON)                                            
C        -----------------------------------                            
C---           X0=AL0*SIN((XLON1-SLON1)*CK)                             
C---           Y0=AL0*COS((XLON1-SLON1)*CK)                             
         RLAT1=-RLAT*RAD                                                
         RLON1= RLON*RAD                                                
         ALI=R0*(TAN(PI4-RLAT1*0.5))**CK                                
C        ---------------------------                                    
         RXLON = RLON-XLON                                              
         RXLON = MOD (RXLON+900., 360.) - 180.                          
         RXLON = RXLON * RAD * CK                                       
         DI=ALI*SIN(RXLON + XSLON)                                      
         DJ=ALI*COS(RXLON + XSLON)                                      
C        ---------------------------                                    
C----         DI=ALI*SIN((RLON1-SLON1)*CK)                              
C----         DJ=ALI*COS((RLON1-SLON1)*CK)                              
         GI=(DI-X0)/DELS+XI                                             
         GJ=(Y0-DJ)/DELS+XJ                                             
      END IF                                                            
C---------------------------------------------------------------------- 
      RETURN ; END                                                      
