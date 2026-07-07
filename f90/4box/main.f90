      PROGRAM MAIN

      USE CHEMEQ
      USE O2
      USE CO2

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Ocean 4-box model
!  Calculation of phosphate, DIC, alkalinity, CO2, and PCO2
!
!  References: Toggweiler et al. (1999, PA)
!
! ---------------------------------------------------------------------

      INTEGER :: I
      INTEGER, PARAMETER :: NBOX=4
!     1:H, 2:L, 3:N, 4:D

      REAL(8) :: DT, DELTA
      REAL(8) :: CV1, CV2, CV3, CV4, CV5

      REAL(8) :: VOCN, AOCN, VATM
      REAL(8) :: ZOCNH, DEH, AOCNH, FAOCNH, VOCNH
      REAL(8) :: ZOCNL, DEL, AOCNL, FAOCNL, VOCNL
      REAL(8) :: ZOCNN, DEN, AOCNN, FAOCNN, VOCNN
      REAL(8) :: ZOCND,                     VOCND

      REAL(8) :: TRAN, FHD, FLD, FND, FHL, FLN
      REAL(8) :: PVH, PVL, PVN
      REAL(8) :: FAH, FAL, FAN
      REAL(8) :: R, LF, HSC
      REAL(8) :: RRC, RCP, RNP, RO2P, CEPH
      REAL(8) :: BT(NBOX)
      REAL(8) :: K0(NBOX), K1(NBOX), K2(NBOX)
      REAL(8) :: KB(NBOX), KW(NBOX), FH(NBOX)

      REAL(8) :: TEMH,  SALH,   TEMHX!  SALHX
      REAL(8) :: PO4H,  ALKH,   DICH,   DO2H
      REAL(8) :: PO4HX, ALKHX,  DICHX,  DO2HX
      REAL(8) :: CO2H,  HCO3H,  CO32H,  PCO2H, O2SATH, AOUH, EPH

      REAL(8) :: TEML,  SALL,   TEMLX!  SALLX
      REAL(8) :: PO4L,  ALKL,   DICL,   DO2L
      REAL(8) :: PO4LX, ALKLX,  DICLX,  DO2LX
      REAL(8) :: CO2L,  HCO3L,  CO32L,  PCO2L, O2SATL, AOUL, EPL

      REAL(8) :: TEMN,  SALN,   TEMNX!  SALNX
      REAL(8) :: PO4N,  ALKN,   DICN,   DO2N
      REAL(8) :: PO4NX, ALKNX,  DICNX,  DO2NX
      REAL(8) :: CO2N,  HCO3N,  CO32N,  PCO2N, O2SATN, AOUN, EPN

      REAL(8) :: TEMD,  SALD,   TEMDX!  SALDX
      REAL(8) :: PO4D,  ALKD,   DICD,   DO2D
      REAL(8) :: PO4DX, ALKDX,  DICDX,  DO2DX
      REAL(8) :: CO2D,  HCO3D,  CO32D,  PCO2D, O2SATD, AOUD

      REAL(8) :: PCO2A, PCO2AX
      REAL(8) :: TOCINI, TOCFIN
      REAL(8) :: OLD_PO4H

      DATA DT    / 8.64D4 /
      DATA DELTA / 1.0D-6 /

      DATA CV1 / 1.0250D+3 / !! [mol/kg] to [mol/m^3]
      DATA CV2 / 9.7561D-4 / !! [mol/m^3] to [mol/kg]
      DATA CV3 / 1.0000D+6 / !! [atm] to [ppmv]
      DATA CV4 / 3.1536D+7 / !! [yr] to [s]
      DATA CV5 / 8.6400D+4 / !! [day] to [s]

!--   Temperature [K]
!--   Salinity [PSU]
      DATA TEMH, TEML, TEMN, TEMD &
     & / 1.0, 21.5, 3.0, 1.75 / !! Table 1,2 in Toggweiler+ (1999)
      DATA SALH, SALL, SALN, SALD &
     & / 34.7, 34.7, 34.7, 34.7 / !! Table 1,2 in Toggweiler+ (1999)

!--   Ocean volume [m^3]
!--   Ocean area [m^2]
!--   Atmospheric molers [mol] = [mol/atm]
      DATA VOCN, AOCN, VATM / 1.292D18, 3.49D14, 1.773D20 / !! Table 1 in Toggweiler+ (1999)

!--   Ocean depth [m]
      DATA ZOCNH, ZOCNL, ZOCNN, ZOCND &
     & / 250.D0, 100.D0, 250.D0, 4000.D0 / !! Table 1,2 in Toggweiler+ (1999)

!--   Ocean depth of euphotic zone [m]
      DATA DEH, DEL, DEN &
     & / 100.D0, 100.D0, 100.D0 /

!--   Ocean area [%]
      DATA FAOCNH, FAOCNL, FAOCNN &
     & / 0.15D0, 0.75D0, 0.10D0 / !! Table 1 in Toggweiler+ (1999)

!--   Transport [m^3/s]
!--   Diffusivity [m^3/s]
      DATA TRAN &
     & / 20.0D6 /
      DATA FHD, FLD, FND, FHL, FLN &
     & / 60.0D6, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / !! 3-300 Sv / Table 1,2 in Toggweiler+ (1999)

!--   Piston velocity [m/day]
      DATA PVH, PVL, PVN / 3.0D0, 3.0D0, 3.0D0 / !! Table 1,2 in Toggweiler+ (1999)

!--   Bio-production efficiency [/yr]
      DATA R / 1.0D0 /

!--   Proportional to an annual averaged solar radiation
      DATA LF / 5.0D-1 /

!--   Half-saturation constant [mol/kg]
      DATA HSC / 2.0D-8 /

!--   Rain ratio, C/P, N/P, O2/P
      DATA RRC  / 0.250D0 / !! Table 1 in Toggweiler+ (1999)
     !DATA RCP  / 106.0D0 /
     !DATA RNP  /  16.0D0 /
     !DATA RO2P / 172.0D0 /
      DATA RCP  / 162.5D0 / !! Table 1 in Toggweiler+ (1999)
      DATA RNP  /  15.0D0 / !! Table 1 in Toggweiler+ (1999)
      DATA RO2P / 169.0D0 / !! Table 1 in Toggweiler+ (1999)
 
!--   Biological sinking flux [molC/m^2/yr]
     !DATA CEPH / 0.75D0 / !! 0.075-7.5 molC/m^2/yr / Table 1 in Toggweiler+ (1999)
      DATA CEPH / 1.00D0 /

!--   Initial atomospheric PCO2 [ppmv]
      DATA PCO2A / 280.D0 /

!-- Initial condition (etc)

!--   Ocean area [m^2]
      AOCNH = AOCN * FAOCNH
      AOCNL = AOCN * FAOCNL
      AOCNN = AOCN * FAOCNN

!--   Ocean volume [m^3]
      VOCNH = AOCNH * ZOCNH
      VOCNL = AOCNL * ZOCNL
      VOCNN = AOCNN * ZOCNN
      VOCND = VOCN - VOCNH - VOCNL - VOCNN

!--   Biological sinking flux [molP/yr]
      EPH = (CEPH / RCP) * AOCNH

!--   Gas exchange [m^3/day]
      FAH = PVH * AOCNH
      FAL = PVL * AOCNL
      FAN = PVN * AOCNN

!--  Initial condition (biogeochemical parameter)

!--   Phosphate [mol/kg]
      DATA PO4H, PO4L, PO4N, PO4D &
     & / 2.09D-6, 2.09D-6, 2.09D-6, 2.09D-6 / !! Toggweiler and Sarmiento (1985)
!    & / 2.884D-6, 0.0D-6, 0.0D-6, 2.142D-6 /

!--   Alkalinity [eq/kg]
      DATA ALKH, ALKL, ALKN, ALKD &
     & / 2.371D-3, 2.371D-3, 2.371D-3, 2.371D-3 / !! Toggweiler and Sarmiento (1985)
!    & / 2.3626D-3, 2.277D-3, 2.277D-3, 2.3738D-3 /

!--   Dissolved inorganic carbon [mol/kg]
      DATA DICH, DICL, DICN, DICD &
     & / 2.258D-3, 2.258D-3, 2.258D-3, 2.258D-3 / !! Toggweiler and Sarmiento (1985)
!    & / 2.1316D-3, 2.070D-3, 1.933D-3, 2.2561D-3 /

!--   Dissolved oxygen [mol/kg]
      DATA DO2H, DO2L, DO2N, DO2D &
     & / 1.60D-4, 1.60D-4, 1.60D-4, 1.60D-4 /
!    & / 3.50D-4, 3.00D-4, 3.00D-4, 2.50D-4 /

!-- Unit conversion

      EPH   = EPH   / CV4 !! [/yr]    to [/s]
      FAH   = FAH   / CV5 !! [/day]   to [/s]
      FAL   = FAL   / CV5 !! [/day]   to [/s]
      FAN   = FAN   / CV5 !! [/day]   to [/s]
      HSC   = HSC   * CV1 !! [mol/kg] to [mol/m^3]
      R     = R     / CV4 !! [/yr]    to [/s]
      PO4H  = PO4H  * CV1 !! [mol/kg] to [mol/m^3]
      PO4L  = PO4L  * CV1 !! [mol/kg] to [mol/m^3]
      PO4N  = PO4N  * CV1 !! [mol/kg] to [mol/m^3]
      PO4D  = PO4D  * CV1 !! [mol/kg] to [mol/m^3]
      DICH  = DICH  * CV1 !! [mol/kg] to [mol/m^3]
      DICL  = DICL  * CV1 !! [mol/kg] to [mol/m^3]
      DICN  = DICN  * CV1 !! [mol/kg] to [mol/m^3]
      DICD  = DICD  * CV1 !! [mol/kg] to [mol/m^3]
      ALKH  = ALKH  * CV1 !! [mol/kg] to [mol/m^3]
      ALKL  = ALKL  * CV1 !! [mol/kg] to [mol/m^3]
      ALKN  = ALKN  * CV1 !! [mol/kg] to [mol/m^3]
      ALKD  = ALKD  * CV1 !! [mol/kg] to [mol/m^3]
      DO2H  = DO2H  * CV1 !! [mol/kg] to [mol/m^3]
      DO2L  = DO2L  * CV1 !! [mol/kg] to [mol/m^3]
      DO2N  = DO2N  * CV1 !! [mol/kg] to [mol/m^3]
      DO2D  = DO2D  * CV1 !! [mol/kg] to [mol/m^3]
      PCO2A = PCO2A / CV3 !! [ppmv]   to [atm]

!-- Loop

      DO I = 1, 100000

!-- Calculation of equilibrium constant
!   1:H, 2:L, 3:N, 4:D

         CALL CHEMEQCONST(                                             &
     &                     TEMH, SALH,                                 &
     &                     BT(1), K0(1), K1(1), K2(1),                 &
     &                     KB(1), KW(1), FH(1)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEML, SALL,                                 &
     &                     BT(2), K0(2), K1(2), K2(2),                 &
     &                     KB(2), KW(2), FH(2)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMN, SALN,                                 &
     &                     BT(3), K0(3), K1(3), K2(3),                 &
     &                     KB(3), KW(3), FH(3)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMD, SALD,                                 &
     &                     BT(4), K0(4), K1(4), K2(4),                 &
     &                     KB(4), KW(4), FH(4)                         &
     &                   )

!-- Initial CO2 and PCO2
!-- Initial total carbon

         IF (I == 1) THEN
            CALL CO2_NIBUN(                                            &
     &                      BT(1), K0(1), K1(1), K2(1), KB(1), KW(1),  &
     &                      ALKH, DICH, CO2H, HCO3H, CO32H, PCO2H      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(2), K0(2), K1(2), K2(2), KB(2), KW(2),  &
     &                      ALKL, DICL, CO2L, HCO3L, CO32L, PCO2L      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(3), K0(3), K1(3), K2(3), KB(3), KW(3),  &
     &                      ALKN, DICN, CO2N, HCO3N, CO32N, PCO2N      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(4), K0(4), K1(4), K2(4), KB(4), KW(4),  &
     &                      ALKD, DICD, CO2D, HCO3D, CO32D, PCO2D      &
     &                    )

            CALL O2SAT(TEMH, SALH, O2SATH)
            CALL O2SAT(TEML, SALL, O2SATL)
            CALL O2SAT(TEMN, SALN, O2SATN)
            CALL O2SAT(TEMD, SALD, O2SATD)

            TOCINI = + DICH  * VOCNH                                   &
     &               + DICL  * VOCNL                                   &
     &               + DICN  * VOCNN                                   &
     &               + DICD  * VOCND                                   &
     &               + PCO2A * VATM
         ENDIF

!-- Calculation of biogeochemical tracers

!        EPL = (TRAN + FHL) * (PO4H - PO4L) + FLD * (PO4D - PO4L)
         EPL = R * DEL * LF * PO4L * (PO4L / (HSC + PO4L)) * VOCNL
         EPN = R * DEN * LF * PO4N * (PO4N / (HSC + PO4N)) * VOCNN

         TEMHX = TEMH
         TEMLX = TEML
         TEMNX = TEMN
         TEMDX = TEMD
!        TEMDX = TEMD
!    &         + (                                                     &
!    &             + TRAN * (TEMN - TEMD)                              &
!    &             + FHD * (TEMH - TEMD)                               &
!    &             + FLD * (TEML - TEMD)                               &
!    &           )                                                     &
!    &           * (DT / VOCND)

         PO4HX = PO4H                                                  &
     &         + (                                                     &
     &             + (TRAN + FHD) * (PO4D - PO4H)                      &
     &             + FHL * (PO4L - PO4H)                               &
     &             - EPH                                               &
     &           )                                                     &
     &           * (DT / VOCNH)
         PO4LX = PO4L                                                  &
     &         + (                                                     &
     &             + (TRAN + FHL) * (PO4H - PO4L)                      &
     &             + FLN * (PO4N - PO4L)                               &
     &             + FLD * (PO4D - PO4L)                               &
     &             - EPL                                               &
     &           )                                                     &
     &           * (DT / VOCNL)
         PO4NX = PO4N                                                  &
     &         + (                                                     &
     &             + (TRAN + FLN) * (PO4L - PO4N)                      &
     &             + FND * (PO4D - PO4N)                               &
     &             - EPN                                               &
     &           )                                                     &
     &           * (DT / VOCNN)
         PO4DX = PO4D                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (PO4N - PO4D)                      &
     &             + FHD * (PO4H - PO4D)                               &
     &             + FLD * (PO4L - PO4D)                               &
     &             + (EPH + EPL + EPN)                                 &
     &           )                                                     &
     &           * (DT / VOCND)

         ALKHX = ALKH                                                  &
     &         + (                                                     &
     &             + (TRAN + FHD) * (ALKD - ALKH)                      &
     &             + FHL * (ALKL - ALKH)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPH                   &
     &           )                                                     &
     &           * (DT / VOCNH)
         ALKLX = ALKL                                                  &
     &         + (                                                     &
     &             + (TRAN + FHL) * (ALKH - ALKL)                      &
     &             + FLN * (ALKN - ALKL)                               &
     &             + FLD * (ALKD - ALKL)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPL                   &
     &           )                                                     &
     &           * (DT / VOCNL)
         ALKNX = ALKN                                                  &
     &         + (                                                     &
     &             + (TRAN + FLN) * (ALKL - ALKN)                      &
     &             + FND * (ALKD - ALKN)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPN                   &
     &           )                                                     &
     &           * (DT / VOCNN)
         ALKDX = ALKD                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (ALKN - ALKD)                      &
     &             + FHD * (ALKH - ALKD)                               &
     &             + FLD * (ALKL - ALKD)                               &
     &             + (2.0D0 * RRC * RCP - RNP) * (EPH + EPL + EPN)     &
     &           )                                                     &
     &           * (DT / VOCND)

         DICHX = DICH                                                  &
     &         + (                                                     &
     &             + (TRAN + FHD) * (DICD - DICH)                      &
     &             + FHL * (DICL - DICH)                               &
     &             - (1.0D0 + RRC) * RCP * EPH                         &
     &             + FAH * CV1 * K0(1) * (PCO2A - PCO2H)               &
     &           )                                                     &
     &           * (DT / VOCNH)
         DICLX = DICL                                                  &
     &         + (                                                     &
     &             + (TRAN + FHL) * (DICH - DICL)                      &
     &             + FLN * (DICN - DICL)                               &
     &             + FLD * (DICD - DICL)                               &
     &             - (1.0D0 + RRC) * RCP * EPL                         &
     &             + FAL * CV1 * K0(2) * (PCO2A - PCO2L)               &
     &           )                                                     &
     &           * (DT / VOCNL)
         DICNX = DICN                                                  &
     &         + (                                                     &
     &             + (TRAN + FLN) * (DICL - DICN)                      &
     &             + FND * (DICD - DICN)                               &
     &             - (1.0D0 + RRC) * RCP * EPN                         &
     &             + FAN * CV1 * K0(3) * (PCO2A - PCO2N)               &
     &           )                                                     &
     &           * (DT / VOCNN)
         DICDX = DICD                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (DICN - DICD)                      &
     &             + FHD * (DICH - DICD)                               &
     &             + FLD * (DICL - DICD)                               &
     &             + (1.0D0 + RRC) * RCP * (EPH + EPL + EPN)           &
     &           )                                                     &
     &           * (DT / VOCND)

         DO2HX = DO2H                                                  &
     &         + (                                                     &
     &             + (TRAN + FHD) * (DO2D - DO2H)                      &
     &             + FHL * (O2SATL - DO2H)                             &
     &             + RO2P * EPH                                        &
     &             + FAH * (O2SATH - DO2H)                             &
     &           )                                                     &
     &           * (DT / VOCNH)
         DO2LX = DO2L                                                  &
     &         + (                                                     &
     &             + (TRAN + FHL) * (DO2H - DO2L)                      &
     &             + FLN * (DO2N - DO2L)                               &
     &             + FLD * (DO2D - DO2L)                               &
     &             + RO2P * EPL                                        &
     &             + FAL * (O2SATL - DO2L)                             &
     &           )                                                     &
     &           * (DT / VOCNL)
         DO2NX = DO2N                                                  &
     &         + (                                                     &
     &             + (TRAN + FLN) * (DO2L - DO2N)                      &
     &             + FND * (DO2D - DO2N)                               &
     &             + RO2P * EPN                                        &
     &             + FAN * (O2SATN - DO2N)                             &
     &           )                                                     &
     &           * (DT / VOCNN)
         DO2DX = DO2D                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (DO2N - DO2D)                      &
     &             + FHD * (DO2H - DO2D)                               &
     &             + FLD * (DO2L - DO2D)                               &
     &             - RO2P * (EPH + EPL + EPN)                          &
     &           )                                                     &
     &           * (DT / VOCND)

!-- Calculation of carbonate system

         CALL CO2_NIBUN(                                               &
     &                   BT(1), K0(1), K1(1), K2(1), KB(1), KW(1),     &
     &                   ALKHX, DICHX, CO2H, HCO3H, CO32H, PCO2H       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(2), K0(2), K1(2), K2(2), KB(2), KW(2),     &
     &                   ALKLX, DICLX, CO2L, HCO3L, CO32L, PCO2L       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(3), K0(3), K1(3), K2(3), KB(3), KW(3),     &
     &                   ALKNX, DICNX, CO2N, HCO3N, CO32N, PCO2N       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(4), K0(4), K1(4), K2(4), KB(4), KW(4),     &
     &                   ALKDX, DICDX, CO2D, HCO3D, CO32D, PCO2D       &
     &                 )

!-- Calculation of atmospheric PCO2

         PCO2AX = PCO2A                                                &
     &          + (                                                    &
     &              + FAH * CV1 * K0(1) * (PCO2H - PCO2A)              &
     &              + FAL * CV1 * K0(2) * (PCO2L - PCO2A)              &
     &              + FAN * CV1 * K0(3) * (PCO2N - PCO2A)              &
     &            )                                                    &
     &            * (DT / VATM)

!-- Next loop

         OLD_PO4H = PO4H
         TEMH  = TEMHX
         TEML  = TEMLX
         TEMN  = TEMNX
         TEMD  = TEMDX
         PO4H  = PO4HX
         PO4L  = PO4LX
         PO4N  = PO4NX
         PO4D  = PO4DX
         ALKH  = ALKHX
         ALKL  = ALKLX
         ALKN  = ALKNX
         ALKD  = ALKDX
         DICH  = DICHX
         DICL  = DICLX
         DICN  = DICNX
         DICD  = DICDX
         DO2H  = DO2HX
         DO2L  = DO2LX
         DO2N  = DO2NX
         DO2D  = DO2DX
         PCO2A = PCO2AX

!-- Equilibrium condition

         IF (ABS((PO4H / OLD_PO4H) - 1.0D0) <= DELTA) THEN
           WRITE(*,*) 'TIMESTEP =', I
           EXIT
         ENDIF

      ENDDO

!-- Final total carbon

       TOCFIN = + DICH  * VOCNH                                        &
     &          + DICL  * VOCNL                                        &
     &          + DICN  * VOCNN                                        &
     &          + DICD  * VOCND                                        &
     &          + PCO2A * VATM

!-- Oxygen saturation

      CALL O2SAT(TEMH, SALH, O2SATH)
      CALL O2SAT(TEML, SALL, O2SATL)
      CALL O2SAT(TEMN, SALN, O2SATN)
      CALL O2SAT(TEMD, SALD, O2SATD)

      AOUH = (O2SATH - DO2H)
      AOUL = (O2SATL - DO2L)
      AOUN = (O2SATN - DO2N)
      AOUD = (O2SATD - DO2D)

!-- Unit conversion

      EPH   = CV4 * EPH * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      EPL   = CV4 * EPL * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      EPN   = CV4 * EPN * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      PO4H  = CV2 * 1.D+6 * PO4H              !! [mol/m^3] to [umol/kg]
      PO4L  = CV2 * 1.D+6 * PO4L              !! [mol/m^3] to [umol/kg]
      PO4N  = CV2 * 1.D+6 * PO4N              !! [mol/m^3] to [umol/kg]
      PO4D  = CV2 * 1.D+6 * PO4D              !! [mol/m^3] to [umol/kg]
      DICH  = CV2 * 1.D+6 * DICH              !! [mol/m^3] to [umol/kg]
      DICL  = CV2 * 1.D+6 * DICL              !! [mol/m^3] to [umol/kg]
      DICN  = CV2 * 1.D+6 * DICN              !! [mol/m^3] to [umol/kg]
      DICD  = CV2 * 1.D+6 * DICD              !! [mol/m^3] to [umol/kg]
      ALKH  = CV2 * 1.D+6 * ALKH              !! [mol/m^3] to [umol/kg]
      ALKL  = CV2 * 1.D+6 * ALKL              !! [mol/m^3] to [umol/kg]
      ALKN  = CV2 * 1.D+6 * ALKN              !! [mol/m^3] to [umol/kg]
      ALKD  = CV2 * 1.D+6 * ALKD              !! [mol/m^3] to [umol/kg]
      DO2H  = CV2 * 1.D+6 * DO2H              !! [mol/m^3] to [umol/kg]
      DO2L  = CV2 * 1.D+6 * DO2L              !! [mol/m^3] to [umol/kg]
      DO2N  = CV2 * 1.D+6 * DO2N              !! [mol/m^3] to [umol/kg]
      DO2D  = CV2 * 1.D+6 * DO2D              !! [mol/m^3] to [umol/kg]
      AOUH  = CV2 * 1.D+6 * AOUH              !! [mol/m^3] to [umol/kg]
      AOUL  = CV2 * 1.D+6 * AOUL              !! [mol/m^3] to [umol/kg]
      AOUN  = CV2 * 1.D+6 * AOUN              !! [mol/m^3] to [umol/kg]
      AOUD  = CV2 * 1.D+6 * AOUD              !! [mol/m^3] to [umol/kg]
      CO2H  = CV2 * 1.D+6 * CO2H              !! [mol/m^3] to [umol/kg]
      CO2L  = CV2 * 1.D+6 * CO2L              !! [mol/m^3] to [umol/kg]
      CO2N  = CV2 * 1.D+6 * CO2N              !! [mol/m^3] to [umol/kg]
      CO2D  = CV2 * 1.D+6 * CO2D              !! [mol/m^3] to [umol/kg]
      HCO3H = CV2 * 1.D+6 * HCO3H             !! [mol/m^3] to [umol/kg]
      HCO3L = CV2 * 1.D+6 * HCO3L             !! [mol/m^3] to [umol/kg]
      HCO3N = CV2 * 1.D+6 * HCO3N             !! [mol/m^3] to [umol/kg]
      HCO3D = CV2 * 1.D+6 * HCO3D             !! [mol/m^3] to [umol/kg]
      CO32H = CV2 * 1.D+6 * CO32H             !! [mol/m^3] to [umol/kg]
      CO32L = CV2 * 1.D+6 * CO32L             !! [mol/m^3] to [umol/kg]
      CO32N = CV2 * 1.D+6 * CO32N             !! [mol/m^3] to [umol/kg]
      CO32D = CV2 * 1.D+6 * CO32D             !! [mol/m^3] to [umol/kg]
      PCO2H = CV3 * PCO2H                     !! [atm]     to [ppmv]
      PCO2L = CV3 * PCO2L                     !! [atm]     to [ppmv]
      PCO2N = CV3 * PCO2N                     !! [atm]     to [ppmv]
      PCO2D = CV3 * PCO2D                     !! [atm]     to [ppmv]
      PCO2A = CV3 * PCO2A                     !! [atm]     to [ppmv]
      TOCINI = TOCINI * 12.D-15               !! [molC]    to [PgC]
      TOCFIN = TOCFIN * 12.D-15               !! [molC]    to [PgC]

!-- Output

      WRITE(*,*) 'Temperature (H)          :', TEMH, 'K'
      WRITE(*,*) 'Temperature (L)          :', TEML, 'K'
      WRITE(*,*) 'Temperature (N)          :', TEMN, 'K'
      WRITE(*,*) 'Temperature (D)          :', TEMD, 'K'
      WRITE(*,*) ''
      WRITE(*,*) 'Salinity (H)             :', SALH, 'PSU'
      WRITE(*,*) 'Salinity (L)             :', SALL, 'PSU'
      WRITE(*,*) 'Salinity (N)             :', SALN, 'PSU'
      WRITE(*,*) 'Salinity (D)             :', SALD, 'PSU'
      WRITE(*,*) ''
      WRITE(*,*) 'EP (H)                   :', EPH, 'PgC/yr'
      WRITE(*,*) 'EP (L)                   :', EPL, 'PgC/yr'
      WRITE(*,*) 'EP (N)                   :', EPN, 'PgC/yr'
      WRITE(*,*) ''
      WRITE(*,*) 'PO4 (H)                  :', PO4H, 'umol/kg'
      WRITE(*,*) 'PO4 (L)                  :', PO4L, 'umol/kg'
      WRITE(*,*) 'PO4 (N)                  :', PO4N, 'umol/kg'
      WRITE(*,*) 'PO4 (D)                  :', PO4D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DIC (H)                  :', DICH, 'umol/kg'
      WRITE(*,*) 'DIC (L)                  :', DICL, 'umol/kg'
      WRITE(*,*) 'DIC (N)                  :', DICN, 'umol/kg'
      WRITE(*,*) 'DIC (D)                  :', DICD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'ALK (H)                  :', ALKH, 'ueq/kg'
      WRITE(*,*) 'ALK (L)                  :', ALKL, 'ueq/kg'
      WRITE(*,*) 'ALK (N)                  :', ALKN, 'ueq/kg'
      WRITE(*,*) 'ALK (D)                  :', ALKD, 'ueq/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DO2 (H)                  :', DO2H, 'umol/kg'
      WRITE(*,*) 'DO2 (L)                  :', DO2L, 'umol/kg'
      WRITE(*,*) 'DO2 (N)                  :', DO2N, 'umol/kg'
      WRITE(*,*) 'DO2 (D)                  :', DO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'AOU (H)                  :', AOUH, 'umol/kg'
      WRITE(*,*) 'AOU (L)                  :', AOUL, 'umol/kg'
      WRITE(*,*) 'AOU (N)                  :', AOUN, 'umol/kg'
      WRITE(*,*) 'AOU (D)                  :', AOUD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO2 (H)                  :', CO2H, 'umol/kg'
      WRITE(*,*) 'CO2 (L)                  :', CO2L, 'umol/kg'
      WRITE(*,*) 'CO2 (N)                  :', CO2N, 'umol/kg'
      WRITE(*,*) 'CO2 (D)                  :', CO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'HCO3 (H)                 :', HCO3H, 'umol/kg'
      WRITE(*,*) 'HCO3 (L)                 :', HCO3L, 'umol/kg'
      WRITE(*,*) 'HCO3 (N)                 :', HCO3N, 'umol/kg'
      WRITE(*,*) 'HCO3 (D)                 :', HCO3D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO32 (H)                 :', CO32H, 'umol/kg'
      WRITE(*,*) 'CO32 (L)                 :', CO32L, 'umol/kg'
      WRITE(*,*) 'CO32 (N)                 :', CO32N, 'umol/kg'
      WRITE(*,*) 'CO32 (D)                 :', CO32D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'PCO2 (H)                 :', PCO2H, 'ppmv'
      WRITE(*,*) 'PCO2 (L)                 :', PCO2L, 'ppmv'
      WRITE(*,*) 'PCO2 (N)                 :', PCO2N, 'ppmv'
      WRITE(*,*) 'PCO2 (D)                 :', PCO2D, 'ppmv'
      WRITE(*,*) 'PCO2 (ATM)               :', PCO2A, 'ppmv'
      WRITE(*,*) ''
      WRITE(*,*) 'Total Carbon (Initial)   :', TOCINI, 'PgC'
      WRITE(*,*) 'Total Carbon (Final)     :', TOCFIN, 'PgC'
      WRITE(*,*) ''

      END PROGRAM
