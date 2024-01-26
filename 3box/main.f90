      PROGRAM MAIN

      USE CHEMEQ
      USE O2
      USE CO2

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Ocean 3-box model
!  Calculation of phosphate, DIC, alkalinity, CO2, and PCO2
!
!  References: Toggweiler et al. (1999, PA)
!
! ---------------------------------------------------------------------

      INTEGER :: I, ISTEP
      INTEGER, PARAMETER :: NBOX=3
!     1:H, 2:L, 3:D

      REAL(8) :: DT, DELTA
      REAL(8) :: CV1, CV2, CV3, CV4, CV5

      REAL(8) :: VOCN, AOCN, VATM
      REAL(8) :: ZOCNH, DEH, AOCNH, VOCNH
      REAL(8) :: ZOCNL, DEL, AOCNL, VOCNL
      REAL(8) :: ZOCND, VOCND

      REAL(8) :: T, FHD, FLD, FLH, PVH, PVL, FAH, FAL
      REAL(8) :: R, LF, HSC
      REAL(8) :: RRC, RCP, RNP, RO2P, CEPH
      REAL(8) :: BT(NBOX)
      REAL(8) :: K0(NBOX), K1(NBOX), K2(NBOX)
      REAL(8) :: KB(NBOX), KW(NBOX), FH(NBOX)

      REAL(8) :: TEMH, SALH,   TEMHX, SALHX
      REAL(8) :: PO4H,  ALKH,   DICH,   DO2H
      REAL(8) :: PO4HX, ALKHX,  DICHX,  DO2HX
      REAL(8) :: CO2H,  HCO3H,  CO32H,  PCO2H, O2SATH, AOUH, EPH
      REAL(8) :: CO2HX, HCO3HX, CO32HX, PCO2HX

      REAL(8) :: TEML, SALL,   TEMLX, SALLX
      REAL(8) :: PO4L,  ALKL,   DICL,   DO2L
      REAL(8) :: PO4LX, ALKLX,  DICLX,  DO2LX
      REAL(8) :: CO2L,  HCO3L,  CO32L,  PCO2L, O2SATL, AOUL, EPL
      REAL(8) :: CO2LX, HCO3LX, CO32LX, PCO2LX

      REAL(8) :: TEMD, SALD,   TEMDX, SALDX
      REAL(8) :: PO4D,  ALKD,   DICD,   DO2D
      REAL(8) :: PO4DX, ALKDX,  DICDX,  DO2DX
      REAL(8) :: CO2D,  HCO3D,  CO32D,  PCO2D, O2SATD, AOUD, EPD
      REAL(8) :: CO2DX, HCO3DX, CO32DX, PCO2DX

      REAL(8) :: PCO2A, PCO2AX
      REAL(8) :: TOCINI, TOCFIN

      DATA DT    / 8.64D4 /
      DATA DELTA / 1.0D-6 /

      DATA CV1 / 1.0250D+3 / !! [mol/kg] to [mol/m^3]
      DATA CV2 / 9.7561D-4 / !! [mol/m^3] to [mol/kg]
      DATA CV3 / 1.0000D+6 / !! [atm] to [ppmv]
      DATA CV4 / 3.1536D+7 / !! [yr] to [s]
      DATA CV5 / 8.6400D+4 / !! [day] to [s]

!--   Temperature [K]
!--   Salinity [PSU]
      DATA TEMH, TEML, TEMD &
     & / 2.0, 21.5, 1.75 /
      DATA SALH, SALL, SALD &
     & / 34.7, 34.7, 34.7 /

!--   Ocean volume [m^3]
!--   Ocean area [m^2]
!--   Atmospheric molers [mol] = [mol/atm]
      DATA VOCN, AOCN, VATM &
     & / 1.292D18, 3.49D14, 1.773D20 /

!--   Ocean depth [m]
      DATA ZOCNH, ZOCNL, ZOCND &
     & / 2.5D2, 1.0D2, 4.0D3 /

!--   Ocean depth of euphotic zone [m]
      DATA DEH, DEL / 2.5D2, 1.0D2 /

!--   Transport [m^3/s]
!--   Diffusivity [m^3/s]
      DATA T &
     & / 2.0D7 /
      DATA FHD, FLD, FLH &
     & / 6.0D7, 0.0D0, 0.0D0 /

!--   Piston velocity [m/day]
      DATA PVH, PVL / 3.0D0, 3.0D0 /

!--   Bio-production efficiency [/yr]
      R = 1.0D0 / CV4

!--   Proportional to an annual averaged solar radiation
      DATA LF / 5.0D-1 /

!--   Half-saturation constant [mol/kg]
      DATA HSC / 2.0D-8 /

!--   Rain ratio, C/P, N/P, O2/P
      DATA RRC  / 0.25D0 /
      DATA RCP  / 106.D0 /
      DATA RNP  /  16.D0 /
      DATA RO2P / 172.D0 /

!--   Biological sinking flux [molC/m^2/yr]
      DATA CEPH / 7.5D-1 /

!--   Initial atomospheric PCO2 [ppmv]
      DATA PCO2A / 280.D0 /

!-- Initial condition (etc)

!--   Ocean area [m^2]
!--   H, L
!--   15% , 85%
      AOCNH = AOCN * 1.5D-1
      AOCNL = AOCN * 8.5D-1

!--   Ocean volume [m^3]
      VOCNH = AOCNH * ZOCNH
      VOCNL = AOCNL * ZOCNL
      VOCND = VOCN - VOCNH - VOCNL

!--   Biological sinking flux [molP/yr]
      EPH = (CEPH / RCP) * AOCNH

!--   Gas exchange [m^3/day]
      FAH = PVH * AOCNH
      FAL = PVL * AOCNL

!--  Initial condition (biogeochemical parameter)

!--   Phosphate [mol/kg]
      DATA PO4H, PO4L, PO4D &
     & / 2.10D-6, 2.10D-6, 2.10D-6 /
!    & / 1.421D-6, 0.2D-6, 2.142D-6 /

!--   Alkalinity [eq/kg]
      DATA ALKH, ALKL, ALKD &
     & / 2.374D-3, 2.374D-3, 2.374D-3 /
!    & / 2.3378D-3, 2.2767D-3, 2.3738D-3 /

!--   Dissolved inorganic carbon [mol/kg]
      DATA DICH, DICL, DICD &
     & / 2.235D-3, 2.235D-3, 2.235D-3 /
!    & / 2.1391D-3, 1.9327D-3, 2.2561D-3 /

!--   Dissolved oxygen [mol/kg]
      DATA DO2H, DO2D &
     & / 3.25D-4, 1.70D-4 /

!-- Unit conversion

      EPH   = EPH   / CV4 !! [/yr]    to [/s]
      FAH   = FAH   / CV5 !! [/day]   to [/s]
      FAL   = FAL   / CV5 !! [/day]   to [/s]
      PO4H  = PO4H  * CV1 !! [mol/kg] to [mol/m^3]
      PO4L  = PO4L  * CV1 !! [mol/kg] to [mol/m^3]
      PO4D  = PO4D  * CV1 !! [mol/kg] to [mol/m^3]
      DICH  = DICH  * CV1 !! [mol/kg] to [mol/m^3]
      DICL  = DICL  * CV1 !! [mol/kg] to [mol/m^3]
      DICD  = DICD  * CV1 !! [mol/kg] to [mol/m^3]
      ALKH  = ALKH  * CV1 !! [mol/kg] to [mol/m^3]
      ALKL  = ALKL  * CV1 !! [mol/kg] to [mol/m^3]
      ALKD  = ALKD  * CV1 !! [mol/kg] to [mol/m^3]
      DO2H  = DO2H  * CV1 !! [mol/kg] to [mol/m^3]
      DO2D  = DO2D  * CV1 !! [mol/kg] to [mol/m^3]
      PCO2A = PCO2A / CV3 !! [ppmv] to [atm]

!-- Loop

      DO I = 1, 100000

!-- Calculation of equilibrium constant
!   1:H, 2:L, 3:D

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
     &                     TEMD, SALD,                                 &
     &                     BT(3), K0(3), K1(3), K2(3),                 &
     &                     KB(3), KW(3), FH(3)                         &
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
     &                      ALKD, DICD, CO2D, HCO3D, CO32D, PCO2D      &
     &                    )
            TOCINI = + DICH  * VOCNH                                   &
     &               + DICL  * VOCNL                                   &
     &               + DICD  * VOCND                                   &
     &               + PCO2A * VATM
         ENDIF

!-- Calculation of biogeochemical tracers
!2345&78901234567890123456789012345678901234567890123456789012345678901&


!         EPL = T * PO4D
!!--      EP_H = R * DE_H * LF * PO_H * (PO_H / (HSC + PO_H))
!!--      EP_L = R * DE_L * LF * PO_L * (PO_L / (HSC + PO_L))
!
!         TEMDX = TEMD &
!     &   + ((T + FHD) * (TEMH - TEMD) + FLD * (TEML - TEMD)) &
!     &   * (DT / VOCND)
!
!         PO4HX = PO4H &
!     &   + ((T + FLH) * (PO4L - PO4H) + FHD * (PO4D - PO4H) - EPH)
!     &   * (DT / VOCNH)
!         PO4LX = PO4L &
!     &   + ((T + FLD) * (PO4D - PO4L) + FLH * (PO4H - PO4L) - EPL)
!     &   * (DT / VOCNL)
!         PO4DX = PO4D &
!     &   + ((T + FHD) * (PO4H - PO4D) + FLD * (PO4L - PO4D) + (EPH + EPL))
!     &   * (DT / VOCND)
!
!         ALKHX = ALKH &
!     &   + ((T + FLH) * (ALKL - ALKH) + FHD * (ALKD - ALKH) - (2.0D0 * RRC * RCP - RNP) * EPH)
!     &   * (DT / VOCNH)
!         ALKLX = ALKL &
!     &   + ((T + FLD) * (ALKD - ALKL) + FLH * (ALKH - ALKL) - (2.0D0 * RRC * RCP - RNP) * EPL)
!     &   * (DT / VOCNL)
!         ALKDX = ALKD &
!     &   + ((T + FHD) * (ALKH - ALKD) + FLD * (ALKL - ALKD) + (2.0D0 * RRC * RCP - RNP) * (EPH + EPL))
!     &   * (DT / VOCND)
!
!         DICHX = DICH &
!     &   + ((T + FLH) * (DICL - DICH) + FHD * (DICD - DICH) &
!         &     - (1.0D0 + RRC) * RCP * EPH + FAH * CV1 * K0(1) * (PCO2A - PCO2H))
!     & * (DT / VOCNH)
!         DICLX = DICL &
!     &   + ((T + FLD) * (DICD - DICL) + FLH * (DICH - DICL) &
!         &     - (1.0D0 + RRC) * RCP * EPL + FAL * CV1 * K0(2) * (PCO2A - PCO2L))
!     &   * (DT / VOCNL)
!         DICDX = DICD &
!     &   + ((T + FHD) * (DICH - DICD) + FLD * (DICL - DICD) + (1.0D0 + RRC) * RCP * (EPH + EPL))
!     & * (DT / VOCND)
!
!         DO2HX = DO2H &
!     &   + ((T + FLH) * (O2SATL - DO2H) + FHD * (O2SATH - DO2H) + RDO2P * EPH + FAH * (O2SATH - DO2H))
!     &   * (DT / VOCNH)
!!        DO2LX = DO2L + ((T + FLD) * (DO2D - DO2L) + FLH * (DO2H - DO2L) + RDO2P * EPL)
!!    &* (DT / VOCNL)
!         DO2DX = DO2D &
!     &   + ((T + FHD) * (DO2H - DO2D) + FLD * (O2SATL - DO2D) - RDO2P * (EPH + EPL))
!     &* (DT / VOCND)

!-- Calculation of carbonate system

         CALL CO2_NIBUN(                                               &
     &                   BT(1), K0(1), K1(1), K2(1), KB(1), KW(1),     &
     &                   ALKHX, DICHX, CO2HX, HCO3H, CO32H, PCO2H      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(2), K0(2), K1(2), K2(2), KB(2), KW(2),     &
     &                   ALKLX, DICLX, CO2LX, HCO3L, CO32L, PCO2L      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(3), K0(3), K1(3), K2(3), KB(3), KW(3),     &
     &                   ALKDX, DICDX, CO2DX, HCO3D, CO32D, PCO2D      &
     &                 )

!-- Calculation of atmospheric PCO2

         PCO2AX = PCO2A                                                &
     &    + (                                                          &
     &        + K0(1) * FAH * (PCO2H - PCO2A)                          &
     &        + K0(2) * FAL * (PCO2L - PCO2A)                          &
     &      )                                                          &
     &      * (DT / VATM)

!-- Equilibrium condition

         IF (ABS((PO4HX / PO4H) - 1.0D0) <= DELTA) THEN
           WRITE(*,*) 'TIMESTEP =', I
           EXIT
         ENDIF

!-- Next loop

!        TEMD  = TEMDX
         PO4H  = PO4HX
         PO4L  = PO4LX
         PO4D  = PO4DX
         ALKH  = ALKHX
         ALKL  = ALKLX
         ALKD  = ALKDX
         DICH  = DICHX
         DICL  = DICLX
         DICD  = DICDX
         DO2H  = DO2HX
         DO2D  = DO2DX
         PCO2A = PCO2AX

      ENDDO

!-- Final total carbon

      TOCFIN = + DICH  * VOCNH                                         &
     &         + DICL  * VOCNL                                         &
     &         + DICD  * VOCND                                         &
     &         + PCO2A * VATM

!-- Oxygen utilization

      CALL O2SAT(TEMH, SALH, O2SATH)
      CALL O2SAT(TEML, SALL, O2SATL)
      CALL O2SAT(TEMD, SALD, O2SATD)

      AOUH = (O2SATH - DO2H)
      AOUL = (O2SATL - DO2L)
      AOUD = (O2SATD - DO2D)

!-- Unit conversion

      EPH   = CV4 * EPH * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      EPL   = CV4 * EPL * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      PO4H  = CV2 * 1.D+6 * PO4H              !! [mol/m^3] to [umol/kg]
      PO4L  = CV2 * 1.D+6 * PO4L              !! [mol/m^3] to [umol/kg]
      PO4D  = CV2 * 1.D+6 * PO4D              !! [mol/m^3] to [umol/kg]
      DICH  = CV2 * 1.D+6 * DICH              !! [mol/m^3] to [umol/kg]
      DICL  = CV2 * 1.D+6 * DICL              !! [mol/m^3] to [umol/kg]
      DICD  = CV2 * 1.D+6 * DICD              !! [mol/m^3] to [umol/kg]
      ALKH  = CV2 * 1.D+6 * ALKH              !! [mol/m^3] to [umol/kg]
      ALKL  = CV2 * 1.D+6 * ALKL              !! [mol/m^3] to [umol/kg]
      ALKD  = CV2 * 1.D+6 * ALKD              !! [mol/m^3] to [umol/kg]
      DO2H  = CV2 * 1.D+6 * DO2H              !! [mol/m^3] to [umol/kg]
      DO2L  = CV2 * 1.D+6 * DO2L              !! [mol/m^3] to [umol/kg]
      DO2D  = CV2 * 1.D+6 * DO2D              !! [mol/m^3] to [umol/kg]
      AOUH  = CV2 * 1.D+6 * AOUH              !! [mol/m^3] to [umol/kg]
      AOUL  = CV2 * 1.D+6 * AOUL              !! [mol/m^3] to [umol/kg]
      AOUD  = CV2 * 1.D+6 * AOUD              !! [mol/m^3] to [umol/kg]
      CO2H  = CV2 * 1.D+6 * CO2H              !! [mol/m^3] to [umol/kg]
      CO2L  = CV2 * 1.D+6 * CO2L              !! [mol/m^3] to [umol/kg]
      CO2D  = CV2 * 1.D+6 * CO2D              !! [mol/m^3] to [umol/kg]
      HCO3H = CV2 * 1.D+6 * HCO3H             !! [mol/m^3] to [umol/kg]
      HCO3L = CV2 * 1.D+6 * HCO3L             !! [mol/m^3] to [umol/kg]
      HCO3D = CV2 * 1.D+6 * HCO3D             !! [mol/m^3] to [umol/kg]
      CO32H = CV2 * 1.D+6 * CO32H             !! [mol/m^3] to [umol/kg]
      CO32L = CV2 * 1.D+6 * CO32L             !! [mol/m^3] to [umol/kg]
      CO32D = CV2 * 1.D+6 * CO32D             !! [mol/m^3] to [umol/kg]
      PCO2H = CV3 * PCO2H                     !! [atm] to [ppmv]
      PCO2L = CV3 * PCO2L                     !! [atm] to [ppmv]
      PCO2D = CV3 * PCO2D                     !! [atm] to [ppmv]
      PCO2A = CV3 * PCO2A                     !! [atm] to [ppmv]

!-- Output

      WRITE(*,*) 'Temperature (H)          :', TEMH, 'K'
      WRITE(*,*) 'Temperature (L)          :', TEML, 'K'
      WRITE(*,*) 'Temperature (D)          :', TEMD, 'K'
      WRITE(*,*) ''
      WRITE(*,*) 'Salinity (H)             :', SALH, 'PSU'
      WRITE(*,*) 'Salinity (L)             :', SALL, 'PSU'
      WRITE(*,*) 'Salinity (D)             :', SALD, 'PSU'
      WRITE(*,*) ''
      WRITE(*,*) 'EP (H)                   :', EPH, 'PgC/yr'
      WRITE(*,*) ''
      WRITE(*,*) 'PO4 (H)                  :', PO4H, 'umol/kg'
      WRITE(*,*) 'PO4 (L)                  :', PO4L, 'umol/kg'
      WRITE(*,*) 'PO4 (D)                  :', PO4D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DIC (H)                  :', DICH, 'umol/kg'
      WRITE(*,*) 'DIC (L)                  :', DICL, 'umol/kg'
      WRITE(*,*) 'DIC (D)                  :', DICD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'ALK (H)                  :', ALKH, 'ueq/kg'
      WRITE(*,*) 'ALK (L)                  :', ALKL, 'ueq/kg'
      WRITE(*,*) 'ALK (D)                  :', ALKD, 'ueq/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DO2 (H)                  :', DO2H, 'umol/kg'
      WRITE(*,*) 'DO2 (L)                  :', DO2L, 'umol/kg'
      WRITE(*,*) 'DO2 (D)                  :', DO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'AOU (H)                  :', AOUH, 'umol/kg'
      WRITE(*,*) 'AOU (L)                  :', AOUL, 'umol/kg'
      WRITE(*,*) 'AOU (D)                  :', AOUD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO2 (H)                  :', CO2H, 'umol/kg'
      WRITE(*,*) 'CO2 (L)                  :', CO2L, 'umol/kg'
      WRITE(*,*) 'CO2 (D)                  :', CO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'HCO3 (H)                 :', HCO3H, 'umol/kg'
      WRITE(*,*) 'HCO3 (L)                 :', HCO3L, 'umol/kg'
      WRITE(*,*) 'HCO3 (D)                 :', HCO3D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO32 (H)                 :', CO32H, 'umol/kg'
      WRITE(*,*) 'CO32 (L)                 :', CO32L, 'umol/kg'
      WRITE(*,*) 'CO32 (D)                 :', CO32D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'PCO2 (H)                 :', PCO2H, 'ppmv'
      WRITE(*,*) 'PCO2 (L)                 :', PCO2L, 'ppmv'
      WRITE(*,*) 'PCO2 (D)                 :', PCO2D, 'ppmv'
      WRITE(*,*) 'PCO2 (ATM)               :', PCO2A, 'ppmv'
      WRITE(*,*) ''
      WRITE(*,*) 'Total Carbon (Initial)   :', TOCINI, 'mol'
      WRITE(*,*) 'Total Carbon (Final)     :', TOCFIN, 'mol'
      WRITE(*,*) ''

      END PROGRAM
