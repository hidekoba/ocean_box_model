      PROGRAM MAIN

      USE CHEMEQ
      USE O2
      USE CO2

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Ocean 2-box model
!  Calculation of phosphate, DIC, alkalinity, CO2, and PCO2
!
!  References: Toggweiler et al. (1999, PA)
!              Emerson and Hedges. (2008)
!
! ---------------------------------------------------------------------

      INTEGER :: I
      INTEGER, PARAMETER :: NBOX=2

      REAL(8) :: DT, DELTA
      REAL(8) :: CV1, CV2, CV3, CV4, CV5

      REAL(8) :: VOCN, AOCN, VATM
      REAL(8) :: ZOCNS, AOCNS, VOCNS
      REAL(8) :: ZOCND, VOCND

      REAL(8) :: TRAN, FSD, PVS, FAS
      REAL(8) :: R, LF, HSC
      REAL(8) :: RRC, RCP, RNP, RO2P, CEPS
      REAL(8) :: BTS, K0S, K1S, K2S, KBS, KWS, FHS
      REAL(8) :: BTD, K0D, K1D, K2D, KBD, KWD, FHD
      REAL(8) :: TEMS, SALS, TEMSX! SALSX
      REAL(8) :: TEMD, SALD, TEMDX! SALDX

      REAL(8) :: PO4S,  ALKS,  DICS,  DO2S
      REAL(8) :: PO4SX, ALKSX, DICSX, DO2SX
      REAL(8) :: CO2S,  HCO3S, CO32S, PCO2S, O2SATS, AOUS, EPS
      REAL(8) :: PO4D,  ALKD,  DICD,  DO2D
      REAL(8) :: PO4DX, ALKDX, DICDX, DO2DX
      REAL(8) :: CO2D,  HCO3D, CO32D, PCO2D, O2SATD, AOUD
      REAL(8) :: PCO2A, PCO2AX
      REAL(8) :: TOCINI, TOCFIN

      LOGICAL :: AOGE
      LOGICAL :: OFIXDB

      DATA DT    / 8.64D4 /
      DATA DELTA / 1.0D-6 /

      DATA CV1 / 1.0250D+3 / !! [mol/kg] to [mol/m^3]
      DATA CV2 / 9.7561D-4 / !! [mol/m^3] to [mol/kg]
      DATA CV3 / 1.0000D+6 / !! [atm] to [ppmv]
      DATA CV4 / 3.1536D+7 / !! [yr] to [s]
      DATA CV5 / 8.6400D+4 / !! [day] to [s]

!--   Temperature [K]
!--   Salinity [PSU]
!     DATA TEMS, TEMD / 2.0, 2.0 /
!     DATA TEMS, TEMD / 15.0, 2.0 /
      DATA TEMS, TEMD / 20.0, 2.0 /
!     DATA TEMS, TEMD / 25.0, 2.0 /
      DATA SALS,  SALD  / 34.7, 34.7 /

!--   Ocean volume [m^3]
!--   Ocean area [m^2]
!--   Atmospheric molers [mol] = [mol/atm]
      DATA VOCN, AOCN, VATM / 1.292D18, 3.49D14, 1.773D20 /

!--   Ocean depth [m]
      DATA ZOCNS, ZOCND / 250.D0, 4000.D0 /

!--   Transport [m^3/s]
!--   Diffusivity [m^3/s]
!     DATA TRAN, FSD / 0.0D0, 6.0D7 /
!     DATA TRAN, FSD / 0.0D0, 0.0D7 /
!     DATA TRAN, FSD / 0.0D0, 8.193D7 / !! tau_mix =  500 yr
      DATA TRAN, FSD / 0.0D0, 4.096D7 / !! tau_mix = 1000 yr
!     DATA TRAN, FSD / 0.0D0, 2.731D7 / !! tau_mix = 1500 yr
!     1 / 1000              [1/yr]
!     1 / 1000 * VOCN       [m^3/yr]
!     1 / 1000 * VOCN / CV4 [m^3/s]

!--   Piston velocity [m/day]
      DATA PVS / 3.0D0 /

!--   Bio-production efficiency [/yr]
      R = 1.0D0 / CV4

!--   Proportional to an annual averaged solar radiation
      DATA LF / 5.0D-1 /

!--   Half-saturation constant [mol/kg]
      DATA HSC / 2.0D-8 /

!--   Rain ratio, C/P, N/P, O2/P
!     DATA RRC  / 0.25D0 /
!     DATA RRC  / 0.100D0 / !! POC:CaCO3 = 10.0:1
      DATA RRC  / 0.288D0 / !! POC:CaCO3 =  3.5:1
!     DATA RRC  / 0.667D0 / !! POC:CaCO3 =  1.5:1
      DATA RCP  / 106.D0 /
      DATA RNP  /  16.D0 /
      DATA RO2P / 172.D0 /

!--   Biological sinking flux [molC/m^2/yr]
      DATA CEPS / 0.611D0 / !! standard
!     DATA CEPS / 0.000D0 / !! minimum
!     DATA CEPS / 0.900D0 / !! maximum
!     DATA CEPS / 1.000D0 / !!  4.188 [PgC/yr]
!     DATA CEPS / 2.000D0 / !!  8.376 [PgC/yr]
!     DATA CEPS / 2.388D0 / !! 10.000 [PgC/yr]
!     10                      [PgC/yr]
!     10 / 12 /               [PmolC/yr]
!     10 / 12 / AOCN          [PmolC/m^2/yr]
!     10 / 12 / AOCN * 1.D+15 [molC/m^2/yr]

!--   Initial atomospheric PCO2 [ppmv]
      DATA PCO2A / 280.D0 /

!--   Air-sea gas exchange ?
!     DATA AOGE / .TRUE. /
      DATA AOGE / .FALSE. /

!--   fix Deep box ?
      DATA OFIXDB / .TRUE. /
!     DATA OFIXDB / .FALSE. /

!-- Initial condition (etc)

!--   Ocean area [m^2]
      AOCNS = AOCN * 1.0D0

!--   Ocean volume [m^3]
      VOCNS = AOCNS * ZOCNS
      VOCND = VOCN - VOCNS

!--   Biological sinking flux [molP/yr]
      EPS = (CEPS / RCP) * AOCNS

!--   Gas exchange [m^3/day]
      FAS = PVS * AOCNS

!--  Initial condition (biogeochemical parameter)

!--   Phosphate [mol/kg]
      DATA PO4S, PO4D / 2.20D-6, 2.20D-6 /

!--   Alkalinity [eq/kg]
      DATA ALKS, ALKD / 2.371D-3, 2.371D-3 /

!--   Dissolved inorganic carbon [mol/kg]
      DATA DICS, DICD / 2.258D-3, 2.258D-3 /

!--   Dissolved oxygen [mol/kg]
      DATA DO2S, DO2D / 1.70D-4, 1.70D-4 /

!-- Unit conversion

      EPS   = EPS   / CV4 !! [/yr]    to [/s]
      FAS   = FAS   / CV5 !! [/day]   to [/s]
      PO4S  = PO4S  * CV1 !! [mol/kg] to [mol/m^3]
      PO4D  = PO4D  * CV1 !! [mol/kg] to [mol/m^3]
      DICS  = DICS  * CV1 !! [mol/kg] to [mol/m^3]
      DICD  = DICD  * CV1 !! [mol/kg] to [mol/m^3]
      ALKS  = ALKS  * CV1 !! [mol/kg] to [mol/m^3]
      ALKD  = ALKD  * CV1 !! [mol/kg] to [mol/m^3]
      DO2S  = DO2S  * CV1 !! [mol/kg] to [mol/m^3]
      DO2D  = DO2D  * CV1 !! [mol/kg] to [mol/m^3]
      PCO2A = PCO2A / CV3 !! [ppmv]   to [atm]

!-- Loop

      DO I = 1, 100000

!-- Calculation of equilibrium constant

         CALL CHEMEQCONST(
                           TEMS, SALS,                                 &
                           BTS, K0S, K1S, K2S, KBS, KWS, FHS           &
                         )
         CALL CHEMEQCONST(
                           TEMD, SALD,                                 &
                           BTD, K0D, K1D, K2D, KBD, KWD, FHD           &
                         )

!-- Initial CO2 and PCO2
!-- Initial total carbon

         IF ( I == 1 ) THEN
            CALL CO2_NIBUN(                                            &
     &                      BTS, K0S, K1S, K2S, KBS, KWS,              &
     &                      ALKS, DICS, CO2S, HCO3S, CO32S, PCO2S      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BTS, K0S, K1S, K2S, KBS, KWS,              &
     &                      ALKD, DICD, CO2D, HCO3D, CO32D, PCO2D      &
     &                    )
            TOCINI = DICS * VOCNS + DICD * VOCND + PCO2A * VATM
         ENDIF

!-- Calculation of biogeochemical tracers

!        EPS = R * DE_S * LF * PO_S * (PO_S / (HSC + PO_S))

         TEMSX = TEMS
         TEMDX = TEMD
!        TEMSX = TEMS                                                  &
!    &          + (                                                    &
!    &              + (TRAN + FSD) * (TEMD - TEMS)                     &
!    &            )                                                    &
!    &            * (DT / VOCNS)
!        TEMDX = TEMD                                                  &
!    &          + (                                                    &
!    &              + (TRAN + FSD) * (TEMS - TEMD)                     &
!    &            )                                                    &
!    &            * (DT / VOCND)

         PO4SX = PO4S                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (PO4D - PO4S)                      &
     &             + (- EPS)                                           &
     &           )                                                     &
     &           * (DT / VOCNS)
         IF (.NOT.OFIXDB) THEN
         PO4DX = PO4D                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (PO4S - PO4D)                      &
     &             + (+ EPS)                                           &
     &           )                                                     &
     &           * (DT / VOCND)
         ELSEIF (OFIXDB) THEN
         PO4DX = PO4D
         ENDIF

         ALKSX = ALKS                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (ALKD - ALKS)                      &
     &             + ( + 2.0D0 * RRC * RCP - RNP) * (- EPS)            &
     &           )                                                     &
     &           * (DT / VOCNS)
         IF (.NOT.OFIXDB) THEN
         ALKDX = ALKD                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (ALKS - ALKD)                      &
     &             + ( + 2.0D0 * RRC * RCP - RNP) * (+ EPS)            &
     &           )                                                     &
     &           * (DT / VOCND)
         ELSEIF (OFIXDB) THEN
         ALKDX = ALKD
         ENDIF

         IF (AOGE) THEN
         DICSX = DICS                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (DICD - DICS)                      &
     &             + (1.0D0 + RRC) * RCP * (- EPS)                     &
     &             + FAS * CV1 * K0S * (PCO2A - PCO2S)                 &
     &           )                                                     &
     &           * (DT / VOCNS)
         ELSEIF (.NOT.AOGE) THEN
         DICSX = DICS                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (DICD - DICS)                      &
     &             + (1.0D0 + RRC) * RCP * (- EPS)                     &
     &           )                                                     &
     &           * (DT / VOCNS)
         ENDIF
         IF (.NOT.OFIXDB) THEN
         DICDX = DICD                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (DICS - DICD)                      &
     &             + (1.0D0 + RRC) * RCP * (+ EPS)                     &
     &           )                                                     &
     &           * (DT / VOCND)
         ELSEIF (OFIXDB) THEN
         DICDX = DICD
         ENDIF

         DO2SX = DO2S                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (DO2D - DO2S)                      &
     &             - RO2P * (- EPS)                                    &
     &           )                                                     &
     &           * (DT / VOCNS)
         IF (.NOT.OFIXDB) THEN
         DO2DX = DO2D                                                  &
     &         + (                                                     &
     &             + (TRAN + FSD) * (DO2S - DO2D)                      &
     &             - RO2P * (+ EPS)                                    &
     &           )                                                     &
     &           * (DT / VOCND)
         ELSEIF (OFIXDB) THEN
         DO2DX = DO2D
         ENDIF

!-- Calculation of carbonate system

         CALL CO2_NIBUN(                                               &
     &                   BTS, K0S, K1S, K2S, KBS, KWS,                 &
     &                   ALKSX, DICSX, CO2S, HCO3S, CO32S, PCO2S       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BTD, K0D, K1D, K2D, KBD, KWD,                 &
     &                   ALKDX, DICDX, CO2D, HCO3D, CO32D, PCO2D       &
     &                 )

!-- Calculation of atmospheric PCO2

         IF (AOGE) THEN
         PCO2AX = PCO2A                                                &
     &          + (                                                    &
     &              + K0S * FAS * CV1* (PCO2S - PCO2A)                 &
     &            )                                                    &
     &            * (DT / VATM)
         ELSEIF (.NOT.AOGE) THEN
         PCO2AX = PCO2A
         ENDIF


!-- Equilibrium condition

!        WRITE(*,*) I, ABS((PO4SX/PO4S) - 1.0D0)
!        IF (ABS((PO4SX / PO4S) - 1.0D0) <= DELTA) THEN
!          WRITE(*,*) 'TIMESTEP =', I
!          EXIT
!        ENDIF

!-- Next loop

         TEMS  = TEMSX
         TEMD  = TEMDX
         PO4S  = PO4SX
         PO4D  = PO4DX
         DICS  = DICSX
         DICD  = DICDX
         ALKS  = ALKSX
         ALKD  = ALKDX
         DO2S  = DO2SX
         DO2D  = DO2DX
         PCO2A = PCO2AX

      ENDDO

!-- Final total carbon

      TOCFIN = DICS * VOCNS + DICD * VOCND + PCO2A * VATM

!-- Oxygen utilization

      CALL O2SAT(TEMS, SALS, O2SATS)
      CALL O2SAT(TEMD, SALD, O2SATD)

      AOUS = (O2SATS - DO2S)
      AOUD = (O2SATD - DO2D)

!-- Unit conversion

      EPS   = CV4 * EPS * RCP * 1.D-15 * 12.0  !! [molP/s]  to [PgC/yr]
      PO4S  = CV2 * 1.D+6 * PO4S               !! [mol/m^3] to [umol/kg]
      PO4D  = CV2 * 1.D+6 * PO4D               !! [mol/m^3] to [umol/kg]
      DICS  = CV2 * 1.D+6 * DICS               !! [mol/m^3] to [umol/kg]
      DICD  = CV2 * 1.D+6 * DICD               !! [mol/m^3] to [umol/kg]
      ALKS  = CV2 * 1.D+6 * ALKS               !! [mol/m^3] to [umol/kg]
      ALKD  = CV2 * 1.D+6 * ALKD               !! [mol/m^3] to [umol/kg]
      DO2S  = CV2 * 1.D+6 * DO2S               !! [mol/m^3] to [umol/kg]
      DO2D  = CV2 * 1.D+6 * DO2D               !! [mol/m^3] to [umol/kg]
      AOUS  = CV2 * 1.D+6 * AOUS               !! [mol/m^3] to [umol/kg]
      AOUD  = CV2 * 1.D+6 * AOUD               !! [mol/m^3] to [umol/kg]
      CO2S  = CV2 * 1.D+6 * CO2S               !! [mol/m^3] to [umol/kg]
      CO2D  = CV2 * 1.D+6 * CO2D               !! [mol/m^3] to [umol/kg]
      HCO3S = CV2 * 1.D+6 * HCO3S              !! [mol/m^3] to [umol/kg]
      HCO3D = CV2 * 1.D+6 * HCO3D              !! [mol/m^3] to [umol/kg]
      CO32S = CV2 * 1.D+6 * CO32S              !! [mol/m^3] to [umol/kg]
      CO32D = CV2 * 1.D+6 * CO32D              !! [mol/m^3] to [umol/kg]
      PCO2S = CV3 * PCO2S                      !! [atm]     to [ppmv]
      PCO2D = CV3 * PCO2D                      !! [atm]     to [ppmv]
      PCO2A = CV3 * PCO2A                      !! [atm]     to [ppmv]

!-- Output

      WRITE(*,*) 'Temperature (S)          :', TEMS, 'K'
      WRITE(*,*) 'Temperature (D)          :', TEMD, 'K'
      WRITE(*,*) ''
      WRITE(*,*) 'Salinity (S)             :', SALS, 'PSU'
      WRITE(*,*) 'Salinity (D)             :', SALD, 'PSU'
      WRITE(*,*) ''
      WRITE(*,*) 'EP (S)                   :', EPS, 'PgC/yr'
      WRITE(*,*) ''
      WRITE(*,*) 'PO4 (S)                  :', PO4S, 'umol/kg'
      WRITE(*,*) 'PO4 (D)                  :', PO4D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DIC (S)                  :', DICS, 'umol/kg'
      WRITE(*,*) 'DIC (D)                  :', DICD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'ALK (S)                  :', ALKS, 'ueq/kg'
      WRITE(*,*) 'ALK (D)                  :', ALKD, 'ueq/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DO2 (S)                  :', DO2S, 'umol/kg'
      WRITE(*,*) 'DO2 (D)                  :', DO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'AOU (S)                  :', AOUS, 'umol/kg'
      WRITE(*,*) 'AOU (D)                  :', AOUD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO2 (S)                  :', CO2S, 'umol/kg'
      WRITE(*,*) 'CO2 (D)                  :', CO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'HCO3 (S)                 :', HCO3S, 'umol/kg'
      WRITE(*,*) 'HCO3 (D)                 :', HCO3D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO32 (S)                 :', CO32S, 'umol/kg'
      WRITE(*,*) 'CO32 (D)                 :', CO32D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'PCO2 (S)                 :', PCO2S, 'ppmv'
      WRITE(*,*) 'PCO2 (D)                 :', PCO2D, 'ppmv'
      WRITE(*,*) 'PCO2 (ATM)               :', PCO2A, 'ppmv'
      WRITE(*,*) ''
!     WRITE(*,*) 'BTS',BTS
!     WRITE(*,*) 'K0S',K0S
!     WRITE(*,*) 'K1S',K1S
!     WRITE(*,*) 'K2S',K2S
!     WRITE(*,*) 'KBS',KBS
!     WRITE(*,*) 'KWS',KWS
!     WRITE(*,*) ''
!     WRITE(*,*) 'BTD',BTD
!     WRITE(*,*) 'K0D',K0D
!     WRITE(*,*) 'K1D',K1D
!     WRITE(*,*) 'K2D',K2D
!     WRITE(*,*) 'KBD',KBD
!     WRITE(*,*) 'KWD',KWD
!     WRITE(*,*) ''
      WRITE(*,*) 'Total Carbon (Initial)   :', TOCINI, 'mol'
      WRITE(*,*) 'Total Carbon (Final)     :', TOCFIN, 'mol'
      WRITE(*,*) ''

      END PROGRAM
