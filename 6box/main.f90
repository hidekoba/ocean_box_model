      PROGRAM MAIN

      USE CHEMEQ
      USE O2
      USE CO2

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Ocean 6-box model
!  Calculation of phosphate, DIC, alkalinity, CO2, and PCO2
!
!  References: Toggweiler et al. (1999, PA)
!
! ---------------------------------------------------------------------

      INTEGER :: I
      INTEGER, PARAMETER :: NBOX=6
!     1:P, 2:S, 3:L, 4:N, 5:M, 6:D

      REAL(8) :: DT, DELTA
      REAL(8) :: CV1, CV2, CV3, CV4, CV5

      REAL(8) :: VOCN, AOCN, VATM
      REAL(8) :: ZOCNP, DEP, AOCNP, FAOCNP, VOCNP
      REAL(8) :: ZOCNS, DES, AOCNS, FAOCNS, VOCNS
      REAL(8) :: ZOCNL, DEL, AOCNL, FAOCNL, VOCNL
      REAL(8) :: ZOCNN, DEN, AOCNN, FAOCNN, VOCNN
      REAL(8) :: ZOCNM,      AOCNM, FAOCNM, VOCNM
      REAL(8) :: ZOCND,                     VOCND

      REAL(8) :: TRAN, FPD, FND, FSM, FLM, FPS, FSL, FLN
      REAL(8) :: PVP, PVS, PVL, PVN
      REAL(8) :: FAP, FAS, FAL, FAN
      REAL(8) :: CEPP, CEPS, CEPN
      REAL(8) :: R, LF, HSC
      REAL(8) :: RRC, RCP, RNP, RO2P, GM, RGM
      REAL(8) :: BT(NBOX)
      REAL(8) :: K0(NBOX), K1(NBOX), K2(NBOX)
      REAL(8) :: KB(NBOX), KW(NBOX), FH(NBOX)

      REAL(8) :: TEMP,  SALP,   TEMPX!  SALPX
      REAL(8) :: PO4P,  ALKP,   DICP,   DO2P
      REAL(8) :: PO4PX, ALKPX,  DICPX,  DO2PX
      REAL(8) :: CO2P,  HCO3P,  CO32P,  PCO2P, O2SATP, AOUP, EPP

      REAL(8) :: TEMS,  SALS,   TEMSX!  SALSX
      REAL(8) :: PO4S,  ALKS,   DICS,   DO2S
      REAL(8) :: PO4SX, ALKSX,  DICSX,  DO2SX
      REAL(8) :: CO2S,  HCO3S,  CO32S,  PCO2S, O2SATS, AOUS, EPS

      REAL(8) :: TEML,  SALL,   TEMLX!  SALLX
      REAL(8) :: PO4L,  ALKL,   DICL,   DO2L
      REAL(8) :: PO4LX, ALKLX,  DICLX,  DO2LX
      REAL(8) :: CO2L,  HCO3L,  CO32L,  PCO2L, O2SATL, AOUL, EPL

      REAL(8) :: TEMN,  SALN,   TEMNX!  SALNX
      REAL(8) :: PO4N,  ALKN,   DICN,   DO2N
      REAL(8) :: PO4NX, ALKNX,  DICNX,  DO2NX
      REAL(8) :: CO2N,  HCO3N,  CO32N,  PCO2N, O2SATN, AOUN, EPN

      REAL(8) :: TEMM,  SALM,   TEMMX!  SALMX
      REAL(8) :: PO4M,  ALKM,   DICM,   DO2M
      REAL(8) :: PO4MX, ALKMX,  DICMX,  DO2MX
      REAL(8) :: CO2M,  HCO3M,  CO32M,  PCO2M, O2SATM, AOUM

      REAL(8) :: TEMD,  SALD,   TEMDX!  SALDX
      REAL(8) :: PO4D,  ALKD,   DICD,   DO2D
      REAL(8) :: PO4DX, ALKDX,  DICDX,  DO2DX
      REAL(8) :: CO2D,  HCO3D,  CO32D,  PCO2D, O2SATD, AOUD

      REAL(8) :: PCO2A, PCO2AX
      REAL(8) :: TOCINI, TOCFIN
      REAL(8) :: OLD_PO4P

      DATA DT    / 8.64D4 /
      DATA DELTA / 1.0D-6 /

      DATA CV1 / 1.0250D+3 / !! [mol/kg] to [mol/m^3]
      DATA CV2 / 9.7561D-4 / !! [mol/m^3] to [mol/kg]
      DATA CV3 / 1.0000D+6 / !! [atm] to [ppmv]
      DATA CV4 / 3.1536D+7 / !! [yr] to [s]
      DATA CV5 / 8.6400D+4 / !! [day] to [s]

!--   Temperature [K]
!--   Salinity [PSU]
      DATA TEMP, TEMS, TEML, TEMN, TEMM, TEMD &
     & / 1.0, 8.0, 21.5, 3.0, 10.0, 1.75 / !! Table 1,2,3 in Toggweiler+ (1999)
      DATA SALP, SALS, SALL, SALN, SALM, SALD &
     & / 34.7, 34.7, 34.7, 34.7, 34.7, 34.7 / !! Table 1,2,3 in Toggweiler+ (1999)

!--   Ocean volume [m^3]
!--   Ocean area [m^2]
!--   Atmospheric molers [mol] = [mol/atm]
      DATA VOCN, AOCN, VATM &
     & / 1.292D18, 3.49D14, 1.773D20 / !! Table 1 in Toggweiler+ (1999)

!--   Ocean depth [m]
      DATA ZOCNP, ZOCNS, ZOCNL, ZOCNN, ZOCNM, ZOCND &
     & / 250.D0, 250.D0, 100.D0, 250.D0, 900.0D0, 4000.D0 / !! Table 1,2,3 in Toggweiler+ (1999)

!--   Ocean depth of euphotic zone [m]
      DATA DEP, DES, DEL, DEN &
     & / 100.D0, 100.D0, 100.D0, 100.D0 /

!--   Ocean area [%]
      DATA FAOCNP, FAOCNS, FAOCNL, FAOCNN, FAOCNM &
     & / 0.05D0, 0.10D0, 0.75D0, 0.10D0, 0.85D0 / !! Table 1,2,3 in Toggweiler et al. (1999)

!--   Transport [m^3/s]
!--   Diffusivity [m^3/s]
      DATA TRAN &
     & / 20.0D6 /
      DATA FPD, FND, FSM, FLM, FLN &
     & / 60.0D6, 0.0D0, 0.0D0, 40.0D6, 0.0D0 / !! 3-300 Sv / Table 1,2,3 in Toggweiler et al. (1999)
      DATA FPS, FSL, FLN &
     & / 0.0D0, 0.0D0, 0.0D0 / !! Table 1,2,3 in Toggweiler et al. (1999)

!--   Piston velocity [m/day]
      DATA PVP, PVS, PVL, PVN &
     & / 3.0D0, 3.0D0, 3.0D0, 3.0D0 / !! Table 1,2,3 in Toggweiler et al. (1999)

!--   Bio-production efficiency [/yr]
      DATA R / 1.0D0 /

!--   Proportional to an annual averaged solar radiation
      DATA LF / 5.0D-1 /

!--   Half-saturation constant [mol/kg]
      DATA HSC / 2.0D-8 /

!--   Rain ratio, C/P, N/P, O2/P
      DATA RRC  / 0.25D0 /
     !DATA RCP  / 106.0D0 /
     !DATA RNP  /  16.0D0 /
     !DATA RO2P / 172.0D0 /
      DATA RCP  / 162.5D0 / !! Table 1 in Toggweiler+ (1999)
      DATA RNP  /  15.0D0 / !! Table 1 in Toggweiler+ (1999)
      DATA RO2P / 169.0D0 / !! Table 1 in Toggweiler+ (1999)

!--   Remineralization fraction
      DATA GM / 0.80D0 /
      RGM = 1.0D0 - GM

!--   Biological sinking flux [molC/m^2/yr]
      DATA CEPP, CEPS, CEPN &
!    & / 1.00D0, 1.50D0, 3.00D0 / !! 0.15-15.0 molC/m^2/yr / Table 3 in Toggweiler+ (1999)
     & / 1.00D0, 3.00D0, 3.00D0 /

!--   Initial atomospheric PCO2 [ppmv]
      DATA PCO2A / 280.D0 /

!-- Initial condition (etc)

!--   Ocean area [m^2]
      AOCNP = AOCN * FAOCNP
      AOCNS = AOCN * FAOCNS
      AOCNL = AOCN * FAOCNL
      AOCNN = AOCN * FAOCNN
      AOCNM = AOCN * FAOCNM

!--   Ocean volume [m^3]
      VOCNP = AOCNP * ZOCNP
      VOCNS = AOCNS * ZOCNS
      VOCNL = AOCNL * ZOCNL
      VOCNN = AOCNN * ZOCNN
      VOCNM = AOCNM * ZOCNM - VOCNS - VOCNL
      VOCND = VOCN - VOCNP - VOCNS - VOCNL - VOCNN - VOCNM

!--   Biological sinking flux [molP/yr]
      EPP = (CEPP / RCP) * AOCNP
      EPS = (CEPS / RCP) * AOCNS
      EPN = (CEPN / RCP) * AOCNN

!--   Gas exchange [m^3/day]
      FAP = PVP * AOCNP
      FAS = PVS * AOCNS
      FAL = PVL * AOCNL
      FAN = PVN * AOCNN

!--  Initial condition (biogeochemical parameter)

!--   Phosphate [mol/kg]
      DATA PO4P, PO4S, PO4L, PO4N, PO4M, PO4D &
     & / 2.09D-6, 2.09D-6, 2.09D-6, 2.09D-6, 2.09D-6, 2.09D-6 /
!    & / 2.00D-6, 2.33D-6, 0.0D-6, 1.0D-6, 2.142D-6, 2.142D-6 /

!--   Alkalinity [eq/kg]
      DATA ALKP, ALKS, ALKL, ALKN, ALKM, ALKD &
     & / 2.371D-3, 2.371D-3, 2.371D-3, 2.371D-3, 2.371D-3, 2.371D-3 /
!    & / 2.360D-3, 2.325D-3, 2.277D-3, 2.324D-3, 2.3738D-3, 2.3738D-3 /

!--   Dissolved inorganic carbon [mol/kg]
      DATA DICP, DICS, DICL, DICN, DICM, DICD &
     & / 2.258D-3, 2.258D-3, 2.258D-3, 2.258D-3, 2.258D-3, 2.258D-3 /
!    & / 2.150D-3, 2.107D-3, 1.933D-3, 2.106D-3, 2.2561D-3, 2.2561D-3 /

!--   Dissolved oxygen [mol/kg]
      DATA DO2P, DO2S, DO2L, DO2N, DO2M, DO2D &
     & / 1.60D-4, 1.60D-4, 1.60D-4, 1.60D-4, 1.60D-4, 1.60D-4 /
!    & / 3.50D-4, 3.00D-4, 3.00D-4, 3.00D-4, 1.50D-4, 2.50D-4 /

!-- Unit conversion

      EPP   = EPP   / CV4 !! [/yr]    to [/s]
      EPS   = EPS   / CV4 !! [/yr]    to [/s]
      EPN   = EPN   / CV4 !! [/yr]    to [/s]
      FAP   = FAP   / CV5 !! [/day]   to [/s]
      FAS   = FAS   / CV5 !! [/day]   to [/s]
      FAL   = FAL   / CV5 !! [/day]   to [/s]
      FAN   = FAN   / CV5 !! [/day]   to [/s]
      HSC   = HSC   * CV1 !! [mol/kg] to [mol/m^3]
      R     = R     / CV4 !! [/yr]    to [/s]
      PO4P  = PO4P  * CV1 !! [mol/kg] to [mol/m^3]
      PO4S  = PO4S  * CV1 !! [mol/kg] to [mol/m^3]
      PO4L  = PO4L  * CV1 !! [mol/kg] to [mol/m^3]
      PO4N  = PO4N  * CV1 !! [mol/kg] to [mol/m^3]
      PO4M  = PO4M  * CV1 !! [mol/kg] to [mol/m^3]
      PO4D  = PO4D  * CV1 !! [mol/kg] to [mol/m^3]
      DICP  = DICP  * CV1 !! [mol/kg] to [mol/m^3]
      DICS  = DICS  * CV1 !! [mol/kg] to [mol/m^3]
      DICL  = DICL  * CV1 !! [mol/kg] to [mol/m^3]
      DICN  = DICN  * CV1 !! [mol/kg] to [mol/m^3]
      DICM  = DICM  * CV1 !! [mol/kg] to [mol/m^3]
      DICD  = DICD  * CV1 !! [mol/kg] to [mol/m^3]
      ALKP  = ALKP  * CV1 !! [mol/kg] to [mol/m^3]
      ALKS  = ALKS  * CV1 !! [mol/kg] to [mol/m^3]
      ALKL  = ALKL  * CV1 !! [mol/kg] to [mol/m^3]
      ALKN  = ALKN  * CV1 !! [mol/kg] to [mol/m^3]
      ALKM  = ALKM  * CV1 !! [mol/kg] to [mol/m^3]
      ALKD  = ALKD  * CV1 !! [mol/kg] to [mol/m^3]
      DO2P  = DO2P  * CV1 !! [mol/kg] to [mol/m^3]
      DO2S  = DO2S  * CV1 !! [mol/kg] to [mol/m^3]
      DO2L  = DO2L  * CV1 !! [mol/kg] to [mol/m^3]
      DO2N  = DO2N  * CV1 !! [mol/kg] to [mol/m^3]
      DO2M  = DO2M  * CV1 !! [mol/kg] to [mol/m^3]
      DO2D  = DO2D  * CV1 !! [mol/kg] to [mol/m^3]
      PCO2A = PCO2A / CV3 !! [ppmv]   to [atm]

!-- Loop

      DO I = 1, 100000

!-- Calculation of equilibrium constant
!   1:P, 2:S, 3:L, 4:N, 5:M, 6:D

         CALL CHEMEQCONST(                                             &
     &                     TEMP, SALP,                                 &
     &                     BT(1), K0(1), K1(1), K2(1),                 &
     &                     KB(1), KW(1), FH(1)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMS, SALS,                                 &
     &                     BT(2), K0(2), K1(2), K2(2),                 &
     &                     KB(2), KW(2), FH(2)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEML, SALL,                                 &
     &                     BT(3), K0(3), K1(3), K2(3),                 &
     &                     KB(3), KW(3), FH(3)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMN, SALN,                                 &
     &                     BT(4), K0(4), K1(4), K2(4),                 &
     &                     KB(4), KW(4), FH(4)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMM, SALM,                                 &
     &                     BT(5), K0(5), K1(5), K2(5),                 &
     &                     KB(5), KW(5), FH(5)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMD, SALD,                                 &
     &                     BT(6), K0(6), K1(6), K2(6),                 &
     &                     KB(6), KW(6), FH(6)                         &
     &                   )

!-- Initial CO2 and PCO2
!-- Initial total carbon

         IF (I == 1) THEN
            CALL CO2_NIBUN(                                            &
     &                      BT(1), K0(1), K1(1), K2(1), KB(1), KW(1),  &
     &                      ALKP, DICP, CO2P, HCO3P, CO32P, PCO2P      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(2), K0(2), K1(2), K2(2), KB(2), KW(2),  &
     &                      ALKS, DICS, CO2S, HCO3S, CO32S, PCO2S      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(3), K0(3), K1(3), K2(3), KB(3), KW(3),  &
     &                      ALKL, DICL, CO2L, HCO3L, CO32L, PCO2L      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(4), K0(4), K1(4), K2(4), KB(4), KW(4),  &
     &                      ALKN, DICN, CO2N, HCO3N, CO32N, PCO2N      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(5), K0(5), K1(5), K2(5), KB(5), KW(5),  &
     &                      ALKM, DICM, CO2M, HCO3M, CO32M, PCO2M      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(6), K0(6), K1(6), K2(6), KB(6), KW(6),  &
     &                      ALKD, DICD, CO2D, HCO3D, CO32D, PCO2D      &
     &                    )

            CALL O2SAT(TEMP, SALP, O2SATP)
            CALL O2SAT(TEMS, SALS, O2SATS)
            CALL O2SAT(TEML, SALL, O2SATL)
            CALL O2SAT(TEMN, SALN, O2SATN)
            CALL O2SAT(TEMM, SALM, O2SATM)
            CALL O2SAT(TEMD, SALD, O2SATD)

            TOCINI = + DICP  * VOCNP                                   &
     &               + DICS  * VOCNS                                   &
     &               + DICL  * VOCNL                                   &
     &               + DICN  * VOCNN                                   &
     &               + DICM  * VOCNM                                   &
     &               + DICD  * VOCND                                   &
     &               + PCO2A * VATM
         ENDIF

!-- Calculation of biogeochemical tracers

!        EPL = + FSL * (PO4S - PO4L)                                   &
!    &         + FLN * (PO4N - PO4L)                                   &
!    &         + FLM * (PO4M - PO4L)
         EPL = R * DEL * LF * PO4L * (PO4L / (HSC + PO4L)) * VOCNL

         TEMPX = TEMP
         TEMSX = TEMS
         TEMLX = TEML
         TEMNX = TEMN
         TEMMX = TEMM
         TEMDX = TEMD

         PO4PX = PO4P                                                  &
     &         + (                                                     &
     &             + (TRAN + FPD) * (PO4D - PO4P)                      &
     &             + FPS * (PO4S - PO4P)                               &
     &             - EPP                                               &
     &           )                                                     &
     &           * (DT / VOCNP)
         PO4SX = PO4S                                                  &
     &         + (                                                     &
     &             + (TRAN + FPS) * (PO4P - PO4S)                      &
     &             + FSL * (PO4L - PO4S)                               &
     &             - EPS                                               &
     &           )                                                     &
     &           * (DT / VOCNS)
         PO4LX = PO4L                                                  &
     &         + (                                                     &
     &             + FSL * (PO4S - PO4L)                               &
     &             + FLN * (PO4N - PO4L)                               &
     &             + FLM * (PO4M - PO4L)                               &
     &             - EPL                                               &
     &           )                                                     &
     &           * (DT / VOCNL)
         PO4NX = PO4N                                                  &
     &         + (                                                     &
     &             + TRAN * (PO4M - PO4N)                              &
     &             + FLN * (PO4L - PO4N)                               &
     &             + FND * (PO4D - PO4N)                               &
     &             - EPN                                               &
     &           )                                                     &
     &           * (DT / VOCNN)
         PO4MX = PO4M                                                  &
     &         + (                                                     &
     &             + (TRAN + FSM) * (PO4S - PO4M)                      &
     &             + FLM * (PO4L - PO4M)                               &
     &             + (GM * (EPS + EPL))                                &
     &           )                                                     &
     &           * (DT / VOCNM)
         PO4DX = PO4D                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (PO4N - PO4D)                      &
     &             + FPD * (PO4P - PO4D)                               &
     &             + (RGM * (EPS + EPL) + (EPP + EPN))                 &
     &           )                                                     &
     &           * (DT / VOCND)

         ALKPX = ALKP                                                  &
     &         + (                                                     &
     &             + (TRAN + FPD) * (ALKD - ALKP)                      &
     &             + FPS * (ALKS - ALKP)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPP                   &
     &           )                                                     &
     &           * (DT / VOCNP)
         ALKSX = ALKS                                                  &
     &         + (                                                     &
     &             + (TRAN + FPS) * (ALKP - ALKS)                      &
     &             + FSL * (ALKL - ALKS)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPS                   &
     &           )                                                     &
     &           * (DT / VOCNS)
         ALKLX = ALKL                                                  &
     &         + (                                                     &
     &             + FSL * (ALKS - ALKL)                               &
     &             + FLN * (ALKN - ALKL)                               &
     &             + FLM * (ALKM - ALKL)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPL                   &
     &           )                                                     &
     &           * (DT / VOCNL)
         ALKNX = ALKN                                                  &
     &         + (                                                     &
     &             + TRAN * (ALKM - ALKN)                              &
     &             + FLN * (ALKL - ALKN)                               &
     &             + FND * (ALKD - ALKN)                               &
     &             - (2.0D0 * RRC * RCP - RNP) * EPN                   &
     &           )                                                     &
     &           * (DT / VOCNN)
         ALKMX = ALKM                                                  &
     &         + (                                                     &
     &             + (TRAN + FSM) * (ALKS - ALKM)                      &
     &             + FLM * (ALKL - ALKM)                               &
     &             + (2.0D0 * RRC * RCP - RNP) * (GM * (EPS + EPL))    &
     &           )                                                     &
     &           * (DT / VOCNM)
         ALKDX = ALKD                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (ALKN - ALKD)                      &
     &             + FPD * (ALKP - ALKD)                               &
     &             + (2.0D0 * RRC * RCP - RNP)                         &
     &               * (RGM * (EPS + EPL) + (EPP + EPN))               &
     &           )                                                     &
     &           * (DT / VOCND)

         DICPX = DICP                                                  &
     &         + (                                                     &
     &             + (TRAN + FPD) * (DICD - DICP)                      &
     &             + FPS * (DICS - DICP)                               &
     &             - (1.0D0 + RRC) * RCP * EPP                         &
     &           )                                                     &
     &           * (DT / VOCNP)
         DICSX = DICS                                                  &
     &         + (                                                     &
     &             + (TRAN + FPS) * (DICP - DICS)                      &
     &             + FSL * (DICL - DICS)                               &
     &             - (1.0D0 + RRC) * RCP * EPS                         &
     &           )                                                     &
     &           * (DT / VOCNS)
         DICLX = DICL                                                  &
     &         + (                                                     &
     &             + FSL * (DICS - DICL)                               &
     &             + FLN * (DICN - DICL)                               &
     &             + FLM * (DICM - DICL)                               &
     &             - (1.0D0 + RRC) * RCP * EPL                         &
     &           )                                                     &
     &           * (DT / VOCNL)
         DICNX = DICN                                                  &
     &         + (                                                     &
     &             + TRAN * (DICM - DICN)                              &
     &             + FLN * (DICL - DICN)                               &
     &             + FND * (DICD - DICN)                               &
     &             - (1.0D0 + RRC) * RCP * EPN                         &
     &           )                                                     &
     &           * (DT / VOCNN)
         DICMX = DICM                                                  &
     &         + (                                                     &
     &             + (TRAN + FSM) * (DICS - DICM)                      &
     &             + FLM * (DICL - DICM)                               &
     &             + (1.0D0 + RRC) * RCP * (GM * (EPS + EPL))          &
     &           )                                                     &
     &           * (DT / VOCNM)
         DICDX = DICD                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (DICN - DICD)                      &
     &             + FPD * (DICP - DICD)                               &
     &             + (1.0D0 + RRC) * RCP                               &
     &               * (RGM * (EPS + EPL) + (EPP + EPN))               &
     &           )                                                     &
     &           * (DT / VOCND)

         DO2PX = DO2P                                                  &
     &         + (                                                     &
     &             + (TRAN + FPD) * (DO2D - DO2P)                      &
     &             + FPS * (DO2S - DO2P)                               &
     &             + RO2P * EPP                                        &
     &             + FAP * (O2SATP - DO2P)                             &
     &           )                                                     &
     &           * (DT / VOCNP)
         DO2SX = DO2S                                                  &
     &         + (                                                     &
     &             + (TRAN + FPS) * (DO2P - DO2S)                      &
     &             + FSL * (DO2L - DO2S)                               &
     &             + RO2P * EPS                                        &
     &             + FAS * (O2SATS - DO2S)                             &
     &           )                                                     &
     &           * (DT / VOCNS)
         DO2LX = DO2L                                                  &
     &         + (                                                     &
     &             + FSL * (DO2S - DO2L)                               &
     &             + FLN * (DO2N - DO2L)                               &
     &             + FLM * (DO2M - DO2L)                               &
     &             + RO2P * EPL                                        &
     &             + FAL * (O2SATL - DO2L)                             &
     &           )                                                     &
     &           * (DT / VOCNL)
         DO2NX = DO2N                                                  &
     &         + (                                                     &
     &             + TRAN * (DO2M - DO2N)                              &
     &             + FLN * (DO2L - DO2N)                               &
     &             + FND * (DO2D - DO2N)                               &
     &             + RO2P * EPN                                        &
     &             + FAN * (O2SATN - DO2N)                             &
     &           )                                                     &
     &           * (DT / VOCNN)
         DO2MX = DO2M                                                  &
     &         + (                                                     &
     &             + (TRAN + FSM) * (DO2S - DO2M)                      &
     &             + FLM * (DO2L - DO2M)                               &
     &             - RO2P * (GM * (EPS + EPL))                         &
     &           )                                                     &
     &           * (DT / VOCNM)
         DO2DX = DO2D                                                  &
     &         + (                                                     &
     &             + (TRAN + FND) * (DO2N - DO2D)                      &
     &             + FPD * (DO2P - DO2D)                               &
     &             - RO2P * (RGM * (EPS + EPL) + (EPP + EPN))          &
     &           )                                                     &
     &           * (DT / VOCND)

!-- Calculation of carbonate system

         CALL CO2_NIBUN(                                               &
     &                   BT(1), K0(1), K1(1), K2(1), KB(1), KW(1),     &
     &                   ALKPX, DICPX, CO2P, HCO3P, CO32P, PCO2P       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(2), K0(2), K1(2), K2(2), KB(2), KW(2),     &
     &                   ALKSX, DICSX, CO2S, HCO3S, CO32S, PCO2S       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(3), K0(3), K1(3), K2(3), KB(3), KW(3),     &
     &                   ALKLX, DICLX, CO2L, HCO3L, CO32L, PCO2L       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(4), K0(4), K1(4), K2(4), KB(4), KW(4),     &
     &                   ALKNX, DICNX, CO2N, HCO3N, CO32N, PCO2N       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(5), K0(5), K1(5), K2(5), KB(5), KW(5),     &
     &                   ALKMX, DICMX, CO2M, HCO3M, CO32M, PCO2M       &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(6), K0(6), K1(6), K2(6), KB(6), KW(6),     &
     &                   ALKDX, DICDX, CO2D, HCO3D, CO32D, PCO2D       &
     &                 )

!-- Calculation of atmospheric PCO2

         PCO2AX = PCO2A                                                &
     &          + (                                                    &
     &              + FAP * CV1 * K0(1) * (PCO2P - PCO2A)              &
     &              + FAS * CV1 * K0(2) * (PCO2S - PCO2A)              &
     &              + FAL * CV1 * K0(3) * (PCO2L - PCO2A)              &
     &              + FAN * CV1 * K0(4) * (PCO2N - PCO2A)              &
     &            )                                                    &
     &            * (DT / VATM)

!-- Next loop

         OLD_PO4P = PO4P
         TEMP  = TEMPX
         TEMS  = TEMSX
         TEML  = TEMLX
         TEMN  = TEMNX
         TEMM  = TEMMX
         TEMD  = TEMDX
         PO4P  = PO4PX
         PO4S  = PO4SX
         PO4L  = PO4LX
         PO4N  = PO4NX
         PO4M  = PO4MX
         PO4D  = PO4DX
         ALKP  = ALKPX
         ALKS  = ALKSX
         ALKL  = ALKLX
         ALKN  = ALKNX
         ALKM  = ALKMX
         ALKD  = ALKDX
         DICP  = DICPX
         DICS  = DICSX
         DICL  = DICLX
         DICN  = DICNX
         DICM  = DICMX
         DICD  = DICDX
         DO2P  = DO2PX
         DO2S  = DO2SX
         DO2L  = DO2LX
         DO2N  = DO2NX
         DO2M  = DO2MX
         DO2D  = DO2DX
         PCO2A = PCO2AX

!-- Equilibrium condition

         IF (ABS((PO4P / OLD_PO4P) - 1.0D0) <= DELTA) THEN
            WRITE(*,*) 'TIMESTEP =', I
            EXIT
         ENDIF

      ENDDO

!-- Final total carbon

      TOCFIN = + DICP  * VOCNP                                         &
     &         + DICS  * VOCNS                                         &
     &         + DICL  * VOCNL                                         &
     &         + DICN  * VOCNN                                         &
     &         + DICM  * VOCNM                                         &
     &         + DICD  * VOCND                                         &
     &         + PCO2A * VATM

!-- Oxygen saturation

      CALL O2SAT(TEMP, SALP, O2SATP)
      CALL O2SAT(TEMS, SALS, O2SATS)
      CALL O2SAT(TEML, SALL, O2SATL)
      CALL O2SAT(TEMN, SALN, O2SATN)
      CALL O2SAT(TEMM, SALM, O2SATM)
      CALL O2SAT(TEMD, SALD, O2SATD)

      AOUP = (O2SATP - DO2P)
      AOUS = (O2SATS - DO2S)
      AOUL = (O2SATL - DO2L)
      AOUN = (O2SATN - DO2N)
      AOUM = (O2SATM - DO2M)
      AOUD = (O2SATD - DO2D)

!-- Unit conversion

      EPP   = CV4 * EPP * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      EPS   = CV4 * EPS * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      EPL   = CV4 * EPL * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      EPN   = CV4 * EPN * RCP * 1.D-15 * 12.0 !! [molP/s]  to [PgC/yr]
      PO4P  = CV2 * 1.D+6 * PO4P              !! [mol/m^3] to [umol/kg]
      PO4S  = CV2 * 1.D+6 * PO4S              !! [mol/m^3] to [umol/kg]
      PO4L  = CV2 * 1.D+6 * PO4L              !! [mol/m^3] to [umol/kg]
      PO4N  = CV2 * 1.D+6 * PO4N              !! [mol/m^3] to [umol/kg]
      PO4M  = CV2 * 1.D+6 * PO4M              !! [mol/m^3] to [umol/kg]
      PO4D  = CV2 * 1.D+6 * PO4D              !! [mol/m^3] to [umol/kg]
      DICP  = CV2 * 1.D+6 * DICP              !! [mol/m^3] to [umol/kg]
      DICS  = CV2 * 1.D+6 * DICS              !! [mol/m^3] to [umol/kg]
      DICL  = CV2 * 1.D+6 * DICL              !! [mol/m^3] to [umol/kg]
      DICN  = CV2 * 1.D+6 * DICN              !! [mol/m^3] to [umol/kg]
      DICM  = CV2 * 1.D+6 * DICM              !! [mol/m^3] to [umol/kg]
      DICD  = CV2 * 1.D+6 * DICD              !! [mol/m^3] to [umol/kg]
      ALKP  = CV2 * 1.D+6 * ALKP              !! [mol/m^3] to [umol/kg]
      ALKS  = CV2 * 1.D+6 * ALKS              !! [mol/m^3] to [umol/kg]
      ALKL  = CV2 * 1.D+6 * ALKL              !! [mol/m^3] to [umol/kg]
      ALKN  = CV2 * 1.D+6 * ALKN              !! [mol/m^3] to [umol/kg]
      ALKM  = CV2 * 1.D+6 * ALKM              !! [mol/m^3] to [umol/kg]
      ALKD  = CV2 * 1.D+6 * ALKD              !! [mol/m^3] to [umol/kg]
      DO2P  = CV2 * 1.D+6 * DO2P              !! [mol/m^3] to [umol/kg]
      DO2S  = CV2 * 1.D+6 * DO2S              !! [mol/m^3] to [umol/kg]
      DO2L  = CV2 * 1.D+6 * DO2L              !! [mol/m^3] to [umol/kg]
      DO2N  = CV2 * 1.D+6 * DO2N              !! [mol/m^3] to [umol/kg]
      DO2M  = CV2 * 1.D+6 * DO2M              !! [mol/m^3] to [umol/kg]
      DO2D  = CV2 * 1.D+6 * DO2D              !! [mol/m^3] to [umol/kg]
      AOUP  = CV2 * 1.D+6 * AOUP              !! [mol/m^3] to [umol/kg]
      AOUS  = CV2 * 1.D+6 * AOUS              !! [mol/m^3] to [umol/kg]
      AOUL  = CV2 * 1.D+6 * AOUL              !! [mol/m^3] to [umol/kg]
      AOUN  = CV2 * 1.D+6 * AOUN              !! [mol/m^3] to [umol/kg]
      AOUM  = CV2 * 1.D+6 * AOUM              !! [mol/m^3] to [umol/kg]
      AOUD  = CV2 * 1.D+6 * AOUD              !! [mol/m^3] to [umol/kg]
      CO2P  = CV2 * 1.D+6 * CO2P              !! [mol/m^3] to [umol/kg]
      CO2S  = CV2 * 1.D+6 * CO2S              !! [mol/m^3] to [umol/kg]
      CO2L  = CV2 * 1.D+6 * CO2L              !! [mol/m^3] to [umol/kg]
      CO2N  = CV2 * 1.D+6 * CO2N              !! [mol/m^3] to [umol/kg]
      CO2M  = CV2 * 1.D+6 * CO2M              !! [mol/m^3] to [umol/kg]
      CO2D  = CV2 * 1.D+6 * CO2D              !! [mol/m^3] to [umol/kg]
      HCO3P = CV2 * 1.D+6 * HCO3P             !! [mol/m^3] to [umol/kg]
      HCO3S = CV2 * 1.D+6 * HCO3S             !! [mol/m^3] to [umol/kg]
      HCO3L = CV2 * 1.D+6 * HCO3L             !! [mol/m^3] to [umol/kg]
      HCO3N = CV2 * 1.D+6 * HCO3N             !! [mol/m^3] to [umol/kg]
      HCO3M = CV2 * 1.D+6 * HCO3M             !! [mol/m^3] to [umol/kg]
      HCO3D = CV2 * 1.D+6 * HCO3D             !! [mol/m^3] to [umol/kg]
      CO32P = CV2 * 1.D+6 * CO32P             !! [mol/m^3] to [umol/kg]
      CO32S = CV2 * 1.D+6 * CO32S             !! [mol/m^3] to [umol/kg]
      CO32L = CV2 * 1.D+6 * CO32L             !! [mol/m^3] to [umol/kg]
      CO32N = CV2 * 1.D+6 * CO32N             !! [mol/m^3] to [umol/kg]
      CO32M = CV2 * 1.D+6 * CO32M             !! [mol/m^3] to [umol/kg]
      CO32D = CV2 * 1.D+6 * CO32D             !! [mol/m^3] to [umol/kg]
      PCO2P = CV3 * PCO2P                     !! [atm]     to [ppmv]
      PCO2S = CV3 * PCO2S                     !! [atm]     to [ppmv]
      PCO2L = CV3 * PCO2L                     !! [atm]     to [ppmv]
      PCO2N = CV3 * PCO2N                     !! [atm]     to [ppmv]
      PCO2M = CV3 * PCO2M                     !! [atm]     to [ppmv]
      PCO2D = CV3 * PCO2D                     !! [atm]     to [ppmv]
      PCO2A = CV3 * PCO2A                     !! [atm]     to [ppmv]
      TOCINI = TOCINI * 12.D-15               !! [molC]    to [PgC]
      TOCFIN = TOCFIN * 12.D-15               !! [molC]    to [PgC]

!-- Output

      WRITE(*,*) 'Temperature (P)          :', TEMP, 'K'
      WRITE(*,*) 'Temperature (S)          :', TEMS, 'K'
      WRITE(*,*) 'Temperature (L)          :', TEML, 'K'
      WRITE(*,*) 'Temperature (N)          :', TEMN, 'K'
      WRITE(*,*) 'Temperature (M)          :', TEMM, 'K'
      WRITE(*,*) 'Temperature (D)          :', TEMD, 'K'
      WRITE(*,*) ''
      WRITE(*,*) 'Salinity (P)             :', SALP, 'PSU'
      WRITE(*,*) 'Salinity (S)             :', SALS, 'PSU'
      WRITE(*,*) 'Salinity (L)             :', SALL, 'PSU'
      WRITE(*,*) 'Salinity (N)             :', SALN, 'PSU'
      WRITE(*,*) 'Salinity (M)             :', SALM, 'PSU'
      WRITE(*,*) 'Salinity (D)             :', SALD, 'PSU'
      WRITE(*,*) ''
      WRITE(*,*) 'EP (P)                   :', EPP, 'PgC/yr'
      WRITE(*,*) 'EP (S)                   :', EPS, 'PgC/yr'
      WRITE(*,*) 'EP (L)                   :', EPL, 'PgC/yr'
      WRITE(*,*) 'EP (N)                   :', EPN, 'PgC/yr'
      WRITE(*,*) ''
      WRITE(*,*) 'PO4 (P)                  :', PO4P, 'umol/kg'
      WRITE(*,*) 'PO4 (S)                  :', PO4S, 'umol/kg'
      WRITE(*,*) 'PO4 (L)                  :', PO4L, 'umol/kg'
      WRITE(*,*) 'PO4 (N)                  :', PO4N, 'umol/kg'
      WRITE(*,*) 'PO4 (M)                  :', PO4M, 'umol/kg'
      WRITE(*,*) 'PO4 (D)                  :', PO4D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DIC (P)                  :', DICP, 'umol/kg'
      WRITE(*,*) 'DIC (S)                  :', DICS, 'umol/kg'
      WRITE(*,*) 'DIC (L)                  :', DICL, 'umol/kg'
      WRITE(*,*) 'DIC (N)                  :', DICN, 'umol/kg'
      WRITE(*,*) 'DIC (M)                  :', DICM, 'umol/kg'
      WRITE(*,*) 'DIC (D)                  :', DICD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'ALK (P)                  :', ALKP, 'ueq/kg'
      WRITE(*,*) 'ALK (S)                  :', ALKS, 'ueq/kg'
      WRITE(*,*) 'ALK (L)                  :', ALKL, 'ueq/kg'
      WRITE(*,*) 'ALK (N)                  :', ALKN, 'ueq/kg'
      WRITE(*,*) 'ALK (M)                  :', ALKM, 'ueq/kg'
      WRITE(*,*) 'ALK (D)                  :', ALKD, 'ueq/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DO2 (P)                  :', DO2P, 'umol/kg'
      WRITE(*,*) 'DO2 (S)                  :', DO2S, 'umol/kg'
      WRITE(*,*) 'DO2 (L)                  :', DO2L, 'umol/kg'
      WRITE(*,*) 'DO2 (N)                  :', DO2N, 'umol/kg'
      WRITE(*,*) 'DO2 (M)                  :', DO2M, 'umol/kg'
      WRITE(*,*) 'DO2 (D)                  :', DO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'AOU (P)                  :', AOUP, 'umol/kg'
      WRITE(*,*) 'AOU (S)                  :', AOUS, 'umol/kg'
      WRITE(*,*) 'AOU (L)                  :', AOUL, 'umol/kg'
      WRITE(*,*) 'AOU (N)                  :', AOUN, 'umol/kg'
      WRITE(*,*) 'AOU (M)                  :', AOUM, 'umol/kg'
      WRITE(*,*) 'AOU (D)                  :', AOUD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO2 (P)                  :', CO2P, 'umol/kg'
      WRITE(*,*) 'CO2 (S)                  :', CO2S, 'umol/kg'
      WRITE(*,*) 'CO2 (L)                  :', CO2L, 'umol/kg'
      WRITE(*,*) 'CO2 (N)                  :', CO2N, 'umol/kg'
      WRITE(*,*) 'CO2 (M)                  :', CO2M, 'umol/kg'
      WRITE(*,*) 'CO2 (D)                  :', CO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'HCO3 (P)                 :', HCO3P, 'umol/kg'
      WRITE(*,*) 'HCO3 (S)                 :', HCO3S, 'umol/kg'
      WRITE(*,*) 'HCO3 (L)                 :', HCO3L, 'umol/kg'
      WRITE(*,*) 'HCO3 (N)                 :', HCO3N, 'umol/kg'
      WRITE(*,*) 'HCO3 (M)                 :', HCO3M, 'umol/kg'
      WRITE(*,*) 'HCO3 (D)                 :', HCO3D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO32 (P)                 :', CO32P, 'umol/kg'
      WRITE(*,*) 'CO32 (S)                 :', CO32S, 'umol/kg'
      WRITE(*,*) 'CO32 (L)                 :', CO32L, 'umol/kg'
      WRITE(*,*) 'CO32 (N)                 :', CO32N, 'umol/kg'
      WRITE(*,*) 'CO32 (M)                 :', CO32M, 'umol/kg'
      WRITE(*,*) 'CO32 (D)                 :', CO32D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'PCO2 (P)                 :', PCO2P, 'ppmv'
      WRITE(*,*) 'PCO2 (S)                 :', PCO2S, 'ppmv'
      WRITE(*,*) 'PCO2 (L)                 :', PCO2L, 'ppmv'
      WRITE(*,*) 'PCO2 (N)                 :', PCO2N, 'ppmv'
      WRITE(*,*) 'PCO2 (M)                 :', PCO2M, 'ppmv'
      WRITE(*,*) 'PCO2 (D)                 :', PCO2D, 'ppmv'
      WRITE(*,*) 'PCO2 (ATM)               :', PCO2A, 'ppmv'
      WRITE(*,*) ''
      WRITE(*,*) 'Total Carbon (Initial)   :', TOCINI, 'PgC'
      WRITE(*,*) 'Total Carbon (Final)     :', TOCFIN, 'PgC'
      WRITE(*,*) ''

      END PROGRAM
