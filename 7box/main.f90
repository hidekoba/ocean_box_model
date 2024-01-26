      PROGRAM MAIN

      USE CHEMEQ
      USE O2
      USE CO2

      IMPLICIT NONE

! --- information -----------------------------------------------------
!
!  Ocean 7-box model
!  Calculation of phosphate, DIC, alkalinity, CO2, and PCO2
!
!  References: Toggweiler et al. (1999, PA)
!
! ---------------------------------------------------------------------

      INTEGER :: I
      INTEGER, PARAMETER :: NBOX=7
!     1:P, 2:S, 3:L, 4:N, 5:M, 6:A, 7:D

      REAL(8) :: DT, DELTA
      REAL(8) :: CV1, CV2, CV3, CV4, CV5

      REAL(8) :: VOCN, AOCN, VATM
      REAL(8) :: ZOCNP, DEP, AOCNP, VOCNP
      REAL(8) :: ZOCNS, DES, AOCNS, VOCNS
      REAL(8) :: ZOCNL, DEL, AOCNL, VOCNL
      REAL(8) :: ZOCNN, DEN, AOCNN, VOCNN
      REAL(8) :: ZOCNM,      AOCNM, VOCNM
      REAL(8) :: ZOCNA,      AOCNA, VOCNA
      REAL(8) :: ZOCND,             VOCND

      REAL(8) :: T, FPD, FSM, FLM, FNA, FAD, FPS, FSL, FLN
      REAL(8) :: PVP, PVS, PVL, PVN
      REAL(8) :: FAP, FAS, FAL, FAN
      REAL(8) :: CEPP, CEPS, CEPN
      REAL(8) :: R, LF, HSC
      REAL(8) :: RRC, RCP, RNP, RO2P, GMM, GMA, RGMM, RGMA
      REAL(8) :: BT(NBOX)
      REAL(8) :: K0(NBOX), K1(NBOX), K2(NBOX)
      REAL(8) :: KB(NBOX), KW(NBOX), FH(NBOX)


      REAL(8) :: TEMP, SALP,   TEMPX, SALPX
      REAL(8) :: PO4P,  ALKP,   DICP,   DO2P
      REAL(8) :: PO4PX, ALKPX,  DICPX,  DO2PX
      REAL(8) :: CO2P,  HCO3P,  CO32P,  PCO2P, O2SATP, AOUP, EPP
      REAL(8) :: CO2PX, HCO3PX, CO32PX, PCO2PX

      REAL(8) :: TEMS, SALS,   TEMSX, SALSX
      REAL(8) :: PO4S,  ALKS,   DICS,   DO2S
      REAL(8) :: PO4SX, ALKSX,  DICSX,  DO2SX
      REAL(8) :: CO2S,  HCO3S,  CO32S,  PCO2S, O2SATS, AOUS, EPS
      REAL(8) :: CO2SX, HCO3SX, CO32SX, PCO2SX

      REAL(8) :: TEML, SALL,   TEMLX, SALLX
      REAL(8) :: PO4L,  ALKL,   DICL,   DO2L
      REAL(8) :: PO4LX, ALKLX,  DICLX,  DO2LX
      REAL(8) :: CO2L,  HCO3L,  CO32L,  PCO2L, O2SATL, AOUL, EPL
      REAL(8) :: CO2LX, HCO3LX, CO32LX, PCO2LX

      REAL(8) :: TEMN, SALN,   TEMNX, SALNX
      REAL(8) :: PO4N,  ALKN,   DICN,   DO2N
      REAL(8) :: PO4NX, ALKNX,  DICNX,  DO2NX
      REAL(8) :: CO2N,  HCO3N,  CO32N,  PCO2N, O2SATN, AOUN, EPN
      REAL(8) :: CO2NX, HCO3NX, CO32NX, PCO2NX

      REAL(8) :: TEMM, SALM,   TEMMX, SALMX
      REAL(8) :: PO4M,  ALKM,   DICM,   DO2M
      REAL(8) :: PO4MX, ALKMX,  DICMX,  DO2MX
      REAL(8) :: CO2M,  HCO3M,  CO32M,  PCO2M, O2SATM, AOUM, EPM
      REAL(8) :: CO2MX, HCO3MX, CO32MX, PCO2MX

      REAL(8) :: TEMA, SALA,   TEMAX, SALAX
      REAL(8) :: PO4A,  ALKA,   DICA,   DO2D
      REAL(8) :: PO4AX, ALKAX,  DICAX,  DO2AX
      REAL(8) :: CO2A,  HCO3A,  CO32A,  PCO2A, O2SATA, AOUA, EPA
      REAL(8) :: CO2AX, HCO3AX, CO32AX, PCO2AX

      REAL(8) :: TEMD, SALD,   TEMDX, SALDX
      REAL(8) :: PO4D,  ALKD,   DICD,   DO2D
      REAL(8) :: PO4DX, ALKDX,  DICDX,  DO2DX
      REAL(8) :: CO2D,  HCO3D,  CO32D,  PCO2D, O2SATD, AOUD, EPD
      REAL(8) :: CO2DX, HCO3DX, CO32DX, PCO2DX

      REAL(8) :: PCO2Z, PCO2ZX
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
      DATA TEMP, TEMS, TEML, TEMN, TEMM, TEMA, TEMD &
     & / 1.0, 8.0, 21.5, 3.0, 10.0, 3.0, 1.75 /
      DATA SALP, SALS, SALL, SALN, SALM, SALA, SALD &
     & / 34.7, 34.7, 34.7, 34.7, 34.7, 34.7, 34.7 /

!--   Ocean volume [m^3]
!--   Ocean area [m^2]
!--   Atmospheric molers [mol] = [mol/atm]
      DATA VOCN, AOCN, VATM &
     & / 1.292D18, 3.49D14, 1.773D20 /

!--   Ocean depth [m]
      DATA ZP, ZS, ZL, ZN, ZM, ZA, ZD &
     & / 2.5D2, 2.5D2, 1.0D2, 2.5D2, 9.0D2, 2.0D3, 4.0D3 /

!--   Ocean depth of euphotic zone [m]
      DATA DEP, DES, DEL, DEN &
     & / 2.5D2, 2.5D2, 1.0D2, 2.5D2 /

!--   Transport [m^3/s]
!--   Diffusivity [m^3/s]
      DATA T &
     & / 2.0D7 /
      DATA FPD, FNA, FSM, FLM, FAD &
     & / 6.0D7, 0.0D0, 0.0D0, 4.0D7, 0.0D0 /
      DATA FPS, FSL, FLN &
     & / 0.0D0, 0.0D0, 0.0D0 /

!--   Piston velocity [m/day]
      DATA PVP, PVS, PVL, PVN &
     & / 3.0D0, 3.0D0, 3.0D0, 3.0D0 /

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

!--   Remineralization fraction
      GMM  = 0.80D0
      RGMM = 1.0D0 - GMM
      GMA  = 0.50D0
      RGMA = 1.0D0 - GMA

!--   Biological sinking flux [molC/m^2/yr]
      DATA CEPP, CEPS, CEPN &
     & / 3.0D0, 6.5D-1, 1.0D0 /

!--   Initial atomospheric PCO2 [ppmv]
      DATA PCO2Z / 280.D0 /

!-- Initial condition (etc)

!--   Ocean area [m^2]
!--   P, S, L, N, M, A
!--   5% , 10% , 75% , 10% , 85% ,95%
      AOCNP = AOCN * 5.0D-2
      AOCNS = AOCN * 1.0D-1
      AOCNL = AOCN * 7.5D-1
      AOCNN = AOCN * 1.0D-1
      AOCNM = AOCN * 8.5D-1
      AOCNA = AOCN * 9.5D-1

!--   Ocean volume [m^3]
      VOCNP = AOCNP * ZOCNP
      VOCNS = AOCNS * ZOCNS
      VOCNL = AOCNL * ZOCNL
      VOCNN = AOCNN * ZOCNN
      VOCNM = AOCNM * ZOCNM - VOCNS - VOCNL
      VOCNA = AOCNA * ZOCNA - VOCNS - VOCNL - VOCNN - VOCNM
      VOCND = VOCN - VOCNP - VOCNS - VOCNL - VOCNN - VOCNM - VOCNA

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
      DATA PO4P, PO4S, PO4L, PO4N, PO4M, PO4A, PO4D &
     & / 2.00D-6,2.33D-6,0.0D-6,1.0D-6,2.142D-6,2.142D-6,2.142D-6 /

!--   Alkalinity [eq/kg]
      DATA ALKP, ALKS, ALKL, ALKN, ALKM, ALKA, ALKD &
     & / 2.360D-3,2.325D-3,2.277D-3,2.324D-3,2.3738D-3,2.3738D-3,2.3738D-3 /

!--   Dissolved inorganic carbon [mol/kg]
      DATA DICP, DICS, DICL, DICN, DICM, DICA, DICD &
     & / 2.150D-3,2.107D-3,1.933D-3,2.106D-3,2.2561D-3,2.2561D-3,2.2561D-3 /

!--   Dissolved oxygen [mol/kg]
      DATA DO2P, DO2S, DO2N, DO2M, DO2A, DO2D &
     & / 3.50D-4,3.00D-4,3.00D-4,1.50D-4,2.50D-4,2.50D-4 /

!-- Unit conversion

      EPP   = EPP   / CV4 !! [/yr]    to [/s]
      EPS   = EPS   / CV4 !! [/yr]    to [/s]
      EPN   = EPN   / CV4 !! [/yr]    to [/s]
      FAP   = FAP   / CV5 !! [/day]   to [/s]
      FAS   = FAS   / CV5 !! [/day]   to [/s]
      FAL   = FAL   / CV5 !! [/day]   to [/s]
      FAN   = FAN   / CV5 !! [/day]   to [/s]
      PO4P  = PO4P  * CV1 !! [mol/kg] to [mol/m^3]
      PO4S  = PO4S  * CV1 !! [mol/kg] to [mol/m^3]
      PO4L  = PO4L  * CV1 !! [mol/kg] to [mol/m^3]
      PO4N  = PO4N  * CV1 !! [mol/kg] to [mol/m^3]
      PO4M  = PO4M  * CV1 !! [mol/kg] to [mol/m^3]
      PO4A  = PO4A  * CV1 !! [mol/kg] to [mol/m^3]
      PO4D  = PO4D  * CV1 !! [mol/kg] to [mol/m^3]
      DICP  = DICP  * CV1 !! [mol/kg] to [mol/m^3]
      DICS  = DICS  * CV1 !! [mol/kg] to [mol/m^3]
      DICL  = DICL  * CV1 !! [mol/kg] to [mol/m^3]
      DICN  = DICN  * CV1 !! [mol/kg] to [mol/m^3]
      DICM  = DICM  * CV1 !! [mol/kg] to [mol/m^3]
      DICA  = DICA  * CV1 !! [mol/kg] to [mol/m^3]
      DICD  = DICD  * CV1 !! [mol/kg] to [mol/m^3]
      ALKP  = ALKP  * CV1 !! [mol/kg] to [mol/m^3]
      ALKS  = ALKS  * CV1 !! [mol/kg] to [mol/m^3]
      ALKL  = ALKL  * CV1 !! [mol/kg] to [mol/m^3]
      ALKN  = ALKN  * CV1 !! [mol/kg] to [mol/m^3]
      ALKM  = ALKM  * CV1 !! [mol/kg] to [mol/m^3]
      ALKA  = ALKA  * CV1 !! [mol/kg] to [mol/m^3]
      ALKD  = ALKD  * CV1 !! [mol/kg] to [mol/m^3]
      DO2P  = DO2P  * CV1 !! [mol/kg] to [mol/m^3]
      DO2S  = DO2S  * CV1 !! [mol/kg] to [mol/m^3]
      DO2N  = DO2N  * CV1 !! [mol/kg] to [mol/m^3]
      DO2M  = DO2M  * CV1 !! [mol/kg] to [mol/m^3]
      DO2A  = DO2A  * CV1 !! [mol/kg] to [mol/m^3]
      DO2D  = DO2D  * CV1 !! [mol/kg] to [mol/m^3]
      PCO2Z = PCO2Z / CV3 !! [ppmv]   to [atm]

!-- Loop

      DO I = 1, 100000

!-- Calculation of equilibrium constant
!   1:P, 2:S, 3:L, 4:N, 5:M, 6:A, 7:D

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
     &                     TEMA, SALA,                                 &
     &                     BT(6), K0(6), K1(6), K2(6),                 &
     &                     KB(6), KW(6), FH(6)                         &
     &                   )
         CALL CHEMEQCONST(                                             &
     &                     TEMD, SALD,                                 &
     &                     BT(7), K0(7), K1(7), K2(7),                 &
     &                     KB(7), KW(7), FH(7)                         &
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
     &                      ALKA, DICA, CO2A, HCO3A, CO32A, PCO2A      &
     &                    )
            CALL CO2_NIBUN(                                            &
     &                      BT(7), K0(7), K1(7), K2(7), KB(7), KW(7),  &
     &                      ALKD, DICD, CO2D, HCO3D, CO32D, PCO2D      &
     &                    )
            TOCINI = + DICP  * VOCNP                                   &
     &               + DICS  * VOCNS                                   &
     &               + DICL  * VOCNL                                   &
     &               + DICN  * VOCNN                                   &
     &               + DICM  * VOCNM                                   &
     &               + DICA  * VOCNA                                   &
     &               + DICD  * VOCND                                   &
     &               + PCO2Z * VATM
         ENDIF

!-- Calculation of biogeochemical tracers
!2345&78901234567890123456789012345678901234567890123456789012345678901&

         EPL = FSL * (PO4S - PO4L) + FLN * (PO4N - PO4L) + FLM * (PO4M - PO4L)
!--      EP_H = R * DE_H * LF * PO_H * (PO_H / (HSC + PO_H))
!--      EP_L = R * DE_L * LF * PO_L * (PO_L / (HSC + PO_L))

         !-- 温度
         TEMMX = TEMM + (T * (TEMS - TEMM) + FLM * (TEML - TEMM)) * (DT / VOCNM)
         TEMAX = TEMA + (T * (TEMN - TEMA) + FAD * (TEMD - TEMA)) * (DT / VOCNA)
         TEMDX = TEMD + ((T + FAD) * (TEMN - TEMD) + FPD * (TEMP - TEMD)) * (DT / VOCND)

         !-- 栄養塩
         PO4PX = PO4P + ((T + FPD) * (PO4D - PO4P) + FPS * (PO4S - PO4P) &
         & - EPP) * (DT / VOCNP)
         PO4SX = PO4S + ((T + FPS) * (PO4P - PO4S) + FSL * (PO4L - PO4S) &
         & - EPS) * (DT / VOCNS)
         PO4LX = PO4L + (FSL * (PO4S - PO4L) + FLN * (PO4N - PO4L) + FLM * (PO4M - PO4L) &
         & - EPL) * (DT / VOCNL)
         PO4NX = PO4N + (T * (PO4M - PO4N) + FLN * (PO4L - PO4N) + FNA * (PO4A - PO4N) &
         & - EPN) * (DT / VOCNN)
         PO4MX = PO4M + ((T + FSM) * (PO4S - PO4M) + FLM * (PO4L - PO4M) &
         & + (GMM * (EPS + EPL))) * (DT / VOCNM)
         PO4AX = PO4A + ((T + FNA) * (PO4N - PO4A) + FAD * (PO4D - PO4A) &
         & + (RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
         PO4DX = PO4D + ((T + FAD) * (PO4A - PO4D) + FPD * (PO4P - PO4D) &
         & + (RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

         !-- アルカリ度
         ALKPX = ALKP + ((T + FPD) * (ALKD - ALKP) + FPS * (ALKS - ALKP) &
         & - (2.0D0 * RRC * RCP - RNP) * EPP) * (DT / VOCNP)
         ALKSX = ALKS + ((T + FPS) * (ALKP - ALKS) + FSL * (ALKL - ALKS) &
         & - (2.0D0 * RRC * RCP - RNP) * EPS) * (DT / VOCNS)
         ALKLX = ALKL + (FSL * (ALKS - ALKL) + FLN * (ALKN - ALKL) + FLM * (ALKM - ALKL) &
         & - (2.0D0 * RRC * RCP - RNP) * EPL) * (DT / VOCNL)
         ALKNX = ALKN + (T * (ALKM - ALKN) + FLN * (ALKL - ALKN) + FNA * (ALKA - ALKN) &
         & - (2.0D0 * RRC * RCP - RNP) * EPN) * (DT / VOCNN)
         ALKMX = ALKM + ((T + FSM) * (ALKS - ALKM) + FLM * (ALKL - ALKM) &
         & + (2.0D0 * RRC * RCP - RNP) * (GMM * (EPS + EPL))) * (DT / VOCNM)
         ALKAX = ALKA + ((T + FNA) * (ALKN - ALKA) + FAD * (ALKD - ALKA) &
         & + (2.0D0 * RRC * RCP - RNP) * (RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
         ALKDX = ALKD + ((T + FAD) * (ALKA - ALKD) + FPD * (ALKP - ALKD) &
         & + (2.0D0 * RRC * RCP - RNP) * (RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

         !-- 全炭酸
         DICPX = DICP + ((T + FPD) * (DICD - DICP) + FPS * (DICS - DICP) &
         & - (1.0D0 + RRC) * RCP * EPP) * (DT / VOCNP)
         DICSX = DICS + ((T + FPS) * (DICP - DICS) + FSL * (DICL - DICS) &
         & - (1.0D0 + RRC) * RCP * EPS) * (DT / VOCNS)
         DICLX = DICL + (FSL * (DICS - DICL) + FLN * (DICN - DICL) + FLM * (DICM - DICL) &
         & - (1.0D0 + RRC) * RCP * EPL) * (DT / VOCNL)
         DICNX = DICN + (T * (DICM - DICN) + FLN * (DICL - DICN) + FNA * (DICA - DICN) &
         & - (1.0D0 + RRC) * RCP * EPN) * (DT / VOCNN)
         DICMX = DICM + ((T + FSM) * (DICS - DICM) + FLM * (DICL - DICM) &
         & + (1.0D0 + RRC) * RCP * (GMM * (EPS + EPL))) * (DT / VOCNM)
         DICAX = DICA + ((T + FNA) * (DICN - DICA) + FAD * (DICD - DICA) &
         & + (1.0D0 + RRC) * RCP * (RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
         DICDX = DICD + ((T + FAD) * (DICA - DICD) + FPD * (DICP - DICD) &
         & + (1.0D0 + RRC) * RCP * (RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

         !-- 酸素
         DO2PX = DO2P + ((T + FPD) * (DO2D - DO2P) + FPS * (DO2S - DO2P) &
         & - RDO2P * EPP) * (DT / VOCNP)
         DO2SX = DO2S + ((T + FPS) * (DO2P - DO2S) + FSL * (DO2L - DO2S) &
         & - RDO2P * EPS) * (DT / VOCNS)
         DO2NX = DO2N + (T * (DO2M - DO2N) + FLN * (DO2L - DO2N) + FNA * (DO2A - DO2N) &
         & - RDO2P * EPN) * (DT / VOCNN)
         DO2MX = DO2M + ((T + FSM) * (DO2S - DO2M) + FLM * (DO2L - DO2M) &
         & + RDO2P * (GMM * (EPS + EPL))) * (DT / VOCNM)
         DO2AX = DO2A + ((T + FNA) * (DO2N - DO2A) + FAD * (DO2D - DO2A) &
         & + RDO2P * (RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
         DO2DX = DO2D + ((T + FAD) * (DO2A - DO2D) + FPD * (DO2P - DO2D) &
         & + RDO2P * (RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

!-- Calculation of carbonate system

         CALL CO2_NIBUN(                                               &
     &                   BT(1), K0(1), K1(1), K2(1), KB(1), KW(1),     &
     &                   ALKPX, DICPX, CO2PX, HCO3P, CO32P, PCO2P      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(2), K0(2), K1(2), K2(2), KB(2), KW(2),     &
     &                   ALKSX, DICSX, CO2SX, HCO3S, CO32S, PCO2S      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(3), K0(3), K1(3), K2(3), KB(3), KW(3),     &
     &                   ALKLX, DICLX, CO2LX, HCO3L, CO32L, PCO2L      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(4), K0(4), K1(4), K2(4), KB(4), KW(4),     &
     &                   ALKNX, DICNX, CO2NX, HCO3N, CO32N, PCO2N      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(5), K0(5), K1(5), K2(5), KB(5), KW(5),     &
     &                   ALKMX, DICMX, CO2MX, HCO3M, CO32M, PCO2M      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(6), K0(6), K1(6), K2(6), KB(6), KW(6),     &
     &                   ALKAX, DICAX, CO2AX, HCO3A, CO32A, PCO2A      &
     &                 )
         CALL CO2_NIBUN(                                               &
     &                   BT(7), K0(7), K1(7), K2(7), KB(7), KW(7),     &
     &                   ALKDX, DICDX, CO2DX, HCO3D, CO32D, PCO2D      &
     &                 )

!-- Calculation of atmospheric PCO2

!2345&78901234567890123456789012345678901234567890123456789012345678901&
         PCO2AZ = PCO2Z                                                &
     &          + (                                                    &
     &              + K0(1) * FAP * (PCO2P - PCO2Z)                    &
     &              + K0(2) * FAS * (PCO2S - PCO2Z)                    &
     &              + K0(3) * FAL * (PCO2L - PCO2Z)                    &
     &              + K0(4) * FAN * (PCO2N - PCO2Z)                    &
     &            )                                                    &
     &            * (DT / VOCNATM)

!-- Equilibrium condition

         IF (ABS((PO4PX / PO4P) - 1.0D0) <= DELTA) THEN
           WRITE(*,*) 'TIMESTEP =', I
           EXIT
         END IF

!-- Next loop

!        TEMM  = TEMMX
!        TEMA  = TEMAX
!        TEMD  = TEMDX
         PO4P  = PO4PX
         PO4S  = PO4SX
         PO4L  = PO4LX
         PO4N  = PO4NX
         PO4M  = PO4MX
         PO4A  = PO4AX
         PO4D  = PO4DX
         ALKP  = ALKPX
         ALKS  = ALKSX
         ALKL  = ALKLX
         ALKN  = ALKNX
         ALKM  = ALKMX
         ALKA  = ALKAX
         ALKD  = ALKDX
         DICP  = DICPX
         DICS  = DICSX
         DICL  = DICLX
         DICN  = DICNX
         DICM  = DICMX
         DICA  = DICAX
         DICD  = DICDX
         DO2P  = DO2PX
         DO2S  = DO2SX
         DO2N  = DO2NX
         DO2M  = DO2MX
         DO2A  = DO2AX
         DO2D  = DO2DX
         PCO2Z = PCO2ZX

      ENDDO

!-- Final total carbon

      TOCFIN = + DICP  * VOCNP                                         &
     &         + DICS  * VOCNS                                         &
     &         + DICL  * VOCNL                                         &
     &         + DICN  * VOCNN                                         &
     &         + DICM  * VOCNM                                         &
     &         + DICA  * VOCNA                                         &
     &         + DICD  * VOCND                                         &
     &         + PCO2Z * VATM

!-- Oxygen utilization

      CALL O2SAT(TEMP, SALP, O2SATP)
      CALL O2SAT(TEMS, SALS, O2SATS)
      CALL O2SAT(TEML, SALL, O2SATL)
      CALL O2SAT(TEMN, SALN, O2SATN)
      CALL O2SAT(TEMM, SALM, O2SATM)
      CALL O2SAT(TEMA, SALA, O2SATA)
      CALL O2SAT(TEMD, SALD, O2SATD)

      AOUP = (O2SATP - DO2P)
      AOUS = (O2SATS - DO2S)
      AOUL = (O2SATL - DO2L)
      AOUN = (O2SATN - DO2N)
      AOUM = (O2SATM - DO2M)
      AOUA = (O2SATA - DO2A)
      AOUD = (O2SATD - DO2D)

!-- Unit conversion

      EPP   = CV4 * EPP * RCP * 1.D-15 * 12.0 !! [molP/s] to [PgC/yr]
      EPS   = CV4 * EPS * RCP * 1.D-15 * 12.0 !! [molP/s] to [PgC/yr]
      EPN   = CV4 * EPN * RCP * 1.D-15 * 12.0 !! [molP/s] to [PgC/yr]
      PO4P  = CV2 * 1.D+6 * PO4P              !! [mol/m^3] to [umol/kg]
      PO4S  = CV2 * 1.D+6 * PO4S              !! [mol/m^3] to [umol/kg]
      PO4L  = CV2 * 1.D+6 * PO4L              !! [mol/m^3] to [umol/kg]
      PO4N  = CV2 * 1.D+6 * PO4N              !! [mol/m^3] to [umol/kg]
      PO4M  = CV2 * 1.D+6 * PO4M              !! [mol/m^3] to [umol/kg]
      PO4A  = CV2 * 1.D+6 * PO4A              !! [mol/m^3] to [umol/kg]
      PO4D  = CV2 * 1.D+6 * PO4D              !! [mol/m^3] to [umol/kg]
      DICP  = CV2 * 1.D+6 * DICP              !! [mol/m^3] to [umol/kg]
      DICS  = CV2 * 1.D+6 * DICS              !! [mol/m^3] to [umol/kg]
      DICL  = CV2 * 1.D+6 * DICL              !! [mol/m^3] to [umol/kg]
      DICN  = CV2 * 1.D+6 * DICN              !! [mol/m^3] to [umol/kg]
      DICM  = CV2 * 1.D+6 * DICM              !! [mol/m^3] to [umol/kg]
      DICA  = CV2 * 1.D+6 * DICA              !! [mol/m^3] to [umol/kg]
      DICD  = CV2 * 1.D+6 * DICD              !! [mol/m^3] to [umol/kg]
      ALKP  = CV2 * 1.D+6 * ALKP              !! [mol/m^3] to [umol/kg]
      ALKS  = CV2 * 1.D+6 * ALKS              !! [mol/m^3] to [umol/kg]
      ALKL  = CV2 * 1.D+6 * ALKL              !! [mol/m^3] to [umol/kg]
      ALKN  = CV2 * 1.D+6 * ALKN              !! [mol/m^3] to [umol/kg]
      ALKM  = CV2 * 1.D+6 * ALKM              !! [mol/m^3] to [umol/kg]
      ALKA  = CV2 * 1.D+6 * ALKA              !! [mol/m^3] to [umol/kg]
      ALKD  = CV2 * 1.D+6 * ALKD              !! [mol/m^3] to [umol/kg]
      DO2P  = CV2 * 1.D+6 * DO2P              !! [mol/m^3] to [umol/kg]
      DO2S  = CV2 * 1.D+6 * DO2S              !! [mol/m^3] to [umol/kg]
      DO2L  = CV2 * 1.D+6 * DO2L              !! [mol/m^3] to [umol/kg]
      DO2N  = CV2 * 1.D+6 * DO2N              !! [mol/m^3] to [umol/kg]
      DO2M  = CV2 * 1.D+6 * DO2M              !! [mol/m^3] to [umol/kg]
      DO2A  = CV2 * 1.D+6 * DO2A              !! [mol/m^3] to [umol/kg]
      DO2D  = CV2 * 1.D+6 * DO2D              !! [mol/m^3] to [umol/kg]
      AOUP  = CV2 * 1.D+6 * AOUP              !! [mol/m^3] to [umol/kg]
      AOUS  = CV2 * 1.D+6 * AOUS              !! [mol/m^3] to [umol/kg]
      AOUL  = CV2 * 1.D+6 * AOUL              !! [mol/m^3] to [umol/kg]
      AOUN  = CV2 * 1.D+6 * AOUN              !! [mol/m^3] to [umol/kg]
      AOUM  = CV2 * 1.D+6 * AOUM              !! [mol/m^3] to [umol/kg]
      AOUA  = CV2 * 1.D+6 * AOUA              !! [mol/m^3] to [umol/kg]
      AOUD  = CV2 * 1.D+6 * AOUD              !! [mol/m^3] to [umol/kg]
      CO2P  = CV2 * 1.D+6 * CO2P              !! [mol/m^3] to [umol/kg]
      CO2S  = CV2 * 1.D+6 * CO2S              !! [mol/m^3] to [umol/kg]
      CO2L  = CV2 * 1.D+6 * CO2L              !! [mol/m^3] to [umol/kg]
      CO2N  = CV2 * 1.D+6 * CO2N              !! [mol/m^3] to [umol/kg]
      CO2M  = CV2 * 1.D+6 * CO2M              !! [mol/m^3] to [umol/kg]
      CO2A  = CV2 * 1.D+6 * CO2A              !! [mol/m^3] to [umol/kg]
      CO2D  = CV2 * 1.D+6 * CO2D              !! [mol/m^3] to [umol/kg]
      HCO3P = CV2 * 1.D+6 * HCO3P             !! [mol/m^3] to [umol/kg]
      HCO3S = CV2 * 1.D+6 * HCO3S             !! [mol/m^3] to [umol/kg]
      HCO3L = CV2 * 1.D+6 * HCO3L             !! [mol/m^3] to [umol/kg]
      HCO3N = CV2 * 1.D+6 * HCO3N             !! [mol/m^3] to [umol/kg]
      HCO3M = CV2 * 1.D+6 * HCO3M             !! [mol/m^3] to [umol/kg]
      HCO3A = CV2 * 1.D+6 * HCO3A             !! [mol/m^3] to [umol/kg]
      HCO3D = CV2 * 1.D+6 * HCO3D             !! [mol/m^3] to [umol/kg]
      CO32P = CV2 * 1.D+6 * CO32P             !! [mol/m^3] to [umol/kg]
      CO32S = CV2 * 1.D+6 * CO32S             !! [mol/m^3] to [umol/kg]
      CO32L = CV2 * 1.D+6 * CO32L             !! [mol/m^3] to [umol/kg]
      CO32N = CV2 * 1.D+6 * CO32N             !! [mol/m^3] to [umol/kg]
      CO32M = CV2 * 1.D+6 * CO32M             !! [mol/m^3] to [umol/kg]
      CO32A = CV2 * 1.D+6 * CO32A             !! [mol/m^3] to [umol/kg]
      CO32D = CV2 * 1.D+6 * CO32D             !! [mol/m^3] to [umol/kg]
      PCO2P = CV3 * PCO2P                     !! [atm]     to [ppmv]
      PCO2S = CV3 * PCO2S                     !! [atm]     to [ppmv]
      PCO2L = CV3 * PCO2L                     !! [atm]     to [ppmv]
      PCO2N = CV3 * PCO2N                     !! [atm]     to [ppmv]
      PCO2M = CV3 * PCO2M                     !! [atm]     to [ppmv]
      PCO2A = CV3 * PCO2A                     !! [atm]     to [ppmv]
      PCO2D = CV3 * PCO2D                     !! [atm]     to [ppmv]
      PCO2Z = CV3 * PCO2Z                     !! [atm]     to [ppmv]

!-- Output

      WRITE(*,*) 'Temperature (P)          :', TEMP, 'K'
      WRITE(*,*) 'Temperature (S)          :', TEMS, 'K'
      WRITE(*,*) 'Temperature (L)          :', TEML, 'K'
      WRITE(*,*) 'Temperature (N)          :', TEMN, 'K'
      WRITE(*,*) 'Temperature (M)          :', TEMM, 'K'
      WRITE(*,*) 'Temperature (A)          :', TEMA, 'K'
      WRITE(*,*) 'Temperature (D)          :', TEMD, 'K'
      WRITE(*,*) ''
      WRITE(*,*) 'Salinity (P)             :', SALP, 'PSU'
      WRITE(*,*) 'Salinity (S)             :', SALS, 'PSU'
      WRITE(*,*) 'Salinity (L)             :', SALL, 'PSU'
      WRITE(*,*) 'Salinity (N)             :', SALN, 'PSU'
      WRITE(*,*) 'Salinity (M)             :', SALM, 'PSU'
      WRITE(*,*) 'Salinity (A)             :', SALA, 'PSU'
      WRITE(*,*) 'Salinity (D)             :', SALD, 'PSU'
      WRITE(*,*) ''
      WRITE(*,*) 'EP (P)                   :', EPP, 'PgC/yr'
      WRITE(*,*) 'EP (S)                   :', EPS, 'PgC/yr'
      WRITE(*,*) 'EP (N)                   :', EPN, 'PgC/yr'
      WRITE(*,*) ''
      WRITE(*,*) 'PO4 (P)                  :', PO4P, 'umol/kg'
      WRITE(*,*) 'PO4 (S)                  :', PO4S, 'umol/kg'
      WRITE(*,*) 'PO4 (L)                  :', PO4L, 'umol/kg'
      WRITE(*,*) 'PO4 (N)                  :', PO4N, 'umol/kg'
      WRITE(*,*) 'PO4 (M)                  :', PO4M, 'umol/kg'
      WRITE(*,*) 'PO4 (A)                  :', PO4A, 'umol/kg'
      WRITE(*,*) 'PO4 (D)                  :', PO4D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DIC (P)                  :', DICP, 'umol/kg'
      WRITE(*,*) 'DIC (S)                  :', DICS, 'umol/kg'
      WRITE(*,*) 'DIC (L)                  :', DICL, 'umol/kg'
      WRITE(*,*) 'DIC (N)                  :', DICN, 'umol/kg'
      WRITE(*,*) 'DIC (M)                  :', DICM, 'umol/kg'
      WRITE(*,*) 'DIC (A)                  :', DICA, 'umol/kg'
      WRITE(*,*) 'DIC (D)                  :', DICD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'ALK (P)                  :', ALKP, 'ueq/kg'
      WRITE(*,*) 'ALK (S)                  :', ALKS, 'ueq/kg'
      WRITE(*,*) 'ALK (L)                  :', ALKL, 'ueq/kg'
      WRITE(*,*) 'ALK (N)                  :', ALKN, 'ueq/kg'
      WRITE(*,*) 'ALK (M)                  :', ALKM, 'ueq/kg'
      WRITE(*,*) 'ALK (A)                  :', ALKA, 'ueq/kg'
      WRITE(*,*) 'ALK (D)                  :', ALKD, 'ueq/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'DO2 (P)                  :', DO2P, 'umol/kg'
      WRITE(*,*) 'DO2 (S)                  :', DO2S, 'umol/kg'
      WRITE(*,*) 'DO2 (L)                  :', DO2L, 'umol/kg'
      WRITE(*,*) 'DO2 (N)                  :', DO2N, 'umol/kg'
      WRITE(*,*) 'DO2 (M)                  :', DO2M, 'umol/kg'
      WRITE(*,*) 'DO2 (A)                  :', DO2A, 'umol/kg'
      WRITE(*,*) 'DO2 (D)                  :', DO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'AOU (P)                  :', AOUP, 'umol/kg'
      WRITE(*,*) 'AOU (S)                  :', AOUS, 'umol/kg'
      WRITE(*,*) 'AOU (L)                  :', AOUL, 'umol/kg'
      WRITE(*,*) 'AOU (N)                  :', AOUN, 'umol/kg'
      WRITE(*,*) 'AOU (M)                  :', AOUM, 'umol/kg'
      WRITE(*,*) 'AOU (A)                  :', AOUA, 'umol/kg'
      WRITE(*,*) 'AOU (D)                  :', AOUD, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO2 (P)                  :', CO2P, 'umol/kg'
      WRITE(*,*) 'CO2 (S)                  :', CO2S, 'umol/kg'
      WRITE(*,*) 'CO2 (L)                  :', CO2L, 'umol/kg'
      WRITE(*,*) 'CO2 (N)                  :', CO2N, 'umol/kg'
      WRITE(*,*) 'CO2 (M)                  :', CO2M, 'umol/kg'
      WRITE(*,*) 'CO2 (A)                  :', CO2A, 'umol/kg'
      WRITE(*,*) 'CO2 (D)                  :', CO2D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'HCO3 (P)                 :', HCO3P, 'umol/kg'
      WRITE(*,*) 'HCO3 (S)                 :', HCO3S, 'umol/kg'
      WRITE(*,*) 'HCO3 (L)                 :', HCO3L, 'umol/kg'
      WRITE(*,*) 'HCO3 (N)                 :', HCO3N, 'umol/kg'
      WRITE(*,*) 'HCO3 (M)                 :', HCO3M, 'umol/kg'
      WRITE(*,*) 'HCO3 (A)                 :', HCO3A, 'umol/kg'
      WRITE(*,*) 'HCO3 (D)                 :', HCO3D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'CO32 (P)                 :', CO32P, 'umol/kg'
      WRITE(*,*) 'CO32 (S)                 :', CO32S, 'umol/kg'
      WRITE(*,*) 'CO32 (L)                 :', CO32L, 'umol/kg'
      WRITE(*,*) 'CO32 (N)                 :', CO32N, 'umol/kg'
      WRITE(*,*) 'CO32 (M)                 :', CO32M, 'umol/kg'
      WRITE(*,*) 'CO32 (A)                 :', CO32A, 'umol/kg'
      WRITE(*,*) 'CO32 (D)                 :', CO32D, 'umol/kg'
      WRITE(*,*) ''
      WRITE(*,*) 'PCO2 (P)                 :', PCO2P, 'ppmv'
      WRITE(*,*) 'PCO2 (S)                 :', PCO2S, 'ppmv'
      WRITE(*,*) 'PCO2 (L)                 :', PCO2L, 'ppmv'
      WRITE(*,*) 'PCO2 (N)                 :', PCO2N, 'ppmv'
      WRITE(*,*) 'PCO2 (M)                 :', PCO2M, 'ppmv'
      WRITE(*,*) 'PCO2 (A)                 :', PCO2A, 'ppmv'
      WRITE(*,*) 'PCO2 (D)                 :', PCO2D, 'ppmv'
      WRITE(*,*) 'PCO2 (ATM)               :', PCO2Z, 'ppmv'
      WRITE(*,*) ''
      WRITE(*,*) 'Total Carbon (Initial)   :', TOCINI, 'mol'
      WRITE(*,*) 'Total Carbon (Final)     :', TOCFIN, 'mol'
      WRITE(*,*) ''

      END PROGRAM
