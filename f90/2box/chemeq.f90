      MODULE CHEMEQ

      IMPLICIT NONE
      CONTAINS

      SUBROUTINE CHEMEQCONST( &
     &    TEM,                & !! (in)
     &    SAL,                & !! (in)
     &    BT,                 & !! (out)
     &    K0,                 & !! (out)
     &    K1,                 & !! (out)
     &    K2,                 & !! (out)
     &    KB,                 & !! (out)
     &    KW,                 & !! (out)
     &    FH                  & !! (out)
     &   )

! --- information -----------------------------------------------------
!
!  Calculation of equilibrium constant
!  from temperature and salinity.
!
!    BT : Total Boron [mol/kg]
!    K0 : CO2 solubility [mol/kg/atm]
!    K1 : Equilibrium constant [mol/kg]
!    K2 : Equilibrium constant [mol/kg]
!    KB : Equilibrium constant [mol/kg]
!    KW : Equilibrium constant [mol/kg]
!
! ---------------------------------------------------------------------


      REAL(8), INTENT(IN) :: TEM, SAL
      REAL(8) :: ATEM, TK, SK
      REAL(8), INTENT(OUT) :: BT, K0, K1, K2, KB, KW, FH


!-- Calculation of equilibrium constant

      ATEM = TEM + 273.15D0
      TK   = ATEM * 1.D-2
      SK   = 2.3517D-2 + ( - 2.3656D-2 + 4.7036D-3 * TK ) * TK

      K0 =                                                             &
     &  EXP(                                                           &
     &      - 6.02409D1 + 9.34517D1 / TK + 2.33585D1                   &
     &      * LOG(TK) + SK * SAL                                       &
     &     )

      K1 =                                                             &
     &  EXP(                                                           &
     &       LOG(10.D0) *                                              &
     &       (                                                         &
     &         + 1.37201D1 - 3.1334D-2 * ATEM - 3.23576D3 / ATEM       &
     &         - 1.3D-5 * SAL * ATEM + 1.032D-1 * SQRT(SAL)            &
     &       )                                                         &
     &     )

      K2 =                                                             &
     &  EXP(                                                           &
     &       LOG(10.D0) *                                              &
     &       (                                                         &
     &         - 5.3719645D3 - 1.671221D0 * ATEM - 2.2913D-1 * SAL     &
     &         - 1.83802D1 * LOG10(SAL) + 1.2837528D5 / ATEM           &
     &         + 2.1943005D3 * LOG10(ATEM) + 8.0944D-4 * SAL * ATEM    &
     &         + 5.61711D3 * LOG10(SAL) / ATEM-2.136D0 * SAL / ATEM    &
     &       )                                                         &
     &     )

      KB =                                                             &
     &  EXP(                                                           &
     &       LOG(10.D0) *                                              &
     &       (                                                         &
     &         - 9.26D0 + 8.86D-3 * SAL + 1.D-3 * TEM                  &
     &       )                                                         &
     &     )

      KW =                                                             &
     &  EXP(                                                           &
     &       + 1.489802D2                                              &
     &       - 1.384726D4 / ATEM                                       &
     &       - 2.36521D1 * LOG(ATEM)                                   &
     &       + (                                                       &
     &           - 7.92447D1 + 3.29872D3 / ATEM                        &
     &           + 1.20408D1 * LOG(ATEM)                               &
     &         )                                                       &
     &       * SQRT(SAL) - 1.9813D-2 * SAL                             &
     &     )

      BT = 4.106D-4 * SAL / 3.5D1

      FH =                                                             &
     &  + 1.29D0-2.4D-3 * ATEM + 4.61D-4 * SAL**2.D0                   &
     &  - 1.48D-6 * ATEM * SAL**2.D0

      END SUBROUTINE CHEMEQCONST

      END MODULE CHEMEQ
