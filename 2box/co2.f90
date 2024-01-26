      MODULE CO2

      USE CHEMEQ

      IMPLICIT NONE
      CONTAINS

      SUBROUTINE CO2_NIBUN( &
    &     BT,               & !! (in)
    &     K0,               & !! (in)
    &     K1,               & !! (in)
    &     K2,               & !! (in)
    &     KB,               & !! (in)
    &     KW,               & !! (in)
    &     AT,               & !! (in)
    &     CT,               & !! (in)
    &     CO2,              & !! (out)
    &     HCO3,             & !! (out)
    &     CO32,             & !! (out)
    &     PCO2              & !! (out)
    & )

! --- information -----------------------------------------------------
!
!  Calculation of Carbonate system from alkalinity and DIC.
!
! ---------------------------------------------------------------------

      INTEGER :: I

      REAL(8) :: HMIN, HMAX, DELTA
      REAL(8) :: CV1, CV2

      REAL(8), INTENT(IN) :: BT, K0, K1, K2, KW, KB
      REAL(8), INTENT(IN) :: AT, CT
      REAL(8) :: ATX, CTX

      REAL(8) :: H, ATXX, HX
      REAL(8), INTENT(OUT) :: CO2, HCO3, CO32, PCO2

!-- Constant

      HMIN  = 1.0D-14
      HMAX  = 1.0D0
      DELTA = 1.0D-15

!-- Unit conversion:

      DATA CV1 / 1.0250D+3 / !! [mol/kg] to [mol/m^3]
      DATA CV2 / 9.7561D-4 / !! [mol/m^3] to [mol/kg]


!-- Equilibrium constant (if prescribed)

!     BT = 4.16D-4  !! Total Boron : 4.16 * 10^{-4} [mol/kg]
!     K0 = 3.35D-2  !! EC          : 3.35*10^{-2}   [mol/kg atm]
!     K1 = 1.26D-6  !! EC          : 1.26*10^{-6}   [mol/kg]
!     K2 = 8.55D-10 !! EC          : 8.55*10^{-10}  [mol/kg]
!     KW = 3.18D-14 !! EC          : 3.18*10^{-14}  [mol/kg]
!     KB = 2.09D-9  !! EC          : 2.09*10^{-9}   [mol/kg]

!-- Unit conversion

      ATX = CV2 * AT !! [mol/m^3] to [mol/kg]
      CTX = CV2 * CT !! [mol/m^3] to [mol/kg]

!-- CO2, PCO2

      HX = 0.0D0
      DO I = 1, 100000
        H = (HMIN + HMAX) / 2.0D0
        ATXX =                                                          &
     &       + ( ( 2.0D0 * K1 * K2 * CTX )                              &
     &         / ( H**2.0D0 + K1 * H + K1 * K2 ) )                      &
     &       + ( ( H * K1 * CTX ) / ( H**2.0D0 + K1 * H + K1 * K2 ) )   &
     &       + ( ( KB * BT ) / ( H + KB ) ) + ( KW / H ) - H
        IF ( ABS(ATX - ATXX) <= DELTA ) THEN
          HX = H
          EXIT
        ELSE IF ( ATXX < ATX ) THEN
          HMAX = H
        ELSE
          HMIN = H
        ENDIF
      ENDDO
      PCO2 = (HX**2.0D0 * CTX) / (HX**2.0D0 + K1 * HX + K1 * K2) / K0
      CO2  = (HX**2.0D0 * CTX) / (HX**2.0D0 + K1 * HX + K1 * K2) / CV2
      HCO3 = K1 * CO2 / HX
      CO32 = K2 * HCO3 / HX

      END SUBROUTINE CO2_NIBUN

      END MODULE CO2
