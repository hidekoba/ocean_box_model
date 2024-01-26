      MODULE O2

      IMPLICIT NONE
      CONTAINS

      SUBROUTINE O2SAT( &
    &     TEM,          & !! (in)
    &     SAL,          & !! (in)
    &     O2S           & !! (out)
    & )

! --- information -----------------------------------------------------
!
!  Calculation of oxygen saturation conc.
!  from temperature and salinity.
!
! ---------------------------------------------------------------------

      REAL(8), INTENT(IN) :: TEM, SAL
      REAL(8) :: ATEM
      REAL(8), INTENT(OUT) :: O2S
      REAL(8) :: N1, N2, N3, N4, N5, N6, N7

!-- Constant

      DATA N1/ -1.734292D+2 /
      DATA N2/ +2.496339D+2 /
      DATA N3/ +1.433483D+2 /
      DATA N4/ -2.184920D+1 /
      DATA N5/ -3.309600D-2 /
      DATA N6/ +1.425900D-2 /
      DATA N7/ -1.700000D-3 /

!-- Calculation of oxygen saturation conc.

      ATEM = TEM + 273.15D0

      O2S =                                                             &
     &  EXP(                                                            &
     &       + N1                                                       &
     &       + N2 * 1.0D2 / ATEM                                        &
     &       + N3 * LOG(ATEM / 1.0D2)                                   &
     &       + N4 * ATEM / 1.0D2                                        &
     &       + SAL *                                                    &
     &           (                                                      &
     &             + N5                                                 &
     &             + N6 * ATEM / 1.0D2                                  &
     &             + N7 * (ATEM / 1.0D2)**2.0D0                         &
     &           )                                                      &
     &     )
      O2S = O2S * 4.35D1 * 1.025D-3

      END SUBROUTINE O2SAT

      END MODULE O2
