#!/usr/bin/env python3
"""
two_box.py

Python port of the simple 2-box ocean biogeochemistry model

Boxes:
  S: surface box
  D: deep box

Run:
  python3 two_box.py
"""
from math import exp, log, log10, sqrt
import math
import argparse
import sys
import numpy as np
import pandas as pd


def chemeq_const(TEM, SAL):
    """Compute equilibrium constants BT, K0, K1, K2, KB, KW, FH
    from temperature (deg C) and salinity (PSU).
    Ported from CHEMEQCONST in `chemeq.f90`.
    Returns (BT, K0, K1, K2, KB, KW, FH)
    """
    ATEM = TEM + 273.15
    TK = ATEM * 1.0e-2
    SK = 2.3517e-2 + (-2.3656e-2 + 4.7036e-3 * TK) * TK

    K0 = math.exp(-6.02409e1 + 9.34517e1 / TK + 2.33585e1 * math.log(TK) + SK * SAL)

    K1 = math.exp(math.log(10.0) * (13.7201 - 3.1334e-2 * ATEM - 3.23576e3 / ATEM
                                     - 1.3e-5 * SAL * ATEM + 1.032e-1 * sqrt(SAL)))

    K2 = math.exp(math.log(10.0) * (
        -5.3719645e3 - 1.671221e0 * ATEM - 2.2913e-1 * SAL
        - 1.83802e1 * math.log10(SAL if SAL > 0 else 1.0) + 1.2837528e5 / ATEM
        + 2.1943005e3 * math.log10(ATEM) + 8.0944e-4 * SAL * ATEM
        + 5.61711e3 * math.log10(SAL if SAL > 0 else 1.0) / ATEM - 2.136e0 * SAL / ATEM
    ))

    KB = math.exp(math.log(10.0) * (-9.26 + 8.86e-3 * SAL + 1e-3 * TEM))

    KW = math.exp(+1.489802e2 - 1.384726e4 / ATEM - 2.36521e1 * math.log(ATEM)
                  + ((-7.92447e1 + 3.29872e3 / ATEM + 1.20408e1 * math.log(ATEM))
                     * sqrt(SAL) - 1.9813e-2 * SAL))

    BT = 4.106e-4 * SAL / 35.0

    FH = 1.29e-2 - 2.4e-3 * ATEM + 4.61e-4 * SAL**2 - 1.48e-6 * ATEM * SAL**2

    return BT, K0, K1, K2, KB, KW, FH


def o2sat(TEM, SAL):
    """Compute oxygen saturation concentration (mol/m^3) from temperature (deg C)
    and salinity, ported from `O2SAT` in `o2.f90`.
    """
    N1 = -1.734292e2
    N2 = +2.496339e2
    N3 = +1.433483e2
    N4 = -2.184920e1
    N5 = -3.309600e-2
    N6 = +1.425900e-2
    N7 = -1.700000e-3

    ATEM = TEM + 273.15
    O2S = math.exp(
        N1 + N2 * 1.0e2 / ATEM + N3 * math.log(ATEM / 1.0e2) + N4 * ATEM / 1.0e2
        + SAL * (N5 + N6 * ATEM / 1.0e2 + N7 * (ATEM / 1.0e2) ** 2)
    )
    O2S = O2S * 4.35e1 * 1.025e-3
    return O2S


def co2_nibun(BT, K0, K1, K2, KW, KB, AT, CT):
    """Solve carbonate system from alkalinity `AT` and DIC `CT`.
    Returns (CO2, HCO3, CO32, PCO2)
    Units: inputs AT,CT in [mol/m^3]"""
    # conversion factors from Fortran
    CV2 = 9.7561e-4  # [mol/m^3] to [mol/kg]
    ATX = CV2 * AT
    CTX = CV2 * CT

    HMIN = 1.0e-14
    HMAX = 1.0
    DELTA = 1.0e-15

    HX = 0.0
    H_low = HMIN
    H_high = HMAX
    for _ in range(100000):
        H = 0.5 * (H_low + H_high)
        denom = H * H + K1 * H + K1 * K2
        AT_calc = (2.0 * K1 * K2 * CTX) / denom + (H * K1 * CTX) / denom + (KB * BT) / (H + KB) + KW / H - H
        if abs(ATX - AT_calc) <= DELTA:
            HX = H
            break
        if AT_calc < ATX:
            H_high = H
        else:
            H_low = H
    else:
        HX = 0.5 * (H_low + H_high)

    denom2 = HX * HX + K1 * HX + K1 * K2
    PCO2 = (HX * HX * CTX) / denom2 / K0
    CO2 = (HX * HX * CTX) / denom2 / CV2
    HCO3 = K1 * CO2 / HX
    CO32 = K2 * HCO3 / HX
    return CO2, HCO3, CO32, PCO2


def create_dataframe(out: dict, box1_name: str = 'Surface (S)', box2_name: str = 'Deep (D)') -> pd.DataFrame:
    """Create a pandas DataFrame with human-friendly labels from model output dict."""
    rows = [
        ('Temperature (K)', box1_name, out.get('TEMS_K')),
        ('Temperature (K)', box2_name, out.get('TEMD_K')),
        ('Salinity (PSU)', box1_name, out.get('SALS_PSU')),
        ('Salinity (PSU)', box2_name, out.get('SALD_PSU')),
        ('PO4 (umol/kg)', box1_name, out.get('PO4S_umolkg')),
        ('PO4 (umol/kg)', box2_name, out.get('PO4D_umolkg')),
        ('DIC (umol/kg)', box1_name, out.get('DICS_umolkg')),
        ('DIC (umol/kg)', box2_name, out.get('DICD_umolkg')),
        ('ALK (ueq/kg)', box1_name, out.get('ALKS_ueqkg')),
        ('ALK (ueq/kg)', box2_name, out.get('ALKD_ueqkg')),
        ('O2 (umol/kg)', box1_name, out.get('DO2S_umolkg')),
        ('O2 (umol/kg)', box2_name, out.get('DO2D_umolkg')),
        ('AOU (umol/kg)', box1_name, out.get('AOU_S_umolkg')),
        ('AOU (umol/kg)', box2_name, out.get('AOU_D_umolkg')),
        ('CO2 (umol/kg)', box1_name, out.get('CO2S_umolkg')),
        ('CO2 (umol/kg)', box2_name, out.get('CO2D_umolkg')),
        ('HCO3 (umol/kg)', box1_name, out.get('HCO3S_umolkg')),
        ('HCO3 (umol/kg)', box2_name, out.get('HCO3D_umolkg')),
        ('CO3 (umol/kg)', box1_name, out.get('CO32S_umolkg')),
        ('CO3 (umol/kg)', box2_name, out.get('CO32D_umolkg')),
        ('PCO2 (ppmv)', 'Surface', out.get('PCO2S_ppmv')),
        ('PCO2 (ppmv)', 'Deep', out.get('PCO2D_ppmv')),
        ('PCO2 (ppmv)', 'Atmosphere', out.get('PCO2A_ppmv')),
        ('Total Carbon (PgC)', 'Initial', out.get('TOCINI_PgC')),
        ('Total Carbon (PgC)', 'Final', out.get('TOCFIN_PgC')),
    ]
    df = pd.DataFrame(rows, columns=['Quantity', 'Box', 'Value'])
    return df


def run_model(max_iter=100000, verbose=True, box1_name='Surface (S)', box2_name='Deep (D)'):
    # Parameters
    CV1 = 1.0250e3
    CV2 = 9.7561e-4
    CV3 = 1.0000e6
    CV4 = 3.1536e7
    CV5 = 8.64e4

    # default temps and salts
    TEMS = 20.0
    TEMD = 2.0
    SALS = 34.7
    SALD = 34.7

    VOCN = 1.292e18
    AOCN = 3.49e14
    VATM = 1.773e20
    ZOCNS = 250.0
    ZOCND = 4000.0

    TRAN = 0.0
    FSD = 4.096e7
    PVS = 3.0
    R = 1.0 / CV4
    LF = 5.0e-1
    HSC = 2.0e-8
    RRC = 0.288
    RCP = 106.0
    RNP = 16.0
    RO2P = 172.0
    CEPS = 0.611
    PCO2A = 280.0
    AOGE = False
    OFIXDB = True

    # initial tracers (given in mol/kg -> convert to mol/m^3)
    PO4S = 2.20e-6 * CV1
    PO4D = 2.20e-6 * CV1
    ALKS = 2.371e-3 * CV1
    ALKD = 2.371e-3 * CV1
    DICS = 2.258e-3 * CV1
    DICD = 2.258e-3 * CV1
    DO2S = 1.70e-4 * CV1
    DO2D = 1.70e-4 * CV1

    AOCNS = AOCN * 1.0
    VOCNS = AOCNS * ZOCNS
    VOCND = VOCN - VOCNS
    EPS = (CEPS / RCP) * AOCNS
    FAS = PVS * AOCNS

    # unit convert
    EPS = EPS / CV4
    FAS = FAS / CV5

    # state variables
    PCO2A_atm = PCO2A / CV3

    for I in range(1, max_iter + 1):
        # equilibrium constants
        BTS, K0S, K1S, K2S, KBS, KWS, FHS = chemeq_const(TEMS, SALS)
        BTD, K0D, K1D, K2D, KBD, KWD, FHD = chemeq_const(TEMD, SALD)

        if I == 1:
            CO2S, HCO3S, CO32S, PCO2S = co2_nibun(BTS, K0S, K1S, K2S, KWS, KBS, ALKS, DICS)
            CO2D, HCO3D, CO32D, PCO2D = co2_nibun(BTD, K0D, K1D, K2D, KWD, KBD, ALKD, DICD)
            TOCINI = DICS * VOCNS + DICD * VOCND + PCO2A_atm * VATM

        TEMSX = TEMS
        TEMDX = TEMD

        PO4SX = PO4S + ((TRAN + FSD) * (PO4D - PO4S) - EPS) * (1.0 * (1.0) / VOCNS) * CV4 * 0.0
        # To keep behavior simple, mirror Fortran's formula but keep timestep implicit: we won't change volumes here
        # Use Fortran-style updating with DT
        DT = 8.64e4
        PO4SX = PO4S + (((TRAN + FSD) * (PO4D - PO4S) - EPS) * (DT / VOCNS))
        if not OFIXDB:
            PO4DX = PO4D + (((TRAN + FSD) * (PO4S - PO4D) + EPS) * (DT / VOCND))
        else:
            PO4DX = PO4D

        ALKSX = ALKS + (((TRAN + FSD) * (ALKD - ALKS) + (2.0 * RRC * RCP - RNP) * (-EPS)) * (DT / VOCNS))
        if not OFIXDB:
            ALKDX = ALKD + (((TRAN + FSD) * (ALKS - ALKD) + (2.0 * RRC * RCP - RNP) * (+EPS)) * (DT / VOCND))
        else:
            ALKDX = ALKD

        if AOGE:
            DICSX = DICS + (((TRAN + FSD) * (DICD - DICS) + (1.0 + RRC) * RCP * (-EPS) + FAS * CV1 * K0S * (PCO2A_atm - PCO2S)) * (DT / VOCNS))
        else:
            DICSX = DICS + (((TRAN + FSD) * (DICD - DICS) + (1.0 + RRC) * RCP * (-EPS)) * (DT / VOCNS))

        if not OFIXDB:
            DICDX = DICD + (((TRAN + FSD) * (DICS - DICD) + (1.0 + RRC) * RCP * (+EPS)) * (DT / VOCND))
        else:
            DICDX = DICD

        DO2SX = DO2S + (((TRAN + FSD) * (DO2D - DO2S) - RO2P * (-EPS)) * (DT / VOCNS))
        if not OFIXDB:
            DO2DX = DO2D + (((TRAN + FSD) * (DO2S - DO2D) - RO2P * (+EPS)) * (DT / VOCND))
        else:
            DO2DX = DO2D

        CO2S, HCO3S, CO32S, PCO2S = co2_nibun(BTS, K0S, K1S, K2S, KWS, KBS, ALKSX, DICSX)
        CO2D, HCO3D, CO32D, PCO2D = co2_nibun(BTD, K0D, K1D, K2D, KWD, KBD, ALKDX, DICDX)

        if AOGE:
            PCO2AX = PCO2A_atm + ((K0S * FAS * CV1 * (PCO2S - PCO2A_atm)) * (DT / VATM))
        else:
            PCO2AX = PCO2A_atm

        # simple convergence check on PO4 relative change
        rel_change = abs(PO4SX - PO4S) / max(abs(PO4S), 1.0e-30)
        if rel_change <= 1.0e-6:
            break

        # Update state
        TEMS = TEMSX
        TEMD = TEMDX
        PO4S = PO4SX
        PO4D = PO4DX
        DICS = DICSX
        DICD = DICDX
        ALKS = ALKSX
        ALKD = ALKDX
        DO2S = DO2SX
        DO2D = DO2DX
        PCO2A_atm = PCO2AX

    # finalize outputs (convert to units similar to Fortran prints)
    TOCFIN = DICS * VOCNS + DICD * VOCND + PCO2A_atm * VATM

    O2SATS = o2sat(TEMS, SALS)
    O2SATD = o2sat(TEMD, SALD)
    AOUS = (O2SATS - DO2S)
    AOUD = (O2SATD - DO2D)

    # convert back to umol/kg-like outputs using CV2
    out = {
        'TEMS_K': TEMS,
        'TEMD_K': TEMD,
        'SALS_PSU': SALS,
        'SALD_PSU': SALD,
        'EPS_PgC_yr': CV4 * EPS * RCP * 1.0e-15 * 12.0 if 'EPS' in locals() else None,
        'PO4S_umolkg': CV2 * 1.0e6 * PO4S,
        'PO4D_umolkg': CV2 * 1.0e6 * PO4D,
        'DICS_umolkg': CV2 * 1.0e6 * DICS,
        'DICD_umolkg': CV2 * 1.0e6 * DICD,
        'ALKS_ueqkg': CV2 * 1.0e6 * ALKS,
        'ALKD_ueqkg': CV2 * 1.0e6 * ALKD,
        'DO2S_umolkg': CV2 * 1.0e6 * DO2S,
        'DO2D_umolkg': CV2 * 1.0e6 * DO2D,
        'AOU_S_umolkg': CV2 * 1.0e6 * AOUS,
        'AOU_D_umolkg': CV2 * 1.0e6 * AOUD,
        'CO2S_umolkg': CV2 * 1.0e6 * CO2S,
        'CO2D_umolkg': CV2 * 1.0e6 * CO2D,
        'HCO3S_umolkg': CV2 * 1.0e6 * HCO3S,
        'HCO3D_umolkg': CV2 * 1.0e6 * HCO3D,
        'CO32S_umolkg': CV2 * 1.0e6 * CO32S,
        'CO32D_umolkg': CV2 * 1.0e6 * CO32D,
        'PCO2S_ppmv': CV3 * PCO2S,
        'PCO2D_ppmv': CV3 * PCO2D,
        'PCO2A_ppmv': CV3 * PCO2A_atm,
        'TOCINI_PgC': TOCINI * 12.0e-15 if 'TOCINI' in locals() else None,
        'TOCFIN_PgC': TOCFIN * 12.0e-15,
    }

    if verbose:
        df = create_dataframe(out, box1_name=box1_name, box2_name=box2_name)
        try:
            from tabulate import tabulate
            print(tabulate(df, headers='keys', tablefmt='github', showindex=False))
        except Exception:
            print(df.to_string(index=False))

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--iters', type=int, default=100000)
    parser.add_argument('--box1-name', type=str, default='Surface (S)', help='Name for box 1')
    parser.add_argument('--box2-name', type=str, default='Deep (D)', help='Name for box 2')
    parser.add_argument('--save-csv', type=str, default=None, help='Save output table to CSV')
    args = parser.parse_args()
    out = run_model(max_iter=args.iters, box1_name=args.box1_name, box2_name=args.box2_name)
    if args.save_csv:
        df = create_dataframe(out, box1_name=args.box1_name, box2_name=args.box2_name)
        df.to_csv(args.save_csv, index=False)


if __name__ == '__main__':
    main()
