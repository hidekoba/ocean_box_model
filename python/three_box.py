#!/usr/bin/env python3
"""
three_box.py

Python port of the simple 3-box ocean biogeochemistry model

Boxes:
  H: high-latitude surface box
  L: low-latitude surface box
  D: deep box

Run:
  python3 three_box.py
  python3 three_box.py --save-csv output.csv
"""
import argparse
import math
from math import sqrt
import numpy as np
import pandas as pd


def chemeq_const(TEM, SAL):
    """Compute carbonate equilibrium constants from temperature (deg C)
    and salinity (PSU), ported from the two-box Python/Fortran helper.
    Returns (BT, K0, K1, K2, KB, KW, FH).
    """
    ATEM = TEM + 273.15
    TK = ATEM * 1.0e-2
    SK = 2.3517e-2 + (-2.3656e-2 + 4.7036e-3 * TK) * TK

    K0 = math.exp(-6.02409e1 + 9.34517e1 / TK + 2.33585e1 * math.log(TK) + SK * SAL)

    K1 = math.exp(math.log(10.0) * (
        13.7201 - 3.1334e-2 * ATEM - 3.23576e3 / ATEM
        - 1.3e-5 * SAL * ATEM + 1.032e-1 * sqrt(SAL)
    ))

    safe_sal = SAL if SAL > 0 else 1.0
    K2 = math.exp(math.log(10.0) * (
        -5.3719645e3 - 1.671221e0 * ATEM - 2.2913e-1 * SAL
        - 1.83802e1 * math.log10(safe_sal) + 1.2837528e5 / ATEM
        + 2.1943005e3 * math.log10(ATEM) + 8.0944e-4 * SAL * ATEM
        + 5.61711e3 * math.log10(safe_sal) / ATEM - 2.136e0 * SAL / ATEM
    ))

    KB = math.exp(math.log(10.0) * (-9.26 + 8.86e-3 * SAL + 1e-3 * TEM))

    KW = math.exp(
        +1.489802e2 - 1.384726e4 / ATEM - 2.36521e1 * math.log(ATEM)
        + ((-7.92447e1 + 3.29872e3 / ATEM + 1.20408e1 * math.log(ATEM))
           * sqrt(SAL) - 1.9813e-2 * SAL)
    )

    BT = 4.106e-4 * SAL / 35.0

    FH = 1.29e-2 - 2.4e-3 * ATEM + 4.61e-4 * SAL**2 - 1.48e-6 * ATEM * SAL**2

    return BT, K0, K1, K2, KB, KW, FH


def o2sat(TEM, SAL):
    """Oxygen saturation concentration [mol/m^3] from temperature (deg C)
    and salinity, ported from O2SAT.
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


def co2_nibun(BT, K0, K1, K2, KB, KW, AT, CT):
    """Solve carbonate chemistry from alkalinity AT and DIC CT.

    Inputs:
      AT, CT: mol/m^3
    Returns:
      CO2, HCO3, CO32: mol/m^3
      PCO2: atm
    """
    CV2 = 9.7561e-4  # mol/m^3 -> mol/kg
    ATX = CV2 * AT
    CTX = CV2 * CT

    HMIN = 1.0e-14
    HMAX = 1.0
    DELTA = 1.0e-15

    H_low = HMIN
    H_high = HMAX
    HX = 0.5 * (H_low + H_high)

    for _ in range(100000):
        H = 0.5 * (H_low + H_high)
        denom = H * H + K1 * H + K1 * K2
        AT_calc = (
            (2.0 * K1 * K2 * CTX) / denom
            + (H * K1 * CTX) / denom
            + (KB * BT) / (H + KB)
            + KW / H
            - H
        )

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


def default_params():
    """Parameters"""
    CV1 = 1.0250e3
    CV2 = 9.7561e-4
    CV3 = 1.0000e6
    CV4 = 3.1536e7
    CV5 = 8.64e4

    params = {
        "DT": 8.64e4,
        "DELTA": 1.0e-6,

        "CV1": CV1,
        "CV2": CV2,
        "CV3": CV3,
        "CV4": CV4,
        "CV5": CV5,

        # Temperature [deg C] in the original comments say K, but values are deg C.
        "TEMH": 2.0,
        "TEML": 21.5,
        "TEMD": 1.75,
        "SALH": 34.7,
        "SALL": 34.7,
        "SALD": 34.7,

        # Volumes/areas
        "VOCN": 1.292e18,
        "AOCN": 3.49e14,
        "VATM": 1.773e20,
        "ZOCNH": 2.5e2,
        "ZOCNL": 1.0e2,
        "ZOCND": 4.0e3,
        "DEH": 2.5e2,
        "DEL": 1.0e2,

        # Transport/mixing [m^3/s]
        "T": 2.0e7,
        "FHD": 6.0e7,
        "FLD": 0.0,
        "FLH": 0.0,

        # Piston velocity [m/day]
        "PVH": 3.0,
        "PVL": 3.0,

        # Biological parameters
        "R": 1.0 / CV4,
        "LF": 5.0e-1,
        "HSC": 2.0e-8 * CV1,
        "RRC": 0.20,
        "RCP": 106.0,
        "RNP": 16.0,
        "RO2P": 172.0,
        "CEPH": 7.5e-1,  # mol C m^-2 yr^-1

        # Initial atmospheric PCO2 [ppmv]
        "PCO2A": 280.0,

        # Initial tracers [mol/kg or eq/kg]
        "PO4H0": 2.10e-6,
        "PO4L0": 2.10e-6,
        "PO4D0": 2.10e-6,
        "ALKH0": 2.374e-3,
        "ALKL0": 2.374e-3,
        "ALKD0": 2.374e-3,
        "DICH0": 2.235e-3,
        "DICL0": 2.235e-3,
        "DICD0": 2.235e-3,
        "DO2H0": 1.60e-4,
        "DO2L0": 1.60e-4,
        "DO2D0": 1.60e-4,
    }
    return params


def initialize(params):
    """Initialize volumes, fluxes, and state variables in model units."""
    CV1 = params["CV1"]
    CV3 = params["CV3"]
    CV4 = params["CV4"]
    CV5 = params["CV5"]

    # Areas: H 15%, L 85%
    AOCNH = params["AOCN"] * 1.5e-1
    AOCNL = params["AOCN"] * 8.5e-1
    VOCNH = AOCNH * params["ZOCNH"]
    VOCNL = AOCNL * params["ZOCNL"]
    VOCND = params["VOCN"] - VOCNH - VOCNL

    # Biological sinking flux [mol P/s]
    EPH = (params["CEPH"] / params["RCP"]) * AOCNH / CV4

    # Low-lat biological export
    EPL = None

    # Gas exchange [m^3/s]
    FAH = params["PVH"] * AOCNH / CV5
    FAL = params["PVL"] * AOCNL / CV5

    state = {
        "PO4H": params["PO4H0"] * CV1,
        "PO4L": params["PO4L0"] * CV1,
        "PO4D": params["PO4D0"] * CV1,
        "ALKH": params["ALKH0"] * CV1,
        "ALKL": params["ALKL0"] * CV1,
        "ALKD": params["ALKD0"] * CV1,
        "DICH": params["DICH0"] * CV1,
        "DICL": params["DICL0"] * CV1,
        "DICD": params["DICD0"] * CV1,
        "DO2H": params["DO2H0"] * CV1,
        "DO2L": params["DO2L0"] * CV1,
        "DO2D": params["DO2D0"] * CV1,
        "PCO2A": params["PCO2A"] / CV3,
    }

    aux = {
        "AOCNH": AOCNH,
        "AOCNL": AOCNL,
        "VOCNH": VOCNH,
        "VOCNL": VOCNL,
        "VOCND": VOCND,
        "EPH": EPH,
        "EPL": EPL,
        "FAH": FAH,
        "FAL": FAL,
    }
    return state, aux


def carbonate_all(state, params):
    """Compute carbonate system for H, L, D boxes."""
    constants = {}
    outputs = {}
    for label, temp_key, sal_key, alk_key, dic_key in [
        ("H", "TEMH", "SALH", "ALKH", "DICH"),
        ("L", "TEML", "SALL", "ALKL", "DICL"),
        ("D", "TEMD", "SALD", "ALKD", "DICD"),
    ]:
        BT, K0, K1, K2, KB, KW, FH = chemeq_const(params[temp_key], params[sal_key])
        constants[label] = {"BT": BT, "K0": K0, "K1": K1, "K2": K2, "KB": KB, "KW": KW, "FH": FH}
        CO2, HCO3, CO32, PCO2 = co2_nibun(BT, K0, K1, K2, KB, KW, state[alk_key], state[dic_key])
        outputs[label] = {"CO2": CO2, "HCO3": HCO3, "CO32": CO32, "PCO2": PCO2}
    return constants, outputs


def step(state, params, aux):
    """One forward-Euler step using the intended equations"""
    DT = params["DT"]
    T = params["T"]
    FHD = params["FHD"]
    FLD = params["FLD"]
    FLH = params["FLH"]
    RRC = params["RRC"]
    RCP = params["RCP"]
    RNP = params["RNP"]
    RO2P = params["RO2P"]

    VOCNH = aux["VOCNH"]
    VOCNL = aux["VOCNL"]
    VOCND = aux["VOCND"]
    EPH = aux["EPH"]
    FAH = aux["FAH"]
    FAL = aux["FAL"]

    constants, carb = carbonate_all(state, params)
    PCO2H = carb["H"]["PCO2"]
    PCO2L = carb["L"]["PCO2"]

    O2SATH = o2sat(params["TEMH"], params["SALH"])
    O2SATL = o2sat(params["TEML"], params["SALL"])
    O2SATD = o2sat(params["TEMD"], params["SALD"])

    # Low-lat biological export from the commented Fortran line:
    # EPL = T * PO4D
    # EPL = T * state["PO4D"]
    EPL = (
        params["R"]
        * params["DEL"]
        * params["LF"]
        * state["PO4L"]
        * (state["PO4L"] / (params["HSC"] + state["PO4L"]))
        * VOCNL
    )

    # PO4
    PO4HX = state["PO4H"] + (
        (T + FLH) * (state["PO4L"] - state["PO4H"])
        + FHD * (state["PO4D"] - state["PO4H"])
        - EPH
    ) * (DT / VOCNH)

    PO4LX = state["PO4L"] + (
        (T + FLD) * (state["PO4D"] - state["PO4L"])
        + FLH * (state["PO4H"] - state["PO4L"])
        - EPL
    ) * (DT / VOCNL)

    PO4DX = state["PO4D"] + (
        (T + FHD) * (state["PO4H"] - state["PO4D"])
        + FLD * (state["PO4L"] - state["PO4D"])
        + (EPH + EPL)
    ) * (DT / VOCND)

    alk_factor = 2.0 * RRC * RCP - RNP
    dic_factor = (1.0 + RRC) * RCP

    # Alkalinity
    ALKHX = state["ALKH"] + (
        (T + FLH) * (state["ALKL"] - state["ALKH"])
        + FHD * (state["ALKD"] - state["ALKH"])
        - alk_factor * EPH
    ) * (DT / VOCNH)

    ALKLX = state["ALKL"] + (
        (T + FLD) * (state["ALKD"] - state["ALKL"])
        + FLH * (state["ALKH"] - state["ALKL"])
        - alk_factor * EPL
    ) * (DT / VOCNL)

    ALKDX = state["ALKD"] + (
        (T + FHD) * (state["ALKH"] - state["ALKD"])
        + FLD * (state["ALKL"] - state["ALKD"])
        + alk_factor * (EPH + EPL)
    ) * (DT / VOCND)

    # DIC with air-sea gas exchange
    DICHX = state["DICH"] + (
        (T + FLH) * (state["DICL"] - state["DICH"])
        + FHD * (state["DICD"] - state["DICH"])
        - dic_factor * EPH
        + FAH * params["CV1"] * constants["H"]["K0"] * (state["PCO2A"] - PCO2H)
    ) * (DT / VOCNH)

    DICLX = state["DICL"] + (
        (T + FLD) * (state["DICD"] - state["DICL"])
        + FLH * (state["DICH"] - state["DICL"])
        - dic_factor * EPL
        + FAL * params["CV1"] * constants["L"]["K0"] * (state["PCO2A"] - PCO2L)
    ) * (DT / VOCNL)

    DICDX = state["DICD"] + (
        (T + FHD) * (state["DICH"] - state["DICD"])
        + FLD * (state["DICL"] - state["DICD"])
        + dic_factor * (EPH + EPL)
    ) * (DT / VOCND)

    # Oxygen.
    # The H equation is copied from the commented Fortran, replacing RDO2P with RO2P.
    DO2HX = state["DO2H"] + (
        (T + FLH) * (O2SATL - state["DO2H"])
        + FHD * (O2SATH - state["DO2H"])
        + RO2P * EPH
        + FAH * (O2SATH - state["DO2H"])
    ) * (DT / VOCNH)

    DO2LX = state["DO2L"] + (
        (T + FLD) * (state["DO2D"] - state["DO2L"])
        + FLH * (state["DO2H"] - state["DO2L"])
        + RO2P * EPL
        + FAL * (O2SATL - state["DO2L"])
    ) * (DT / VOCNL)

    DO2DX = state["DO2D"] + (
        (T + FHD) * (state["DO2H"] - state["DO2D"])
        + FLD * (O2SATL - state["DO2D"])
        - RO2P * (EPH + EPL)
    ) * (DT / VOCND)

    # Carbonate system after tracer update, needed for atmospheric PCO2 update.
    tmp_state = dict(state)
    tmp_state.update({
        "ALKH": ALKHX, "ALKL": ALKLX, "ALKD": ALKDX,
        "DICH": DICHX, "DICL": DICLX, "DICD": DICDX,
    })
    constants_x, carb_x = carbonate_all(tmp_state, params)

    PCO2AX = state["PCO2A"] + (
        params["CV1"] * constants_x["H"]["K0"] * FAH * (carb_x["H"]["PCO2"] - state["PCO2A"])
        + params["CV1"] * constants_x["L"]["K0"] * FAL * (carb_x["L"]["PCO2"] - state["PCO2A"])
    ) * (DT / params["VATM"])

    new_state = {
        "PO4H": PO4HX, "PO4L": PO4LX, "PO4D": PO4DX,
        "ALKH": ALKHX, "ALKL": ALKLX, "ALKD": ALKDX,
        "DICH": DICHX, "DICL": DICLX, "DICD": DICDX,
        "DO2H": DO2HX, "DO2L": DO2LX, "DO2D": DO2DX,
        "PCO2A": PCO2AX,
    }

    aux = dict(aux)
    aux["EPL"] = EPL
    return new_state, aux


def make_output(state, params, aux, tocini=None, nstep=None):
    """Convert final state to a dictionary with readable units."""
    CV2 = params["CV2"]
    CV3 = params["CV3"]
    CV4 = params["CV4"]

    constants, carb = carbonate_all(state, params)

    O2SATH = o2sat(params["TEMH"], params["SALH"])
    O2SATL = o2sat(params["TEML"], params["SALL"])
    O2SATD = o2sat(params["TEMD"], params["SALD"])

    AOUH = O2SATH - state["DO2H"]
    AOUL = O2SATL - state["DO2L"]
    AOUD = O2SATD - state["DO2D"]

    to_umolkg = lambda x: CV2 * 1.0e6 * x
    to_ppmv = lambda x: CV3 * x

    tocfin = (
        state["DICH"] * aux["VOCNH"]
        + state["DICL"] * aux["VOCNL"]
        + state["DICD"] * aux["VOCND"]
        + state["PCO2A"] * params["VATM"]
    )

    out = {
        "TIMESTEP": nstep,

        "Temperature_H_degC": params["TEMH"],
        "Temperature_L_degC": params["TEML"],
        "Temperature_D_degC": params["TEMD"],
        "Salinity_H_PSU": params["SALH"],
        "Salinity_L_PSU": params["SALL"],
        "Salinity_D_PSU": params["SALD"],

        "EP_H_PgC_yr": CV4 * aux["EPH"] * params["RCP"] * 1.0e-15 * 12.0,
        "EP_L_PgC_yr": CV4 * (aux["EPL"] if aux["EPL"] is not None else 0.0) * params["RCP"] * 1.0e-15 * 12.0,

        "PO4_H_umolkg": to_umolkg(state["PO4H"]),
        "PO4_L_umolkg": to_umolkg(state["PO4L"]),
        "PO4_D_umolkg": to_umolkg(state["PO4D"]),

        "DIC_H_umolkg": to_umolkg(state["DICH"]),
        "DIC_L_umolkg": to_umolkg(state["DICL"]),
        "DIC_D_umolkg": to_umolkg(state["DICD"]),

        "ALK_H_ueqkg": to_umolkg(state["ALKH"]),
        "ALK_L_ueqkg": to_umolkg(state["ALKL"]),
        "ALK_D_ueqkg": to_umolkg(state["ALKD"]),

        "DO2_H_umolkg": to_umolkg(state["DO2H"]),
        "DO2_L_umolkg": to_umolkg(state["DO2L"]),
        "DO2_D_umolkg": to_umolkg(state["DO2D"]),

        "AOU_H_umolkg": to_umolkg(AOUH),
        "AOU_L_umolkg": to_umolkg(AOUL),
        "AOU_D_umolkg": to_umolkg(AOUD),

        "CO2_H_umolkg": to_umolkg(carb["H"]["CO2"]),
        "CO2_L_umolkg": to_umolkg(carb["L"]["CO2"]),
        "CO2_D_umolkg": to_umolkg(carb["D"]["CO2"]),

        "HCO3_H_umolkg": to_umolkg(carb["H"]["HCO3"]),
        "HCO3_L_umolkg": to_umolkg(carb["L"]["HCO3"]),
        "HCO3_D_umolkg": to_umolkg(carb["D"]["HCO3"]),

        "CO32_H_umolkg": to_umolkg(carb["H"]["CO32"]),
        "CO32_L_umolkg": to_umolkg(carb["L"]["CO32"]),
        "CO32_D_umolkg": to_umolkg(carb["D"]["CO32"]),

        "PCO2_H_ppmv": to_ppmv(carb["H"]["PCO2"]),
        "PCO2_L_ppmv": to_ppmv(carb["L"]["PCO2"]),
        "PCO2_D_ppmv": to_ppmv(carb["D"]["PCO2"]),
        "PCO2_ATM_ppmv": to_ppmv(state["PCO2A"]),

        "Total_Carbon_Initial_PgC": tocini * 12.0e-15,
        "Total_Carbon_Final_PgC": tocfin * 12.0e-15,

    }
    return out


def create_dataframe(out):
    rows = []
    order = [
        ("Temperature (deg C)", "H", "Temperature_H_degC"),
        ("Temperature (deg C)", "L", "Temperature_L_degC"),
        ("Temperature (deg C)", "D", "Temperature_D_degC"),
        ("Salinity (PSU)", "H", "Salinity_H_PSU"),
        ("Salinity (PSU)", "L", "Salinity_L_PSU"),
        ("Salinity (PSU)", "D", "Salinity_D_PSU"),
        ("Export production (PgC/yr)", "H", "EP_H_PgC_yr"),
        ("Export production (PgC/yr)", "L", "EP_L_PgC_yr"),
        ("PO4 (umol/kg)", "H", "PO4_H_umolkg"),
        ("PO4 (umol/kg)", "L", "PO4_L_umolkg"),
        ("PO4 (umol/kg)", "D", "PO4_D_umolkg"),
        ("DIC (umol/kg)", "H", "DIC_H_umolkg"),
        ("DIC (umol/kg)", "L", "DIC_L_umolkg"),
        ("DIC (umol/kg)", "D", "DIC_D_umolkg"),
        ("ALK (ueq/kg)", "H", "ALK_H_ueqkg"),
        ("ALK (ueq/kg)", "L", "ALK_L_ueqkg"),
        ("ALK (ueq/kg)", "D", "ALK_D_ueqkg"),
        ("O2 (umol/kg)", "H", "DO2_H_umolkg"),
        ("O2 (umol/kg)", "L", "DO2_L_umolkg"),
        ("O2 (umol/kg)", "D", "DO2_D_umolkg"),
        ("AOU (umol/kg)", "H", "AOU_H_umolkg"),
        ("AOU (umol/kg)", "L", "AOU_L_umolkg"),
        ("AOU (umol/kg)", "D", "AOU_D_umolkg"),
        ("CO2 (umol/kg)", "H", "CO2_H_umolkg"),
        ("CO2 (umol/kg)", "L", "CO2_L_umolkg"),
        ("CO2 (umol/kg)", "D", "CO2_D_umolkg"),
        ("HCO3 (umol/kg)", "H", "HCO3_H_umolkg"),
        ("HCO3 (umol/kg)", "L", "HCO3_L_umolkg"),
        ("HCO3 (umol/kg)", "D", "HCO3_D_umolkg"),
        ("CO3 (umol/kg)", "H", "CO32_H_umolkg"),
        ("CO3 (umol/kg)", "L", "CO32_L_umolkg"),
        ("CO3 (umol/kg)", "D", "CO32_D_umolkg"),
        ("PCO2 (ppmv)", "H", "PCO2_H_ppmv"),
        ("PCO2 (ppmv)", "L", "PCO2_L_ppmv"),
        ("PCO2 (ppmv)", "D", "PCO2_D_ppmv"),
        ("PCO2 (ppmv)", "Atmosphere", "PCO2_ATM_ppmv"),
        ("Total Carbon (PgC)", "Initial", "Total_Carbon_Initial_PgC"),
        ("Total Carbon (PgC)", "Final", "Total_Carbon_Final_PgC"),
    ]
    for quantity, box, key in order:
        rows.append((quantity, box, out.get(key)))

    df = pd.DataFrame(rows, columns=["Quantity", "Box", "Value"])

    df["Value"] = df["Value"].map(
        lambda x: f"{x:.3f}" if isinstance(x, (int, float)) else x
    )

   #return pd.DataFrame(rows, columns=["Quantity", "Box", "Value"])
    return df


def run_model(max_iter=100000, verbose=True, params=None, return_history=False):
    """Run the 3-box model to approximate steady state."""
    if params is None:
        params = default_params()

    state, aux = initialize(params)

    tocini = (
        state["DICH"] * aux["VOCNH"]
        + state["DICL"] * aux["VOCNL"]
        + state["DICD"] * aux["VOCND"]
        + state["PCO2A"] * params["VATM"]
    )

    history = []
    nstep = max_iter
    for i in range(1, max_iter + 1):
        old_po4h = state["PO4H"]
        state_new, aux_new = step(state, params, aux)

        if return_history and (i == 1 or i % 100 == 0):
            hist_out = make_output(state_new, params, aux_new, tocini=tocini, nstep=i)
            history.append({
                "step": i,
                "PCO2_ATM_ppmv": hist_out["PCO2_ATM_ppmv"],
                "PO4_H_umolkg": hist_out["PO4_H_umolkg"],
                "PO4_L_umolkg": hist_out["PO4_L_umolkg"],
                "PO4_D_umolkg": hist_out["PO4_D_umolkg"],
                "DIC_H_umolkg": hist_out["DIC_H_umolkg"],
                "DIC_L_umolkg": hist_out["DIC_L_umolkg"],
                "DIC_D_umolkg": hist_out["DIC_D_umolkg"],
            })

        rel = abs((state_new["PO4H"] / old_po4h) - 1.0) if old_po4h != 0 else abs(state_new["PO4H"])
        state, aux = state_new, aux_new

        if rel <= params["DELTA"]:
            nstep = i
            break

    out = make_output(state, params, aux, tocini=tocini, nstep=nstep)

    if verbose:
        df = create_dataframe(out)
        try:
            from tabulate import tabulate
            print(f"TIMESTEP = {nstep}")
            print(tabulate(df, headers="keys", tablefmt="github", showindex=False, floatfmt=".6g"))
        except Exception:
            print(f"TIMESTEP = {nstep}")
            print(df.to_string(index=False))

    if return_history:
        return out, pd.DataFrame(history)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iters", type=int, default=100000)
    parser.add_argument("--save-csv", type=str, default=None)
    parser.add_argument("--save-history", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    if args.save_history:
        out, hist = run_model(max_iter=args.iters, verbose=not args.quiet, return_history=True)
        hist.to_csv(args.save_history, index=False)
    else:
        out = run_model(max_iter=args.iters, verbose=not args.quiet)

    if args.save_csv:
        create_dataframe(out).to_csv(args.save_csv, index=False)


if __name__ == "__main__":
    main()
