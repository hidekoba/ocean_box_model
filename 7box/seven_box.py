#!/usr/bin/env python3
"""
seven_box.py

Python port of the simple 7-box ocean biogeochemistry model from main(7).f90.

Boxes:
  P: polar / high-latitude surface box
  S: southern surface box
  L: low-latitude surface box
  N: northern surface box
  M: mid-depth box
  A: abyssal / lower deep box
  D: deepest box

Run:
  python3 seven_box.py
  python3 seven_box.py --save-csv seven_box_output.csv
  python3 seven_box.py --save-history seven_box_history.csv

Notes:
  - The Fortran code calculates DO2LX, but the state-update block in main(7).f90
    appears to omit `DO2L = DO2LX`. This Python version updates DO2L, treating
    the omission as a typo. Use --fortran-compat-do2l-omit to reproduce that
    exact omission.
  - The DIC equations follow main(7).f90: air-sea CO2 exchange affects atmospheric
    PCO2, but the opposite DIC flux is not applied to the surface boxes.
"""
import argparse
import math
from math import sqrt
import pandas as pd


# ---------------------------------------------------------------------
# Chemistry helper routines
# ---------------------------------------------------------------------
def chemeq_const(TEM, SAL):
    """Compute carbonate equilibrium constants from temperature [deg C]
    and salinity [PSU]. Returns (BT, K0, K1, K2, KB, KW, FH).
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
    """Oxygen saturation concentration [mol/m3]."""
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
      AT, CT: mol/m3
    Returns:
      CO2, HCO3, CO32: mol/m3
      PCO2: atm
    """
    CV2 = 9.7561e-4  # mol/m3 -> mol/kg
    ATX = CV2 * AT
    CTX = CV2 * CT

    HMIN = 1.0e-14
    HMAX = 1.0
    DELTA = 1.0e-15

    h_low = HMIN
    h_high = HMAX
    hx = 0.5 * (h_low + h_high)

    for _ in range(100000):
        h = 0.5 * (h_low + h_high)
        denom = h * h + K1 * h + K1 * K2
        at_calc = (
            (2.0 * K1 * K2 * CTX) / denom
            + (h * K1 * CTX) / denom
            + (KB * BT) / (h + KB)
            + KW / h
            - h
        )

        if abs(ATX - at_calc) <= DELTA:
            hx = h
            break

        if at_calc < ATX:
            h_high = h
        else:
            h_low = h
    else:
        hx = 0.5 * (h_low + h_high)

    denom2 = hx * hx + K1 * hx + K1 * K2
    pco2 = (hx * hx * CTX) / denom2 / K0
    co2 = (hx * hx * CTX) / denom2 / CV2
    hco3 = K1 * co2 / hx
    co32 = K2 * hco3 / hx
    return co2, hco3, co32, pco2


# ---------------------------------------------------------------------
# Parameters and initialization
# ---------------------------------------------------------------------
def default_params():
    """Parameters copied from main(7).f90."""
    CV1 = 1.0250e3
    CV2 = 9.7561e-4
    CV3 = 1.0000e6
    CV4 = 3.1536e7
    CV5 = 8.64e4

    return {
        "DT": 8.64e4,
        "DELTA": 1.0e-6,
        "CV1": CV1,
        "CV2": CV2,
        "CV3": CV3,
        "CV4": CV4,
        "CV5": CV5,

        # Temperatures are labelled K in main(7).f90 comments, but values are deg C.
        "TEMP": 1.0,
        "TEMS": 8.0,
        "TEML": 21.5,
        "TEMN": 3.0,
        "TEMM": 10.0,
        "TEMA": 3.0,
        "TEMD": 1.75,
        "SALP": 34.7,
        "SALS": 34.7,
        "SALL": 34.7,
        "SALN": 34.7,
        "SALM": 34.7,
        "SALA": 34.7,
        "SALD": 34.7,

        # Volumes/areas
        "VOCN": 1.292e18,
        "AOCN": 3.49e14,
        "VATM": 1.773e20,
        "ZOCNP": 250.0,
        "ZOCNS": 250.0,
        "ZOCNL": 100.0,
        "ZOCNN": 250.0,
        "ZOCNM": 900.0,
        "ZOCNA": 2000.0,
        "ZOCND": 4000.0,
        "DEP": 100.0,
        "DES": 100.0,
        "DEL": 100.0,
        "DEN": 100.0,
        "FAOCNP": 0.05,
        "FAOCNS": 0.10,
        "FAOCNL": 0.75,
        "FAOCNN": 0.10,
        "FAOCNM": 0.85,
        "FAOCNA": 0.95,

        # Transport/mixing [m3/s]
        "TRAN": 20.0e6,
        "FPD": 60.0e6,
        "FNA": 0.0,
        "FSM": 0.0,
        "FLM": 40.0e6,
        "FAD": 0.0,
        "FPS": 0.0,
        "FSL": 0.0,
        "FLN": 0.0,

        # Piston velocity [m/day]
        "PVP": 3.0,
        "PVS": 3.0,
        "PVL": 3.0,
        "PVN": 3.0,

        # Biology
        "R": 1.0 / CV4,
        "LF": 5.0e-1,
        "HSC": 2.0e-8 * CV1,  # mol/kg -> mol/m3, as in Fortran unit conversion
        "RRC": 0.25,
        "RCP": 162.5,
        "RNP": 15.0,
        "RO2P": 169.0,
        "GMM": 0.80,
        "GMA": 0.50,
        "CEPP": 1.00,
        "CEPS": 3.00,
        "CEPN": 3.00,

        "PCO2A": 280.0,

        # Initial tracers [mol/kg or eq/kg]
        "PO4P0": 2.09e-6,
        "PO4S0": 2.09e-6,
        "PO4L0": 2.09e-6,
        "PO4N0": 2.09e-6,
        "PO4M0": 2.09e-6,
        "PO4A0": 2.09e-6,
        "PO4D0": 2.09e-6,
        "ALKP0": 2.371e-3,
        "ALKS0": 2.371e-3,
        "ALKL0": 2.371e-3,
        "ALKN0": 2.371e-3,
        "ALKM0": 2.371e-3,
        "ALKA0": 2.371e-3,
        "ALKD0": 2.371e-3,
        "DICP0": 2.258e-3,
        "DICS0": 2.258e-3,
        "DICL0": 2.258e-3,
        "DICN0": 2.258e-3,
        "DICM0": 2.258e-3,
        "DICA0": 2.258e-3,
        "DICD0": 2.258e-3,
        "DO2P0": 1.60e-4,
        "DO2S0": 1.60e-4,
        "DO2L0": 1.60e-4,
        "DO2N0": 1.60e-4,
        "DO2M0": 1.60e-4,
        "DO2A0": 1.60e-4,
        "DO2D0": 1.60e-4,
    }


def initialize(params):
    CV1 = params["CV1"]
    CV3 = params["CV3"]
    CV4 = params["CV4"]
    CV5 = params["CV5"]

    # Areas
    AOCNP = params["AOCN"] * params["FAOCNP"]
    AOCNS = params["AOCN"] * params["FAOCNS"]
    AOCNL = params["AOCN"] * params["FAOCNL"]
    AOCNN = params["AOCN"] * params["FAOCNN"]
    AOCNM = params["AOCN"] * params["FAOCNM"]
    AOCNA = params["AOCN"] * params["FAOCNA"]

    # Volumes
    VOCNP = AOCNP * params["ZOCNP"]
    VOCNS = AOCNS * params["ZOCNS"]
    VOCNL = AOCNL * params["ZOCNL"]
    VOCNN = AOCNN * params["ZOCNN"]
    VOCNM = AOCNM * params["ZOCNM"] - VOCNS - VOCNL
    VOCNA = AOCNA * params["ZOCNA"] - VOCNS - VOCNL - VOCNN - VOCNM
    VOCND = params["VOCN"] - VOCNP - VOCNS - VOCNL - VOCNN - VOCNM - VOCNA

    # Biological sinking fluxes [mol P/s]
    EPP = (params["CEPP"] / params["RCP"]) * AOCNP / CV4
    EPS = (params["CEPS"] / params["RCP"]) * AOCNS / CV4
    EPN = (params["CEPN"] / params["RCP"]) * AOCNN / CV4

    # Gas exchange [m3/s]
    FAP = params["PVP"] * AOCNP / CV5
    FAS = params["PVS"] * AOCNS / CV5
    FAL = params["PVL"] * AOCNL / CV5
    FAN = params["PVN"] * AOCNN / CV5

    state = {}
    for tracer in ["PO4", "ALK", "DIC", "DO2"]:
        for box in BOXES:
            state[f"{tracer}{box}"] = params[f"{tracer}{box}0"] * CV1
    state["PCO2A"] = params["PCO2A"] / CV3

    aux = {
        "AOCNP": AOCNP, "AOCNS": AOCNS, "AOCNL": AOCNL, "AOCNN": AOCNN, "AOCNM": AOCNM, "AOCNA": AOCNA,
        "VOCNP": VOCNP, "VOCNS": VOCNS, "VOCNL": VOCNL, "VOCNN": VOCNN, "VOCNM": VOCNM, "VOCNA": VOCNA, "VOCND": VOCND,
        "EPP": EPP, "EPS": EPS, "EPN": EPN, "EPL": None,
        "FAP": FAP, "FAS": FAS, "FAL": FAL, "FAN": FAN,
    }
    return state, aux


# ---------------------------------------------------------------------
# Model core
# ---------------------------------------------------------------------
BOXES = ["P", "S", "L", "N", "M", "A", "D"]
SURFACE_BOXES = ["P", "S", "L", "N"]


def carbonate_all(state, params):
    constants = {}
    outputs = {}
    specs = {
        "P": ("TEMP", "SALP", "ALKP", "DICP"),
        "S": ("TEMS", "SALS", "ALKS", "DICS"),
        "L": ("TEML", "SALL", "ALKL", "DICL"),
        "N": ("TEMN", "SALN", "ALKN", "DICN"),
        "M": ("TEMM", "SALM", "ALKM", "DICM"),
        "A": ("TEMA", "SALA", "ALKA", "DICA"),
        "D": ("TEMD", "SALD", "ALKD", "DICD"),
    }
    for box, (temp_key, sal_key, alk_key, dic_key) in specs.items():
        BT, K0, K1, K2, KB, KW, FH = chemeq_const(params[temp_key], params[sal_key])
        constants[box] = {"BT": BT, "K0": K0, "K1": K1, "K2": K2, "KB": KB, "KW": KW, "FH": FH}
        CO2, HCO3, CO32, PCO2 = co2_nibun(BT, K0, K1, K2, KB, KW, state[alk_key], state[dic_key])
        outputs[box] = {"CO2": CO2, "HCO3": HCO3, "CO32": CO32, "PCO2": PCO2}
    return constants, outputs


def o2sats_all(params):
    return {
        "P": o2sat(params["TEMP"], params["SALP"]),
        "S": o2sat(params["TEMS"], params["SALS"]),
        "L": o2sat(params["TEML"], params["SALL"]),
        "N": o2sat(params["TEMN"], params["SALN"]),
        "M": o2sat(params["TEMM"], params["SALM"]),
        "A": o2sat(params["TEMA"], params["SALA"]),
        "D": o2sat(params["TEMD"], params["SALD"]),
    }


def step(state, params, aux, *, fortran_compat_do2l_omit=False):
    DT = params["DT"]
    TRAN = params["TRAN"]
    FPD = params["FPD"]
    FNA = params["FNA"]
    FSM = params["FSM"]
    FLM = params["FLM"]
    FAD = params["FAD"]
    FPS = params["FPS"]
    FSL = params["FSL"]
    FLN = params["FLN"]
    RRC = params["RRC"]
    RCP = params["RCP"]
    RNP = params["RNP"]
    RO2P = params["RO2P"]
    GMM = params["GMM"]
    GMA = params["GMA"]
    RGMM = 1.0 - GMM
    RGMA = 1.0 - GMA

    EPP = aux["EPP"]
    EPS = aux["EPS"]
    EPN = aux["EPN"]
    FAP = aux["FAP"]
    FAS = aux["FAS"]
    FAL = aux["FAL"]
    FAN = aux["FAN"]
    VOCNP = aux["VOCNP"]
    VOCNS = aux["VOCNS"]
    VOCNL = aux["VOCNL"]
    VOCNN = aux["VOCNN"]
    VOCNM = aux["VOCNM"]
    VOCNA = aux["VOCNA"]
    VOCND = aux["VOCND"]

    constants, carb = carbonate_all(state, params)
    o2satv = o2sats_all(params)

    # Low-latitude export production [mol P/s], updated every time step.
    EPL = (
        params["R"]
        * params["DEL"]
        * params["LF"]
        * state["PO4L"]
        * (state["PO4L"] / (params["HSC"] + state["PO4L"]))
        * VOCNL
    )

    # PO4
    PO4PX = state["PO4P"] + ((TRAN + FPD) * (state["PO4D"] - state["PO4P"]) + FPS * (state["PO4S"] - state["PO4P"]) - EPP) * (DT / VOCNP)
    PO4SX = state["PO4S"] + ((TRAN + FPS) * (state["PO4P"] - state["PO4S"]) + FSL * (state["PO4L"] - state["PO4S"]) - EPS) * (DT / VOCNS)
    PO4LX = state["PO4L"] + (FSL * (state["PO4S"] - state["PO4L"]) + FLN * (state["PO4N"] - state["PO4L"]) + FLM * (state["PO4M"] - state["PO4L"]) - EPL) * (DT / VOCNL)
    PO4NX = state["PO4N"] + (TRAN * (state["PO4M"] - state["PO4N"]) + FLN * (state["PO4L"] - state["PO4N"]) + FNA * (state["PO4A"] - state["PO4N"]) - EPN) * (DT / VOCNN)
    PO4MX = state["PO4M"] + ((TRAN + FSM) * (state["PO4S"] - state["PO4M"]) + FLM * (state["PO4L"] - state["PO4M"]) + GMM * (EPS + EPL)) * (DT / VOCNM)
    PO4AX = state["PO4A"] + ((TRAN + FNA) * (state["PO4N"] - state["PO4A"]) + FAD * (state["PO4D"] - state["PO4A"]) + (0.5 * RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
    PO4DX = state["PO4D"] + ((TRAN + FAD) * (state["PO4A"] - state["PO4D"]) + FPD * (state["PO4P"] - state["PO4D"]) + (0.5 * RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

    alk_factor = 2.0 * RRC * RCP - RNP
    dic_factor = (1.0 + RRC) * RCP

    # ALK
    ALKPX = state["ALKP"] + ((TRAN + FPD) * (state["ALKD"] - state["ALKP"]) + FPS * (state["ALKS"] - state["ALKP"]) - alk_factor * EPP) * (DT / VOCNP)
    ALKSX = state["ALKS"] + ((TRAN + FPS) * (state["ALKP"] - state["ALKS"]) + FSL * (state["ALKL"] - state["ALKS"]) - alk_factor * EPS) * (DT / VOCNS)
    ALKLX = state["ALKL"] + (FSL * (state["ALKS"] - state["ALKL"]) + FLN * (state["ALKN"] - state["ALKL"]) + FLM * (state["ALKM"] - state["ALKL"]) - alk_factor * EPL) * (DT / VOCNL)
    ALKNX = state["ALKN"] + (TRAN * (state["ALKM"] - state["ALKN"]) + FLN * (state["ALKL"] - state["ALKN"]) + FNA * (state["ALKA"] - state["ALKN"]) - alk_factor * EPN) * (DT / VOCNN)
    ALKMX = state["ALKM"] + ((TRAN + FSM) * (state["ALKS"] - state["ALKM"]) + FLM * (state["ALKL"] - state["ALKM"]) + alk_factor * (GMM * (EPS + EPL))) * (DT / VOCNM)
    ALKAX = state["ALKA"] + ((TRAN + FNA) * (state["ALKN"] - state["ALKA"]) + FAD * (state["ALKD"] - state["ALKA"]) + alk_factor * (0.5 * RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
    ALKDX = state["ALKD"] + ((TRAN + FAD) * (state["ALKA"] - state["ALKD"]) + FPD * (state["ALKP"] - state["ALKD"]) + alk_factor * (0.5 * RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

    # DIC
    DICPX = state["DICP"] + ((TRAN + FPD) * (state["DICD"] - state["DICP"]) + FPS * (state["DICS"] - state["DICP"]) - dic_factor * EPP) * (DT / VOCNP)
    DICSX = state["DICS"] + ((TRAN + FPS) * (state["DICP"] - state["DICS"]) + FSL * (state["DICL"] - state["DICS"]) - dic_factor * EPS) * (DT / VOCNS)
    DICLX = state["DICL"] + (FSL * (state["DICS"] - state["DICL"]) + FLN * (state["DICN"] - state["DICL"]) + FLM * (state["DICM"] - state["DICL"]) - dic_factor * EPL) * (DT / VOCNL)
    DICNX = state["DICN"] + (TRAN * (state["DICM"] - state["DICN"]) + FLN * (state["DICL"] - state["DICN"]) + FNA * (state["DICA"] - state["DICN"]) - dic_factor * EPN) * (DT / VOCNN)
    DICMX = state["DICM"] + ((TRAN + FSM) * (state["DICS"] - state["DICM"]) + FLM * (state["DICL"] - state["DICM"]) + dic_factor * (GMM * (EPS + EPL))) * (DT / VOCNM)
    DICAX = state["DICA"] + ((TRAN + FNA) * (state["DICN"] - state["DICA"]) + FAD * (state["DICD"] - state["DICA"]) + dic_factor * (0.5 * RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
    DICDX = state["DICD"] + ((TRAN + FAD) * (state["DICA"] - state["DICD"]) + FPD * (state["DICP"] - state["DICD"]) + dic_factor * (0.5 * RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

    # O2: follows main(7).f90 equations; no air-sea O2 exchange terms are included there.
    DO2PX = state["DO2P"] + ((TRAN + FPD) * (state["DO2D"] - state["DO2P"]) + FPS * (state["DO2S"] - state["DO2P"]) + RO2P * EPP + FAP * (o2satv["P"] - state["DO2P"])) * (DT / VOCNP)
    DO2SX = state["DO2S"] + ((TRAN + FPS) * (state["DO2P"] - state["DO2S"]) + FSL * (state["DO2L"] - state["DO2S"]) + RO2P * EPS + FAS * (o2satv["S"] - state["DO2S"])) * (DT / VOCNS)
    DO2LX = state["DO2L"] + (FSL * (state["DO2S"] - state["DO2L"]) + FLN * (state["DO2N"] - state["DO2L"]) + FLM * (state["DO2M"] - state["DO2L"]) + RO2P * EPL + FAL * (o2satv["L"] - state["DO2L"])) * (DT / VOCNL)
    DO2NX = state["DO2N"] + (TRAN * (state["DO2M"] - state["DO2N"]) + FLN * (state["DO2L"] - state["DO2N"]) + FNA * (state["DO2A"] - state["DO2N"]) + RO2P * EPN + FAN * (o2satv["N"] - state["DO2N"])) * (DT / VOCNN)
    DO2MX = state["DO2M"] + ((TRAN + FSM) * (state["DO2S"] - state["DO2M"]) + FLM * (state["DO2L"] - state["DO2M"]) - RO2P * (GMM * (EPS + EPL))) * (DT / VOCNM)
    DO2AX = state["DO2A"] + ((TRAN + FNA) * (state["DO2N"] - state["DO2A"]) + FAD * (state["DO2D"] - state["DO2A"]) - RO2P * (0.5 * RGMM * (EPS + EPL) + GMA * EPN)) * (DT / VOCNA)
    DO2DX = state["DO2D"] + ((TRAN + FAD) * (state["DO2A"] - state["DO2D"]) + FPD * (state["DO2P"] - state["DO2D"]) - RO2P * (0.5 * RGMM * (EPS + EPL) + RGMA * EPN + EPP)) * (DT / VOCND)

    # Carbonate system after tracer update
    tmp_state = dict(state)
    tmp_state.update({
        "ALKP": ALKPX, "ALKS": ALKSX, "ALKL": ALKLX, "ALKN": ALKNX, "ALKM": ALKMX, "ALKA": ALKAX, "ALKD": ALKDX,
        "DICP": DICPX, "DICS": DICSX, "DICL": DICLX, "DICN": DICNX, "DICM": DICMX, "DICA": DICAX, "DICD": DICDX,
    })
    constants_x, carb_x = carbonate_all(tmp_state, params)

    PCO2AX = state["PCO2A"] + (
        aux["FAP"] * params["CV1"] * constants_x["P"]["K0"] * (carb_x["P"]["PCO2"] - state["PCO2A"])
        + aux["FAS"] * params["CV1"] * constants_x["S"]["K0"] * (carb_x["S"]["PCO2"] - state["PCO2A"])
        + aux["FAL"] * params["CV1"] * constants_x["L"]["K0"] * (carb_x["L"]["PCO2"] - state["PCO2A"])
        + aux["FAN"] * params["CV1"] * constants_x["N"]["K0"] * (carb_x["N"]["PCO2"] - state["PCO2A"])
    ) * (DT / params["VATM"])

    new_state = {
        "PO4P": PO4PX, "PO4S": PO4SX, "PO4L": PO4LX, "PO4N": PO4NX, "PO4M": PO4MX, "PO4A": PO4AX, "PO4D": PO4DX,
        "ALKP": ALKPX, "ALKS": ALKSX, "ALKL": ALKLX, "ALKN": ALKNX, "ALKM": ALKMX, "ALKA": ALKAX, "ALKD": ALKDX,
        "DICP": DICPX, "DICS": DICSX, "DICL": DICLX, "DICN": DICNX, "DICM": DICMX, "DICA": DICAX, "DICD": DICDX,
        "DO2P": DO2PX, "DO2S": DO2SX, "DO2L": DO2LX, "DO2N": DO2NX, "DO2M": DO2MX, "DO2A": DO2AX, "DO2D": DO2DX,
        "PCO2A": PCO2AX,
    }

    if fortran_compat_do2l_omit:
        # Reproduce the apparent omission in main(7).f90's update block.
        new_state["DO2L"] = state["DO2L"]

    aux = dict(aux)
    aux["EPL"] = EPL
    return new_state, aux


# ---------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------
def total_carbon_mol(state, aux, params):
    return (
        state["DICP"] * aux["VOCNP"]
        + state["DICS"] * aux["VOCNS"]
        + state["DICL"] * aux["VOCNL"]
        + state["DICN"] * aux["VOCNN"]
        + state["DICM"] * aux["VOCNM"]
        + state["DICA"] * aux["VOCNA"]
        + state["DICD"] * aux["VOCND"]
        + state["PCO2A"] * params["VATM"]
    )


def make_output(state, params, aux, tocini=None, nstep=None):
    CV2 = params["CV2"]
    CV3 = params["CV3"]
    CV4 = params["CV4"]
    constants, carb = carbonate_all(state, params)
    o2satv = o2sats_all(params)

    to_umolkg = lambda x: CV2 * 1.0e6 * x
    to_ppmv = lambda x: CV3 * x
    tocfin = total_carbon_mol(state, aux, params)

    out = {"TIMESTEP": nstep}

    for box, temp_key, sal_key in [
        ("P", "TEMP", "SALP"), ("S", "TEMS", "SALS"), ("L", "TEML", "SALL"),
        ("N", "TEMN", "SALN"), ("M", "TEMM", "SALM"), ("A", "TEMA", "SALA"), ("D", "TEMD", "SALD"),
    ]:
        out[f"Temperature_{box}_degC"] = params[temp_key]
        out[f"Salinity_{box}_PSU"] = params[sal_key]

    out["EP_P_PgC_yr"] = CV4 * aux["EPP"] * params["RCP"] * 1.0e-15 * 12.0
    out["EP_S_PgC_yr"] = CV4 * aux["EPS"] * params["RCP"] * 1.0e-15 * 12.0
    out["EP_L_PgC_yr"] = CV4 * (aux["EPL"] if aux["EPL"] is not None else 0.0) * params["RCP"] * 1.0e-15 * 12.0
    out["EP_N_PgC_yr"] = CV4 * aux["EPN"] * params["RCP"] * 1.0e-15 * 12.0

    for box in BOXES:
        out[f"PO4_{box}_umolkg"] = to_umolkg(state[f"PO4{box}"])
        out[f"DIC_{box}_umolkg"] = to_umolkg(state[f"DIC{box}"])
        out[f"ALK_{box}_ueqkg"] = to_umolkg(state[f"ALK{box}"])
        out[f"DO2_{box}_umolkg"] = to_umolkg(state[f"DO2{box}"])
        out[f"AOU_{box}_umolkg"] = to_umolkg(o2satv[box] - state[f"DO2{box}"])
        out[f"CO2_{box}_umolkg"] = to_umolkg(carb[box]["CO2"])
        out[f"HCO3_{box}_umolkg"] = to_umolkg(carb[box]["HCO3"])
        out[f"CO32_{box}_umolkg"] = to_umolkg(carb[box]["CO32"])
        out[f"PCO2_{box}_ppmv"] = to_ppmv(carb[box]["PCO2"])

    out["PCO2_ATM_ppmv"] = to_ppmv(state["PCO2A"])
    out["Total_Carbon_Initial_PgC"] = tocini * 12.0e-15 if tocini is not None else None
    out["Total_Carbon_Final_PgC"] = tocfin * 12.0e-15
    return out


def create_dataframe(out):
    rows = []
    for box in BOXES:
        rows.append(("Temperature (deg C)", box, out.get(f"Temperature_{box}_degC")))
    for box in BOXES:
        rows.append(("Salinity (PSU)", box, out.get(f"Salinity_{box}_PSU")))
    for box in SURFACE_BOXES:
        rows.append(("Export production (PgC/yr)", box, out.get(f"EP_{box}_PgC_yr")))
    for q, key_suffix in [
        ("PO4 (umol/kg)", "PO4_{box}_umolkg"),
        ("DIC (umol/kg)", "DIC_{box}_umolkg"),
        ("ALK (ueq/kg)", "ALK_{box}_ueqkg"),
        ("O2 (umol/kg)", "DO2_{box}_umolkg"),
        ("AOU (umol/kg)", "AOU_{box}_umolkg"),
        ("CO2 (umol/kg)", "CO2_{box}_umolkg"),
        ("HCO3 (umol/kg)", "HCO3_{box}_umolkg"),
        ("CO3 (umol/kg)", "CO32_{box}_umolkg"),
        ("PCO2 (ppmv)", "PCO2_{box}_ppmv"),
    ]:
        for box in BOXES:
            rows.append((q, box, out.get(key_suffix.format(box=box))))
    rows.append(("PCO2 (ppmv)", "Atmosphere", out.get("PCO2_ATM_ppmv")))
    rows.append(("Total Carbon (PgC)", "Initial", out.get("Total_Carbon_Initial_PgC")))
    rows.append(("Total Carbon (PgC)", "Final", out.get("Total_Carbon_Final_PgC")))

    df = pd.DataFrame(rows, columns=["Quantity", "Box", "Value"])
    df["Value"] = df["Value"].map(lambda x: f"{x:.3f}" if isinstance(x, (int, float)) else x)
    return df


def run_model(max_iter=100000, verbose=True, params=None, return_history=False, *, fortran_compat_do2l_omit=False):
    if params is None:
        params = default_params()
    state, aux = initialize(params)
    tocini = total_carbon_mol(state, aux, params)

    history = []
    nstep = max_iter
    for i in range(1, max_iter + 1):
        old_po4p = state["PO4P"]
        state_new, aux_new = step(state, params, aux, fortran_compat_do2l_omit=fortran_compat_do2l_omit)
        rel = abs((state_new["PO4P"] / old_po4p) - 1.0) if old_po4p != 0 else abs(state_new["PO4P"])

        state, aux = state_new, aux_new

        if return_history and (i == 1 or i % 100 == 0):
            hist_out = make_output(state, params, aux, tocini=tocini, nstep=i)
            history.append({
                "step": i,
                "PCO2_ATM_ppmv": hist_out["PCO2_ATM_ppmv"],
                **{f"PO4_{box}_umolkg": hist_out[f"PO4_{box}_umolkg"] for box in BOXES},
            })

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
    parser.add_argument(
        "--fortran-compat-do2l-omit",
        action="store_true",
        help="Reproduce the apparent omission of DO2L update in main(7).f90.",
    )
    args = parser.parse_args()

    if args.save_history:
        out, hist = run_model(
            max_iter=args.iters,
            verbose=not args.quiet,
            return_history=True,
            fortran_compat_do2l_omit=args.fortran_compat_do2l_omit,
        )
        hist.to_csv(args.save_history, index=False)
    else:
        out = run_model(
            max_iter=args.iters,
            verbose=not args.quiet,
            fortran_compat_do2l_omit=args.fortran_compat_do2l_omit,
        )

    if args.save_csv:
        create_dataframe(out).to_csv(args.save_csv, index=False)


if __name__ == "__main__":
    main()
