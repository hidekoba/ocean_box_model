#!/usr/bin/env python3
"""
four_box.py

Python port of the simple 4-box ocean biogeochemistry model

Boxes:
  H: high-latitude surface box
  L: low-latitude surface box
  N: northern / intermediate surface-like box
  D: deep box

Run:
  python3 four_box.py
  python3 four_box.py --save-csv output.csv
  python3 four_box.py --save-history history.csv
"""
import argparse
import math
from math import sqrt
import pandas as pd


def chemeq_const(TEM, SAL):
    """Compute carbonate equilibrium constants from temperature [deg C]
    and salinity [PSU]. Ported from the existing two/three-box helpers.
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
    """Oxygen saturation concentration [mol/m^3]."""
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


def default_params():
    """Parameters"""
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

        # Temperatures are numerically in deg C, although the Fortran output labels K.
        "TEMH": 1.0,
        "TEML": 21.5,
        "TEMN": 3.0,
        "TEMD": 1.75,
        "SALH": 34.7,
        "SALL": 34.7,
        "SALN": 34.7,
        "SALD": 34.7,

        # Volumes / areas
        "VOCN": 1.292e18,
        "AOCN": 3.49e14,
        "VATM": 1.773e20,
        "ZOCNH": 250.0,
        "ZOCNL": 100.0,
        "ZOCNN": 250.0,
        "ZOCND": 4000.0,
        "DEH": 100.0,
        "DEL": 100.0,
        "DEN": 100.0,
        "FAOCNH": 0.15,
        "FAOCNL": 0.75,
        "FAOCNN": 0.10,

        # Transport / mixing [m^3/s]
        "TRAN": 20.0e6,
        "FHD": 60.0e6,
        "FLD": 0.0,
        "FND": 0.0,
        "FHL": 0.0,
        "FLN": 0.0,

        # Piston velocity [m/day]
        "PVH": 3.0,
        "PVL": 3.0,
        "PVN": 3.0,

        # Biological parameters
        "R": 1.0 / CV4,       # [/s]
        "LF": 5.0e-1,
        "HSC": 2.0e-8 * CV1,  # [mol/m^3]
        "RRC": 0.250,
        "RCP": 162.5,
        "RNP": 15.0,
        "RO2P": 169.0,
        "CEPH": 1.00,         # [mol C m^-2 yr^-1]

        # Initial atmospheric PCO2 [ppmv]
        "PCO2A": 280.0,

        # Initial tracers [mol/kg or eq/kg]
        "PO4H0": 2.09e-6,
        "PO4L0": 2.09e-6,
        "PO4N0": 2.09e-6,
        "PO4D0": 2.09e-6,
        "ALKH0": 2.371e-3,
        "ALKL0": 2.371e-3,
        "ALKN0": 2.371e-3,
        "ALKD0": 2.371e-3,
        "DICH0": 2.258e-3,
        "DICL0": 2.258e-3,
        "DICN0": 2.258e-3,
        "DICD0": 2.258e-3,
        "DO2H0": 1.60e-4,
        "DO2L0": 1.60e-4,
        "DO2N0": 1.60e-4,
        "DO2D0": 1.60e-4,
    }


def initialize(params):
    """Initialize volumes, fluxes, and state variables in model units."""
    CV1 = params["CV1"]
    CV3 = params["CV3"]
    CV4 = params["CV4"]
    CV5 = params["CV5"]

    # Areas
    AOCNH = params["AOCN"] * params["FAOCNH"]
    AOCNL = params["AOCN"] * params["FAOCNL"]
    AOCNN = params["AOCN"] * params["FAOCNN"]

    # Volumes
    VOCNH = AOCNH * params["ZOCNH"]
    VOCNL = AOCNL * params["ZOCNL"]
    VOCNN = AOCNN * params["ZOCNN"]
    VOCND = params["VOCN"] - VOCNH - VOCNL - VOCNN

    # High-lat export [mol P/s]; low-lat and N exports are diagnosed each step.
    EPH = (params["CEPH"] / params["RCP"]) * AOCNH / CV4

    # Gas exchange [m^3/s]
    FAH = params["PVH"] * AOCNH / CV5
    FAL = params["PVL"] * AOCNL / CV5
    FAN = params["PVN"] * AOCNN / CV5

    state = {
        "PO4H": params["PO4H0"] * CV1,
        "PO4L": params["PO4L0"] * CV1,
        "PO4N": params["PO4N0"] * CV1,
        "PO4D": params["PO4D0"] * CV1,
        "ALKH": params["ALKH0"] * CV1,
        "ALKL": params["ALKL0"] * CV1,
        "ALKN": params["ALKN0"] * CV1,
        "ALKD": params["ALKD0"] * CV1,
        "DICH": params["DICH0"] * CV1,
        "DICL": params["DICL0"] * CV1,
        "DICN": params["DICN0"] * CV1,
        "DICD": params["DICD0"] * CV1,
        "DO2H": params["DO2H0"] * CV1,
        "DO2L": params["DO2L0"] * CV1,
        "DO2N": params["DO2N0"] * CV1,
        "DO2D": params["DO2D0"] * CV1,
        "PCO2A": params["PCO2A"] / CV3,
    }

    aux = {
        "AOCNH": AOCNH,
        "AOCNL": AOCNL,
        "AOCNN": AOCNN,
        "VOCNH": VOCNH,
        "VOCNL": VOCNL,
        "VOCNN": VOCNN,
        "VOCND": VOCND,
        "EPH": EPH,
        "EPL": None,
        "EPN": None,
        "FAH": FAH,
        "FAL": FAL,
        "FAN": FAN,
    }
    return state, aux


BOXES = {
    "H": {"temp": "TEMH", "sal": "SALH", "alk": "ALKH", "dic": "DICH"},
    "L": {"temp": "TEML", "sal": "SALL", "alk": "ALKL", "dic": "DICL"},
    "N": {"temp": "TEMN", "sal": "SALN", "alk": "ALKN", "dic": "DICN"},
    "D": {"temp": "TEMD", "sal": "SALD", "alk": "ALKD", "dic": "DICD"},
}


def carbonate_all(state, params):
    """Compute carbonate system for all boxes."""
    constants = {}
    outputs = {}
    for label, keys in BOXES.items():
        BT, K0, K1, K2, KB, KW, FH = chemeq_const(params[keys["temp"]], params[keys["sal"]])
        constants[label] = {"BT": BT, "K0": K0, "K1": K1, "K2": K2, "KB": KB, "KW": KW, "FH": FH}
        CO2, HCO3, CO32, PCO2 = co2_nibun(
            BT, K0, K1, K2, KB, KW, state[keys["alk"]], state[keys["dic"]]
        )
        outputs[label] = {"CO2": CO2, "HCO3": HCO3, "CO32": CO32, "PCO2": PCO2}
    return constants, outputs


def oxygen_saturation_all(params):
    return {
        "H": o2sat(params["TEMH"], params["SALH"]),
        "L": o2sat(params["TEML"], params["SALL"]),
        "N": o2sat(params["TEMN"], params["SALN"]),
        "D": o2sat(params["TEMD"], params["SALD"]),
    }


def step(state, params, aux):
    """One forward-Euler step using main(5).f90 equations."""
    DT = params["DT"]
    TRAN = params["TRAN"]
    FHD = params["FHD"]
    FLD = params["FLD"]
    FND = params["FND"]
    FHL = params["FHL"]
    FLN = params["FLN"]
    RRC = params["RRC"]
    RCP = params["RCP"]
    RNP = params["RNP"]
    RO2P = params["RO2P"]

    VOCNH = aux["VOCNH"]
    VOCNL = aux["VOCNL"]
    VOCNN = aux["VOCNN"]
    VOCND = aux["VOCND"]
    EPH = aux["EPH"]
    FAH = aux["FAH"]
    FAL = aux["FAL"]
    FAN = aux["FAN"]

    constants, carb = carbonate_all(state, params)
    o2sat_all = oxygen_saturation_all(params)

    EPL = (
        params["R"]
        * params["DEL"]
        * params["LF"]
        * state["PO4L"]
        * (state["PO4L"] / (params["HSC"] + state["PO4L"]))
        * VOCNL
    )
    EPN = (
        params["R"]
        * params["DEN"]
        * params["LF"]
        * state["PO4N"]
        * (state["PO4N"] / (params["HSC"] + state["PO4N"]))
        * VOCNN
    )

    # PO4
    PO4HX = state["PO4H"] + (
        (TRAN + FHD) * (state["PO4D"] - state["PO4H"])
        + FHL * (state["PO4L"] - state["PO4H"])
        - EPH
    ) * (DT / VOCNH)
    PO4LX = state["PO4L"] + (
        (TRAN + FHL) * (state["PO4H"] - state["PO4L"])
        + FLN * (state["PO4N"] - state["PO4L"])
        + FLD * (state["PO4D"] - state["PO4L"])
        - EPL
    ) * (DT / VOCNL)
    PO4NX = state["PO4N"] + (
        (TRAN + FLN) * (state["PO4L"] - state["PO4N"])
        + FND * (state["PO4D"] - state["PO4N"])
        - EPN
    ) * (DT / VOCNN)
    PO4DX = state["PO4D"] + (
        (TRAN + FND) * (state["PO4N"] - state["PO4D"])
        + FHD * (state["PO4H"] - state["PO4D"])
        + FLD * (state["PO4L"] - state["PO4D"])
        + (EPH + EPL + EPN)
    ) * (DT / VOCND)

    alk_factor = 2.0 * RRC * RCP - RNP
    dic_factor = (1.0 + RRC) * RCP

    # ALK
    ALKHX = state["ALKH"] + (
        (TRAN + FHD) * (state["ALKD"] - state["ALKH"])
        + FHL * (state["ALKL"] - state["ALKH"])
        - alk_factor * EPH
    ) * (DT / VOCNH)
    ALKLX = state["ALKL"] + (
        (TRAN + FHL) * (state["ALKH"] - state["ALKL"])
        + FLN * (state["ALKN"] - state["ALKL"])
        + FLD * (state["ALKD"] - state["ALKL"])
        - alk_factor * EPL
    ) * (DT / VOCNL)
    ALKNX = state["ALKN"] + (
        (TRAN + FLN) * (state["ALKL"] - state["ALKN"])
        + FND * (state["ALKD"] - state["ALKN"])
        - alk_factor * EPN
    ) * (DT / VOCNN)
    ALKDX = state["ALKD"] + (
        (TRAN + FND) * (state["ALKN"] - state["ALKD"])
        + FHD * (state["ALKH"] - state["ALKD"])
        + FLD * (state["ALKL"] - state["ALKD"])
        + alk_factor * (EPH + EPL + EPN)
    ) * (DT / VOCND)

    # DIC
    DICHX = state["DICH"] + (
        (TRAN + FHD) * (state["DICD"] - state["DICH"])
        + FHL * (state["DICL"] - state["DICH"])
        - dic_factor * EPH
        + FAH * params["CV1"] * constants["H"]["K0"] * (state["PCO2A"] - carb["H"]["PCO2"])
    ) * (DT / VOCNH)
    DICLX = state["DICL"] + (
        (TRAN + FHL) * (state["DICH"] - state["DICL"])
        + FLN * (state["DICN"] - state["DICL"])
        + FLD * (state["DICD"] - state["DICL"])
        - dic_factor * EPL
        + FAL * params["CV1"] * constants["L"]["K0"] * (state["PCO2A"] - carb["L"]["PCO2"])
    ) * (DT / VOCNL)
    DICNX = state["DICN"] + (
        (TRAN + FLN) * (state["DICL"] - state["DICN"])
        + FND * (state["DICD"] - state["DICN"])
        - dic_factor * EPN
        + FAN * params["CV1"] * constants["N"]["K0"] * (state["PCO2A"] - carb["N"]["PCO2"])
    ) * (DT / VOCNN)
    DICDX = state["DICD"] + (
        (TRAN + FND) * (state["DICN"] - state["DICD"])
        + FHD * (state["DICH"] - state["DICD"])
        + FLD * (state["DICL"] - state["DICD"])
        + dic_factor * (EPH + EPL + EPN)
    ) * (DT / VOCND)

    # O2
    DO2HX = state["DO2H"] + (
        (TRAN + FHD) * (state["DO2D"] - state["DO2H"])
        + FHL * (o2sat_all["L"] - state["DO2H"])
        + RO2P * EPH
        + FAH * (o2sat_all["H"] - state["DO2H"])
    ) * (DT / VOCNH)
    DO2LX = state["DO2L"] + (
        (TRAN + FHL) * (state["DO2H"] - state["DO2L"])
        + FLN * (state["DO2N"] - state["DO2L"])
        + FLD * (state["DO2D"] - state["DO2L"])
        + RO2P * EPL
        + FAL * (o2sat_all["L"] - state["DO2L"])
    ) * (DT / VOCNL)
    DO2NX = state["DO2N"] + (
        (TRAN + FLN) * (state["DO2L"] - state["DO2N"])
        + FND * (state["DO2D"] - state["DO2N"])
        + RO2P * EPN
        + FAN * (o2sat_all["N"] - state["DO2N"])
    ) * (DT / VOCNN)
    DO2DX = state["DO2D"] + (
        (TRAN + FND) * (state["DO2N"] - state["DO2D"])
        + FHD * (state["DO2H"] - state["DO2D"])
        + FLD * (state["DO2L"] - state["DO2D"])
        - RO2P * (EPH + EPL + EPN)
    ) * (DT / VOCND)

    tmp_state = dict(state)
    tmp_state.update({
        "ALKH": ALKHX, "ALKL": ALKLX, "ALKN": ALKNX, "ALKD": ALKDX,
        "DICH": DICHX, "DICL": DICLX, "DICN": DICNX, "DICD": DICDX,
    })
    constants_x, carb_x = carbonate_all(tmp_state, params)

    PCO2AX = state["PCO2A"] + (
        aux["FAH"] * params["CV1"] * constants_x["H"]["K0"] * (carb_x["H"]["PCO2"] - state["PCO2A"])
        + aux["FAL"] * params["CV1"] * constants_x["L"]["K0"] * (carb_x["L"]["PCO2"] - state["PCO2A"])
        + aux["FAN"] * params["CV1"] * constants_x["N"]["K0"] * (carb_x["N"]["PCO2"] - state["PCO2A"])
    ) * (DT / params["VATM"])

    new_state = {
        "PO4H": PO4HX, "PO4L": PO4LX, "PO4N": PO4NX, "PO4D": PO4DX,
        "ALKH": ALKHX, "ALKL": ALKLX, "ALKN": ALKNX, "ALKD": ALKDX,
        "DICH": DICHX, "DICL": DICLX, "DICN": DICNX, "DICD": DICDX,
        "DO2H": DO2HX, "DO2L": DO2LX, "DO2N": DO2NX, "DO2D": DO2DX,
        "PCO2A": PCO2AX,
    }

    aux = dict(aux)
    aux["EPL"] = EPL
    aux["EPN"] = EPN
    return new_state, aux


def total_carbon(state, params, aux):
    return (
        state["DICH"] * aux["VOCNH"]
        + state["DICL"] * aux["VOCNL"]
        + state["DICN"] * aux["VOCNN"]
        + state["DICD"] * aux["VOCND"]
        + state["PCO2A"] * params["VATM"]
    )


def make_output(state, params, aux, tocini=None, nstep=None):
    CV2 = params["CV2"]
    CV3 = params["CV3"]
    CV4 = params["CV4"]

    constants, carb = carbonate_all(state, params)
    o2sat_all = oxygen_saturation_all(params)

    to_umolkg = lambda x: CV2 * 1.0e6 * x
    to_ppmv = lambda x: CV3 * x

    tocfin = total_carbon(state, params, aux)
    tocini = total_carbon(state, params, aux) if tocini is None else tocini

    out = {
        "TIMESTEP": nstep,
    }

    for box in ["H", "L", "N", "D"]:
        out[f"Temperature_{box}_degC"] = params[f"TEM{box}"]
        out[f"Salinity_{box}_PSU"] = params[f"SAL{box}"]

    out["EP_H_PgC_yr"] = CV4 * aux["EPH"] * params["RCP"] * 1.0e-15 * 12.0
    out["EP_L_PgC_yr"] = CV4 * (aux["EPL"] if aux["EPL"] is not None else 0.0) * params["RCP"] * 1.0e-15 * 12.0
    out["EP_N_PgC_yr"] = CV4 * (aux["EPN"] if aux["EPN"] is not None else 0.0) * params["RCP"] * 1.0e-15 * 12.0

    tracer_map = {
        "PO4": "PO4",
        "DIC": "DIC",
        "ALK": "ALK",
        "DO2": "DO2",
    }
    for box in ["H", "L", "N", "D"]:
        for short, prefix in tracer_map.items():
            key = f"{prefix}{box}"
            label = "DIC" if short == "DIC" else ("ALK" if short == "ALK" else short)
            unit_suffix = "ueqkg" if short == "ALK" else "umolkg"
            out[f"{label}_{box}_{unit_suffix}"] = to_umolkg(state[key])

        out[f"AOU_{box}_umolkg"] = to_umolkg(o2sat_all[box] - state[f"DO2{box}"])
        out[f"CO2_{box}_umolkg"] = to_umolkg(carb[box]["CO2"])
        out[f"HCO3_{box}_umolkg"] = to_umolkg(carb[box]["HCO3"])
        out[f"CO32_{box}_umolkg"] = to_umolkg(carb[box]["CO32"])
        out[f"PCO2_{box}_ppmv"] = to_ppmv(carb[box]["PCO2"])

    out["PCO2_ATM_ppmv"] = to_ppmv(state["PCO2A"])
    out["Total_Carbon_Initial_PgC"] = tocini * 12.0e-15
    out["Total_Carbon_Final_PgC"] = tocfin * 12.0e-15

    return out


def create_dataframe(out):
    rows = []

    for box in ["H", "L", "N", "D"]:
        rows.append(("Temperature (deg C)", box, out.get(f"Temperature_{box}_degC")))
    for box in ["H", "L", "N", "D"]:
        rows.append(("Salinity (PSU)", box, out.get(f"Salinity_{box}_PSU")))

    for box in ["H", "L", "N"]:
        rows.append(("Export production (PgC/yr)", box, out.get(f"EP_{box}_PgC_yr")))

    quantities = [
        ("PO4 (umol/kg)", "PO4", "umolkg"),
        ("DIC (umol/kg)", "DIC", "umolkg"),
        ("ALK (ueq/kg)", "ALK", "ueqkg"),
        ("O2 (umol/kg)", "DO2", "umolkg"),
        ("AOU (umol/kg)", "AOU", "umolkg"),
        ("CO2 (umol/kg)", "CO2", "umolkg"),
        ("HCO3 (umol/kg)", "HCO3", "umolkg"),
        ("CO3 (umol/kg)", "CO32", "umolkg"),
        ("PCO2 (ppmv)", "PCO2", "ppmv"),
    ]
    for quantity, key, unit in quantities:
        for box in ["H", "L", "N", "D"]:
            rows.append((quantity, box, out.get(f"{key}_{box}_{unit}")))

    rows.append(("PCO2 (ppmv)", "Atmosphere", out.get("PCO2_ATM_ppmv")))
    rows.append(("Total Carbon (PgC)", "Initial", out.get("Total_Carbon_Initial_PgC")))
    rows.append(("Total Carbon (PgC)", "Final", out.get("Total_Carbon_Final_PgC")))

    df = pd.DataFrame(rows, columns=["Quantity", "Box", "Value"])
    df["Value"] = df["Value"].map(
        lambda x: f"{x:.3f}" if isinstance(x, (int, float)) else x
    )
    return df


def run_model(max_iter=100000, verbose=True, params=None, return_history=False):
    if params is None:
        params = default_params()

    state, aux = initialize(params)
    tocini = total_carbon(state, params, aux)

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
                "PO4_N_umolkg": hist_out["PO4_N_umolkg"],
                "PO4_D_umolkg": hist_out["PO4_D_umolkg"],
                "DIC_H_umolkg": hist_out["DIC_H_umolkg"],
                "DIC_L_umolkg": hist_out["DIC_L_umolkg"],
                "DIC_N_umolkg": hist_out["DIC_N_umolkg"],
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
