"""Rigorous equilibrium model for eluate / 1 M acetic acid / 1 M Tris mixtures (the simulator's 'truth').

Charge balance with Davies activity coefficients, iterated on ionic strength.
Constants (sources in the chat): acetic acid pKa 4.756 at 25 C; Tris pKa 8.07 at 25 C with -0.028 per C
(i.e. 8.21 at 20 C); water pKw 14.00 at 25 C (14.17 at 20 C).
A liquid is a dict of total concentrations (M): Ac (acetate+acetic acid), Tr (Tris), Na (strong cations,
incl. Lys/Arg), Cl (strong anions), Pc (protein carboxyls, pKa ~4), Ph (protein histidines, pKa ~6.3).
"""
import numpy as np
from scipy.optimize import brentq

KEYS = ("Ac", "Tr", "Na", "Cl", "Pc", "Ph")
ACID = dict(Ac=1.0)          # 1 M acetic acid
TRIS = dict(Tr=1.0)          # 1 M Tris base

def mix(parts):
    """parts: list of (volume, liquid) -> mixed liquid (volume-weighted)."""
    V = sum(v for v, _ in parts)
    return {k: sum(v * liq.get(k, 0.0) for v, liq in parts) / V for k in KEYS}

def davies(I):
    s = np.sqrt(I)
    return 10 ** (-0.51 * (s / (1 + s) - 0.3 * I))          # |z| = 1

def ph_of(liq, T=21.0, pKa_c=4.0, pKa_h=6.3):
    c = {k: liq.get(k, 0.0) for k in KEYS}
    pKa_ac, pKa_tr = 4.756, 8.07 - 0.028 * (T - 25.0)
    pKw = 14.00 + 0.034 * (25.0 - T)
    I = 0.1
    for _ in range(30):
        g = davies(I)
        def charges(pH):
            aH = 10 ** -pH
            H = aH / g
            OH = 10 ** -pKw / aH / g
            ac_m = c["Ac"] / (1 + aH * g / 10 ** -pKa_ac)          # [A-] = Ka [HA] / (aH g)
            pc_m = c["Pc"] / (1 + aH * g / 10 ** -pKa_c)
            tr_p = c["Tr"] / (1 + 10 ** -pKa_tr * g / aH)          # [BH+] = aH [B] / (Ka g)
            ph_p = c["Ph"] / (1 + 10 ** -pKa_h * g / aH)
            return H, OH, ac_m, pc_m, tr_p, ph_p
        def balance(pH):
            H, OH, ac_m, pc_m, tr_p, ph_p = charges(pH)
            return H + c["Na"] + tr_p + ph_p - OH - c["Cl"] - ac_m - pc_m
        pH = brentq(balance, -1.0, 15.0, xtol=1e-10)
        H, OH, ac_m, pc_m, tr_p, ph_p = charges(pH)
        I_new = 0.5 * (H + OH + c["Na"] + c["Cl"] + ac_m + pc_m + tr_p + ph_p)
        if abs(I_new - I) < 1e-7:
            break
        I = 0.5 * (I + I_new)
    return pH

def acid_well(eluate, r, T=21.0):
    """pH after mixing 1 mL eluate with r mL 1 M acetic acid."""
    return ph_of(mix([(1.0, eluate), (r, ACID)]), T)

def pool_liquid(eluate, Rp):
    return mix([(1.0, eluate), (Rp, ACID)])

def base_well(pool, t, T=21.0):
    """pH after mixing 1 mL pool with t mL 1 M Tris."""
    return ph_of(mix([(1.0, pool), (t, TRIS)]), T)
