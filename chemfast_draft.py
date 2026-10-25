"""Draft of the vectorised chemistry model that goes into EDP_tecan_gp_master.py (v14). Same physics as chem.py:
charge balance with Davies activity coefficients, acetic acid pKa 4.756, Tris pKa 8.07 at 25 C (-0.028 per C),
protein carboxyls pKa 4.0, histidine-like groups pKa 6.3, water pKw 14.00 at 25 C."""
import numpy as np

KEYS = ("Ac", "Tr", "Na", "Cl", "Pc", "Ph")


def ph_batch(c, T, it=200):
    """pH of many mixtures at once. c: dict of arrays (M) for KEYS; T: temperature (C). Newton steps on the charge
    balance (limited to 1 pH unit), with the Davies coefficients recomputed from the ionic strength as it converges
    (damped), until both pH and ionic strength settle."""
    c = {k: np.asarray(c.get(k, 0.0), float) for k in KEYS}
    n = np.broadcast(*c.values()).shape
    c = {k: np.broadcast_to(v, n).astype(float) for k, v in c.items()}
    pka_ac, pka_tr, pka_c, pka_h = 4.756, 8.07 - 0.028 * (T - 25.0), 4.0, 6.3
    kw = 10 ** -(14.00 + 0.034 * (25.0 - T))
    L = np.log(10)
    ph = np.full(n, 7.0)
    I = np.full(n, 0.1)
    for _ in range(it):
        s = np.sqrt(I)
        g = 10 ** (-0.51 * (s / (1 + s) - 0.3 * I))
        aH = 10 ** -ph
        H, OH = aH / g, kw / aH / g
        x_ac = 1 / (1 + aH * g / 10 ** -pka_ac)           # fraction dissociated
        x_c = 1 / (1 + aH * g / 10 ** -pka_c)
        x_tr = 1 / (1 + 10 ** -pka_tr * g / aH)           # fraction protonated
        x_h = 1 / (1 + 10 ** -pka_h * g / aH)
        f = H + c["Na"] + c["Tr"] * x_tr + c["Ph"] * x_h - OH - c["Cl"] - c["Ac"] * x_ac - c["Pc"] * x_c
        df = -L * (H + OH + c["Ac"] * x_ac * (1 - x_ac) + c["Pc"] * x_c * (1 - x_c) + c["Tr"] * x_tr * (1 - x_tr)
                   + c["Ph"] * x_h * (1 - x_h))                    # d f / d pH, always < 0
        step = np.clip(-f / df, -1.0, 1.0)
        I_new = 0.5 * (H + OH + c["Na"] + c["Cl"] + c["Ac"] * x_ac + c["Pc"] * x_c + c["Tr"] * x_tr + c["Ph"] * x_h)
        ph = np.clip(ph + step, -1.0, 15.0)
        dI = I_new - I
        I = I + 0.5 * dI
        if np.all(np.abs(step) < 1e-9) and np.all(np.abs(dI) < 1e-9):
            break
    return ph


def mix(parts):
    """parts: list of (volume array, dict of concentrations) -> mixed concentrations (dict of arrays)."""
    V = sum(np.asarray(v, float) for v, _ in parts)
    return {k: sum(np.asarray(v, float) * liq.get(k, 0.0) for v, liq in parts) / V for k in KEYS}


if __name__ == "__main__":
    import sys, os, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import chem
    rng = np.random.default_rng(1)
    worst = 0.0
    t0 = time.time()
    for trial in range(300):
        el = dict(Ac=rng.uniform(0, 0.2), Tr=rng.uniform(0, 0.15) * (rng.random() < 0.3), Na=rng.uniform(0, 0.2),
                  Cl=rng.uniform(0, 0.05), Pc=rng.uniform(0, 0.1), Ph=rng.uniform(0, 0.02))
        T = rng.uniform(15, 28)
        A, B = rng.uniform(0, 4) * (rng.random() < 0.7), rng.uniform(0, 2) * (rng.random() < 0.6)
        liq = chem.mix([(1.0, el), (A, chem.ACID), (B, chem.TRIS)])
        ref = chem.ph_of(liq, T)
        got = ph_batch(liq, T)
        worst = max(worst, abs(float(got) - ref))
    print(f"300 random mixtures: worst difference from chem.py {worst:.2e} pH")
    # speed on a 40-point batch
    el = dict(Ac=0.1, Na=0.1, Pc=0.02, Ph=0.005)
    A = np.r_[np.linspace(0, 5, 16), np.zeros(12), np.linspace(0.5, 3.5, 12)]
    B = np.r_[np.zeros(16), np.linspace(0.005, 0.3, 12), np.linspace(0.5, 2.0, 12)]
    liq = mix([(np.ones_like(A), el), (A, {"Ac": 1.0}), (B, {"Tr": 1.0})])
    t0 = time.time()
    for _ in range(200):
        ph_batch(liq, 22.0)
    print(f"40-mixture batch: {1000 * (time.time() - t0) / 200:.2f} ms per call")
    ref = [chem.ph_of({k: float(v[i]) for k, v in liq.items()}, 22.0) for i in range(len(A))]
    print("batch vs chem.py worst:", f"{np.max(np.abs(ph_batch(liq, 22.0) - np.array(ref))):.2e}")
