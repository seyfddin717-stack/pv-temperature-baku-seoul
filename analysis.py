"""Temperature-related energy loss of silicon PV modules: Baku vs Seoul.

Hourly irradiance on the tilted plane, ambient temperature and wind speed come
from PVGIS (PVGIS-ERA5 for both cities, so the comparison uses one database).
"""
import os
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "results")
os.makedirs(DATA, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

# city: (latitude, longitude, tilt in degrees = approx. latitude, facing south)
CITIES = {"Baku": (40.41, 49.87, 40), "Seoul": (37.57, 126.98, 37)}
START_YEAR, END_YEAR = 2015, 2020
RAD_DB = "PVGIS-ERA5"
API = "https://re.jrc.ec.europa.eu/api/v5_2/seriescalc"

G_MIN = 50                 # W/m2, hours with lower irradiance are ignored
GAMMAS = [-0.0030, -0.0035, -0.0040, -0.0045]   # temperature coefficients (1/C)
MAIN_GAMMA = -0.0040
NOCT = 45                  # C
U0, U1 = 25.0, 6.84        # Faiman model constants (W/m2K, W s/m3K)


def fetch(city):
    """Download hourly data once and cache it as CSV in data/."""
    path = os.path.join(DATA, f"{city}_{START_YEAR}_{END_YEAR}.csv")
    if os.path.exists(path):
        return pd.read_csv(path, index_col="time", parse_dates=True)
    lat, lon, tilt = CITIES[city]
    params = dict(lat=lat, lon=lon, startyear=START_YEAR, endyear=END_YEAR,
                  angle=tilt, aspect=0, raddatabase=RAD_DB, outputformat="json")
    r = requests.get(API, params=params, timeout=180)
    r.raise_for_status()
    df = pd.DataFrame(r.json()["outputs"]["hourly"])
    df["time"] = pd.to_datetime(df["time"], format="%Y%m%d:%H%M")
    df = df.set_index("time")[["G(i)", "T2m", "WS10m"]]
    df.to_csv(path)
    return df


def cell_temperature(df, model):
    g, ta, ws = df["G(i)"], df["T2m"], df["WS10m"]
    if model == "NOCT":                       # no wind
        return ta + (NOCT - 20) / 800 * g
    if model == "Faiman":                     # includes wind speed
        return ta + g / (U0 + U1 * ws)
    raise ValueError(model)


def energy_loss_pct(df, model, gamma):
    """Irradiance-weighted relative energy loss versus 25 C (percent)."""
    d = df[df["G(i)"] > G_MIN]
    tc = cell_temperature(d, model)
    factor = 1 + gamma * (tc - 25)
    return (1 - (d["G(i)"] * factor).sum() / d["G(i)"].sum()) * 100


data = {c: fetch(c) for c in CITIES}
years = END_YEAR - START_YEAR + 1

# ---- summary table ----
rows = []
for c, df in data.items():
    d = df[df["G(i)"] > G_MIN]
    for model in ("NOCT", "Faiman"):
        tc = cell_temperature(d, model)
        w = d["G(i)"]
        row = {"city": c, "model": model,
               "irradiation_kWh_m2_yr": round(df["G(i)"].sum() / 1000 / years, 0),
               "mean_T_ambient_C": round(d["T2m"].mean(), 1),
                              "mean_wind_10m_ms": round(d["WS10m"].mean(), 2),
               "irr_weighted_T_cell_C": round((tc * w).sum() / w.sum(), 1),
               "max_T_cell_C": round(tc.max(), 1)}
        for g in GAMMAS:
            row[f"loss_pct_gamma{g*100:.2f}"] = round(energy_loss_pct(df, model, g), 2)
        rows.append(row)
summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(OUT, "summary.csv"), index=False)
print(summary.T.to_string(header=False))
# ---- sensitivity to a possible irradiance bias (G scaled), gamma = -0.40 %/C ----
rows2 = []
for c, df in data.items():
    for scale in (1.0, 0.9, 0.8, 0.7):
        scaled = df.copy()
        scaled["G(i)"] = scaled["G(i)"] * scale
        for model in ("NOCT", "Faiman"):
            rows2.append({"city": c, "G_scale": scale, "model": model,
                          "loss_pct": round(energy_loss_pct(scaled, model, MAIN_GAMMA), 2)})
irr_sens = pd.DataFrame(rows2).pivot(index=["city", "G_scale"], columns="model",
                                     values="loss_pct")
irr_sens.to_csv(os.path.join(OUT, "irradiance_sensitivity.csv"))
print(irr_sens)
# ---- monthly loss, main scenario (Faiman, gamma = -0.40 %/C) ----
monthly = {}
for c, df in data.items():
    d = df[df["G(i)"] > G_MIN].copy()
    d["tc"] = cell_temperature(d, "Faiman")
    d["f"] = 1 + MAIN_GAMMA * (d["tc"] - 25)
    d["gf"] = d["G(i)"] * d["f"]
    m = d.groupby(d.index.month)
    monthly[c] = (1 - m["gf"].sum() / m["G(i)"].sum()) * 100
monthly = pd.DataFrame(monthly).round(2)
monthly.index.name = "month"
monthly.to_csv(os.path.join(OUT, "monthly_loss_pct.csv"))
print(monthly)

# ---- figures ----
fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))

x = np.arange(1, 13)
width = 0.38
axs[0].bar(x - width / 2, monthly["Baku"], width, label="Baku")
axs[0].bar(x + width / 2, monthly["Seoul"], width, label="Seoul")
axs[0].axhline(0, color="k", lw=0.8)
axs[0].set_xticks(x)
axs[0].set_xlabel("Month")
axs[0].set_ylabel("Energy loss vs 25 C (%)")
axs[0].set_title("Monthly temperature loss (Faiman, -0.40 %/C)")
axs[0].legend()
axs[0].grid(True, axis="y")

for c, df in data.items():
    d = df[df["G(i)"] > G_MIN]
    tc = cell_temperature(d, "Faiman")
    axs[1].hist(tc, bins=np.arange(-20, 71, 2), weights=d["G(i)"] / d["G(i)"].sum() * 100,
                alpha=0.55, label=c)
axs[1].axvline(25, color="k", ls="--", lw=1, label="25 C (STC)")
axs[1].set_xlabel("Cell temperature (C)")
axs[1].set_ylabel("Share of irradiation (%)")
axs[1].set_title("Cell temperature distribution (irradiance-weighted)")
axs[1].legend()
axs[1].grid(True)

labels = [f"{g*100:.2f}" for g in GAMMAS]
xx = np.arange(len(GAMMAS))
bw = 0.2
k = 0
for c in CITIES:
    for model in ("NOCT", "Faiman"):
        r = summary[(summary.city == c) & (summary.model == model)].iloc[0]
        vals = [r[f"loss_pct_gamma{g*100:.2f}"] for g in GAMMAS]
        axs[2].bar(xx + (k - 1.5) * bw, vals, bw, label=f"{c}, {model}")
        k += 1
axs[2].set_xticks(xx)
axs[2].set_xticklabels(labels)
axs[2].set_xlabel("Temperature coefficient (%/C)")
axs[2].set_ylabel("Annual energy loss (%)")
axs[2].set_title("Sensitivity to model and coefficient")
axs[2].legend(fontsize=8)
axs[2].grid(True, axis="y")

plt.tight_layout()
plt.savefig(os.path.join(OUT, "baku_vs_seoul.png"), dpi=150)
