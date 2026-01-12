
import numpy as np
import pandas as pd
from dataclasses import dataclass

_LAMBDA_17 = 0.528  # ln-based slope for 17O

@dataclass
class SnowIsoOptions:
    # Whether to include seasonal cycle for δ18O via sin/cos of day-of-year
    include_seasonal: bool = True
    # Clip bounds for δ18O predictions (‰). Use None to disable.
    clip_d18O_min: float | None = -40.0
    clip_d18O_max: float | None = 0.0
    # If True, enforce δ2H = m*δ18O + c using LMWL from rainfall
    force_lmwl_for_d2H: bool = True

def _design_matrix_d18O(doy: np.ndarray, T: np.ndarray, include_seasonal: bool):
    X = [np.ones_like(T), T]
    if include_seasonal:
        ang = 2 * np.pi * doy / 365.25
        X += [np.sin(ang), np.cos(ang)]
    return np.vstack(X).T

def _ols(X, y):
    # Ordinary least squares with pseudo-inverse
    beta = np.linalg.pinv(X) @ y
    return beta

def _predict(beta, X):
    return X @ beta

def _compute_lmwl(d18, d2H):
    # δ2H = m*δ18O + c
    X = np.vstack([d18, np.ones_like(d18)]).T
    m, c = np.linalg.pinv(X) @ d2H
    return float(m), float(c)

def _ln1p_permil(x):
    # ln(1 + x/1000)
    return np.log1p(x / 1000.0)

def _exp_permil(x):
    # (exp(x) - 1) * 1000
    return (np.exp(x) - 1.0) * 1000.0

def _mean_17O_excess(d18, d17):
    # 17O-excess = ln(1+δ17/1000) - λ * ln(1+δ18/1000)
    return np.nanmean(_ln1p_permil(d17) - _LAMBDA_17 * _ln1p_permil(d18))

def estimate_snow_isotopes(df: pd.DataFrame, opts: SnowIsoOptions | None = None) -> pd.DataFrame:
    """
    Estimate snowfall δ18O, δ2H, δ17O from rainfall isotope data and temperature.
    Input DataFrame must include:
      - date (datetime-like index or column)
      - T_C (air temperature, °C)
      - rainfall_mm (>=0)
      - snowfall_swe_mm (>=0)
      - d18O_rain_permil, d2H_rain_permil (‰) for rainfall events
      - d17O_rain_permil (‰) is optional but recommended for 17O-excess carryover

    Returns original df with appended columns on dates where snowfall occurs (snowfall_swe_mm>0):
      - d18O_snow_est_permil, d2H_snow_est_permil, d17O_snow_est_permil
    """
    if opts is None:
        opts = SnowIsoOptions()

    df = df.copy()
    if "date" in df.columns:
        df["date"] = pd.to_datetime(df["date"])
        df = df.set_index("date")
    df = df.sort_index()

    # Prepare training set from rainfall-only rows with valid isotopes
    rain_mask = (df.get("rainfall_mm", 0).fillna(0) > 0)
    valid_rain_iso = rain_mask & df["d18O_rain_permil"].notna() & df["d2H_rain_permil"].notna()
    train = df.loc[valid_rain_iso].copy()
    if train.empty:
        raise ValueError("No rainfall rows with isotopes found for training.")

    # δ18O ~ a + b*T + seasonal terms
    doy = train.index.dayofyear.values.astype(float)
    T = train["T_C"].values.astype(float)
    X = _design_matrix_d18O(doy, T, include_seasonal=opts.include_seasonal)
    y18 = train["d18O_rain_permil"].values.astype(float)
    beta18 = _ols(X, y18)

    # LMWL from rainfall for δ2H
    m_lmwl, c_lmwl = _compute_lmwl(train["d18O_rain_permil"].values.astype(float),
                                   train["d2H_rain_permil"].values.astype(float))

    # 17O-excess (ln definition) from rainfall; fallback to 0 if no d17O provided
    if "d17O_rain_permil" in train.columns and train["d17O_rain_permil"].notna().any():
        mean_17Oex = _mean_17O_excess(train["d18O_rain_permil"].values.astype(float),
                                      train["d17O_rain_permil"].values.astype(float))
    else:
        mean_17Oex = 0.0  # neutral fallback

    # Predict δ18O for snowfall dates
    snow_mask = (df.get("snowfall_swe_mm", 0).fillna(0) > 0)
    if not snow_mask.any():
        # nothing to predict
        df["d18O_snow_est_permil"] = np.nan
        df["d2H_snow_est_permil"] = np.nan
        df["d17O_snow_est_permil"] = np.nan
        return df

    snow = df.loc[snow_mask].copy()
    Xs = _design_matrix_d18O(snow.index.dayofyear.values.astype(float),
                             snow["T_C"].values.astype(float),
                             include_seasonal=opts.include_seasonal)
    d18_snow = _predict(beta18, Xs)

    # Optional clipping
    if opts.clip_d18O_min is not None or opts.clip_d18O_max is not None:
        d18_snow = np.clip(d18_snow,
                           opts.clip_d18O_min if opts.clip_d18O_min is not None else -np.inf,
                           opts.clip_d18O_max if opts.clip_d18O_max is not None else np.inf)

    # δ2H via LMWL
    if opts.force_lmwl_for_d2H:
        d2H_snow = m_lmwl * d18_snow + c_lmwl
    else:
        # If not forcing LMWL, build separate δ2H model (T + season)
        y2H = train["d2H_rain_permil"].values.astype(float)
        beta2H = _ols(X, y2H)
        d2H_snow = _predict(beta2H, Xs)

    # δ17O from mean 17O-excess carryover
    d17_snow = _exp_permil(_LAMBDA_17 * np.log1p(d18_snow/1000.0) + mean_17Oex)

    # Write back
    df.loc[snow_mask, "d18O_snow_est_permil"] = d18_snow
    df.loc[snow_mask, "d2H_snow_est_permil"] = d2H_snow
    df.loc[snow_mask, "d17O_snow_est_permil"] = d17_snow

    # Store fitted parameters for transparency
    df.attrs["delta18O_beta"] = beta18.tolist()
    df.attrs["lmwl_m"] = m_lmwl
    df.attrs["lmwl_c"] = c_lmwl
    df.attrs["mean_17Oexcess_ln"] = mean_17Oex
    return df

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Estimate snowfall isotopes from rainfall isotopes + temperature.")
    p.add_argument("--in", dest="infile", required=True, help="Input CSV with columns: date,T_C,rainfall_mm,snowfall_swe_mm,d18O_rain_permil,d2H_rain_permil[,d17O_rain_permil]")
    p.add_argument("--out", dest="outfile", required=True, help="Output CSV path")
    p.add_argument("--no-season", action="store_true", help="Disable seasonal terms for δ18O model")
    args = p.parse_args()

    df = pd.read_csv(args.infile, parse_dates=["date"])
    opts = SnowIsoOptions(include_seasonal=not args.no_season)
    out = estimate_snow_isotopes(df, opts)
    out.to_csv(args.out)
    print(f"Wrote {args.out}")
