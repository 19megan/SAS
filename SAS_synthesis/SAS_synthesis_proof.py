# =============================================================================
# PROOF / VALIDATION for SAS-derived metric #2 ("bottom vs top" plot)
# in SAS_synthesis.py
#
# Goal: build three synthetic catchments with a KNOWN storage effect, run each
# through the mesas SAS model, then check that the "bottom (dryness) vs
# top (young water fraction)" plot correctly identifies which effect is present.
#
#   1. Inverse Storage Effect (ISE): prefer YOUNG water when storage is HIGH
#   2. Direct  Storage Effect (DSE): prefer OLD   water when storage is HIGH
#   3. No      Storage Effect (NSE): preference independent of storage
#
# Mechanism: the SAS function is Kumaraswamy with b=1, so its CDF is
#   Omega(x) = x^a   on x = S_T / S  in [0, 1]   (S_T = age-ranked storage)
# and its pdf is  omega(x) = a * x^(a-1).
#   a < 1  -> weight on x~0 (YOUNG water)  -> young-water preference
#   a > 1  -> weight on x~1 (OLD   water)  -> old-water preference
#   a = 1  -> uniform sampling (no age preference)
#
# We make `a` (column "k") a function of normalized storage w in [0,1]:
#   ISE: a = a_young + (1-w)*(a_old-a_young)   -> high storage => a_young (<1)
#   DSE: a = a_young +    w *(a_old-a_young)   -> high storage => a_old   (>1)
#   NSE: a = 1.0 (constant)
#
# Date: 2026-06-30
# =============================================================================

#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mesas.sas.model import Model

rng = np.random.default_rng(42)

# -----------------------------------------------------------------------------
# 1. Build a synthetic catchment (daily linear-reservoir bucket model)
# -----------------------------------------------------------------------------
N_YEARS   = 8
N         = N_YEARS * 365          # days
SPINUP    = 3 * 365                # days discarded from analysis
MAX_AGE   = 1460                   # days tracked (TTD fully decays well before this)
TAU       = 250.0                  # linear reservoir residence time (days)
S0        = 500.0                 # initial / mean storage (mm)
A_YOUNG   = 0.15                   # Kumaraswamy a for strong young-water pref
A_OLD     = 5.0                    # Kumaraswamy a for strong old-water pref
YW_AGE    = 90                     # young-water age threshold (days)

dates = pd.date_range('2000-01-01', periods=N, freq='D')
doy   = dates.dayofyear.to_numpy()
season = 2 * np.pi * doy / 365.0

# Seasonal potential ET (mm/d): low in winter, high in summer
PET = 0.25 - .025 * np.cos(season)          # range ~[0.5, 4.5]

# Stochastic precip: wetter in winter, gamma-distributed intensities
# p_wet = np.clip(0.30 + 0.15 * np.cos(season), 0.05, 0.95)
# wet   = rng.random(N) < p_wet
# J     = wet * rng.gamma(shape=0.8, scale=11.0, size=N)   # mm/d
def simulate_daily_rainfall_rate(
    n_days=30*365,         # number of days to simulate
    prob_rain_per_day=0.3,   # probability of rain each day
    shape=2.0,            # shape parameter for gamma distribution (controls skew)
    scale=1.5,            # scale parameter for gamma distribution  (controls mean intensity)
    hours_per_day=24
):
    """
    Simulate daily rainfall rate (mm/hr) using a Bernoulli process for rain occurrence
    and a gamma distribution for intensity on rainy days.
    Returns a pandas DataFrame with hourly rainfall rate.
    """

    # Rain days (1 if rain, 0 if no rain)
    rain_days = rng.random(n_days) < prob_rain_per_day

    # Rain intensity
    rain_mm_day = np.zeros(n_days)
    rain_mm_day[rain_days] = rng.gamma(shape, scale, size=rain_days.sum())
    
    return rain_mm_day


prpd = [0.1, 0.5, .9]
shape = [.25, .5, 1, 1.5]
scale = [10, 15, 20]
b_ISE = []
b_DSE = []
b_NSE = []
t_ISE = []
t_DSE = []
t_NSE = []

for i in prpd:
    for j in shape:
        for k in scale:
            J = simulate_daily_rainfall_rate(n_days=N, prob_rain_per_day=i, shape=j, scale=k)



            # Integrate the bucket: Q = S/TAU (linear reservoir), AET limited by storage
            S  = np.zeros(N)
            Q  = np.zeros(N)
            ET = np.zeros(N)
            s  = S0
            for t in range(N):
                S[t]  = s                       # storage at start of day = SAS scale
                q     = s / TAU
                aet   = PET[t] * s / (s + 50.0) # Budyko-style storage limitation
                Q[t]  = q
                ET[t] = aet
                s     = max(s + J[t] - q - aet, 1.0)

            # Normalized storage as a PERCENTILE/rank (0 = driest day, 1 = wettest day).
            # Rank-based (not min-max) so the wettest/driest days actually exercise the
            # full young/old preference range instead of being compressed by outliers.
            from scipy.stats import rankdata
            w = (rankdata(S) - 1) / (N - 1)
            w = (S-S.min()) / (S.max()-S.min())  # alternative: min-max normalization
            # -----------------------------------------------------------------------------
            # 2. Encode the three storage effects as the Kumaraswamy `a` (= column "k")
            # -----------------------------------------------------------------------------
            k_cases = {
                'Inverse SE (young when wet)': A_YOUNG + (1 - w) * (A_OLD - A_YOUNG),
                'Direct SE (old when wet)'   : A_YOUNG +      w  * (A_OLD - A_YOUNG),
                'No SE (constant a=1)'       : np.full(N, 1.0),
            }

            # Seed initial age-ranked storage as exponential w/ total ~ S0 (steady state)
            age = np.arange(MAX_AGE)
            sT_init = (S0 / TAU) * np.exp(-age / TAU)

            def build_model(k):
                df = pd.DataFrame({'J': J, 'Q': Q, 'ET': ET, 'S': S, 'k': k}, index=dates)
                sas_specs = {
                    'Q':  {'Q SAS':  {'func': 'kumaraswamy',
                                    'args': {'loc': 0.0, 'scale': 'S', 'a': 'k', 'b': 1.0}}},
                    'ET': {'ET SAS': {'func': 'kumaraswamy',
                                    'args': {'loc': 0.0, 'scale': 'S', 'a': 0.5, 'b': 1.0}}},
                }
                return Model(data_df=df, sas_specs=sas_specs, influx='J',
                            dt=1, max_age=MAX_AGE, sT_init=sT_init,
                            verbose=False, record_state=True)

            # -----------------------------------------------------------------------------
            # 3. Run each case, compute top (young water fraction) & bottom (dryness)
            # -----------------------------------------------------------------------------
            # bottom (dryness) is the SAME across cases: it only depends on Q
            bottom_full = (Q.max() - Q) / (Q.max() - Q.min())

            results = {}
            for label, k in k_cases.items():
                print(f'Running: {label} ...')
                m = build_model(k)
                m.run()
                pq = m.get_pQ('Q')                              # (max_age, N)
                PQ = np.cumsum(pq, axis=0) * m.options['dt']    # cumulative TTD
                top = PQ[YW_AGE, :]                             # young water fraction
                sl = slice(SPINUP, N)                           # drop spin-up
                b, tp = bottom_full[sl], top[sl]
                slope, intercept = np.polyfit(b, tp, 1)
                corr = np.corrcoef(b, tp)[0, 1]
                # diagnostic: young water fraction in the wettest 20% vs driest 20% of days
                wet_mask = b < np.quantile(b, 0.20)     # low dryness
                dry_mask = b > np.quantile(b, 0.80)     # high dryness
                yw_wet, yw_dry = tp[wet_mask].mean(), tp[dry_mask].mean()
                results[label] = dict(bottom=b, top=tp, slope=slope, intercept=intercept,
                                    corr=corr, k=k[sl], yw_wet=yw_wet, yw_dry=yw_dry,
                                    yw_mean=tp.mean())
                #save results to csv
                if label[0]=='I':
                    b_ISE.append(b)
                    t_ISE.append(tp)
                elif label[0]=='D':
                    b_DSE.append(b)
                    t_DSE.append(tp)
                elif label[0]=='N':
                    b_NSE.append(b)
                    t_NSE.append(tp)


            

colors = {'Inverse SE (young when wet)': 'tab:green',
          'Direct SE (old when wet)'   : 'tab:blue',
          'No SE (constant a=1)'       : 'tab:orange'}

fig, ax = plt.subplots(1, 3, figsize=[16, 5])

for label, k in k_cases.items():
    if label[0]=='I':
        for i in range(len(b_ISE)):
            ax[0].plot(b_ISE[i], t_ISE[i], '.', ms=3, alpha=0.4, color=colors[label])
    elif label[0]=='D':
        for j in range(len(b_DSE)):
            ax[1].plot(b_DSE[j], t_DSE[j], '.', ms=3, alpha=0.4, color=colors[label])
    elif label[0]=='N':
        for k in range(len(b_NSE)):
            ax[2].plot(b_NSE[k], t_NSE[k], '.', ms=3, alpha=0.4, color=colors[label])
    
ax[0].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[0].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[0].set_title('Inverse storage effect Synthetic proof: bottom vs top')
ax[1].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[1].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[1].set_title('Direct storage effect Synthetic proof: bottom vs top')
ax[2].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[2].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[2].set_title('No storage effect Synthetic proof: bottom vs top')
plt.tight_layout()

#%%
# -----------------------------------------------------------------------------
# 4. Verdict: slope sign should identify the effect
# -----------------------------------------------------------------------------
# print('\n' + '=' * 70)
# print('VERDICT  (bottom = dryness, top = young water fraction)')
# print('  Inverse SE  => NEGATIVE slope (young water dominates when wet)')
# print('  Direct  SE  => POSITIVE slope (old   water dominates when wet)')
# print('  No      SE  => ~FLAT   slope (~0)')
# print('=' * 70)
# for label, r in results.items():
#     sign = 'NEGATIVE' if r['slope'] < -0.02 else ('POSITIVE' if r['slope'] > 0.02 else 'FLAT (~0)')
#     print(f'{label:32s} slope={r["slope"]:+.3f}  corr={r["corr"]:+.3f}  -> {sign}')

# print('\nDIAGNOSTIC: young water fraction, wettest 20% vs driest 20% of days')
# print(f'{"case":32s} {"YWF_wet":>8s} {"YWF_dry":>8s} {"wet-dry":>8s} {"YWF_mean":>9s}')
# for label, r in results.items():
#     print(f'{label:32s} {r["yw_wet"]:8.3f} {r["yw_dry"]:8.3f} '
#           f'{r["yw_wet"]-r["yw_dry"]:+8.3f} {r["yw_mean"]:9.3f}')
# # Compare each case to the No-SE baseline (removes the flux-driven confound)
# base = results['No SE (constant a=1)']
# print('\nDIFFERENCE FROM NO-SE BASELINE  (isolates the storage effect):')
# for label, r in results.items():
#     d_wet = r['yw_wet'] - base['yw_wet']
#     d_dry = r['yw_dry'] - base['yw_dry']
#     print(f'{label:32s} dYWF_wet={d_wet:+.3f}  dYWF_dry={d_dry:+.3f}')

# # ---- Automated PASS/FAIL: the metric must ORDER the three effects correctly --
# s_inv = results['Inverse SE (young when wet)']['slope']
# s_dir = results['Direct SE (old when wet)']['slope']
# s_non = results['No SE (constant a=1)']['slope']
# print('\n' + '=' * 70)
# ordered = s_inv < s_non < s_dir
# print(f'ORDERING TEST  slope(Inverse) < slope(None) < slope(Direct):')
# print(f'   {s_inv:+.3f}  <  {s_non:+.3f}  <  {s_dir:+.3f}   ->  '
#       f'{"PASS - plot correctly ranks the storage effects" if ordered else "FAIL"}')
# print('CAVEATS revealed by this proof:')
# print('  * The No-SE null is NOT flat and NOT the 1:1 line: a flux confound')
# print('    (fresh rain dominates high flow) makes young-water rise when wet,')
# print(f'    giving a weakly NEGATIVE null slope ({s_non:+.3f}), not 0 and not +1.')
# print('  * => the right null reference is the a=1 run, not bottom-vs-bottom.')
# print('  * A weakly-negative slope is ambiguous (weak inverse vs pure confound);')
# print('    only clearly-positive (Direct) or steeply-negative (Inverse) are safe.')
# print('=' * 70)

# -----------------------------------------------------------------------------
# 5. Plot: replicate the "bottom vs top" plot for all three synthetic cases
# -----------------------------------------------------------------------------
colors = {'Inverse SE (young when wet)': 'tab:green',
          'Direct SE (old when wet)'   : 'tab:blue',
          'No SE (constant a=1)'       : 'tab:orange'}

fig, ax = plt.subplots(1, 2, figsize=[14, 5])

for label, r in results.items():
    ax[0].plot(r['bottom'], r['top'], '.', ms=3, alpha=0.4, color=colors[label],
               label=f'{label} (slope={r["slope"]:+.2f})')
    xs = np.array([0, 1])
    ax[0].plot(xs, r['intercept'] + r['slope'] * xs, '-', color=colors[label], lw=2)

ax[0].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[0].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[0].set_title('Synthetic proof: bottom vs top')
ax[0].legend(fontsize=8, loc='best')

# Fit lines only, with the EMPIRICAL no-SE run as the correct null reference.
# b_ref = np.linspace(0, 1, 50)
# for label, r in results.items():
#     ax[1].plot(b_ref, r['intercept'] + r['slope'] * b_ref, '-',
#                color=colors[label], lw=2.5,
#                label=f'{label} fit (slope={r["slope"]:+.2f})')

# ax[1].set_xlabel('Dryness  (0 = wet,  1 = dry)')
# ax[1].set_ylabel(f'Young water fraction (Age <= {YW_AGE} d)')
# ax[1].set_title('Fitted slopes\n'
#                 '(true null = orange a=1 run, weakly negative)')
# ax[1].legend(fontsize=8, loc='best')


from scipy.stats import gamma
x=np.linspace(gamma.ppf(0.0001, shape, scale=scale),
                  gamma.ppf(0.999, shape, scale=scale), 100)
y=gamma.pdf(x, shape, scale=scale)
ax[1].plot(x, y, label=f'Gamma distribution of rainfall with reservoir storage {TAU}')
ax[1].set_xlabel('rainfall amount (mm)')
ax[1].set_ylabel('gamma distribution')
ax[1].set_title('Precipitation distribution')
ax[1].legend()

plt.tight_layout()
plt.savefig('SAS_synthesis_proof.png', dpi=130)
print('\nSaved figure -> SAS_synthesis_proof.png')
plt.show()
# %%
function = 'exp fit'  # 'log fit', 'poly fit', 'exp fit'

from scipy.optimize import curve_fit
# fit functions to scatter plots

def exponential_model(x, a, b, c):
    """
    a: scale factor / initial value
    b: growth/decay rate constant
    c: vertical offset (shift)
    """
    return a * np.exp(b * x) + c



fig, ax = plt.subplots(1, 3, figsize=[16, 5])

for label in results.keys():
    if label[0]=='I':
        b = b_ISE
        t = t_ISE
        idx = 0
        init_params = [1, -1, 0]
    elif label[0]=='D':
        b = b_DSE
        t = t_DSE
        idx = 1
        init_params = [1, 10, 0]
    elif label[0]=='N':
        b = b_NSE
        t = t_NSE
        idx = 2
        init_params = [1, 0, 0]

    b_all = [item for column in b for item in column]
    t_all = [item for column in t for item in column]

    for i in range(len(b)):
        ax[idx].plot(b[i], t[i], '.', ms=3, alpha=0.4, color=colors[label])

    if function == 'log fit':
        #log fit
        coeff = np.polyfit(b_all, np.log(t_all), 1)
        print(f'Fitted log-linear coefficients for {label}: slope={coeff[0]:.3f}, intercept={coeff[1]:.3f}')
        poly_func = np.polyval(coeff, b_all)
        b_smooth = np.linspace(min(b_all), max(b_all), 100)
        t_smooth = np.polyval(coeff, b_smooth)
        ax[idx].plot(b_smooth, np.exp(t_smooth), color='black', lw=2, label=f'{function}: coefficients a={coeff[0]:.2f}, b={coeff[1]:.2f}')
    elif function == 'poly fit':
        # polynomial fit
        coeff = np.polyfit(b_all, t_all, 2) #[0.9213325, -1.41866932, 0.54897838]
        print(f'Fitted polynomial coefficients for {label}:\na = {coeff[0]:.3f}\nb = {coeff[1]:.3f}\nc = {coeff[2]:.3f}')
        poly_func = np.polyval(coeff, b_all)
        b_smooth = np.linspace(min(b_all), max(b_all), 100)
        t_smooth = np.polyval(coeff, b_smooth)
        ax[idx].plot(b_smooth, t_smooth, color='black', lw=2, label=f'{function}: coefficients a={coeff[0]:.2f}, b={coeff[1]:.2f}, c={coeff[2]:.2f}')
    elif function == 'exp fit':
        #exponential fit
        popt, pcov = curve_fit(exponential_model, b_all, t_all, p0=init_params)
        # Extract optimized parameters
        a_opt, b_opt, c_opt = popt
        print(f"Fitted exponential parameters for {label}:\na = {a_opt:.3f}\nb = {b_opt:.3f}\nc = {c_opt:.3f}")
        b_smooth = np.linspace(min(b_all), max(b_all), 100)
        t_smooth = exponential_model(b_smooth, a_opt, b_opt, c_opt)
        ax[idx].plot(b_smooth, t_smooth, color='black', lw=2, label=f'{function}: a={a_opt:.2f}, b={b_opt:.2f}, c={c_opt:.2f}')

    
    ax[idx].legend(fontsize=8, loc='best')


ax[0].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[0].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[0].set_title('Inverse storage effect Synthetic proof: bottom vs top')
ax[1].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[1].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[1].set_title('Direct storage effect Synthetic proof: bottom vs top')
ax[2].set_xlabel('Dryness  (0 = wet / high storage,  1 = dry / low storage)')
ax[2].set_ylabel(f'Young water fraction  (Age <= {YW_AGE} d)')
ax[2].set_title('No storage effect Synthetic proof: bottom vs top')
plt.tight_layout()



# %%
