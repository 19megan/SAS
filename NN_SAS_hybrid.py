# Hybrid NN / SAS model for ORPB 18O
# Date created by Claude code: 04/07/2026

# ----------------------------------------------------------------------------
# Idea: keep the SAS *structure* from the gamma-split model (separate quickflow
# and baseflow TTDs convolved against precip-18O history, then flux-weighted)
# but let a small neural network learn the mapping
#     hydrologic_state(t) -> TTD shape parameters(t)
# for BOTH the quickflow and baseflow gamma distributions.
#
# Gamma-split structure (mirrors make_gamma_split_model_from in ORPB_SAS.py):
#     Quickflow: Gamma(a_qf, scale_qf)  — young water, short travel time
#     Baseflow:  Gamma(a_bf, scale_bf)   — old water, storage-dependent travel time
#     C_pred(t) = qf_weight(t) * C_qf(t) + bf_weight(t) * C_bf(t)
#   where C_qf and C_bf are each computed by convolving their respective TTD
#   against the influx-weighted precip 18O history.
#
# At each timestep t the NN predicts 4 parameters:
#     (a_qf, scale_qf, a_bf, scale_bf) = NN( storage, discharge, qf_weight, doy_sin, doy_cos )
#
# The gamma pdf is differentiable in (a, scale) via torch.lgamma, so the whole
# pipeline trains end-to-end with autograd against observed ORPB 18O.
#
# Why this is expected to do better than the pure-NN model in NN_hydro_similarity.py:
#   - Mass balance, convolution structure, and the qf/bf split are HARD-CODED.
#   - The NN only has to learn how TTD shapes vary with storage state — a few
#     hundred labeled samples are enough for that, whereas learning the entire
#     input->output map from scratch is not.
#
# Why this might do better than pure ORPB_SAS.py (gamma-split model):
#   - Pure SAS uses fixed a_qf, t_qf and a linear S_scale = lamda*(S - S_c) for
#     baseflow. The NN can learn an arbitrary nonlinear mapping for BOTH flow
#     paths, which may matter for the 2015 drought / 2018 wet extremes that a
#     linear S_scale cannot represent well.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
from sklearn.metrics import r2_score

# =============================================================================
# SECTION 1: Load data (same preprocessing as ORPB_SAS.py / NN_hydro_similarity.py)
# =============================================================================
data_df = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0])
data_df['precip 18O']      = data_df['precip 18O'].ffill().bfill()
data_df['influx (mm/hr)']  = data_df[['rainfall (mm/hr)', 'snowmelt (mm/hr)']].sum(axis=1)

# Quickflow / baseflow split — same as ORPB_SAS.py
data_df['quickflow (mm/hr)'] = data_df['discharge (mm/hr)'] - data_df['baseflow 1 (mm/hr)']
data_df['bf1_weight'] = data_df['baseflow 1 (mm/hr)'] / data_df['discharge (mm/hr)']
data_df['qf_weight']  = data_df['quickflow (mm/hr)']   / data_df['discharge (mm/hr)']

issample = data_df['ORPB 18O'].notna().values
obs_idx  = np.where(issample)[0]
N        = len(data_df)

print(f"Total hourly timesteps: {N}")
print(f"Stream 18O observations: {issample.sum()}")

# =============================================================================
# SECTION 2: Build per-observation history windows
# =============================================================================
# For each observed timestep we need the prior MAX_AGE hours of (influx, C_precip)
# so we can convolve them with a per-timestep gamma TTD.
#
# MAX_AGE=2000 hr (~12 weeks) — long enough to capture the slow tail seen in the
# ORPB SAS calibration but short enough to keep the (N_obs, MAX_AGE) tensors small.
# Tweak downward if memory is tight.

MAX_AGE = 2000
DT      = 1.0   # hours

influx    = data_df['influx (mm/hr)'].values.astype(np.float32)
c_precip  = data_df['precip 18O'].values.astype(np.float32)
storage   = data_df['storage (mm)'].values.astype(np.float32)
discharge = data_df['discharge (mm/hr)'].values.astype(np.float32)
bf_weight = data_df['bf1_weight'].fillna(0.5).values.astype(np.float32)
qf_weight = data_df['qf_weight'].fillna(0.5).values.astype(np.float32)

# Pad the front of the record with zeros so windows for early observations work.
pad_zeros_inf = np.zeros(MAX_AGE, dtype=np.float32)
pad_cprecip   = np.full(MAX_AGE, np.nan, dtype=np.float32)  # marks "older than record"
influx_pad    = np.concatenate([pad_zeros_inf, influx])
cprecip_pad   = np.concatenate([pad_cprecip,   c_precip])

# For each observed timestep build the window of the previous MAX_AGE hours.
# Convention: lag k=0 is the current timestep, k=MAX_AGE-1 is MAX_AGE-1 hours ago.
N_obs = len(obs_idx)
infl_win = np.zeros((N_obs, MAX_AGE), dtype=np.float32)
cpre_win = np.zeros((N_obs, MAX_AGE), dtype=np.float32)
old_mask = np.zeros((N_obs, MAX_AGE), dtype=np.float32)   # 1 where the lag falls before t=0

for i, t in enumerate(obs_idx):
    # window indices in the padded arrays: lag 0 = position MAX_AGE+t
    end   = MAX_AGE + t + 1            # exclusive
    start = end - MAX_AGE
    w_inf = influx_pad[start:end][::-1]    # reverse so index k = lag k
    w_cpr = cprecip_pad[start:end][::-1]
    nan_m = np.isnan(w_cpr)
    w_cpr[nan_m] = 0.0                      # will be replaced by C_old contribution
    infl_win[i] = w_inf
    cpre_win[i] = w_cpr
    old_mask[i] = nan_m.astype(np.float32)

# State features for the NN at each observed timestep
doy = data_df.index.dayofyear.values.astype(np.float32)
doy_sin = np.sin(2 * np.pi * doy / 365.25)
doy_cos = np.cos(2 * np.pi * doy / 365.25)

state = np.stack([
    storage[obs_idx],
    discharge[obs_idx],
    qf_weight[obs_idx],
    doy_sin[obs_idx],
    doy_cos[obs_idx],
], axis=1).astype(np.float32)

# z-score the state features (use full observed window — no fold splits here yet)
state_mean = state.mean(axis=0, keepdims=True)
state_std  = state.std(axis=0, keepdims=True) + 1e-9
state_z    = (state - state_mean) / state_std

target = data_df['ORPB 18O'].values[obs_idx].astype(np.float32)
times  = data_df.index[obs_idx]

print(f"History window tensor: {infl_win.shape}  (~{infl_win.nbytes/1e6:.1f} MB)")

# =============================================================================
# SECTION 3: Hybrid model — NN predicts gamma TTD params, convolution is fixed
# =============================================================================
# The NN outputs unconstrained (log_a, log_scale); we exponentiate so that both
# stay positive. C_old is a single learned scalar (initialised at -7.6 like SAS).

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Precompute lag axis (k=0..MAX_AGE-1) on the device
lag_k = torch.arange(MAX_AGE, dtype=torch.float32, device=device) * DT + 0.5 * DT
log_lag_k = torch.log(lag_k)

def gamma_log_pdf(log_k, k, log_a, log_scale):
    """
    log of Gamma(a, scale) pdf evaluated at lag k (broadcasts over batch).
    log_a, log_scale: (batch, 1)
    log_k, k:         (1, MAX_AGE)
    Returns:          (batch, MAX_AGE)
    """
    a     = torch.exp(log_a)
    scale = torch.exp(log_scale)
    # log pdf = (a-1)*log(k) - k/scale - a*log(scale) - lgamma(a)
    return (a - 1.0) * log_k - k / scale - a * log_scale - torch.lgamma(a)

def _convolve_ttd(log_a, log_scale, infl_win, cpre_win, old_mask, c_old):
    """Convolve a single gamma TTD against the influx-weighted precip history.
    Returns predicted concentration (B,)."""
    log_w = gamma_log_pdf(log_lag_k.unsqueeze(0), lag_k.unsqueeze(0),
                          log_a, log_scale)
    log_w = log_w - log_w.max(dim=1, keepdim=True).values
    w     = torch.exp(log_w)

    flux_w      = w * infl_win
    in_rec_num  = (flux_w * (1.0 - old_mask) * cpre_win).sum(dim=1)
    old_rec_num = (flux_w * old_mask).sum(dim=1) * c_old
    denom       = flux_w.sum(dim=1) + 1e-9
    return (in_rec_num + old_rec_num) / denom

class HybridSplitSASNet(nn.Module):
    """
    Tiny MLP that maps catchment state ->
        (log_a_qf, log_scale_qf, log_a_bf, log_scale_bf)
    i.e. 4 gamma TTD parameters: 2 for quickflow and 2 for baseflow.
    Stream concentration = qf_weight * C_qf + bf_weight * C_bf
    mirroring the gamma-split model in ORPB_SAS.py.
    """
    def __init__(self, n_state, hidden=16):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_state, hidden),
            nn.Tanh(),
            nn.Linear(hidden, hidden),
            nn.Tanh(),
            nn.Linear(hidden, 4),   # 4 outputs: (log_a_qf, log_scale_qf, log_a_bf, log_scale_bf)
        )
        # Initialise near the ORPB SAS calibration values:
        #   quickflow: a_qf ~ 0.5 (young, peaked), scale_qf ~ 24 hr (fast)
        #   baseflow:  a_bf ~ 1.26, scale_bf ~ 200 hr (slow, storage-dependent)
        with torch.no_grad():
            self.net[-1].bias.copy_(torch.tensor([
                np.log(0.5),    # log_a_qf
                np.log(24.0),   # log_scale_qf
                np.log(1.26),   # log_a_bf
                np.log(200.0),  # log_scale_bf
            ]))
            self.net[-1].weight.mul_(0.01)

        # Learnable C_old, initialised to the SAS calibration value
        self.c_old = nn.Parameter(torch.tensor(-7.6))

    def forward(self, state, infl_win, cpre_win, old_mask, bf_wt):
        """
        state    : (B, n_state)
        infl_win : (B, MAX_AGE)  influx at lag k
        cpre_win : (B, MAX_AGE)  precip 18O at lag k (zero where unknown)
        old_mask : (B, MAX_AGE)  1 where the lag predates the record
        bf_wt    : (B,)          baseflow fraction of discharge at each obs
        """
        h = self.net(state)
        log_a_qf, log_s_qf = h[:, 0:1], h[:, 1:2]
        log_a_bf, log_s_bf = h[:, 2:3], h[:, 3:4]

        # Convolve each TTD separately against the same precip history
        c_qf = _convolve_ttd(log_a_qf, log_s_qf, infl_win, cpre_win,
                             old_mask, self.c_old)
        c_bf = _convolve_ttd(log_a_bf, log_s_bf, infl_win, cpre_win,
                             old_mask, self.c_old)

        # Flux-weighted combination (same as ORPB_SAS.py gamma-split)
        return (1.0 - bf_wt) * c_qf + bf_wt * c_bf

# =============================================================================
# SECTION 4: K-fold training (same contiguous-block scheme as the pure-NN script)
# =============================================================================
K_FOLDS = 5
fold_size = N_obs // K_FOLDS
fold_indices = [(i*fold_size, (i+1)*fold_size if i < K_FOLDS-1 else N_obs)
                for i in range(K_FOLDS)]

print(f"\nK-fold setup: {K_FOLDS} folds")
for i, (s, e) in enumerate(fold_indices):
    print(f"  Fold {i+1}: samples {s}-{e}  ({times[s].strftime('%Y-%m')} -> {times[e-1].strftime('%Y-%m')})")

# Move full data tensors once
state_t   = torch.tensor(state_z,  device=device)
infl_t    = torch.tensor(infl_win, device=device)
cpre_t    = torch.tensor(cpre_win, device=device)
oldm_t    = torch.tensor(old_mask, device=device)
target_t  = torch.tensor(target,   device=device)
bfwt_t    = torch.tensor(bf_weight[obs_idx], device=device)

N_EPOCHS = 400
PATIENCE = 60
LR       = 5e-3
BATCH    = 128
fold_colors = ['steelblue', 'seagreen', 'darkorange', 'mediumpurple', 'firebrick']

all_t_test, all_y_true, all_y_pred = [], [], []
fold_r2s, fold_rmses = [], []
all_train_losses, all_val_losses = [], []

for fold_idx, (test_start, test_end) in enumerate(fold_indices):
    print(f"\nFold {fold_idx+1}/{K_FOLDS}  test block: "
          f"{times[test_start].strftime('%Y-%m')} -> {times[test_end-1].strftime('%Y-%m')}")

    test_mask = np.zeros(N_obs, dtype=bool)
    test_mask[test_start:test_end] = True
    train_mask = ~test_mask
    train_idx = np.where(train_mask)[0]
    test_idx  = np.where(test_mask)[0]

    train_idx_t = torch.tensor(train_idx, device=device)
    test_idx_t  = torch.tensor(test_idx,  device=device)

    model = HybridSplitSASNet(n_state=state_z.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                           factor=0.5, patience=15)

    best_val = np.inf
    best_state = None
    patience_counter = 0
    train_losses, val_losses = [], []

    for epoch in range(1, N_EPOCHS + 1):
        model.train()
        perm = train_idx_t[torch.randperm(len(train_idx_t), device=device)]
        batch_losses = []
        for b0 in range(0, len(perm), BATCH):
            idx = perm[b0:b0+BATCH]
            pred = model(state_t[idx], infl_t[idx], cpre_t[idx], oldm_t[idx], bfwt_t[idx])
            loss = torch.mean((pred - target_t[idx])**2)
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            batch_losses.append(loss.item())
        train_losses.append(float(np.mean(batch_losses)))

        model.eval()
        with torch.no_grad():
            val_pred = model(state_t[test_idx_t], infl_t[test_idx_t],
                             cpre_t[test_idx_t], oldm_t[test_idx_t],
                             bfwt_t[test_idx_t])
            val_loss = torch.mean((val_pred - target_t[test_idx_t])**2).item()
        val_losses.append(val_loss)
        scheduler.step(val_loss)

        if val_loss < best_val:
            best_val = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
        if patience_counter >= PATIENCE:
            print(f"  Early stop @ epoch {epoch}  best val MSE = {best_val:.4f}")
            break

    all_train_losses.append(train_losses)
    all_val_losses.append(val_losses)

    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        y_pred = model(state_t[test_idx_t], infl_t[test_idx_t],
                       cpre_t[test_idx_t], oldm_t[test_idx_t],
                       bfwt_t[test_idx_t]).cpu().numpy()
    y_true = target[test_idx]

    r2   = r2_score(y_true, y_pred)
    rmse = float(np.sqrt(np.mean((y_pred - y_true)**2)))
    fold_r2s.append(r2)
    fold_rmses.append(rmse)
    print(f"  Fold {fold_idx+1} R^2 = {r2:.3f}  RMSE = {rmse:.3f} per mil  "
          f"C_old learned = {model.c_old.item():.3f}")

    all_t_test.extend(times[test_idx])
    all_y_true.extend(y_true)
    all_y_pred.extend(y_pred)

all_y_true = np.array(all_y_true)
all_y_pred = np.array(all_y_pred)
all_t_test = np.array(all_t_test)
overall_r2   = r2_score(all_y_true, all_y_pred)
overall_rmse = float(np.sqrt(np.mean((all_y_pred - all_y_true)**2)))

print("\n" + "="*60)
print(f"Hybrid NN/SAS k-fold results ({K_FOLDS} folds):")
for i, (r2, rmse) in enumerate(zip(fold_r2s, fold_rmses)):
    print(f"  Fold {i+1}: R^2={r2:.3f}  RMSE={rmse:.3f} per mil")
print(f"  Mean R^2 = {np.mean(fold_r2s):.3f} +/- {np.std(fold_r2s):.3f}")
print(f"  Overall R^2 (all OOF) = {overall_r2:.3f}  RMSE = {overall_rmse:.3f} per mil")

# =============================================================================
# SECTION 5: Plots
# =============================================================================
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

ax = axes[0]
for i, (tl, vl) in enumerate(zip(all_train_losses, all_val_losses)):
    ax.plot(tl, color=fold_colors[i], linestyle='-',  alpha=0.7, label=f'Fold {i+1} train')
    ax.plot(vl, color=fold_colors[i], linestyle='--', alpha=0.7, label=f'Fold {i+1} val')
ax.set_yscale('log'); ax.set_xlabel('Epoch'); ax.set_ylabel('MSE')
ax.set_title('Hybrid NN/SAS — train/val loss per fold')
ax.legend(fontsize=7, ncol=2)

sort_idx = np.argsort(all_t_test)
ax = axes[1]
ax.plot(all_t_test[sort_idx], all_y_true[sort_idx], 'o', color='black',
        markersize=3, alpha=0.7, label='Observed')
ax.plot(all_t_test[sort_idx], all_y_pred[sort_idx], '-', color='firebrick',
        linewidth=1, alpha=0.9, label='Predicted (OOF)')
ax.set_ylabel('ORPB $^{18}$O (per mil)')
ax.set_title(f'Out-of-fold predictions — overall R^2={overall_r2:.3f}, '
             f'RMSE={overall_rmse:.3f} per mil')
ax.legend(fontsize=8)

ax = axes[2]
fold_labels = np.empty(N_obs, dtype=int)
for i, (s, e) in enumerate(fold_indices):
    fold_labels[s:e] = i
lims = [min(all_y_true.min(), all_y_pred.min()) - 0.2,
        max(all_y_true.max(), all_y_pred.max()) + 0.2]
ax.plot(lims, lims, 'k--', linewidth=1)
for i in range(K_FOLDS):
    m = fold_labels == i
    ax.scatter(all_y_true[m], all_y_pred[m], s=15, color=fold_colors[i],
               alpha=0.7, label=f'Fold {i+1} (R^2={fold_r2s[i]:.2f})')
ax.set_xlabel('Observed ORPB $^{18}$O (per mil)')
ax.set_ylabel('Predicted ORPB $^{18}$O (per mil)')
ax.set_title('1:1 predicted vs observed — all folds')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('NN_SAS_hybrid_results.png', dpi=150, bbox_inches='tight')
plt.show()
print("\nFigure saved to NN_SAS_hybrid_results.png")
print("Done.")
