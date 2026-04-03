# Neural Network for testing hydrologic similarity on catchments
# by training with precip isotope concentrations and predicting outflow isotope concentrations.
# Can a complex, abstract NN learn underlying hydrologic similarity between catchments?

# First train NN on one catchment precip data and outflow data, then test on another catchment's precip data to see if it can predict the outflow data.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score

# =============================================================================
# SECTION 1: Load and preprocess data
# =============================================================================
# WHY: The raw CSV has hourly timestamps with precip 18O (continuous, some NaNs)
# and ORPB 18O (sparse — only measured at discrete grab samples).
# We need to:
#   (a) fill NaNs in the input (precip 18O) so the NN always has something to learn from
#   (b) keep track of which timesteps have observed stream 18O for training/evaluation

data_df = pd.read_csv('ORPB_isotope_data.csv', index_col=0, parse_dates=[0])

# Fill NaNs in precip 18O using forward/backward fill (same as ORPB_SAS.py)
# WHY ffill/bfill instead of mean: preserves the temporal isotope signal rather
# than injecting an artificial "average event" into the sequence
data_df['precip 18O'] = data_df['precip 18O'].ffill().bfill()

# Compute influx (rainfall + snowmelt) — same as SAS model
data_df['influx (mm/hr)'] = data_df[['rainfall (mm/hr)', 'snowmelt (mm/hr)']].sum(axis=1)

# Boolean mask: True only at timesteps where stream 18O was actually measured
issample = data_df['ORPB 18O'].notna()

print(f"Total hourly timesteps: {len(data_df)}")
print(f"Stream 18O observations: {issample.sum()}")
print(f"Date range: {data_df.index[0]} to {data_df.index[-1]}")


# =============================================================================
# SECTION 2: Build physics-informed features
# =============================================================================
# WHY pre-computed mixing features instead of raw sequences:
#   The LSTM approach requires the model to discover, from scratch, how far back in time
#   precip 18O matters. With sparse grab-sample observations (~few hundred points), there
#   is not enough labeled data to learn this from raw hourly sequences.
#
#   Instead, we encode the mixing physics directly: compute exponentially-weighted
#   averages (EWAs) of precip 18O at multiple timescales (halflifes). Each EWA is like
#   asking "what was the flux-weighted mean precip 18O over the last N hours?", which
#   mirrors what the SAS TTD does analytically. The NN then learns the right blend of
#   these timescales — a much simpler task requiring far fewer samples.
#
# Halflifes chosen to span event (days) → seasonal (months) timescales.
# influx-weighted: weight precip 18O by influx at each timestep before averaging,
#   so large rain events count more than dry periods with filled NaNs.

TARGET_COL = 'ORPB 18O'

halflifes_hr = [24, 72, 168, 500, 1000, 2000]  # 1d, 3d, 1wk, ~3wk, ~6wk, ~12wk

# Flux-weighted precip 18O: weight isotope signal by how much water fell
# WHY: a tiny drizzle with depleted 18O should move the stream less than a heavy storm
influx = data_df['influx (mm/hr)'].values
precip18O = data_df['precip 18O'].values
flux_weighted = influx * precip18O   # zero when no rain, signal when rain falls

feature_dict = {}
for hl in halflifes_hr:
    # alpha = decay rate per timestep for this halflife
    # EWA formula: y_t = alpha * x_t + (1-alpha) * y_{t-1}
    alpha = 1.0 - np.exp(-np.log(2) / hl)

    # Exponentially weighted average of flux-weighted precip 18O
    ewa_flux   = pd.Series(flux_weighted).ewm(alpha=alpha, adjust=False).mean().values
    ewa_influx = pd.Series(influx).ewm(alpha=alpha, adjust=False).mean().values

    # Normalize by accumulated influx to get a concentration (not a flux)
    # Add small epsilon to avoid division by zero during dry periods
    feature_dict[f'precip18O_ewa_{hl}h'] = ewa_flux / (ewa_influx + 1e-9)

# After computing feature_dict, before building features_df
# Interaction: storage state × isotope signal at each timescale
# WHY: precip 18O only influences stream when new water is moving through
#      (high storage). This mirrors the storage-dependent SAS function.
storage_norm = (data_df['storage (mm)'] - data_df['storage (mm)'].mean()) / data_df['storage (mm)'].std()
for hl in halflifes_hr:
    col = f'precip18O_ewa_{hl}h'
    feature_dict[f'storage_x_ewa_{hl}h'] = storage_norm.values * feature_dict[col]


# Add current discharge and storage as hydrologic state variables
# WHY: these tell the NN how "full" the catchment is, which affects how much
# old vs. new water is in the stream — same role as S_scale in the SAS model
feature_dict['discharge (mm/hr)'] = data_df['discharge (mm/hr)'].values
feature_dict['storage (mm)']      = data_df['storage (mm)'].values

features_df = pd.DataFrame(feature_dict, index=data_df.index)
FEATURE_COLS = list(features_df.columns)
print(f"Features ({len(FEATURE_COLS)}): {FEATURE_COLS}")

target = data_df[TARGET_COL].values



# =============================================================================
# SECTION 3: Normalize features and target
# =============================================================================
# WHY normalize: features span very different scales; normalization prevents any
# one feature from dominating gradient updates during training.

scaler_X = StandardScaler()
scaler_y = StandardScaler()

features_scaled = scaler_X.fit_transform(features_df.values)

# Fit target scaler only on observed stream values (not NaN rows)
obs_vals = target[issample].reshape(-1, 1)
scaler_y.fit(obs_vals)

target_scaled = np.full_like(target, np.nan)
target_scaled[issample] = scaler_y.transform(obs_vals).ravel()



# =============================================================================
# SECTION 4: Build sample arrays at observation timesteps
# =============================================================================
# WHY no sliding window needed now:
#   The temporal memory is already encoded in the EWA features — each row of
#   features_scaled already summarizes the relevant history at that timestep.
#   So we just extract rows where stream 18O was measured.

obs_idx = np.where(issample)[0]   # integer indices of observed timesteps

seq_X = features_scaled[obs_idx].astype(np.float32)   # (N_obs, N_features)
seq_y = target_scaled[obs_idx].astype(np.float32)      # (N_obs,)
seq_t = data_df.index[obs_idx]

print(f"\nSequences available for training/testing: {len(seq_X)}")

# --- Feature-target correlation check ---
import pandas as pd
obs_features = features_df.iloc[obs_idx]   # features only at observation times
obs_target   = data_df[TARGET_COL].iloc[obs_idx]
corr = obs_features.corrwith(obs_target)
print("\nPearson correlation of each feature with ORPB 18O:")
print(corr.sort_values().to_string())

# =============================================================================
# SECTION 5: K-fold cross-validation setup
# =============================================================================
# WHY k-fold instead of a single chronological split:
#   A single 80/20 split happened to put the anomalous 2018 wet year entirely in the
#   test set, so the model was trained on drought dynamics and evaluated on wet dynamics.
#   K-fold rotates the test window through the full record, so every fold's training set
#   contains samples from both wet and dry periods, giving a fairer estimate of
#   generalisation performance.
#
# WHY contiguous time blocks (not random folds):
#   Random assignment would leak future information — a sample from 2019 in the training
#   set would help predict a nearby 2019 sample in the test set, inflating R².
#   Contiguous blocks preserve temporal realism.
#
# K_FOLDS = 5: divides ~2758 samples into 5 blocks of ~550 each (~14 months per fold).
#   Each iteration trains on 4 blocks (~2200 samples) and tests on 1 (~550).

K_FOLDS   = 5
fold_size = len(seq_X) // K_FOLDS
fold_indices = [(i * fold_size, (i + 1) * fold_size if i < K_FOLDS - 1 else len(seq_X))
                for i in range(K_FOLDS)]

print(f"\nK-fold setup: {K_FOLDS} folds")
for i, (start, end) in enumerate(fold_indices):
    print(f"  Fold {i+1}: samples {start}–{end}  ({seq_t[start].strftime('%Y-%m')} to {seq_t[end-1].strftime('%Y-%m')})")


# =============================================================================
# SECTION 6: Dataset class and hyperparameters
# =============================================================================

class IsotopeDataset(Dataset):
    """Wraps numpy arrays into a PyTorch-compatible dataset."""
    def __init__(self, X, y):
        self.X = torch.tensor(X)
        self.y = torch.tensor(y)

    def __len__(self):
        return len(self.y)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

BATCH_SIZE  = 32
N_FEATURES  = len(FEATURE_COLS)
HIDDEN_SIZE = 32
N_EPOCHS    = 200
PATIENCE    = 40
criterion   = nn.MSELoss()
device      = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"\nUsing device: {device}")


# =============================================================================
# SECTION 7: Define the feedforward Neural Network
# =============================================================================
# WHY feedforward (not LSTM):
#   The temporal memory is already encoded in the EWA features — each row is a
#   flat vector summarising the relevant history. A feedforward network just needs
#   to learn the right weighted combination of those pre-computed mixing signals.
#
# Architecture: Input → Linear(32) → ReLU → Dropout(0.3) → Linear(16) → ReLU → Linear(1)
# Dropout(0.3): randomly zeros 30% of neurons during training to prevent overfitting.

class FeedforwardIsotopeNet(nn.Module):
    def __init__(self, n_features, hidden_size=32, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, 16),
            nn.ReLU(),
            nn.Linear(16, 1)
        )

    def forward(self, x):
        return self.net(x).squeeze(-1)   # (batch,)


# =============================================================================
# SECTION 8-10: K-fold training and evaluation loop
# =============================================================================
# Each fold: hold out one contiguous time block as the test set, train on the rest.
# WHY reinitialise the model each fold: we want to measure how well the architecture
#   generalises, not how well a pre-trained model fine-tunes. Fresh weights each fold
#   gives an honest estimate.
# Collect all out-of-fold predictions so we can evaluate across the full record.

fold_colors = ['steelblue', 'seagreen', 'darkorange', 'mediumpurple', 'firebrick']

# Arrays to accumulate out-of-fold results across all folds
all_t_test    = []
all_y_true    = []
all_y_pred    = []
fold_r2s      = []
fold_rmses    = []
all_train_losses = []   # loss curves from every fold (for plotting)
all_val_losses   = []
fold_models   = []      # best weights from each fold, used for ensemble transfer

print("\n" + "="*60)
for fold_idx, (test_start, test_end) in enumerate(fold_indices):
    print(f"\nFold {fold_idx+1}/{K_FOLDS}  —  test block: {seq_t[test_start].strftime('%Y-%m')} to {seq_t[test_end-1].strftime('%Y-%m')}")

    # ---- Split this fold ----
    test_mask  = np.zeros(len(seq_X), dtype=bool)
    test_mask[test_start:test_end] = True
    train_mask = ~test_mask

    X_train_f = seq_X[train_mask]
    y_train_f = seq_y[train_mask]
    X_test_f  = seq_X[test_mask]
    y_test_f  = seq_y[test_mask]
    t_test_f  = seq_t[test_mask]

    train_loader = DataLoader(IsotopeDataset(X_train_f, y_train_f), batch_size=BATCH_SIZE, shuffle=True)
    test_loader  = DataLoader(IsotopeDataset(X_test_f,  y_test_f),  batch_size=BATCH_SIZE, shuffle=False)

    # ---- Fresh model for each fold ----
    model     = FeedforwardIsotopeNet(N_FEATURES, hidden_size=HIDDEN_SIZE).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=10)

    train_losses  = []
    val_losses    = []
    best_val_loss = np.inf
    patience_counter = 0
    best_model_state = None

    # ---- Training loop ----
    for epoch in range(1, N_EPOCHS + 1):
        model.train()
        batch_losses = []
        for X_batch, y_batch in train_loader:
            optimizer.zero_grad()
            pred = model(X_batch.to(device))
            loss = criterion(pred, y_batch.to(device))
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            batch_losses.append(loss.item())

        train_loss = np.mean(batch_losses)
        train_losses.append(train_loss)

        model.eval()
        with torch.no_grad():
            val_preds, val_targets = [], []
            for X_batch, y_batch in test_loader:
                val_preds.extend(model(X_batch.to(device)).cpu().numpy())
                val_targets.extend(y_batch.numpy())
            val_loss = np.mean((np.array(val_preds) - np.array(val_targets))**2)
        val_losses.append(val_loss)

        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = {k: v.clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1

        if patience_counter >= PATIENCE:
            print(f"  Early stopping at epoch {epoch}  (best val loss: {best_val_loss:.4f})")
            break

    all_train_losses.append(train_losses)
    all_val_losses.append(val_losses)

    # ---- Evaluate this fold ----
    model.load_state_dict(best_model_state)
    model.eval()
    with torch.no_grad():
        y_pred_scaled = model(torch.tensor(X_test_f).to(device)).cpu().numpy()

    y_pred = scaler_y.inverse_transform(y_pred_scaled.reshape(-1, 1)).ravel()
    y_true = scaler_y.inverse_transform(y_test_f.reshape(-1, 1)).ravel()

    fold_r2   = r2_score(y_true, y_pred)
    fold_rmse = np.sqrt(np.mean((y_pred - y_true)**2))
    fold_r2s.append(fold_r2)
    fold_rmses.append(fold_rmse)
    print(f"  Fold {fold_idx+1} R² = {fold_r2:.3f}  |  RMSE = {fold_rmse:.3f} ‰")

    all_t_test.extend(t_test_f)
    all_y_true.extend(y_true)
    all_y_pred.extend(y_pred)

    # Save this fold's trained model for use in the ensemble
    fold_models.append({k: v.clone() for k, v in best_model_state.items()})

# Overall metrics across all out-of-fold predictions
all_y_true = np.array(all_y_true)
all_y_pred = np.array(all_y_pred)
all_t_test = np.array(all_t_test)
overall_r2   = r2_score(all_y_true, all_y_pred)
overall_rmse = np.sqrt(np.mean((all_y_pred - all_y_true)**2))

print(f"\n{'='*60}")
print(f"K-fold results ({K_FOLDS} folds):")
for i, (r2, rmse) in enumerate(zip(fold_r2s, fold_rmses)):
    print(f"  Fold {i+1}: R²={r2:.3f}  RMSE={rmse:.3f} ‰")
print(f"  Mean R² = {np.mean(fold_r2s):.3f} ± {np.std(fold_r2s):.3f}")
print(f"  Overall R² (all OOF predictions) = {overall_r2:.3f}  RMSE = {overall_rmse:.3f} ‰")


# --- Linear regression baseline (same k-fold splits) ---
# WHY: if linear regression also scores ~0, the features genuinely don't predict
# stream 18O and no model will fix that. If linear scores higher than the NN,
# the NN is underfitting (dropout too strong, architecture too simple).
from sklearn.linear_model import LinearRegression

lr = LinearRegression()
lr_preds = np.empty(len(seq_X))
for test_start, test_end in fold_indices:
    test_mask = np.zeros(len(seq_X), dtype=bool)
    test_mask[test_start:test_end] = True
    lr.fit(seq_X[~test_mask], seq_y[~test_mask])
    lr_preds[test_mask] = lr.predict(seq_X[test_mask])

lr_true = scaler_y.inverse_transform(seq_y.reshape(-1, 1)).ravel()
lr_pred = scaler_y.inverse_transform(lr_preds.reshape(-1, 1)).ravel()
print(f"\nLinear regression k-fold R² = {r2_score(lr_true, lr_pred):.3f}")


# =============================================================================
# SECTION 11: Plot results
# =============================================================================
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# --- Panel 1: Loss curves for each fold ---
ax = axes[0]
for i, (tl, vl) in enumerate(zip(all_train_losses, all_val_losses)):
    ax.plot(tl, color=fold_colors[i], linestyle='-',  alpha=0.7, label=f'Fold {i+1} train')
    ax.plot(vl, color=fold_colors[i], linestyle='--', alpha=0.7, label=f'Fold {i+1} val')
ax.set_xlabel('Epoch')
ax.set_ylabel('MSE (normalized)')
ax.set_title('Training and Validation Loss — all folds')
ax.legend(fontsize=7, ncol=2)
ax.set_yscale('log')

# --- Panel 2: Full timeseries of out-of-fold predictions ---
# Sort by time so the line plot connects chronologically
sort_idx = np.argsort(all_t_test)
t_sorted    = all_t_test[sort_idx]
true_sorted = all_y_true[sort_idx]
pred_sorted = all_y_pred[sort_idx]

ax = axes[1]
ax.plot(t_sorted, true_sorted, 'o', color='black',     markersize=3, label='Observed',   alpha=0.7)
ax.plot(t_sorted, pred_sorted, '-', color='firebrick',  linewidth=1,  label='Predicted (OOF)', alpha=0.9)
ax.set_ylabel('ORPB $^{18}$O (‰)')
ax.set_title(f'Out-of-fold predictions across full record  (Overall R²={overall_r2:.3f}, RMSE={overall_rmse:.3f} ‰)')
ax.legend(fontsize=8)

# --- Panel 3: 1:1 scatter coloured by fold ---
ax = axes[2]
sort_idx_orig = np.argsort(all_t_test)   # reuse time sort to assign fold colours
# Recover which fold each point belongs to
fold_labels = np.empty(len(all_t_test), dtype=int)
for i, (start, end) in enumerate(fold_indices):
    fold_labels[start:end] = i
fold_labels_sorted = fold_labels[sort_idx]

lims = [min(all_y_true.min(), all_y_pred.min()) - 0.2,
        max(all_y_true.max(), all_y_pred.max()) + 0.2]
ax.plot(lims, lims, 'k--', linewidth=1)
for i in range(K_FOLDS):
    mask = fold_labels == i
    ax.scatter(all_y_true[mask], all_y_pred[mask], s=15,
               color=fold_colors[i], alpha=0.7, label=f'Fold {i+1} (R²={fold_r2s[i]:.2f})')
ax.set_xlabel('Observed ORPB $^{18}$O (‰)')
ax.set_ylabel('Predicted ORPB $^{18}$O (‰)')
ax.set_title(f'1:1 Predicted vs. Observed — all folds (Overall R²={overall_r2:.3f})')
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig('NN_ORPB_results.png', dpi=150, bbox_inches='tight')
plt.show()
print("Figure saved to NN_ORPB_results.png")


# =============================================================================
# SECTION 12: Save ensemble of fold models for catchment transfer
# =============================================================================
# WHY save all 5 fold models instead of retraining on all data:
#   Each fold model was trained on a different 80% slice of the record, so together
#   they have all seen all hydrologic regimes (drought 2015, wet 2018, etc.).
#   To predict on a new catchment, run all 5 models and average their outputs —
#   this ensemble is more robust than any single model because disagreement between
#   folds flags uncertain predictions.
#
# To apply to another catchment:
#   1. Load this file and reconstruct 5 FeedforwardIsotopeNet instances
#   2. Build EWA features the same way (same halflifes_hr) for the new catchment
#   3. Normalise with scaler_X_mean / scaler_X_scale
#   4. Run each model in eval() + no_grad(), average the 5 outputs
#   5. Inverse-transform with scaler_y_mean / scaler_y_scale

torch.save({
    'fold_model_states': fold_models,       # list of 5 state_dicts
    'scaler_X_mean':     scaler_X.mean_,
    'scaler_X_scale':    scaler_X.scale_,
    'scaler_y_mean':     scaler_y.mean_,
    'scaler_y_scale':    scaler_y.scale_,
    'feature_cols':      FEATURE_COLS,
    'halflifes_hr':      halflifes_hr,
    'hidden_size':       HIDDEN_SIZE,
    'n_features':        N_FEATURES,
    'kfold_mean_r2':     float(np.mean(fold_r2s)),
    'kfold_overall_r2':  float(overall_r2),
    'kfold_overall_rmse': float(overall_rmse),
}, 'NN_ORPB_model.pt')
print("Ensemble of 5 fold models saved to NN_ORPB_model.pt")
print("\nDone.")
