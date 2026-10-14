"""
vmbdiet.cst — vaginal community state type (CST) assignment from 16S relative abundance.

The field standard is VALENCIA (France et al. 2020), which assigns each sample to a CST by the
Yue-Clayton similarity of its taxon vector to a set of reference centroids. This module implements:

  * `yue_clayton` — the same similarity VALENCIA uses, so a user-supplied VALENCIA centroid table
    can be scored directly (`assign_cst_valencia`);
  * `assign_cst` — a transparent, dependency-free dominant-taxon classifier (the widely used
    "dominant Lactobacillus > threshold" scheme) for when centroids aren't on hand.

CST scheme: I = L. crispatus, II = L. gasseri, III = L. iners, V = L. jensenii,
IV = diverse / anaerobic (Gardnerella, Prevotella, etc.) — the BV-associated, Lactobacillus-depleted state.
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd

_CENTROID_PATH = os.path.join(os.path.dirname(__file__), "data", "CST_centroids.csv")


def load_valencia_centroids() -> pd.DataFrame:
    """Load the packaged VALENCIA reference centroids (sub_CST × taxa). Source: ravel-lab/VALENCIA (MIT)."""
    return pd.read_csv(_CENTROID_PATH, index_col=0)

DOM_THRESHOLD = 0.50
# species keyword -> CST (dominant-taxon rule)
LACTO_CST = {"crispatus": "I", "gasseri": "II", "iners": "III", "jensenii": "V"}


def _relabund(abund: pd.DataFrame) -> pd.DataFrame:
    """Row-normalize counts/abundances to fractions (rows = samples, cols = taxa)."""
    tot = abund.sum(axis=1).replace(0, np.nan)
    return abund.div(tot, axis=0).fillna(0.0)


def yue_clayton(p: np.ndarray, q: np.ndarray) -> float:
    """Yue & Clayton theta similarity between two relative-abundance vectors (VALENCIA's metric)."""
    p = np.asarray(p, float); q = np.asarray(q, float)
    dot = float((p * q).sum())
    denom = float((p * p).sum() + (q * q).sum() - dot)
    return dot / denom if denom > 0 else 0.0


def lactobacillus_fraction(abund: pd.DataFrame) -> pd.Series:
    """Total Lactobacillus relative abundance per sample (by column-name match)."""
    rel = _relabund(abund)
    lac = [c for c in rel.columns if "lactobacillus" in c.lower()
           or any(k in c.lower() for k in LACTO_CST)]
    return rel[lac].sum(axis=1) if lac else pd.Series(0.0, index=rel.index)


def assign_cst(abund: pd.DataFrame, dom_threshold: float = DOM_THRESHOLD) -> pd.DataFrame:
    """Dominant-taxon CST per sample.

    abund: DataFrame (samples × taxa), counts or relative abundances.
    Returns: DataFrame indexed by sample with columns
      cst, subcst, dominant_taxon, dominant_frac, lacto_frac, lacto_dominant.
    """
    rel = _relabund(abund)
    lacfrac = lactobacillus_fraction(abund)
    rows = []
    for s in rel.index:
        v = rel.loc[s]
        dom = v.idxmax(); domf = float(v.max())
        name = str(dom).lower()
        cst = "IV"; sub = "IV"
        if domf >= dom_threshold:
            hit = next((LACTO_CST[k] for k in LACTO_CST if k in name), None)
            if hit:
                cst = hit; sub = {"I": "I", "II": "II", "III": "III", "V": "V"}[hit]
            elif "lactobacillus" in name:
                cst = "III"; sub = "III-like"   # unresolved Lactobacillus dominance
        rows.append(dict(sample=s, cst=cst, subcst=sub, dominant_taxon=str(dom),
                         dominant_frac=round(domf, 4), lacto_frac=round(float(lacfrac[s]), 4),
                         lacto_dominant=bool(lacfrac[s] >= dom_threshold)))
    return pd.DataFrame(rows).set_index("sample")


def assign_cst_valencia(abund: pd.DataFrame, centroids: pd.DataFrame) -> pd.DataFrame:
    """CST by nearest VALENCIA centroid (max Yue-Clayton).

    centroids: DataFrame (CST label × taxa), the published VALENCIA reference centroids.
    Taxa are aligned on shared column names; both are row-normalized first.
    """
    rel = _relabund(abund)
    cen = _relabund(centroids)
    taxa = [t for t in rel.columns if t in cen.columns]
    if not taxa:
        raise ValueError("no shared taxa between samples and centroids — harmonize taxonomy names first")
    R = rel[taxa].to_numpy(); C = cen[taxa].to_numpy()
    rows = []
    for i, s in enumerate(rel.index):
        sims = [yue_clayton(R[i], C[j]) for j in range(len(cen))]
        j = int(np.argmax(sims))
        sub = str(cen.index[j])
        rows.append(dict(sample=s, subcst=sub, cst=sub.split("-")[0].rstrip("ABl"),
                         score=round(float(sims[j]), 4)))
    return pd.DataFrame(rows).set_index("sample")
