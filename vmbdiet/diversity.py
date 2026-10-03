"""vmbdiet.diversity — beta-diversity distances + ordination for the diet↔community analyses."""
from __future__ import annotations
import numpy as np
import pandas as pd


def _relabund(abund: pd.DataFrame) -> pd.DataFrame:
    tot = abund.sum(axis=1).replace(0, np.nan)
    return abund.div(tot, axis=0).fillna(0.0)


def bray_curtis(abund: pd.DataFrame) -> pd.DataFrame:
    """Pairwise Bray-Curtis dissimilarity (samples × samples)."""
    rel = _relabund(abund); X = rel.to_numpy(); n = len(rel)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(X[i] - X[j]).sum(); den = (X[i] + X[j]).sum()
            d = num / den if den > 0 else 0.0
            D[i, j] = D[j, i] = d
    return pd.DataFrame(D, index=rel.index, columns=rel.index)


def shannon(abund: pd.DataFrame) -> pd.Series:
    """Shannon alpha-diversity per sample."""
    rel = _relabund(abund)
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -(rel * np.log(rel.where(rel > 0))).sum(axis=1)
    return H.rename("shannon")


def pcoa(dist: pd.DataFrame, n_axes: int = 2):
    """Classical MDS (principal coordinates) of a distance matrix.
    Returns (coords DataFrame [sample × PCo1..], proportion_explained array)."""
    D = dist.to_numpy(); n = len(D)
    A = -0.5 * D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = J @ A @ J
    vals, vecs = np.linalg.eigh(B)
    order = np.argsort(vals)[::-1]
    vals = vals[order]; vecs = vecs[:, order]
    pos = vals > 0
    coords = vecs[:, pos] * np.sqrt(vals[pos])
    k = min(n_axes, coords.shape[1])
    cols = [f"PCo{i+1}" for i in range(k)]
    prop = (vals[pos] / vals[pos].sum())[:k]
    return pd.DataFrame(coords[:, :k], index=dist.index, columns=cols), prop
