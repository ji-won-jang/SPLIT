"""Spatial Preference in Local Interaction (SPLIT) from an origin-destination flow table.

Run `python split.py` to compute SPLIT on a small hypothetical dataset.
"""

import numpy as np
import pandas as pd

RADII_KM = [5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200]
NO_ALTERNATIVE = 0.9999  # SPLIT is set to 0 when S_exp is at or above this value


def cosine_similarity_matrix(X):
    """Pairwise cosine similarity between the rows of X (NaN for all-zero rows)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        U = X / np.linalg.norm(X, axis=1, keepdims=True)
    return U @ U.T


def compute_split(flows, tracts, attributes, radii_km=RADII_KM):
    """Outflow share, unadjusted SPLIT (S_obs), structural similarity (S_exp) and SPLIT.

    flows      : DataFrame with columns Origin, Destination, Frequency. Frequency is the
                 annual population-rescaled flow; self-loops (Origin == Destination) are allowed.
    tracts     : DataFrame indexed by GEOID with centroid coordinates x, y (km) and the
                 composition shares listed in `attributes`.
    attributes : {name: [composition columns]}, e.g. {"race": [...], "income": [...]}.
    radii_km   : buffer radii s.

    Returns a DataFrame with one row per tract. Undefined values (tract without outflows)
    are -1 for the outflow share and S_obs, and NaN for S_exp and SPLIT.
    """
    n = len(tracts)
    pos = pd.Series(np.arange(n), index=tracts.index)
    o = pos.loc[flows["Origin"]].to_numpy()
    d = pos.loc[flows["Destination"]].to_numpy()
    w = flows["Frequency"].to_numpy(dtype=float)

    xy = tracts[["x", "y"]].to_numpy(dtype=float)
    dist = np.linalg.norm(xy[:, None, :] - xy[None, :, :], axis=-1)
    similarity = {name: cosine_similarity_matrix(tracts[cols].to_numpy(dtype=float))
                  for name, cols in attributes.items()}

    total_out = np.bincount(o, weights=w, minlength=n)  # all outflows, self-loop included
    cols = {}
    for s in radii_km:
        # B_s(i): other tracts whose centroid lies within s km of tract i.
        in_buffer = dist <= s
        np.fill_diagonal(in_buffer, False)
        flow_in_buffer = in_buffer[o, d]
        out_in_buffer = np.bincount(o, weights=np.where(flow_in_buffer, w, 0.0), minlength=n)

        with np.errstate(invalid="ignore", divide="ignore"):
            cols[f"out_shareP_{s}"] = np.where(total_out > 0, out_in_buffer / total_out, -1)

            for name, S in similarity.items():
                # S_obs: outflow-weighted mean similarity to destinations within the buffer.
                weighted = np.where(flow_in_buffer, w * S[o, d], 0.0)
                s_obs = np.bincount(o, weights=weighted, minlength=n) / out_in_buffer
                s_obs = np.where(out_in_buffer > 0, s_obs, -1)

                # S_exp: unweighted mean similarity to every tract within the buffer.
                s_exp = np.where(in_buffer, S, 0.0).sum(axis=1) / in_buffer.sum(axis=1)

                split = (s_obs - s_exp) / (1 - s_exp)
                split = np.where(s_exp >= NO_ALTERNATIVE, 0.0, split)
                split = np.where(s_obs == -1, np.nan, split)

                cols[f"unadjusted_{name}_{s}"] = s_obs
                cols[f"str_sim_{name}_{s}"] = s_exp
                cols[f"SPLIT_{name.capitalize()}_{s}"] = split

    order = [f"out_shareP_{s}" for s in radii_km]
    for template in ("unadjusted_{n}_{s}", "str_sim_{n}_{s}", "SPLIT_{N}_{s}"):
        for name in attributes:
            order += [template.format(n=name, N=name.capitalize(), s=s) for s in radii_km]
    return pd.DataFrame(cols, index=tracts.index)[order]


if __name__ == "__main__":
    # Hypothetical tracts: centroid coordinates (km) and composition shares.
    tracts = pd.DataFrame(
        {
            "x": [0, 6, 0, 30, 40, 90],
            "y": [0, 0, 9, 0, 28, 0],
            "race_a": [0.8, 0.7, 0.1, 0.2, 0.6, 0.3],
            "race_b": [0.1, 0.2, 0.8, 0.2, 0.3, 0.3],
            "race_c": [0.1, 0.1, 0.1, 0.6, 0.1, 0.4],
            "inc_low": [0.2, 0.3, 0.6, 0.4, 0.1, 0.3],
            "inc_mid": [0.5, 0.5, 0.3, 0.4, 0.4, 0.4],
            "inc_high": [0.3, 0.2, 0.1, 0.2, 0.5, 0.3],
        },
        index=pd.Index(list("ABCDEF"), name="GEOID"),
    )
    attributes = {"race": ["race_a", "race_b", "race_c"],
                  "income": ["inc_low", "inc_mid", "inc_high"]}

    # Hypothetical annual OD flows, expressed as shares of each origin's outflows (F has none).
    flows = pd.DataFrame(
        [("A", "A", 0.50), ("A", "B", 0.30), ("A", "C", 0.05), ("A", "E", 0.15),
         ("B", "B", 0.40), ("B", "A", 0.40), ("B", "D", 0.20),
         ("C", "C", 0.60), ("C", "A", 0.30), ("C", "F", 0.10),
         ("D", "D", 0.50), ("D", "B", 0.25), ("D", "E", 0.25),
         ("E", "E", 0.40), ("E", "A", 0.45), ("E", "D", 0.15)],
        columns=["Origin", "Destination", "Frequency"],
    )

    result = compute_split(flows, tracts, attributes, radii_km=[10, 50, 100])
    with pd.option_context("display.width", 120, "display.max_columns", None):
        print(result.round(3).T)
