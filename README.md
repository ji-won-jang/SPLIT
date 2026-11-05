# High resolution measures of preference-based local segregation for the United States, 2018-2023

This repo contains the code used to develop the Spatial Preference in Local Interaction (SPLIT).

The data is archived and downloadable at Zenodo: https://doi.org/10.5281/zenodo.23671422

## Procedure

`split.py` shows how SPLIT is computed from an origin–destination (OD) flow table. It needs two inputs:

- **Flows**: `Origin`, `Destination`, `Frequency`, the annual flow between two tracts. In the dataset, monthly Advan device visits are rescaled to population (`Frequency × Total Population / devices residing` of the origin tract) and averaged over the months in which each OD pair is observed.
- **Tracts**: `GEOID`, centroid coordinates `x`, `y` (km), and composition shares (racial/ethnic shares and household-income shares across five ranges).

For each tract *i* and buffer radius *s*, with `B_s(i)` the other tracts whose centroid lies within *s* km of *i*:

1. **Outflow share**: outflows from *i* to `B_s(i)`, divided by all outflows of *i* (self-loop included).
2. **Unadjusted SPLIT (S_obs)**: outflow-weighted mean cosine similarity between the composition of *i* and that of its destinations in `B_s(i)`.
3. **Structural similarity (S_exp)**: unweighted mean cosine similarity between *i* and every tract in `B_s(i)`.
4. **SPLIT** = (S_obs − S_exp) / (1 − S_exp), set to 0 when S_exp ≥ 0.9999.

To run the example on a small hypothetical dataset:

```bash
conda env create -f environment.yml
conda activate split-metrics
python split.py
```

[`example/example.ipynb`](example/example.ipynb) walks through each step by hand for one tract and checks the result against `compute_split`.

### Output columns

`{s}` is the buffer radius in km (5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200).

| Column | Metric | Description |
|---|---|---|
| `GEOID` | Identifier | 11-digit FIPS code of the residential (origin) census tract. |
| `out_shareP_{s}` | Outflow share | Rescaled outflows to other tracts whose centroid lies within *s* km, divided by the tract's total outflows (including the self-loop). -1 if undefined. |
| `unadjusted_race_{s}` | Unadjusted SPLIT (race/ethnicity) | S_obs: outflow-weighted mean cosine similarity of racial/ethnic composition between the tract and its destinations within *s* km (self-loop excluded). -1 if undefined. |
| `unadjusted_income_{s}` | Unadjusted SPLIT (income) | S_obs: outflow-weighted mean cosine similarity of household income composition between the tract and its destinations within *s* km (self-loop excluded). -1 if undefined. |
| `str_sim_race_{s}` | Structural similarity (race/ethnicity) | S_exp: unweighted mean cosine similarity of racial/ethnic composition between the tract and all other tracts within *s* km, independent of mobility. NaN if undefined. |
| `str_sim_income_{s}` | Structural similarity (income) | S_exp: unweighted mean cosine similarity of household income composition between the tract and all other tracts within *s* km. NaN if undefined. |
| `SPLIT_Race_{s}` | SPLIT (race/ethnicity) | (S_obs − S_exp) / (1 − S_exp); 0 if S_exp ≥ 0.9999. NaN if S_obs or S_exp is undefined. |
| `SPLIT_Income_{s}` | SPLIT (income) | (S_obs − S_exp) / (1 − S_exp); 0 if S_exp ≥ 0.9999. NaN if S_obs or S_exp is undefined. |
