# American Gut: diet → gut → health mediation prototype (Route B, fully public)

Does diet quality act on BMI **through** the gut microbiome? This is the fully-open prototype of
[`docs/diet-gut-health-landscape.md`](../../docs/diet-gut-health-landscape.md) Route B, run on the one
public cohort that carries detailed diet + microbiome + a health outcome with **no application** —
[American Gut / The Microsetta Initiative](https://journals.asm.org/doi/10.1128/msystems.00031-18).

        diet quality (HEI-2010)  →  gut alpha diversity (Shannon)  →  BMI        (adjusted age + sex)

It runs the counterfactual mediation in `vmbdiet.mediation`, decomposing diet's total effect on BMI
into the part flowing **through** diversity (natural indirect effect, NIE) and the rest (direct, NDE).

## Data — one public file, no login

American Gut's manuscript package ships, per rarefaction depth, a single mapping file
(`ag_map_with_alpha.txt`) that already has the **VioScreen FFQ appended** (HEI-2010 total +
components + nutrients), **alpha diversity** (Shannon/Faith's PD/observed OTUs), and the survey
metadata (BMI, age, sex, a `healthy_subset` flag) — one fecal sample per participant. So the whole
diet+microbiome+health table is one ~33 MB download.

```bash
python fetch_americangut.py                 # -> agp_mediation_input.tsv (depth 10000)
python analyze.py --data agp_mediation_input.tsv
python analyze.py --data agp_mediation_input.tsv --healthy-only   # AGP healthy subset sensitivity
python analyze.py                            # no --data: runs the synthetic demo (make_demo.py)
```

Source: `https://ftp.microbio.me/AmericanGut/manuscript-package/<depth>/` (CC-licensed public data).
The distilled input and the demo data are git-ignored (reproducible via the scripts, not
redistributed); the synthetic demo is validated by `tests/test_american_gut_mediation.py`.

## Real-data run (depth 10000, n=839 with HEI + Shannon + BMI + age + sex) — result

| quantity | estimate (per +1 SD HEI) | 95% CI |
|---|---|---|
| path a (HEI → Shannon) | +0.094 | — |
| path b (Shannon → BMI \| HEI) | −0.382 | — |
| **total** effect (HEI → BMI) | **−0.715** | [−1.058, −0.379] |
| direct (NDE) | −0.680 | [−1.031, −0.342] |
| indirect (NIE, via diversity) | −0.036 | **[−0.089, +0.001]** |
| proportion mediated | 5.0% | [−0.001, +0.146] |

**Conclusion.** Better diet quality is robustly associated with **lower BMI** (total effect clearly
negative) — but **gut alpha diversity mediates essentially none of it** (~5%, NIE 95% CI includes 0).
The diet→BMI association is almost entirely *direct*. That is biologically sensible: a single scalar
diversity index is too coarse to carry the diet→adiposity mechanism, which would more plausibly run
through specific taxa or metabolic **function** (SCFA, fiber fermentation) than through α-diversity.
It also matches the landscape caveat — 16S gives a measured but coarse mediator; the mechanistic
middle wants measured function (shotgun/HUMAnN).

## Honest limits (why this is a prototype, not established causation)

- **Cross-sectional**, one sample per participant → the order diet→diversity→BMI is *assumed*, not
  shown; the estimate is associational (see `vmbdiet.mediation`'s assumptions block).
- **Self-reported** BMI and diet; **16S** → diversity is the mediator, not measured function.
- American Gut is a **WEIRD, self-selected** citizen-science cohort; the VioScreen FFQ is a **subset**
  (~900 of the depth-10000 samples here).

## Natural next iteration

The same manuscript package ships a `picrust/` folder (PICRUSt-inferred KEGG function, collapsed
L1–L3) for each depth. Swapping the mediator from Shannon to an **inferred SCFA / fiber-fermentation
pathway** gives a function-based mediation on the same public data — still inferred (not HUMAnN), but
a more mechanistic middle node than diversity. The measured-function version needs shotgun + diet +
outcome together (MetaCardis/PREDICT, controlled) or FRESH.
