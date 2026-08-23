# Microbiome-Differential-Abundance-Benchmark

Which microbiome differential abundance method should you trust? Three tools, three
different statistical philosophies, tested against simulated data where the answer is
known in advance.

> Learning project. Simulated data, meant to be read and re-run.

---

## The problem

Microbiome differential abundance methods frequently disagree on which taxa are
associated with a phenotype. A comparison across 38 real datasets found substantial
disagreement among 14 commonly used methods (Nearing et al. 2022, *Nature
Communications*). Because the true signal in real data is unknown, there is no way to
score which method was right — disagreement alone doesn't say who is wrong.

Simulation solves that by inverting the problem: generate data where the
differentially-abundant taxa are planted on purpose, then check which methods recover
them.

---

## The three methods, and why these three

| Tool | Normalisation | Transform | Distribution | Compositionally aware |
|---|---|---|---|---|
| DESeq2 | RLE | None | Negative binomial | No |
| MaAsLin2 | TSS | AST (log) | Normal | No |
| ANCOM-BC2 | Bias correction | — | — | **Yes** |

Not three arbitrary popular names — three different answers to the same design
question. DESeq2 models count data directly and serves here as the
**compositionally-naive baseline**. MaAsLin2 normalises to proportions and assumes
approximate normality after transformation — a middle ground. ANCOM-BC2 corrects for
compositional bias explicitly and is the closest thing here to a principled treatment
of the constraint.

**Not included, and why:** the wider comparison (Nearing et al.'s 14 methods, reproduced
in part below) exists to map an entire ecosystem. That isn't this repo's goal, and
running all 14 trades feasibility for comprehensiveness neither this analysis nor one
person's time budget needs. ALDEx2, Corncob and metagenomeSeq are the natural next
additions — CLR-based, beta-binomial, and zero-inflated respectively — if this is
extended later.

---

## Simulation design

Synthetic count tables generated with **SparseDOSSA2**, fit to a public template (Human
Microbiome Project) rather than to any of my own unpublished data, so the noise
structure — sparsity, dispersion, taxon correlation — matches real microbiome data
rather than an arbitrary distribution. A naive per-taxon negative binomial simulation,
generated independently of any real template, would lack this structure and make every
method look better than it would on real data.

**Ground truth** is planted on top of the fitted background: a known subset of taxa are
assigned a genuine effect size in one simulated group, everything else carries no signal
beyond the noise model.

**Grid:**

| Axis | Levels |
|---|---|
| Sample size | small, large |
| Effect size | weak, strong |
| Signal sparsity | few true positives, many true positives |

Every combination is simulated independently and run through all three methods, so the
result is a performance surface rather than a single number per tool.

---

## Evaluation

For each method, on each scenario, against the known ground truth:

- **Sensitivity** — proportion of true positives recovered
- **Precision** — proportion of reported hits that are real
- **False discovery rate** — reported vs. the method's own nominal FDR threshold, which
  is the specific failure Hawinkel et al. (2019, *Briefings in Bioinformatics*) found in
  several DA methods on simulated data
- **Runtime** — wall-clock time per method per scenario

FDR is reported against the *nominal* threshold deliberately: a method can have
reasonable precision while still failing to control FDR at the level it claims to,
which is a different and more specific failure than "too many false positives" in the
abstract.

---

## Layout

```text
Microbiome-Differential-Abundance-Benchmark/
│
├── code/
│   ├── 01_fit_template.R      # SparseDOSSA2 fit to HMP
│   ├── 02_simulate_data.R     # grid generation, ground truth saved alongside each table
│   ├── 03_run_methods.R       # DESeq2, MaAsLin2, ANCOM-BC2 on every simulated dataset
│   ├── 04_evaluate.R          # sensitivity, precision, FDR, runtime vs ground truth
│   └── util_functions.R
│
├── data/                      # simulated tables + ground truth (RDS), per scenario
├── figures/
├── results/                   # per-method, per-scenario metric tables
│
└── README.md
```

---

## Running it

```r
BiocManager::install(c("SparseDOSSA2", "DESeq2", "ANCOMBC", "Maaslin2"))
```

Run the four scripts in `code/` in order from the repo root. `01_fit_template.R` fits
once and caches the model; the remaining scripts read from `data/` and `results/`
relative to the repo root.

---

## What isn't here

- **The other 11 methods** from the full DA ecosystem — see above.
- **Real-data validation.** A method that wins on simulation can still disappoint on
  real data, because no simulator fully captures actual biological structure. This
  repo's conclusions are about simulated data; a companion real-dataset check would be
  the natural next step, not a substitute for one.
- **Covariate adjustment and random effects.** All three tools support them; this
  benchmark runs the simplest two-group comparison each was designed for, so the
  results isolate the core distributional assumption rather than covariate handling.

---

## References

- Nearing JT et al. (2022). Microbiome differential abundance methods produce different
  results across 38 datasets. *Nature Communications* 13:342.
- Hawinkel S et al. (2019). A broken promise: microbiome differential abundance methods
  do not control the false discovery rate. *Briefings in Bioinformatics* 20:210–221.
- Ma S et al. (2022). Population structure discovery in meta-analyzed microbial
  communities and inflammatory bowel disease using SparseDOSSA2. *Genome Biology*
  23:95.
- Love MI, Huber W & Anders S (2014). Moderated estimation of fold change and dispersion
  for RNA-seq data with DESeq2. *Genome Biology* 15:550.
- Mallick H et al. (2021). Multivariable association discovery in population-scale
  meta-omics studies. *PLoS Computational Biology* 17:e1009442.
- Lin H & Peddada SD (2020). Analysis of compositions of microbiomes with bias
  correction. *Nature Communications* 11:3514.

## License

MIT
