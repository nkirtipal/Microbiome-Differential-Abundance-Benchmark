# Microbiome Differential Abundance Benchmark

Which microbiome differential abundance method should you trust? Three popular tools,
three different ways of approaching the same problem, tested on simulated data where
the right answer is already known.

> Learning project. Simulated data, meant to be read and re-run.

---

## The problem

Different microbiome differential abundance methods often give different answers about
which taxa are linked to a phenotype. A study across 38 real datasets found that 14
commonly used methods frequently disagreed with each other (Nearing et al. 2022,
*Nature Communications*).

The trouble is, on real data you don't know which method is right. There's no answer
key to check against — just disagreement, with no way to say who's wrong.

Simulation fixes that. Build the data yourself, decide in advance which taxa are truly
different between groups, and then see which method actually finds them.

---

## Three methods, three approaches

| Tool | Normalisation | Transform | Distribution | Compositionally aware |
|---|---|---|---|---|
| DESeq2 | RLE | None | Negative binomial | No |
| MaAsLin2 | TSS | AST (log) | Normal | No |
| ANCOM-BC2 | Bias correction | — | — | **Yes** |

These aren't just three popular names picked at random — they represent three different
ways of handling the same underlying problem: microbiome data is compositional, meaning
the numbers are proportions of a whole rather than independent counts.

**DESeq2** models the counts directly and ignores that issue entirely. It's the
baseline here — what happens if compositionality is never addressed.

**MaAsLin2** sits in the middle. It normalises the data into proportions first, then
assumes the result behaves roughly like a normal distribution.

**ANCOM-BC2** takes compositionality seriously from the start and corrects for it
directly, rather than working around it.

**Why not more methods?** The paper above tested 14. Running all of them would turn
this into a replication of that study rather than a focused comparison, and most people
building something like this don't have the time or infrastructure for 14 tools. Three
is enough to show that the choice of method matters. ALDEx2, Corncob, and metagenomeSeq
would be natural additions later.

---

## How the data is simulated

Real microbiome data is sparse, noisy, and correlated across taxa in ways that are hard
to fake convincingly. So instead of generating counts from a simple textbook
distribution, this uses **SparseDOSSA2**, fit to a public reference dataset (the Human
Microbiome Project). That keeps the simulated data's noise and sparsity close to what
real data looks like — a simpler simulation would make every method look better than it
actually performs in practice.

On top of that realistic background, a known signal is planted: a chosen set of taxa
are made genuinely different between the two groups, and everything else is left as
pure noise. Because this signal is planted on purpose, it's the answer key every method
is later checked against.

**The grid:**

| What varies | Levels |
|---|---|
| Sample size | small, large |
| Effect size | weak, strong |
| How many taxa carry signal | few, many |

Every combination is simulated separately and tested with all three methods, so the
result is a full picture of where each method holds up and where it breaks down —
not just one number per tool.

---

## What gets measured

For each method, on each scenario, checked against the planted signal:

- **Sensitivity** — how many of the true differences it actually found
- **Precision** — of what it reported, how much was actually real
- **False discovery rate** — compared against what the method itself promises to
  control. A method can look reasonably precise and still be quietly breaking its own
  promise — that's the specific problem Hawinkel et al. (2019, *Briefings in
  Bioinformatics*) found when they ran a similar test.
- **Runtime** — how long each method takes per scenario

---

## How to run it

```r
BiocManager::install(c("SparseDOSSA2", "DESeq2", "ANCOMBC", "Maaslin2"))
```

Run the scripts in `code/` in order from the repo root. The first one fits the
simulation template and caches it; everything after that reads and writes to `data/`
and `results/` relative to the repo root.

---

## What this doesn't cover

- **The other 11 methods** from the wider comparison — see above for why.
- **Real data.** A method that wins here can still disappoint on an actual dataset,
  because no simulation captures every quirk of real biology. This is a test on
  simulated data, not a final verdict — checking the conclusions against a real dataset
  would be the natural next step.
- **Covariates and random effects.** All three tools can handle them; this benchmark
  keeps to the simplest two-group comparison so the result reflects each method's core
  assumptions rather than how it handles extra complexity.

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
