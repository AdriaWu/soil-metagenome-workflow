# Empirical FDR via Label Permutations

A practical guide to permutation-based false-discovery-rate calibration for differential
abundance analysis of metagenomic data. Based on the BIOL3161 Soil Genomics Lab (2026)
workflow.

---

## Why empirical FDR?

When we run a differential abundance test across thousands of contigs (or genes, or MAGs),
we generate thousands of p-values. To control the false discovery rate, the standard tool
is the **Benjamini–Hochberg (BH) correction**: it adjusts p-values assuming the tests are
roughly independent, and lets us call any contig with `BH-FDR < 0.05` significant.

That assumption breaks in two common situations:

1. **Correlated tests.** Contigs from the same organism move together. Co-occurring
   microbes move together. So contig-level tests are not really independent, and BH
   tends to be too permissive when correlation is strong.
2. **Post-hoc design choices.** If we merged or dropped groups *after* looking at the
   data (for example, merging water and chemical fertiliser into "no_compost"), BH has
   no way to know about that decision, so it can't account for the inflation it might
   introduce.

**Empirical FDR** sidesteps both problems by building the null distribution from the data
itself, using random label permutations. Instead of trusting a theoretical formula, we
*ask the data what a null result looks like*.

---

## The core idea, in plain words

1. Run the test on the real data. Record the p-values.
2. Shuffle the treatment labels randomly. Run the test again. Any "significant" hits are
   false positives by definition — the labels are random.
3. Repeat many times.
4. At any p-value cutoff `τ`, compare how many hits we get with real labels versus the
   average across all the shuffled runs. That ratio is the empirical FDR at `τ`.
5. Find the `τ` where empirical FDR first hits our target (e.g. 5% or 10%). Use that
   `τ` as the significance threshold.

Formally:

```
                       mean null discoveries at τ
empirical FDR(τ) = ───────────────────────────────────
                       observed discoveries at τ
```

If BH and empirical FDR agree, great — BH was fine for this dataset. If empirical FDR
gives a *tighter* threshold (smaller `τ`), BH was being too permissive and we should
trust empirical. If empirical never reaches the target FDR at all, the design is
underpowered for that target.

---

## A worked walk-through

Let's make this concrete. Suppose we have **9 river-soil samples**: 3 in `compost`,
6 in `no_compost`. We're testing every contig for differential abundance between the
two groups.

### Step 1: Observed analysis

Run edgeR with the real labels. Suppose we get something like this (toy numbers for
illustration):

| Contig    | p-value | BH-FDR |
|-----------|---------|--------|
| ctg_001   | 0.0001  | 0.02   |
| ctg_002   | 0.0008  | 0.04   |
| ctg_003   | 0.003   | 0.08   |
| ctg_004   | 0.015   | 0.18   |
| ...       | ...     | ...    |

Count how many contigs fall below each candidate `τ`. For example:

- `τ = 0.01` → 12 contigs below (call these `obs_hits(0.01) = 12`)
- `τ = 0.02` → 25 contigs below
- `τ = 0.05` → 70 contigs below

### Step 2: Permute and repeat

Shuffle the `compost` / `no_compost` labels and re-run edgeR. The shuffle preserves how
many of each label exist — there are still 3 compost slots and 6 no_compost slots, just
assigned to different samples. Now any "hit" is a false positive.

For each shuffled run, count hits at each `τ`. Across many shuffles, take the mean:

- Permutation 1: 4 hits at `τ = 0.01`
- Permutation 2: 0 hits at `τ = 0.01`
- Permutation 3: 2 hits at `τ = 0.01`
- ...
- **Mean across all permutations**: ~1.8 hits at `τ = 0.01`

### Step 3: Build the FDR curve

For each `τ`, divide mean null hits by observed hits:

| τ     | obs_hits | mean_null_hits | empirical FDR |
|-------|----------|----------------|---------------|
| 0.005 | 5        | 0.4            | 0.080         |
| 0.010 | 12       | 1.8            | 0.150         |
| 0.014 | 18       | 0.9            | **0.050**     |
| 0.020 | 25       | 2.5            | 0.100         |
| 0.050 | 70       | 12             | 0.171         |

The empirical FDR curve crosses 5% at `τ ≈ 0.014`. So we declare any contig with raw
`PValue ≤ 0.014` significant — keeping roughly 18 contigs, with an expected ~5% being
false positives.

Compare this to BH: `BH-FDR < 0.05` might correspond to `τ ≈ 0.020` and ~25 hits. The
empirical procedure is tighter — BH was being slightly too permissive because of the
post-hoc merge and contig correlation.

Visually, the procedure looks like this — observed and null discovery counts diverge as
`τ` grows, and their ratio (the empirical FDR) crosses the target threshold at a
specific cutoff:

![Empirical FDR curve and tau selection](figures/fig2_fdr_curve.svg)

---

## How many permutations actually exist?

This is where sample size becomes the bottleneck — and where the math gets interesting.

If you have **n** total samples and split them into a **k vs (n−k)** group comparison,
the number of distinct ways to assign labels is:

```
C(n, k) = n! / (k! × (n−k)!)
```

This is **not n!** — we don't care about the order of samples within a group, only which
group each sample belongs to. Note also that `C(n, k) = C(n, n−k)` because choosing
`k` samples for compost is equivalent to choosing the complementary `n−k` for no_compost.

For the BIOL3161 dataset (Saige's compost experiment):

| Soil   | n | k_compost | C(n, k) unique permutations |
|--------|---|-----------|-----------------------------|
| Forest | 7 | 3 (or 4)  | **35**                      |
| River  | 9 | 3 (or 6)  | **84**                      |

That's the total number of distinct label assignments — full stop. You cannot generate
more shuffles than this without repeating yourself.

The figure below enumerates the 35 forest permutations explicitly. Each row is one
assignment (filled = compost, empty = no_compost). The real labels appear as one of
the 35, and the "mirror image" — the maximally different assignment — appears as
another:

![All 35 label permutations for forest soil](figures/fig1_permutation_enumeration.svg)

### Why `sample()` 500 times can be wasteful

The standard implementation does:

```r
for (i in seq_len(500)) {
  meta_perm$treatment_merged <- sample(meta_perm$treatment_merged)
  ...
}
```

`sample()` shuffles the existing labels, preserving the count of each. So every shuffled
vector still has 3 compost / 6 no_compost (for river) — drawn at random from the 84
unique possibilities. With 500 draws from a pool of 84, each unique permutation gets
drawn ~6 times on average. **The information content is fixed at 84 unique p-value
vectors; oversampling buys nothing.**

For forest (35 unique perms), 500 draws gives ~14 redraws per unique permutation.

When `C(n, k)` is small enough to enumerate, **enumerate** — it's faster *and* gives the
exact null instead of a Monte Carlo approximation.

### The real labels are inside the null

Here's the subtle point. Among the `C(n, k)` permutations, **one is the real label
assignment itself**. We're enumerating every possible way to assign `k` samples to
compost — and the truth is one of those ways. When `sample()` happens to return that
exact ordering (or when enumeration reaches it), edgeR runs on the actual data and
returns exactly the observed p-values.

This is unavoidable and has an important consequence (see next section).

### What about the "mirror image"?

A two-sided test cares about the *magnitude* of the test statistic, not its direction.
Among the `C(n, k)` permutations, there's also an "anti-real" assignment — the most
extreme alternative, where the `k` compost slots are filled by samples that were
originally no_compost. That assignment produces logFCs with opposite sign but similar
magnitude, and a similarly extreme p-value distribution.

Strict mathematical statement: only **the real labels** are guaranteed to reproduce the
observed result exactly. Approximate practical statement: by test symmetry, roughly
**2 of the `C(n, k)` permutations** look as extreme as the real data.

---

## Sample size sets a floor on empirical FDR

This is the key result. Two distinct floors are at play.

### Floor 1: Minimum resolvable p-value

A permutation-based null can only resolve probabilities down to `1 / N_unique`:

- Forest (35 perms): minimum resolvable p ≈ **0.029**
- River (84 perms): minimum resolvable p ≈ **0.012**

This isn't what edgeR itself does — edgeR p-values come from a GLM asymptotic test, not
direct permutation. But the *empirical FDR* uses the permuted-label p-value distributions
as its null, and that null can't be more refined than the underlying combinatorics
support.

### Floor 2: The empirical FDR floor

The numerator of empirical FDR is the **mean** number of null hits at `τ`, averaged
across all permutations:

```
mean_null_hits(τ) = (1 / N_perm) × Σᵢ null_hits_i(τ)
```

When we enumerate all `C(n, k)` permutations, **one of them is the real labels**. For
that permutation, `null_hits = obs_hits` — exactly. That single contribution alone gives
a floor: even if every *other* permutation produced zero hits at `τ`, the mean would
still be at least `obs_hits(τ) / C(n, k)`. So:

```
empirical FDR(τ) ≥ obs_hits(τ) / [C(n, k) × obs_hits(τ)] = 1 / C(n, k)
```

If we also count the "anti-real" permutation as contributing roughly the same magnitude
(test symmetry), the floor doubles to `2 / C(n, k)`.

For Saige's dataset:

| Soil   | C(n, k) | Strict floor 1/C(n,k) | Symmetric floor 2/C(n,k) |
|--------|---------|-----------------------|--------------------------|
| Forest | 35      | 2.9%                  | **5.7%**                 |
| River  | 84      | 1.2%                  | **2.4%**                 |

Side by side, the consequence is clear: forest's floor sits above the conventional 5%
target, while river's sits well below it:

![How sample size sets the empirical FDR floor](figures/fig3_sample_size_floor.svg)

**For forest, this is a hard mathematical wall.** Targeting 5% FDR is impossible by
construction — the floor is at least 5.7%. The minimum reasonable target is **10%**.

**For river, 5% is achievable** (floor is 2.4%, comfortably below 5%). But if you want
forest and river to be reported in the same framework (same target FDR across panels of
a figure, for example), using 10% for both is the cleanest call.

---

## How to choose your target FDR

Decision rule, using the symmetric floor:

```
empirical FDR floor ≈ 2 / C(n, k)
```

Pick a target that sits safely above this floor:

| Floor                   | Target to use |
|-------------------------|---------------|
| Floor < 2%              | 5% is fine    |
| Floor between 2% and 5% | 5% is borderline; the curve hits the threshold steeply |
| Floor > 5%              | Use 10%       |

For BIOL3161:

| Soil   | C(n, k) | Symmetric floor | Achievable target          |
|--------|---------|-----------------|----------------------------|
| Forest | 35      | 5.7%            | **10% (mandatory)**        |
| River  | 84      | 2.4%            | 5% works; 10% for consistency |

The choice between 5% and 10% is a precision/recall trade-off, not a statistical truth:

- **5%**: stricter, fewer hits, cleaner list, misses more real signal
- **10%**: more permissive, more hits, dirtier list, recovers more real signal

For a headline result in a paper, 5% is the convention. For an exploratory or
supplementary analysis that feeds further downstream investigation (e.g. "are these
compost-up contigs seen in other experiments?"), 10% is well-defended.

---

## Reproducible workflow

Here's the full pipeline, with the enumeration approach baked in.

### Helper 1: edgeR wrapper

```r
run_edgeR <- function(counts, design, coef_name) {
  y <- DGEList(counts = counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]

  y <- calcNormFactors(y, method = "TMMwsp")
  y <- estimateDisp(y, design, robust = TRUE)

  fit <- glmQLFit(y, design, robust = TRUE)
  coef_idx <- which(colnames(design) == coef_name)
  if (length(coef_idx) == 0) coef_idx <- ncol(design)
  qlf <- glmQLFTest(fit, coef = coef_idx)

  tab <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$contig <- rownames(tab)
  list(table = tab, y = y, design = design)
}
```

### Helper 2: empirical FDR calculation

```r
calculate_empirical_fdr <- function(obs_p, perm_p_list, m_tests = NULL) {
  if (is.null(m_tests)) m_tests <- length(obs_p)

  # Tau grid -- finer near 0 where p-values cluster
  grid <- sort(unique(c(seq(1e-6, 0.2, by = 1e-3), 0.01, 0.05, 0.1)))

  # Observed discoveries at each tau
  obs_hits <- sapply(grid, function(tau) sum(obs_p <= tau, na.rm = TRUE))

  # Mean null discoveries across permutations
  null_hits_mat <- sapply(perm_p_list, function(pv) {
    sapply(grid, function(tau) sum(pv <= tau, na.rm = TRUE))
  })
  null_mean_hits <- rowMeans(null_hits_mat, na.rm = TRUE)

  empirical_fdr <- null_mean_hits / pmax(obs_hits, 1)
  bh_fdr        <- (m_tests * grid) / pmax(obs_hits, 1)

  tibble(
    tau              = grid,
    obs_hits         = obs_hits,
    null_mean_hits   = null_mean_hits,
    expected_uniform = m_tests * grid,
    empirical_fdr    = pmin(empirical_fdr, 1),
    bh_fdr           = pmin(bh_fdr, 1)
  )
}
```

### Helper 3: find τ at the target FDR

```r
find_tau <- function(fdr_curve, target_fdr = 0.05,
                     fdr_col = "empirical_fdr") {
  df  <- as.data.frame(fdr_curve)
  ord <- order(df[[fdr_col]], df$tau)
  df  <- df[ord, , drop = FALSE]
  df  <- df[!duplicated(df[[fdr_col]]), , drop = FALSE]

  approx(
    x    = df[[fdr_col]],
    y    = df$tau,
    xout = target_fdr,
    ties = "ordered"
  )$y
}
```

Returns `NA` if the curve never reaches the target — the test telling you it can't
calibrate at that FDR with the available data.

### The main wrapper: enumeration version

Use this when `C(n, k)` is small enough to enumerate (anything under a few thousand is
fine).

```r
run_compost_contrast_enum <- function(subset_meta, subset_counts,
                                      target_fdr = 0.10, label = "") {
  cat("\n=== ", label, " ===\n", sep = "")

  # Observed analysis
  design <- model.matrix(~ treatment_merged, data = subset_meta)
  res <- run_edgeR(
    counts    = subset_counts,
    design    = design,
    coef_name = "treatment_mergedcompost"
  )

  cat("  BH significant (FDR < 0.05):",
      sum(res$table$FDR < 0.05), "contigs\n")

  # Enumerate every assignment of k samples to compost
  n_samples <- nrow(subset_meta)
  k_compost <- sum(subset_meta$treatment_merged == "compost")
  all_perms <- combn(n_samples, k_compost, simplify = FALSE)
  n_perm    <- length(all_perms)

  cat("  Enumerating all", n_perm, "label permutations (C(",
      n_samples, ",", k_compost, "))\n")

  perm_p <- vector("list", n_perm)
  for (i in seq_along(all_perms)) {
    meta_perm <- subset_meta
    meta_perm$treatment_merged <- "no_compost"
    meta_perm$treatment_merged[all_perms[[i]]] <- "compost"
    meta_perm$treatment_merged <- factor(
      meta_perm$treatment_merged,
      levels = c("no_compost", "compost")
    )
    design_perm <- model.matrix(~ treatment_merged, data = meta_perm)

    perm_res <- tryCatch(
      run_edgeR(subset_counts, design_perm, "treatment_mergedcompost"),
      error = function(e) NULL
    )
    perm_p[[i]] <- if (!is.null(perm_res)) perm_res$table$PValue else numeric(0)

    if (i %% 10 == 0) cat("    Permutation", i, "of", n_perm, "\n")
  }

  # FDR curves and threshold
  m_tests   <- nrow(res$table)
  fdr_curve <- calculate_empirical_fdr(res$table$PValue, perm_p,
                                       m_tests = m_tests)
  tau       <- find_tau(fdr_curve, target_fdr = target_fdr,
                        fdr_col = "empirical_fdr")
  tau_bh    <- find_tau(fdr_curve, target_fdr = target_fdr,
                        fdr_col = "bh_fdr")

  # Floor diagnostics
  floor_emp <- 2 / n_perm
  cat("  Tested universe (m_tests)        :", m_tests, "contigs\n")
  cat("  Empirical FDR floor (2/C(n,k))   :", signif(floor_emp, 3), "\n")
  cat("  Target FDR                       :", target_fdr, "\n")
  cat("  Empirical FDR threshold (tau_emp):", signif(tau,    3), "\n")
  cat("  BH-style FDR threshold  (tau_bh) :", signif(tau_bh, 3), "\n")
  cat("  Discoveries at empirical tau     :",
      sum(res$table$PValue <= tau, na.rm = TRUE), "\n")

  list(res = res, fdr_curve = fdr_curve,
       tau = tau, tau_bh = tau_bh,
       n_perm = n_perm, floor_emp = floor_emp,
       m_tests = m_tests, target_fdr = target_fdr)
}
```

### Running it

```r
forest <- run_compost_contrast_enum(
  subset_meta   = meta_forest,
  subset_counts = counts_forest,
  target_fdr    = 0.10,
  label         = "Forest soil"
)

river <- run_compost_contrast_enum(
  subset_meta   = meta_river,
  subset_counts = counts_river,
  target_fdr    = 0.10,
  label         = "River soil"
)
```

---

## How to read the output

A healthy run looks like this:

```
=== River soil ===
  Retained 3812 of 10000 contigs after filtering
  BH significant (FDR < 0.05): 1487 contigs
  Enumerating all 84 label permutations (C( 9 , 3 ))
  ...
  Tested universe (m_tests)        : 3812 contigs
  Empirical FDR floor (2/C(n,k))   : 0.0238
  Target FDR                       : 0.1
  Empirical FDR threshold (tau_emp): 0.026
  BH-style FDR threshold  (tau_bh) : 0.039
  Discoveries at empirical tau     : 845
```

Things to check:

1. **`tau_emp` is well above the floor.** If `tau_emp ≈ floor_emp`, the target FDR is
   too aggressive for the design — relax it.
2. **`tau_emp < tau_bh`** means the empirical procedure tightened BH (the more common
   case in correlated data).
3. **`tau_emp > tau_bh`** means BH was unusually strict for this dataset — uncommon, but
   not pathological. Report both.
4. **`tau_emp = NA`** means the curve never reaches the target FDR. Two responses:
   either relax the target, or report only BH and note the empirical curve doesn't
   support a stricter call.

---

## Methodological notes for write-up

Two choices worth documenting in any methods section:

1. **Are the real labels included in the null?** This implementation keeps them in
   (they appear as one of the `C(n, k)` enumerated permutations). Some implementations
   strictly exclude them, giving a denominator of `C(n, k) − 1`. Including them is the
   more conservative choice and matches what `sample()` does in expectation. The
   difference is negligible for moderate sample sizes but worth a one-line footnote.

2. **Why enumerate instead of Monte Carlo?** For small `n`, the full enumeration *is* the
   true permutation null — random sampling would just be a noisy estimate of the same
   thing. For large `n`, where `C(n, k)` is unmanageable, Monte Carlo with `sample()` is
   the only practical option.

---

## Quick reference

```
C(n, k)             = number of unique label permutations
1 / C(n, k)         = minimum non-zero p-value resolvable from the null
1 / C(n, k)         = strict empirical FDR floor (real labels in null)
2 / C(n, k)         = symmetric empirical FDR floor (real + mirror)
target_fdr          = the FDR you commit to controlling
tau                 = raw p-value cutoff where empirical FDR(tau) = target_fdr
```

If `target_fdr < 2 / C(n, k)`, the design is underpowered for that target.
Relax `target_fdr`, increase sample size, or report BH only with appropriate caveats.

---

## References

- Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a
  practical and powerful approach to multiple testing. *JRSS-B* 57:289–300.
- Storey, J. D. (2002). A direct approach to false discovery rates. *JRSS-B* 64:479–498.
- Phipson, B. & Smyth, G. K. (2010). Permutation P-values should never be zero:
  calculating exact P-values when permutations are randomly drawn. *SAGMB* 9:Article 39.
- BIOL3161 Soil Genomics Lab Manual (2026), Part C — Empirical FDR via label
  permutations.
