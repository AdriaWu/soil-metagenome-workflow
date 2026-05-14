---

editor_options: 
  markdown: 
    wrap: 72
---

# Empirical FDR via Label Permutations

A practical guide to permutation-based false-discovery-rate calibration for differential abundance analysis of metagenomic data.

------------------------------------------------------------------------

## Why empirical FDR?

When we run a differential abundance test across thousands of contigs (or genes, or MAGs), we generate thousands of p-values. To control the false discovery rate, the standard tool is the **Benjamini–Hochberg (BH) correction**: it adjusts p-values assuming the tests are roughly independent, and lets us call any contig with `BH-FDR < 0.05` significant.

That assumption breaks in two common situations:

1.  **Correlated tests.** Contigs from the same organism move together. Co-occurring microbes move together. So contig-level tests are not really independent, and BH tends to be too permissive when correlation is strong.
2.  **Post-hoc design choices.** If we merged or dropped groups *after* looking at the data (for example, merging water and chemical fertiliser into "no_compost"), BH has no way to know about that decision, so it can't account for the inflation it might introduce.

**Empirical FDR** sidesteps both problems by building the null distribution from the data itself, using random label permutations. Instead of trusting a theoretical formula, we *ask the data what a null result looks like*.

------------------------------------------------------------------------

## The core idea, in plain words

1.  Run the test on the real data. Record the p-values.
2.  Shuffle the treatment labels randomly. Run the test again. Any "significant" hits are false positives by definition — the labels are random.
3.  Repeat many times.
4.  At any p-value cutoff `τ`, compare how many hits we get with real labels versus the average across all the shuffled runs. That ratio is the empirical FDR at `τ`.
5.  Find the `τ` where empirical FDR first hits our target (e.g. 5% or 10%). Use that `τ` as the significance threshold.

Formally:

```         
                       mean null discoveries at τ
empirical FDR(τ) = ───────────────────────────────────
                       observed discoveries at τ
```

If BH and empirical FDR agree, great — BH was fine for this dataset. If empirical FDR gives a *tighter* threshold (smaller `τ`), BH was being too permissive and we should trust empirical. If empirical never reaches the target FDR at all, the design is underpowered for that target.

------------------------------------------------------------------------

## A worked walk-through

Let's make this concrete. Suppose we have **9 river-soil samples**: 3 in `compost`, 6 in `no_compost`. We're testing every contig for differential abundance between the two groups.

### Step 1: Observed analysis

Run edgeR with the real labels. Suppose we get something like this (toy numbers for illustration):

| Contig  | p-value | BH-FDR |
|---------|---------|--------|
| ctg_001 | 0.0001  | 0.02   |
| ctg_002 | 0.0008  | 0.04   |
| ctg_003 | 0.003   | 0.08   |
| ctg_004 | 0.015   | 0.18   |
| ...     | ...     | ...    |

Count how many contigs fall below each candidate `τ`. For example:

- `τ = 0.01` → 12 contigs below (call these `obs_hits(0.01) = 12`)
- `τ = 0.02` → 25 contigs below
- `τ = 0.05` → 70 contigs below

### Step 2: Permute and repeat

Shuffle the `compost` / `no_compost` labels and re-run edgeR. The shuffle preserves how many of each label exist — there are still 3 compost slots and 6 no_compost slots, just assigned to different samples. Now any "hit" is a false positive.

For each shuffled run, count hits at each `τ`. Across many shuffles, take the mean:

- Permutation 1: 4 hits at `τ = 0.01`
- Permutation 2: 0 hits at `τ = 0.01`
- Permutation 3: 2 hits at `τ = 0.01`
- ...
- **Mean across all permutations**: \~1.8 hits at `τ = 0.01`

### Step 3: Build the FDR curve

For each `τ`, divide mean null hits by observed hits:

| τ     | obs_hits | mean_null_hits | empirical FDR |
|-------|----------|----------------|---------------|
| 0.005 | 5        | 0.4            | 0.080         |
| 0.010 | 12       | 1.8            | 0.150         |
| 0.014 | 18       | 0.9            | **0.050**     |
| 0.020 | 25       | 2.5            | 0.100         |
| 0.050 | 70       | 12             | 0.171         |

The empirical FDR curve crosses 5% at `τ ≈ 0.014`. So we declare any contig with raw `PValue ≤ 0.014` significant — keeping roughly 18 contigs, with an expected \~5% being false positives.

Compare this to BH: `BH-FDR < 0.05` might correspond to `τ ≈ 0.020` and \~25 hits. The empirical procedure is tighter — BH was being slightly too permissive because of the post-hoc merge and contig correlation.

![Top: observed and null discovery counts diverge as the p-value threshold τ grows. Bottom: their ratio is the empirical FDR, and the cutoff τ is read off where this curve first meets the target.](data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTAwJSIgdmlld0JveD0iMCAwIDY4MCA0NDAiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+CjxzdHlsZT4KLnQgIHsgZm9udC1mYW1pbHk6IC1hcHBsZS1zeXN0ZW0sIEJsaW5rTWFjU3lzdGVtRm9udCwgIlNlZ29lIFVJIiwgc2Fucy1zZXJpZjsgZm9udC1zaXplOiAxNHB4OyBmaWxsOiAjMkMyQzJBOyB9Ci50cyB7IGZvbnQtZmFtaWx5OiAtYXBwbGUtc3lzdGVtLCBCbGlua01hY1N5c3RlbUZvbnQsICJTZWdvZSBVSSIsIHNhbnMtc2VyaWY7IGZvbnQtc2l6ZTogMTJweDsgZmlsbDogIzVGNUU1QTsgfQoudGggeyBmb250LWZhbWlseTogLWFwcGxlLXN5c3RlbSwgQmxpbmtNYWNTeXN0ZW1Gb250LCAiU2Vnb2UgVUkiLCBzYW5zLXNlcmlmOyBmb250LXNpemU6IDE0cHg7IGZvbnQtd2VpZ2h0OiA1MDA7IGZpbGw6ICMyQzJDMkE7IH0KPC9zdHlsZT4KCjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMzQwIiB5PSIyOCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+RmluZGluZyDPhCB3aGVyZSBlbXBpcmljYWwgRkRSIG1lZXRzIHRoZSB0YXJnZXQ8L3RleHQ+Cgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSg2MCwgNjApIj4KPHRleHQgY2xhc3M9InRoIiB4PSIwIiB5PSIwIj5Ub3A6IGRpc2NvdmVyaWVzIGF0IGVhY2ggcC12YWx1ZSBjdXRvZmY8L3RleHQ+Cgo8bGluZSB4MT0iMCIgeTE9IjIwIiB4Mj0iMCIgeTI9IjE2MCIgc3Ryb2tlPSIjODg4NzgwIiBzdHJva2Utd2lkdGg9IjAuNSIvPgo8bGluZSB4MT0iMCIgeTE9IjE2MCIgeDI9IjU2MCIgeTI9IjE2MCIgc3Ryb2tlPSIjODg4NzgwIiBzdHJva2Utd2lkdGg9IjAuNSIvPgo8dGV4dCBjbGFzcz0idHMiIHg9Ii04IiB5PSIyNCIgdGV4dC1hbmNob3I9ImVuZCI+aGl0czwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIwIiB5PSIxNzgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPjA8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0cyIgeD0iMTQwIiB5PSIxNzgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPjAuMDE8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0cyIgeD0iMjgwIiB5PSIxNzgiIHRleHQtYW5jaG9yPSJtaWRkbGUiPjAuMDI1PC90ZXh0Pgo8dGV4dCBjbGFzcz0idHMiIHg9IjQyMCIgeT0iMTc4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj4wLjA1PC90ZXh0Pgo8dGV4dCBjbGFzcz0idHMiIHg9IjU2MCIgeT0iMTc4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj4wLjEwPC90ZXh0Pgo8dGV4dCBjbGFzcz0idHMiIHg9IjI4MCIgeT0iMTk2IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5wLXZhbHVlIHRocmVzaG9sZCAoz4QpPC90ZXh0PgoKPHBhdGggZD0iTTAgMTYwIEw0MCAxNDUgTDgwIDEyMCBMMTQwIDkwIEwyMDAgNzAgTDI2MCA1NSBMMzIwIDQ1IEwzODAgNDAgTDQ0MCAzNiBMNTAwIDMzIEw1NjAgMzEiIHN0cm9rZT0iIzE4NUZBNSIgc3Ryb2tlLXdpZHRoPSIyIiBmaWxsPSJub25lIi8+CjxwYXRoIGQ9Ik0wIDE2MCBMNDAgMTU4IEw4MCAxNTYgTDE0MCAxNTIgTDIwMCAxNDggTDI2MCAxNDIgTDMyMCAxMzQgTDM4MCAxMjQgTDQ0MCAxMTAgTDUwMCA5NSBMNTYwIDc4IiBzdHJva2U9IiM4ODg3ODAiIHN0cm9rZS13aWR0aD0iMiIgZmlsbD0ibm9uZSIgc3Ryb2tlLWRhc2hhcnJheT0iNCAzIi8+Cgo8Y2lyY2xlIGN4PSIxODAiIGN5PSI4MCIgcj0iMyIgZmlsbD0iIzE4NUZBNSIvPgo8dGV4dCBjbGFzcz0idHMiIHg9IjE4NiIgeT0iNzYiIGZpbGw9IiMwQzQ0N0MiPm9ic19oaXRzPC90ZXh0Pgo8Y2lyY2xlIGN4PSIxODAiIGN5PSIxNTAiIHI9IjMiIGZpbGw9IiM4ODg3ODAiLz4KPHRleHQgY2xhc3M9InRzIiB4PSIxODYiIHk9IjE0NiIgZmlsbD0iIzQ0NDQ0MSI+bnVsbF9tZWFuX2hpdHM8L3RleHQ+CjwvZz4KCjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDYwLCAyODApIj4KPHRleHQgY2xhc3M9InRoIiB4PSIwIiB5PSIwIj5Cb3R0b206IGVtcGlyaWNhbCBGRFIgPSBudWxsIGhpdHMgLyBvYnNlcnZlZCBoaXRzPC90ZXh0PgoKPGxpbmUgeDE9IjAiIHkxPSIyMCIgeDI9IjAiIHkyPSIxNDAiIHN0cm9rZT0iIzg4ODc4MCIgc3Ryb2tlLXdpZHRoPSIwLjUiLz4KPGxpbmUgeDE9IjAiIHkxPSIxNDAiIHgyPSI1NjAiIHkyPSIxNDAiIHN0cm9rZT0iIzg4ODc4MCIgc3Ryb2tlLXdpZHRoPSIwLjUiLz4KPHRleHQgY2xhc3M9InRzIiB4PSItOCIgeT0iMjQiIHRleHQtYW5jaG9yPSJlbmQiPkZEUjwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSItOCIgeT0iNjAiIHRleHQtYW5jaG9yPSJlbmQiPjAuMTA8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0cyIgeD0iLTgiIHk9IjkyIiB0ZXh0LWFuY2hvcj0iZW5kIj4wLjA1PC90ZXh0Pgo8dGV4dCBjbGFzcz0idHMiIHg9Ii04IiB5PSIxNDIiIHRleHQtYW5jaG9yPSJlbmQiPjA8L3RleHQ+Cgo8bGluZSB4MT0iMCIgeTE9IjYwIiB4Mj0iNTYwIiB5Mj0iNjAiIHN0cm9rZT0iI0EzMkQyRCIgc3Ryb2tlLXdpZHRoPSIwLjgiIHN0cm9rZS1kYXNoYXJyYXk9IjQgMyIvPgo8dGV4dCBjbGFzcz0idHMiIHg9IjU1NiIgeT0iNTYiIHRleHQtYW5jaG9yPSJlbmQiIGZpbGw9IiNBMzJEMkQiPnRhcmdldCBGRFIgPSAwLjEwPC90ZXh0PgoKPHBhdGggZD0iTTAgMTQwIEwyMCAxMzggTDQwIDEzNCBMNjAgMTI4IEw4MCAxMTggTDEwMCAxMDUgTDEyMCA4OCBMMTQwIDcwIEwxNjAgNjAgTDE4MCA1NSBMMjIwIDQ3IEwyNjAgNDIgTDMyMCAzOCBMMzgwIDM2IEw0NDAgMzQgTDUwMCAzMiBMNTYwIDMwIiBzdHJva2U9IiMwRjZFNTYiIHN0cm9rZS13aWR0aD0iMiIgZmlsbD0ibm9uZSIvPgoKPGxpbmUgeDE9IjE2MCIgeTE9IjYwIiB4Mj0iMTYwIiB5Mj0iMTQwIiBzdHJva2U9IiMwRjZFNTYiIHN0cm9rZS13aWR0aD0iMC44IiBzdHJva2UtZGFzaGFycmF5PSIzIDMiLz4KPGNpcmNsZSBjeD0iMTYwIiBjeT0iNjAiIHI9IjQiIGZpbGw9IiMwRjZFNTYiLz4KPHRleHQgY2xhc3M9InRoIiB4PSIxNjAiIHk9IjE1OCIgdGV4dC1hbmNob3I9Im1pZGRsZSIgZmlsbD0iIzA0MzQyQyI+z4Qg4omIIDAuMDI1PC90ZXh0PgoKPHRleHQgY2xhc3M9InRzIiB4PSIyODAiIHk9IjE5NiIgdGV4dC1hbmNob3I9Im1pZGRsZSI+cC12YWx1ZSB0aHJlc2hvbGQgKM+EKTwvdGV4dD4KPC9nPgo8L3N2Zz4K){alt="Top: observed and null discovery counts diverge as the p-value threshold τ grows. Bottom: their ratio is the empirical FDR, and the cutoff τ is read off where this curve first meets the target."}

------------------------------------------------------------------------

## How many permutations actually exist?

This is where sample size becomes the bottleneck — and where the math gets interesting.

If you have **n** total samples and split them into a **k vs (n−k)** group comparison, the number of distinct ways to assign labels is:

```         
C(n, k) = n! / (k! × (n−k)!)
```

This is **not n!** — we don't care about the order of samples within a group, only which group each sample belongs to. Note also that `C(n, k) = C(n, n−k)` because choosing `k` samples for compost is equivalent to choosing the complementary `n−k` for no_compost.

For the BIOL3161 dataset (Saige's compost experiment):

| Soil   | n   | k_compost | C(n, k) unique permutations |
|--------|-----|-----------|-----------------------------|
| Forest | 7   | 3 (or 4)  | **35**                      |
| River  | 9   | 3 (or 6)  | **84**                      |

That's the total number of distinct label assignments — full stop. You cannot generate more shuffles than this without repeating yourself.

![All 35 ways to assign 3 of 7 forest samples to compost. The real labels and their mirror image both sit inside the enumeration — the source of the empirical FDR floor.](data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTAwJSIgdmlld0JveD0iMCAwIDY4MCA0NjAiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+CjxzdHlsZT4KLnQgIHsgZm9udC1mYW1pbHk6IC1hcHBsZS1zeXN0ZW0sIEJsaW5rTWFjU3lzdGVtRm9udCwgIlNlZ29lIFVJIiwgc2Fucy1zZXJpZjsgZm9udC1zaXplOiAxNHB4OyBmaWxsOiAjMkMyQzJBOyB9Ci50cyB7IGZvbnQtZmFtaWx5OiAtYXBwbGUtc3lzdGVtLCBCbGlua01hY1N5c3RlbUZvbnQsICJTZWdvZSBVSSIsIHNhbnMtc2VyaWY7IGZvbnQtc2l6ZTogMTJweDsgZmlsbDogIzVGNUU1QTsgfQoudGggeyBmb250LWZhbWlseTogLWFwcGxlLXN5c3RlbSwgQmxpbmtNYWNTeXN0ZW1Gb250LCAiU2Vnb2UgVUkiLCBzYW5zLXNlcmlmOyBmb250LXNpemU6IDE0cHg7IGZvbnQtd2VpZ2h0OiA1MDA7IGZpbGw6ICMyQzJDMkE7IH0KPC9zdHlsZT4KCjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMzQwIiB5PSIzMCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+Qyg3LCAzKSA9IDM1IHVuaXF1ZSBsYWJlbCBhc3NpZ25tZW50czwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIzNDAiIHk9IjQ4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5FYWNoIHJvdyBzaG93cyBvbmUgd2F5IHRvIGFzc2lnbiAzIG9mIDcgc2FtcGxlcyB0byBjb21wb3N0IChmaWxsZWQgPSBjb21wb3N0KTwvdGV4dD4KCjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDYwLCA4MCkiPgo8Zz48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsMTgpIj48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsMzYpIj48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsNTQpIj48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsNzIpIj48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsOTApIj48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsMTA4KSI+PGNpcmNsZSBjeD0iMCIgIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iMTQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMjgiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iNDIiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iNTYiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNzAiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iODQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PC9nPgo8L2c+Cgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgyMDAsIDgwKSI+CjxnPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwxOCkiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwzNikiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCw1NCkiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCw3MikiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCw5MCkiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwxMDgpIj48Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI0MiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI1NiIgY3k9IjAiIHI9IjUiIGZpbGw9IiM3Rjc3REQiLz48Y2lyY2xlIGN4PSI3MCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48Y2lyY2xlIGN4PSI4NCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iIzdGNzdERCIgc3Ryb2tlLXdpZHRoPSIwLjgiLz48L2c+CjwvZz4KCjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDM0MCwgODApIj4KPGc+PGNpcmNsZSBjeD0iMCIgIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMTQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMjgiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNDIiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNTYiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iNzAiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iODQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PC9nPgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgwLDE4KSI+PGNpcmNsZSBjeD0iMCIgIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMTQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMjgiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNDIiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iNTYiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNzAiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iODQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PC9nPgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgwLDM2KSI+PGNpcmNsZSBjeD0iMCIgIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMTQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iMjgiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PGNpcmNsZSBjeD0iNDIiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNTYiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iNzAiIGN5PSIwIiByPSI1IiBmaWxsPSIjN0Y3N0REIi8+PGNpcmNsZSBjeD0iODQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiM3Rjc3REQiIHN0cm9rZS13aWR0aD0iMC44Ii8+PC9nPgoKPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCw3MikiPjxyZWN0IHg9Ii0xMCIgeT0iLTEwIiB3aWR0aD0iMTA1IiBoZWlnaHQ9IjIwIiByeD0iNCIgZmlsbD0iI0VFRURGRSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuNSIvPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNSIgZmlsbD0iIzUzNEFCNyIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzUzNEFCNyIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzUzNEFCNyIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjwvZz4KCjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsMTA4KSI+PHJlY3QgeD0iLTEwIiB5PSItMTAiIHdpZHRoPSIxMDUiIGhlaWdodD0iMjAiIHJ4PSI0IiBmaWxsPSIjRkFFRURBIiBzdHJva2U9IiNCQTc1MTciIHN0cm9rZS13aWR0aD0iMS41Ii8+PGNpcmNsZSBjeD0iMCIgIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiNCQTc1MTciIHN0cm9rZS13aWR0aD0iMS4yIi8+PGNpcmNsZSBjeD0iMTQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiNCQTc1MTciIHN0cm9rZS13aWR0aD0iMS4yIi8+PGNpcmNsZSBjeD0iMjgiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiNCQTc1MTciIHN0cm9rZS13aWR0aD0iMS4yIi8+PGNpcmNsZSBjeD0iNDIiIGN5PSIwIiByPSI1IiBmaWxsPSIjQkE3NTE3Ii8+PGNpcmNsZSBjeD0iNTYiIGN5PSIwIiByPSI1IiBmaWxsPSIjQkE3NTE3Ii8+PGNpcmNsZSBjeD0iNzAiIGN5PSIwIiByPSI1IiBmaWxsPSIjQkE3NTE3Ii8+PGNpcmNsZSBjeD0iODQiIGN5PSIwIiByPSI0IiBmaWxsPSJub25lIiBzdHJva2U9IiNCQTc1MTciIHN0cm9rZS13aWR0aD0iMS4yIi8+PC9nPgo8L2c+Cgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSg0ODAsIDgwKSI+CjxnPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwxOCkiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCwzNikiPjxjaXJjbGUgY3g9IjAiICBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjE0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjI4IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjN0Y3N0REIiBzdHJva2Utd2lkdGg9IjAuOCIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNSIgZmlsbD0iIzdGNzdERCIvPjwvZz4KPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoMCw1NCkiPjx0ZXh0IGNsYXNzPSJ0cyIgeD0iNDIiIHk9IjQiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZpbGw9IiM4ODg3ODAiPuKApjwvdGV4dD48L2c+CjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDAsNzIpIj48dGV4dCBjbGFzcz0idHMiIHg9IjQyIiB5PSI0IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmaWxsPSIjODg4NzgwIj4zNSB0b3RhbDwvdGV4dD48L2c+CjwvZz4KCjxsaW5lIHgxPSI2MCIgeTE9IjIxOCIgeDI9IjYyMCIgeTI9IjIxOCIgc3Ryb2tlPSIjRDNEMUM3IiBzdHJva2Utd2lkdGg9IjAuNSIvPgoKPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoODAsIDI1MCkiPgo8cmVjdCB4PSItMTAiIHk9Ii0xMCIgd2lkdGg9IjEyMCIgaGVpZ2h0PSIyMCIgcng9IjQiIGZpbGw9IiNFRUVERkUiIHN0cm9rZT0iIzUzNEFCNyIgc3Ryb2tlLXdpZHRoPSIxLjUiLz4KPGNpcmNsZSBjeD0iMCIgIGN5PSIwIiByPSI1IiBmaWxsPSIjNTM0QUI3Ii8+PGNpcmNsZSBjeD0iMTQiIGN5PSIwIiByPSI1IiBmaWxsPSIjNTM0QUI3Ii8+PGNpcmNsZSBjeD0iMjgiIGN5PSIwIiByPSI1IiBmaWxsPSIjNTM0QUI3Ii8+CjxjaXJjbGUgY3g9IjQyIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjxjaXJjbGUgY3g9IjU2IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjxjaXJjbGUgY3g9IjcwIiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjNTM0QUI3IiBzdHJva2Utd2lkdGg9IjEuMiIvPgo8dGV4dCBjbGFzcz0idGgiIHg9IjEyMCIgeT0iNCIgZmlsbD0iIzI2MjE1QyI+VGhlIHJlYWwgbGFiZWxzPC90ZXh0Pgo8dGV4dCBjbGFzcz0idHMiIHg9IjAiIHk9IjI4Ij5Qcm9kdWNlcyB0aGUgb2JzZXJ2ZWQgcmVzdWx0IGV4YWN0bHkgd2hlbjwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIwIiB5PSI0MiI+ZW51bWVyYXRpb24gaGl0cyB0aGlzIGFzc2lnbm1lbnQ8L3RleHQ+CjwvZz4KCjxnIHRyYW5zZm9ybT0idHJhbnNsYXRlKDM2MCwgMjUwKSI+CjxyZWN0IHg9Ii0xMCIgeT0iLTEwIiB3aWR0aD0iMTIwIiBoZWlnaHQ9IjIwIiByeD0iNCIgZmlsbD0iI0ZBRUVEQSIgc3Ryb2tlPSIjQkE3NTE3IiBzdHJva2Utd2lkdGg9IjEuNSIvPgo8Y2lyY2xlIGN4PSIwIiAgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iI0JBNzUxNyIgc3Ryb2tlLXdpZHRoPSIxLjIiLz48Y2lyY2xlIGN4PSIxNCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iI0JBNzUxNyIgc3Ryb2tlLXdpZHRoPSIxLjIiLz48Y2lyY2xlIGN4PSIyOCIgY3k9IjAiIHI9IjQiIGZpbGw9Im5vbmUiIHN0cm9rZT0iI0JBNzUxNyIgc3Ryb2tlLXdpZHRoPSIxLjIiLz4KPGNpcmNsZSBjeD0iNDIiIGN5PSIwIiByPSI1IiBmaWxsPSIjQkE3NTE3Ii8+PGNpcmNsZSBjeD0iNTYiIGN5PSIwIiByPSI1IiBmaWxsPSIjQkE3NTE3Ii8+PGNpcmNsZSBjeD0iNzAiIGN5PSIwIiByPSI1IiBmaWxsPSIjQkE3NTE3Ii8+CjxjaXJjbGUgY3g9Ijg0IiBjeT0iMCIgcj0iNCIgZmlsbD0ibm9uZSIgc3Ryb2tlPSIjQkE3NTE3IiBzdHJva2Utd2lkdGg9IjEuMiIvPgo8dGV4dCBjbGFzcz0idGgiIHg9IjEyMCIgeT0iNCIgZmlsbD0iIzQxMjQwMiI+VGhlIG1pcnJvciBpbWFnZTwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIwIiB5PSIyOCI+TWF4aW1hbGx5IGRpZmZlcmVudCBmcm9tIHJlYWwgbGFiZWxzIOKAlDwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIwIiB5PSI0MiI+b3Bwb3NpdGUgc2lnbiwgc2ltaWxhciBtYWduaXR1ZGU8L3RleHQ+CjwvZz4KCjxyZWN0IHg9IjYwIiB5PSIzMzAiIHdpZHRoPSI1NjAiIGhlaWdodD0iMTAwIiByeD0iOCIgZmlsbD0iI0YxRUZFOCIgc3Ryb2tlPSIjQjRCMkE5IiBzdHJva2Utd2lkdGg9IjAuNSIvPgo8dGV4dCBjbGFzcz0idGgiIHg9IjgwIiB5PSIzNTgiIGZpbGw9IiMyQzJDMkEiPldoeSB0aGlzIG1hdHRlcnM8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0cyIgeD0iODAiIHk9IjM4MCI+VGhlIHJlYWwgbGFiZWxzIGFwcGVhciBhcyBvbmUgb2YgdGhlIDM1IHBlcm11dGF0aW9ucyBkdXJpbmcgZW51bWVyYXRpb24uIFdoZW4gdGhlIGxvb3AgaGl0cyB0aGVtLDwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSI4MCIgeT0iMzk2Ij5lZGdlUiBydW5zIG9uIHRoZSBhY3R1YWwgZGF0YSBhbmQgcmV0dXJucyB0aGUgb2JzZXJ2ZWQgcC12YWx1ZXMg4oCUIGNvbnRyaWJ1dGluZyBvYnNfaGl0cy8zNSB0byB0aGUgbnVsbCBtZWFuLjwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSI4MCIgeT0iNDEyIj5CeSB0ZXN0IHN5bW1ldHJ5LCB0aGUgbWlycm9yIGltYWdlIGNvbnRyaWJ1dGVzIGEgc2ltaWxhciBhbW91bnQuIFRoYXQgZ2l2ZXMgYSBoYXJkIGZsb29yIG9uIGVtcGlyaWNhbCBGRFI8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0cyIgeD0iODAiIHk9IjQyOCI+b2Ygcm91Z2hseSAyLzM1IOKJiCA1LjclIOKAlCBpbmRlcGVuZGVudCBvZiBob3cgc3Ryb25nIHRoZSBzaWduYWwgaXMuPC90ZXh0Pgo8L3N2Zz4K){alt="All 35 ways to assign 3 of 7 forest samples to compost. The real labels and their mirror image both sit inside the enumeration — the source of the empirical FDR floor."}

### Why `sample()` 500 times can be wasteful

The standard implementation does:

``` r
for (i in seq_len(500)) {
  meta_perm$treatment_merged <- sample(meta_perm$treatment_merged)
  ...
}
```

`sample()` shuffles the existing labels, preserving the count of each. So every shuffled vector still has 3 compost / 6 no_compost (for river) — drawn at random from the 84 unique possibilities. With 500 draws from a pool of 84, each unique permutation gets drawn \~6 times on average. **The information content is fixed at 84 unique p-value vectors; oversampling buys nothing.**

For forest (35 unique perms), 500 draws gives \~14 redraws per unique permutation.

When `C(n, k)` is small enough to enumerate, **enumerate** — it's faster *and* gives the exact null instead of a Monte Carlo approximation.

### The real labels are inside the null

Here's the subtle point. Among the `C(n, k)` permutations, **one is the real label assignment itself**. We're enumerating every possible way to assign `k` samples to compost — and the truth is one of those ways. When `sample()` happens to return that exact ordering (or when enumeration reaches it), edgeR runs on the actual data and returns exactly the observed p-values.

This is unavoidable and has an important consequence (see next section).

### What about the "mirror image"?

A two-sided test cares about the *magnitude* of the test statistic, not its direction. Among the `C(n, k)` permutations, there's also an "anti-real" assignment — the most extreme alternative, where the `k` compost slots are filled by samples that were originally no_compost. That assignment produces logFCs with opposite sign but similar magnitude, and a similarly extreme p-value distribution.

Strict mathematical statement: only **the real labels** are guaranteed to reproduce the observed result exactly. Approximate practical statement: by test symmetry, roughly **2 of the `C(n, k)` permutations** look as extreme as the real data.

------------------------------------------------------------------------

## Sample size sets a floor on empirical FDR

This is the key result. Two distinct floors are at play.

### Floor 1: Minimum resolvable p-value

A permutation-based null can only resolve probabilities down to `1 / N_unique`:

- Forest (35 perms): minimum resolvable p ≈ **0.029**
- River (84 perms): minimum resolvable p ≈ **0.012**

This isn't what edgeR itself does — edgeR p-values come from a GLM asymptotic test, not direct permutation. But the *empirical FDR* uses the permuted-label p-value distributions as its null, and that null can't be more refined than the underlying combinatorics support.

### Floor 2: The empirical FDR floor

The numerator of empirical FDR is the **mean** number of null hits at `τ`, averaged across all permutations:

```         
mean_null_hits(τ) = (1 / N_perm) × Σᵢ null_hits_i(τ)
```

When we enumerate all `C(n, k)` permutations, **one of them is the real labels**. For that permutation, `null_hits = obs_hits` — exactly. That single contribution alone gives a floor: even if every *other* permutation produced zero hits at `τ`, the mean would still be at least `obs_hits(τ) / C(n, k)`. So:

```         
empirical FDR(τ) ≥ obs_hits(τ) / [C(n, k) × obs_hits(τ)] = 1 / C(n, k)
```

If we also count the "anti-real" permutation as contributing roughly the same magnitude (test symmetry), the floor doubles to `2 / C(n, k)`.

For Saige's dataset:

| Soil   | C(n, k) | Strict floor 1/C(n,k) | Symmetric floor 2/C(n,k) |
|--------|---------|-----------------------|--------------------------|
| Forest | 35      | 2.9%                  | **5.7%**                 |
| River  | 84      | 1.2%                  | **2.4%**                 |

![Forest (n=7) versus river (n=9). The number of unique permutations directly sets the minimum achievable empirical FDR; 5% is unreachable for forest but comfortable for river.](data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTAwJSIgdmlld0JveD0iMCAwIDY4MCAzNjAiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+CjxzdHlsZT4KLnQgIHsgZm9udC1mYW1pbHk6IC1hcHBsZS1zeXN0ZW0sIEJsaW5rTWFjU3lzdGVtRm9udCwgIlNlZ29lIFVJIiwgc2Fucy1zZXJpZjsgZm9udC1zaXplOiAxNHB4OyBmaWxsOiAjMkMyQzJBOyB9Ci50cyB7IGZvbnQtZmFtaWx5OiAtYXBwbGUtc3lzdGVtLCBCbGlua01hY1N5c3RlbUZvbnQsICJTZWdvZSBVSSIsIHNhbnMtc2VyaWY7IGZvbnQtc2l6ZTogMTJweDsgZmlsbDogIzVGNUU1QTsgfQoudGggeyBmb250LWZhbWlseTogLWFwcGxlLXN5c3RlbSwgQmxpbmtNYWNTeXN0ZW1Gb250LCAiU2Vnb2UgVUkiLCBzYW5zLXNlcmlmOyBmb250LXNpemU6IDE0cHg7IGZvbnQtd2VpZ2h0OiA1MDA7IGZpbGw6ICMyQzJDMkE7IH0KPC9zdHlsZT4KCjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMzQwIiB5PSIyOCIgdGV4dC1hbmNob3I9Im1pZGRsZSI+U2FtcGxlIHNpemUgc2V0cyBhIGhhcmQgZmxvb3Igb24gZW1waXJpY2FsIEZEUjwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIzNDAiIHk9IjQ2IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIj5mbG9vciDiiYggMiAvIEMobiwgaykg4oCUIHRoZSBzbWFsbGVyIHRoZSBzYW1wbGUsIHRoZSBoaWdoZXIgdGhlIGZsb29yPC90ZXh0PgoKPGcgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoNjAsIDgwKSI+CjxyZWN0IHg9IjAiIHk9IjAiIHdpZHRoPSIyNjAiIGhlaWdodD0iMjQwIiByeD0iMTIiIGZpbGw9IiNGQUVFREEiIHN0cm9rZT0iI0JBNzUxNyIgc3Ryb2tlLXdpZHRoPSIwLjUiLz4KPHRleHQgY2xhc3M9InRoIiB4PSIxMzAiIHk9IjI4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmaWxsPSIjNDEyNDAyIj5Gb3Jlc3Qgc29pbDwvdGV4dD4KPHRleHQgY2xhc3M9InRzIiB4PSIxMzAiIHk9IjQ2IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmaWxsPSIjNjMzODA2Ij5uID0gNywgayA9IDMgY29tcG9zdDwvdGV4dD4KCjx0ZXh0IGNsYXNzPSJ0cyIgeD0iMjAiIHk9IjgwIiBmaWxsPSIjNjMzODA2Ij5VbmlxdWUgcGVybXV0YXRpb25zPC90ZXh0Pgo8dGV4dCBjbGFzcz0idGgiIHg9IjI0MCIgeT0iODAiIHRleHQtYW5jaG9yPSJlbmQiIGZpbGw9IiM0MTI0MDIiPkMoNywgMykgPSAzNTwvdGV4dD4KCjx0ZXh0IGNsYXNzPSJ0cyIgeD0iMjAiIHk9IjEwOCIgZmlsbD0iIzYzMzgwNiI+U3RyaWN0IGZsb29yICgxL0MpPC90ZXh0Pgo8dGV4dCBjbGFzcz0idGgiIHg9IjI0MCIgeT0iMTA4IiB0ZXh0LWFuY2hvcj0iZW5kIiBmaWxsPSIjNDEyNDAyIj4yLjklPC90ZXh0PgoKPHRleHQgY2xhc3M9InRzIiB4PSIyMCIgeT0iMTM2IiBmaWxsPSIjNjMzODA2Ij5TeW1tZXRyaWMgZmxvb3IgKDIvQyk8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMjQwIiB5PSIxMzYiIHRleHQtYW5jaG9yPSJlbmQiIGZpbGw9IiM0MTI0MDIiPjUuNyU8L3RleHQ+Cgo8bGluZSB4MT0iMjAiIHkxPSIxNTYiIHgyPSIyNDAiIHkyPSIxNTYiIHN0cm9rZT0iI0JBNzUxNyIgc3Ryb2tlLXdpZHRoPSIwLjUiIG9wYWNpdHk9IjAuNCIvPgoKPHRleHQgY2xhc3M9InRzIiB4PSIyMCIgeT0iMTgwIiBmaWxsPSIjNjMzODA2Ij5UYXJnZXQgRkRSID0gNSU8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMjQwIiB5PSIxODAiIHRleHQtYW5jaG9yPSJlbmQiIGZpbGw9IiNBMzJEMkQiPkJlbG93IGZsb29yPC90ZXh0PgoKPHRleHQgY2xhc3M9InRzIiB4PSIyMCIgeT0iMjA4IiBmaWxsPSIjNjMzODA2Ij5UYXJnZXQgRkRSID0gMTAlPC90ZXh0Pgo8dGV4dCBjbGFzcz0idGgiIHg9IjI0MCIgeT0iMjA4IiB0ZXh0LWFuY2hvcj0iZW5kIiBmaWxsPSIjMDQzNDJDIj5BY2hpZXZhYmxlPC90ZXh0Pgo8L2c+Cgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgzNjAsIDgwKSI+CjxyZWN0IHg9IjAiIHk9IjAiIHdpZHRoPSIyNjAiIGhlaWdodD0iMjQwIiByeD0iMTIiIGZpbGw9IiNFMUY1RUUiIHN0cm9rZT0iIzBGNkU1NiIgc3Ryb2tlLXdpZHRoPSIwLjUiLz4KPHRleHQgY2xhc3M9InRoIiB4PSIxMzAiIHk9IjI4IiB0ZXh0LWFuY2hvcj0ibWlkZGxlIiBmaWxsPSIjMDQzNDJDIj5SaXZlciBzb2lsPC90ZXh0Pgo8dGV4dCBjbGFzcz0idHMiIHg9IjEzMCIgeT0iNDYiIHRleHQtYW5jaG9yPSJtaWRkbGUiIGZpbGw9IiMwODUwNDEiPm4gPSA5LCBrID0gMyBjb21wb3N0PC90ZXh0PgoKPHRleHQgY2xhc3M9InRzIiB4PSIyMCIgeT0iODAiIGZpbGw9IiMwODUwNDEiPlVuaXF1ZSBwZXJtdXRhdGlvbnM8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMjQwIiB5PSI4MCIgdGV4dC1hbmNob3I9ImVuZCIgZmlsbD0iIzA0MzQyQyI+Qyg5LCAzKSA9IDg0PC90ZXh0PgoKPHRleHQgY2xhc3M9InRzIiB4PSIyMCIgeT0iMTA4IiBmaWxsPSIjMDg1MDQxIj5TdHJpY3QgZmxvb3IgKDEvQyk8L3RleHQ+Cjx0ZXh0IGNsYXNzPSJ0aCIgeD0iMjQwIiB5PSIxMDgiIHRleHQtYW5jaG9yPSJlbmQiIGZpbGw9IiMwNDM0MkMiPjEuMiU8L3RleHQ+Cgo8dGV4dCBjbGFzcz0idHMiIHg9IjIwIiB5PSIxMzYiIGZpbGw9IiMwODUwNDEiPlN5bW1ldHJpYyBmbG9vciAoMi9DKTwvdGV4dD4KPHRleHQgY2xhc3M9InRoIiB4PSIyNDAiIHk9IjEzNiIgdGV4dC1hbmNob3I9ImVuZCIgZmlsbD0iIzA0MzQyQyI+Mi40JTwvdGV4dD4KCjxsaW5lIHgxPSIyMCIgeTE9IjE1NiIgeDI9IjI0MCIgeTI9IjE1NiIgc3Ryb2tlPSIjMEY2RTU2IiBzdHJva2Utd2lkdGg9IjAuNSIgb3BhY2l0eT0iMC40Ii8+Cgo8dGV4dCBjbGFzcz0idHMiIHg9IjIwIiB5PSIxODAiIGZpbGw9IiMwODUwNDEiPlRhcmdldCBGRFIgPSA1JTwvdGV4dD4KPHRleHQgY2xhc3M9InRoIiB4PSIyNDAiIHk9IjE4MCIgdGV4dC1hbmNob3I9ImVuZCIgZmlsbD0iIzA0MzQyQyI+QWNoaWV2YWJsZTwvdGV4dD4KCjx0ZXh0IGNsYXNzPSJ0cyIgeD0iMjAiIHk9IjIwOCIgZmlsbD0iIzA4NTA0MSI+VGFyZ2V0IEZEUiA9IDEwJTwvdGV4dD4KPHRleHQgY2xhc3M9InRoIiB4PSIyNDAiIHk9IjIwOCIgdGV4dC1hbmNob3I9ImVuZCIgZmlsbD0iIzA0MzQyQyI+Q29tZm9ydGFibGU8L3RleHQ+CjwvZz4KPC9zdmc+Cg==){alt="Forest (n=7) versus river (n=9). The number of unique permutations directly sets the minimum achievable empirical FDR; 5% is unreachable for forest but comfortable for river."}

**For forest, this is a hard mathematical wall.** Targeting 5% FDR is impossible by construction — the floor is at least 5.7%. The minimum reasonable target is **10%**.

**For river, 5% is achievable** (floor is 2.4%, comfortably below 5%). But if you want forest and river to be reported in the same framework (same target FDR across panels of a figure, for example), using 10% for both is the cleanest call.

------------------------------------------------------------------------

## How to choose your target FDR

Decision rule, using the symmetric floor:

```         
empirical FDR floor ≈ 2 / C(n, k)
```

Pick a target that sits safely above this floor:

| Floor | Target to use |
|---------------------------------------------|---------------------------|
| Floor \< 2% | 5% is fine |
| Floor between 2% and 5% | 5% is borderline; the curve hits the threshold steeply |
| Floor \> 5% | Use 10% |

For BIOL3161:

| Soil   | C(n, k) | Symmetric floor | Achievable target             |
|--------|---------|-----------------|-------------------------------|
| Forest | 35      | 5.7%            | **10% (mandatory)**           |
| River  | 84      | 2.4%            | 5% works; 10% for consistency |

The choice between 5% and 10% is a precision/recall trade-off, not a statistical truth:

- **5%**: stricter, fewer hits, cleaner list, misses more real signal
- **10%**: more permissive, more hits, dirtier list, recovers more real signal

For a headline result in a paper, 5% is the convention. For an exploratory or supplementary analysis that feeds further downstream investigation (e.g. "are these compost-up contigs seen in other experiments?"), 10% is well-defended.

------------------------------------------------------------------------

## Reproducible workflow

Here's the full pipeline, with the enumeration approach baked in.

### Helper 1: edgeR wrapper

``` r
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

``` r
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

``` r
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

Returns `NA` if the curve never reaches the target — the test telling you it can't calibrate at that FDR with the available data.

### The main wrapper: enumeration version

Use this when `C(n, k)` is small enough to enumerate (anything under a few thousand is fine).

``` r
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

``` r
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

------------------------------------------------------------------------

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

1.  **`tau_emp` is well above the floor.** If `tau_emp ≈ floor_emp`, the target FDR is too aggressive for the design — relax it.
2.  **`tau_emp < tau_bh`** means the empirical procedure tightened BH (the more common case in correlated data).
3.  **`tau_emp > tau_bh`** means BH was unusually strict for this dataset — uncommon, but not pathological. Report both.
4.  **`tau_emp = NA`** means the curve never reaches the target FDR. Two responses: either relax the target, or report only BH and note the empirical curve doesn't support a stricter call.

------------------------------------------------------------------------

## Methodological notes for write-up

Two choices worth documenting in any methods section:

1.  **Are the real labels included in the null?** This implementation keeps them in (they appear as one of the `C(n, k)` enumerated permutations). Some implementations strictly exclude them, giving a denominator of `C(n, k) − 1`. Including them is the more conservative choice and matches what `sample()` does in expectation. The difference is negligible for moderate sample sizes but worth a one-line footnote.

2.  **Why enumerate instead of Monte Carlo?** For small `n`, the full enumeration *is* the true permutation null — random sampling would just be a noisy estimate of the same thing. For large `n`, where `C(n, k)` is unmanageable, Monte Carlo with `sample()` is the only practical option.

------------------------------------------------------------------------

## Quick reference

```         
C(n, k)             = number of unique label permutations
1 / C(n, k)         = minimum non-zero p-value resolvable from the null
1 / C(n, k)         = strict empirical FDR floor (real labels in null)
2 / C(n, k)         = symmetric empirical FDR floor (real + mirror)
target_fdr          = the FDR you commit to controlling
tau                 = raw p-value cutoff where empirical FDR(tau) = target_fdr
```

If `target_fdr < 2 / C(n, k)`, the design is underpowered for that target. Relax `target_fdr`, increase sample size, or report BH only with appropriate caveats.

------------------------------------------------------------------------

## References

- Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *JRSS-B* 57:289–300.
- Storey, J. D. (2002). A direct approach to false discovery rates. *JRSS-B* 64:479–498.
- Phipson, B. & Smyth, G. K. (2010). Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *SAGMB* 9:Article 39.
- BIOL3161 Soil Genomics Lab Manual (2026), Part C — Empirical FDR via label permutations.
