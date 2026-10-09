# Handoff: Ribosome package, tests first, then document and consolidate

Written 2026-10-09 on obed by a Claude Code session working in
`~/Repositories/polysomeModel` (session 1c32c4d6). Audience: the next
Claude session (and MG) working in `~/Repositories/Ribosome`.

Recommended model: Opus for Part B (deciding what is live, dead, or
duplicate). Sonnet is sufficient for Part A (the tests are specified
below) and Part D (writing docs from the map). Use haiku subagents to
read or grep the two large files; never read them whole into an Opus
context.

## 0. Where things are

- Code: `~/Repositories/Ribosome`, clone of the private GitHub repo
  `rurquidi/Ribosome`, commit f70c00a (2023-07-25, 11 commits). MG has
  read access; do not push to it without asking.
- Review that motivates this work:
  `~/Repositories/polysomeModel/Claude.Code_2026-10-09_v1.1_simulated.editorial.review.md`,
  Section 9 (code verification) and items 3.2 to 3.4.
- Independent reference solver, base R, no limSolve:
  `~/Repositories/Ribosome/Claude.Code_2026-10-09_reference.solver_check.share.R`
  (copied from the review session scratchpad). It reproduces the Fig 10
  curves of the paper to three figures and is the oracle for the tests.
- Manuscript: `~/Repositories/polysomeModel/mRNA_Populations_PLOS_GENETICS.tex`.

## 1. What is known about the code (verified this session)

File roles:

| File | Lines | Role |
|---|---|---|
| `R/Polysome_functions_RAUC.R` | 685 | Solver and summary statistics. Live. |
| `Text/Paper final figures.Rmd` | ~1700 | Produces the paper's figures. Live. |
| `R/Manuscript_Figures.R` | 3954 | Earlier version of the Rmd. Mostly superseded. |
| `R/Model Figures.R` | 1322 | Older exploratory figures. Probably dead. |
| `R/samples.R` | 743 | `*NP` solver variants and plot helpers. Probably dead. |
| `R/Continuous approximation.R` | 392 | PDE/continuum approximation. Unknown. |
| `R/Parameter histograms.R` | 406 | Fig 2 and parameter mining from data. Live but needs data. |
| `R/Reanalysis_20200914.R` | 1835 | Unknown. |
| `R/hello.R`, `polysome_model_results.Rmd` | | RStudio templates. Dead. |

Solver facts (from `Polysome_functions_RAUC.R`, confirmed against the
reference solver):

- Parameter vector order: `kappa, mu, lambda, tau, delta, s2`.
  `tau` is the elongation rate in codons per second. `s2` is unused.
- Termination from class i: `tau * i / (9 * imax)`.
- Initiation into class i+1 from i: `kappa * (1 - i/imax)` when
  `"kappa"` is in `which.interfere`, else `kappa`.
- Capped steady state: tridiagonal solve, `Solve.tridiag` from limSolve.
- Decapped steady state: `m*_i = mu * sum_{j>=i} m_j / tau_i` for i >= 1,
  `m*_0 = mu * sum_j m_j / delta`.
- Three near-duplicate solver sets exist: `CalcUnmarkedClass` /
  `CalcMarkedClass` / `Classes` (lines 41 to 157), the `_newtau` set
  (412 to 530), and the `*NP` set in `samples.R`. The Rmd calls the
  first set only. The `_newtau` set has stray `print(tau)` and
  `print(kappa)` debug calls, and `CalcMarkedClass_newtau` falls back
  to `CalcUnmarkedClass`, not `CalcUnmarkedClass_newtau`.
- Summary statistics (lines 531 to 592): `Ribo_percent`, `Ribo_mean`,
  `Ribo_sd`, `Ribo_mode`, `Ribo_frac`, `Ribo_frac_load`. The vector
  branch of `Ribo_sd` references `mean.tmp` and `square.diff.tmp`, which
  are never defined there. Latent bug; a test will expose it.
- `CHECK_STEADY_STATE` (line 158) exists and can be reused as a
  residual test.

Figure facts (from the Rmd):

- Parameter grid: `CreateParMatrix` uses `expand.grid`, so kappa (5000
  values, `seq(0.0001, 0.5, 0.0001)`) varies fastest and each block of
  5000 rows is one mu. Panels are selected by hard-coded row ranges such
  as `20001:25000`.
- All model figures use `tau = 1`, `lambda = 1`, `delta = 1e5`,
  `imax = 39` (yeast median) or `imax.range = c(4, 12, 39, 98, 194)`.
- Yeast mu list (Presnyak quantiles): `0.01, 0.0057762265, 0.0038508177,
  0.0038508177, 0.0023104906, 0.0016503504, 0.0008251752, 0.0004443251,
  0.0003300701, 0.0002221626`. Note the duplicated third and fourth
  entries.
- Fig 10 curves are blocks 45001:50000 (mu 2.2e-4), 40001:45000
  (mu 3.3e-4) and 20001:25000 (mu 2.3e-3). Peak heights 1.000, 0.692,
  0.143. The reference solver reproduces these.
- Protein production is coded as `colSums(classes) * Ribo_mean(classes)`,
  i.e. total bound ribosomes.
- Data for Fig 2 and Fig 11 are read from `../Data` and `../../Data`,
  which are not in the repo. Until Ricardo supplies them, nothing that
  touches those chunks can be tested.

Environment on obed: R 4.6.1 at `/usr/bin/Rscript`; `data.table`
installed; `limSolve` and `testthat` not installed. `tests/testthat.R`
exists (4 lines) but there is no `tests/testthat/` directory.

## 2. Part A: lock current behaviour with tests (do this first)

Goal: a test suite that passes on commit f70c00a unchanged, so that
every later deletion or merge of duplicate code is checked against it.
Do not fix anything in Part A. A test that fails on the original code
is a finding; record it in a `KNOWN_FAILURES.md` and mark the test
`skip()` with the reason, so the negative is bookkept rather than
hidden.

A.1. Setup

- Install `limSolve` and `testthat` into the user library.
- Create `tests/testthat/` and a `helper-params.R` that defines the
  canonical parameter sets below so every test file shares them.
- Source `R/Polysome_functions_RAUC.R` in the helper (the package does
  not load cleanly; `DESCRIPTION` has `Type: Ribosome`, which is
  invalid, and `NAMESPACE` is a stub). Do not fix DESCRIPTION in Part A.

A.2. Canonical parameter sets (put in the helper)

- `imax`: 39, plus 4 and 194 for the length extremes.
- `tau = 1`, `lambda = 1`, `delta = 1e5`, `s2 = 1`.
- `mu`: the ten yeast values above.
- `kappa`: a 25-point log grid from 1e-4 to 0.5 (not the full 5000;
  keep tests fast).
- Both interference settings used in the paper:
  `which.interfere = c("tau")` and `c("kappa", "tau")`.

A.3. Characterisation (golden) tests

For each parameter set and each of `CalcUnmarkedClass`,
`CalcMarkedClass`, `Classes`, compute the output once on the unmodified
code, save as `tests/testthat/golden/<name>.rds`, and have the test
compare with `expect_equal(tolerance = 1e-10)`. Commit the golden files.
This is the net that protects the refactor.

A.4. Invariant tests (from the mathematics; these are the ones that
would catch a wrong consolidation even if the golden files were
regenerated by mistake)

- Capped total: `sum(unmarked) == lambda / mu`.
- Decapped class 0: `marked[1] == lambda / delta`.
- Decapped recurrence: for i >= 1,
  `tau * i / (9 * imax) * marked[i + 1] == mu * sum(unmarked[(i + 1):(imax + 1)])`.
- Steady-state residual: rebuild the tridiagonal matrix as in the
  reference solver and check `A %*% unmarked + b` is zero to 1e-9.
  `CHECK_STEADY_STATE` may do this already; verify before reusing.
- Probabilities: `sum(Ribo_percent(x)) == 1` for vectors and columns.
- Decapped share of protein output: with
  `pU = sum((0:imax) * unmarked)` and `pM = sum((0:imax) * marked)`,
  `pM / (pU + pM) == 9 * mu * imax / tau / (1 + 9 * mu * imax / tau)`
  for every kappa, to 1e-8. This is the identity derived in the review,
  Section 3.3 and 9.2.
- Decapped monotonicity: `marked[-1]` is non-increasing in i (the
  manuscript asserts this; test it).
- Non-negativity of all classes.

A.5. Equivalence tests between duplicates

- `CalcUnmarkedClass` vs `CalcUnmarkedClass_newtau` vs
  `CalcUnmarkedClassNP` vs the reference solver, same inputs.
- Same for the marked and `Classes` variants.
- Expected: the first set and the reference solver agree. Whether the
  `_newtau` and `NP` sets agree is unknown; record the answer. The
  `_newtau` functions print to stdout; wrap in `capture.output()`.

A.6. Summary statistic tests

- `Ribo_mean`, `Ribo_percent`, `Ribo_mode`, `Ribo_sd` on a vector and
  on a matrix, against hand-computed values on a 4-class toy input.
- `Ribo_sd` vector branch is expected to fail (undefined objects).
  Record and skip.

A.7. Figure regression

- From the Rmd parameter setup, recompute the three Fig 10 curves at
  the 25-point kappa grid and check the normalised peaks are
  1.000, 0.692 and 0.143 to two decimals, and that the decapped share
  is 0.072, 0.104 and 0.448 respectively. This ties the suite to the
  published figure.

A.8. Run `Rscript -e 'testthat::test_dir("tests/testthat")'` and save
the output to a file; report only the summary line and any failures.

## 3. Part B: map the code (Opus)

Produce `Claude.Code_<date>_v1.0_code.map.md` with one row per
function and one row per Rmd chunk: name, file, lines, status (live,
superseded, dead, unknown), what it computes, which figure or table it
feeds, which equation in the manuscript it implements. Evidence for
"dead" must be a grep showing no callers outside its own file. Do not
delete anything in Part B.

Specific questions to settle:

- Is anything in `Manuscript_Figures.R` not reproduced in the Rmd?
- Is `Continuous approximation.R` referenced by the paper (the
  manuscript mentions a truncated Gaussian approximation)?
- What does `Reanalysis_20200914.R` do?
- Which Rmd chunk produced each of the eleven manuscript figures, with
  its row block and therefore its true mu value. The review found Fig 10
  mislabeled; check Figs 3 to 9 the same way.

## 4. Part C: consolidate, one step at a time

Each step is its own commit on a branch, and the full test suite must
pass before and after. Suggested order:

1. Remove the templates (`hello.R`, `polysome_model_results.Rmd`).
2. Remove the `*NP` set in `samples.R` if the equivalence tests show it
   is redundant, else document why it differs.
3. Merge `_newtau` into the base set (or delete it) based on A.5.
4. Replace hard-coded row ranges in the Rmd with lookups by mu value,
   so a panel's parameters are stated, not implied. The A.7 test guards
   this.
5. Fix `Ribo_sd`, `DESCRIPTION` (`Type: Package`, real description,
   license), `NAMESPACE`.
6. Only then consider a fourth solver that takes `mu/tau` as an
   explicit dimensionless parameter (review item 9.4). That is a model
   change and needs MG's sign-off.

## 5. Part D: document

- README: what the model is (three sentences), how to install, how to
  regenerate each figure (chunk name, parameters, data needed), and the
  figure-to-code table from Part B.
- roxygen headers on the live functions, using the parameter
  definitions in Section 1 above so the vocabulary matches the paper.
- Note in the README that the data directory is external and list the
  files expected.

## 6. Constraints

- Private repo belonging to Ricardo: work on a branch, do not push
  without MG's say-so, and do not force anything.
- ASCII only in files. Naming: `Claude.Code_<date>_v<ver>_<topic>.md`
  for session-authored documents.
- Keep large outputs out of context: redirect test runs and file
  listings to files and print tails.
- Do not read the 3954-line or 1700-line files whole; grep, or
  delegate reading to a haiku subagent that returns a table.
- `git commit --only <paths>`, never bare add and commit.
