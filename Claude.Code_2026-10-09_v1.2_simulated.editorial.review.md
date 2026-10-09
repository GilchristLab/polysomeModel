---
title: "Editorial review: Modeling mRNA populations (PLOS Genetics draft)"
author: "Internal pre-submission review, prepared with Claude Code for the coauthors"
date: 2026-10-09
geometry: margin=1in
fontsize: 11pt
colorlinks: true
---

# Scope of this review

This review treats the manuscript as a PLOS Genetics editor would at
initial submission. It covers the current LaTeX source
(`mRNA_Populations_PLOS_GENETICS.tex`, last edited 2025-08-02), the
supplemental text, the figures, and the reference list. The committed PDF
was built on 2024-08-12 from an older draft and does not reflect the
current source, so all comments below refer to the source and to a fresh
compile of it.

**Version 1.1 (2026-10-09).** Adds Section 9, a check of the R code in
Ricardo's `Ribosome` repository (commit f70c00a, 2023-07-25), and
updates items 1.5, 3.2, 3.3 and 3.4, which that check resolves.
Version 1.0 was written from the manuscript alone.
**Version 1.2 (2026-10-09).** Corrects 1.5 and 9.7: the repository is
public, not private.

# Decision

**Return to authors before review.** The model is clean and the central
idea (two coupled polysome populations, capped and decapped) is worth
publishing. In its current form, however, the manuscript is incomplete,
the single validation is weak, the prior-art framing is incorrect, and
there are internal inconsistencies between the equations and the stated
results. An editor would most likely triage it for insufficient strength
of advance rather than send it out.

The code check (Section 9) sharpens this. The headline result, that a
short-lived mRNA produces just under half its protein from decapped
transcripts, is computed at an elongation rate of one codon per second
and at a decapping rate three times smaller than the figure caption
states. At the captioned rate the model gives two thirds, and at
realistic yeast elongation rates it gives one sixth to one third. The
central claim therefore has to be recomputed before the paper goes
anywhere.

Two paths forward:

1. **PLOS Genetics.** Add a genetics-facing result (see Section 4) and a
   second independent validation. Cope, Pak and Gilchrist (2026, PLOS
   Genetics) shows the journal accepts ribosome-scale modeling with a
   genomic application, so the fit is possible.
2. **PLOS Computational Biology.** The more natural home for a theory
   paper with one correlational validation. Chevalier et al. (2023), the
   closest prior model, is published there.

# 1. Completeness and build state

1.1. **The committed PDF is stale.** Built 2024-08-12 with a different
abstract and introduction. Rebuild before anyone reads it.

1.2. **Six citation keys have no bibliography entry** and render as [?]:
MacDonaldEtAl1968, MacDonaldAndGibbs1969, vonHeijneEtAl1987, Chou2003,
GilchristAndWagner2006, ShahAndGilchrist2011.

1.3. **Compile errors in the current source.**
Label `eqn:S` is referenced but the label is `eq:S`.
`\mvechatstart` is an undefined macro (should be `\mvechatstar`).
Two uses of `\MRL^*` in text mode produce "Missing $" and a 1098 pt
overfull line in the Results section on protein production.
The `strip` environment errors in one-column mode.

1.4. **Placeholders visible in the body.** Author Summary reads
"NEED TO WRITE THIS". A bold "[Mike: Ricardo can you revise ...]" note
and an `\mmpar` margin note are still in the Introduction.

1.5. **Code repository is incomplete.**
`github.com/rurquidi/Ribosome` is public (v1.0 and the first draft of
v1.1 called it private; it returned 404 at the time of checking and was
public by the end of the day). Submission needs a Data Availability
Statement pointing at it. The repository lacks the input data: the
scripts read Presnyak 2015, Chan 2018, Weinberg 2016 RPKMs, Dao Duc and
Song 2018 rates, and the yeast and Arabidopsis FASTA files from
`../Data` and `../../Data`, which were never pushed. The package depends
on limSolve, has a placeholder DESCRIPTION (no license, template
description text), and a top-level `polysome_model_results.Rmd` that is
the unedited RStudio template. See Section 9.

1.6. **Bibliography entries never cited:** RN1 (Browning and
Bailey-Serres 2015) and RN37 (Wu and Jaffrey 2016).

# 2. PLOS Genetics format requirements

2.1. Materials and Methods must follow the Discussion unless the cover
letter justifies otherwise. The draft places Methods before Results.

2.2. Author Summary: 150 to 200 words, first person, written for a
non-specialist. Missing.

2.3. Short title (70 characters) is missing. The full title "Modeling
mRNA populations" is too generic to be findable. Say what is coupled and
what is predicted.

2.4. References must be numbered in order of first citation, Vancouver
style, with full DOIs. The hand-written `thebibliography` is alphabetical
and has no DOIs.

2.5. Initial submission should be double spaced with continuous line
numbers. The abstract and summary currently sit outside `\linenumbers`.

2.6. Supplemental text errors. Matrix U: the second diagonal entry lacks
the leading minus sign and parentheses. Matrix R: stray backtick in the
(0,1) entry; the last two diagonal entries use (imax-2)/imax where
(imax-1)/imax is meant.

# 3. Model correctness (verify each against the R code)

3.1. **ODE typos.** `\kappa0` in the general capped row. That row's
coefficient reads kappa_0 (i-1)/imax but should be
kappa_0 (1 - (i-1)/imax). The last decapped row swaps m*_imax and
m_imax.

3.2. **Eq. (decapped_abundance) appears to drop a factor of imax.**
Summing the decapped equations from class k upward gives

    tau (k/imax) m*_k = mu * sum_{j>=k} m_j

so m*_k = (imax/k)(mu/tau) S_k, not (1/k)(mu/tau) S_k. The total
decapped pool, the odds expression, and p*_i all inherit this.
*Resolved in v1.1:* the code divides by tau_i = tau_c i/(9 imax), so it
carries the imax factor. The text is wrong and the figures are right.
Insert imax in Eq. (decapped_abundance), Eq. (decapped_solution), the
total decapped pool, Eq. (odds), and Eq. (decapped_distribution).

3.3. **An exact result the paper misses.** Multiply the identity in 3.2
by tau k/imax and sum over k:

    decapped protein production = mu * sum_j j m_j
    capped protein production   = (tau/imax) * sum_j j m_j

So the decapped share of total protein output is

    (mu imax / tau) / (1 + mu imax / tau)

a one-parameter closed form, independent of kappa. It is the probability
that decapping happens before a ribosome finishes its transit. This is a
cleaner headline than the "41 percent" figure, and it shows the claim
"never more than half" is a statement about the parameter range, not a
theorem: the share exceeds one half whenever the capped half-life is
shorter than ln(2) times the ribosome transit time imax/tau.
*Confirmed numerically in v1.1* (Section 9.2): a re-implementation of
the coded model gives a share that is constant in kappa to four decimal
places and matches the formula.

3.4. **Fig 10 panels B and C are swapped, and the mu values in the
caption and legend are wrong.** *Confirmed in v1.1 from the figure code*
(Section 9.3). The three curves are computed at mu = 2.2e-4, 3.3e-4 and
2.3e-3 per second (half-lives 52, 35 and 5 minutes), not at 2.2e-4,
1.7e-3 and 5.7e-3 (52, 7 and 2 minutes) as captioned. The "41 percent"
in the Results and Discussion is the share at mu = 2.3e-3.

3.5. **The supplement's lambda/mu result is presented as discovered by
"manual exploration".** It follows in one line from summing the capped
ODEs: d/dt (sum m_i) = lambda - mu sum m_i. State it that way.

3.6. **Units.** The initiation-to-elongation ratio kappa' is
dimensionless, not in units of 1/s. Its range is given as 0.001 to 0.1
in Methods, but the figures span 1e-4 to 0.5.

3.7. **The delta assumption is stated two ways:** `delta >> tau` in
Methods, `delta >> lambda` in Discussion. Pick one and justify it.

3.8. **Yeast decapping rate statistics are inconsistent.** Methods: mean
1.3e-3 (SD 1.8e-3). Figures: median 1.7e-3. A right-skewed rate
distribution should have mean above median. The Arabidopsis range
"1.7 x (10^-4 +- 2 x 10^-4)" is garbled.

3.9. **Percentile language is inconsistent.** 2e-4/s is called the
"99th percentile" (of half-life) in Results; 5.7e-3 is called "low
decapping" once in the Discussion. Use one convention, percentiles of
mu.

# 4. Strength of advance and validation

4.1. **The only test against data is Fig 11**, and the text overstates
it ("coincide"). Spearman rho = 0.59 is shown in the figure but not
stated in the text. Observed loads sit well above predicted ones:

| Quantity                 | Predicted | Observed |
|--------------------------|-----------|----------|
| Mean ribosome load range | 1 to 11   | 1.6 to 40|

The absolute scale of the empirical load rests on an ad hoc
"200 / length" correction. Present this as a rank test, explain the
roughly three-fold offset, and add a second independent check. Arava et
al. (2003) polysome-profile data are the obvious one; they also report
that ribosome density falls with ORF length, a direct test of the imax
prediction.

4.2. **The empirical mu distribution drives the headline result, and it
is method dependent.** Chan et al. (2018) report non-invasive yeast
half-lives averaging 4.8 minutes, far shorter than Presnyak et al.
(2015), and conclude that initiation and decay compete. Shorter
half-lives shift mu upward and make the decapped share larger. Discuss
this.

4.3. **Say what a geneticist does with the model.** For example, that
measured half-lives and polysome profiles systematically mix two
populations, and how the model separates them.

# 5. Ribosome drop-off and nonsense errors (to add)

5.1. **The model assumes every initiated ribosome completes the coding
sequence**, in both the capped and decapped states. This assumption is
implicit; it is never stated. Premature termination, drop-off, and
nonsense errors are not mentioned anywhere in the text, the parameter
table, or the Limitations. NMD and no-go decay appear only as excluded
endonucleolytic pathways.

5.2. **Citation.** Cope AL, Pak D, Gilchrist MA. The importance of
nonsense errors: Estimating the rates and implications of ribosome
drop-off during protein synthesis. PLoS Genet. 2026;22(6):e1012162.
doi:10.1371/journal.pgen.1012162. Preprint: bioRxiv 2024.09.05.611510.
Drop-off was already present in the per-transcript models of Gilchrist
(2005, J Theor Biol) and Gilchrist and Wagner (2006), which the draft
cites but has not entered in the bibliography.

5.3. **Where to add it**, in increasing effort.

  a. Introduction, TASEP paragraph: note that drop-off is part of the
     earlier per-transcript models and that Cope et al. provide current
     per-codon estimates from yeast ribosome profiling.
  b. Methods: state the completion assumption explicitly. Limitations:
     say drop-off is omitted, why, and its magnitude per Cope et al.
  c. Model extension: drop-off adds a transition from class i to class
     i-1 at rate i times the per-ribosome drop-off rate, in both states,
     returning a ribosome without producing protein. The system stays
     linear and tridiagonal, so `solve.tridiag` still applies. It
     modifies the exact result in 3.3, because ribosomes on decapped
     mRNAs spend longer in transit and are more exposed to drop-off.

5.4. **Expected size of the effect.** At roughly one error per ten
thousand codons and a median imax of 39 (about 350 codons), the loss per
transit is a few percent. A short paragraph giving the first-order
correction would likely suffice instead of a full re-analysis.

5.5. **Unit conversion.** Cope et al. estimate drop-off per codon; this
model's length unit is one ribosome footprint of nine codons. The
conversion is the same one the paper already uses for tau and should be
stated in the same sentence.

# 6. References: framing and missing prior art

6.1. **"Few have explored the interaction" is not defensible.** Directly
relevant models not cited:

  - Nagar A, Valleriani A, Lipowsky R. Translation by ribosomes with
    mRNA degradation: exclusion processes on aging tracks. J Stat Phys.
    2011;145:1385-1404.
  - Deneke C, Lipowsky R, Valleriani A. Effect of ribosome shielding on
    mRNA stability. Phys Biol. 2013;10:046008.
  - Chevalier C, et al. Physical modeling of ribosomes along messenger
    RNA: estimating kinetic parameters from ribosome profiling
    experiments using a ballistic model. PLoS Comput Biol.
    2023;19(10):e1011522. Includes degradation and predicts monosome
    and polysome distributions.
  - Edri S, Tuller T. Quantifying the effect of ribosomal density on
    mRNA stability. PLoS ONE. 2014;9(7):e102308.
  - Dave P, et al. Single-molecule imaging reveals translation-dependent
    destabilization of mRNAs. Mol Cell. 2023;83(4):589-606. Contains a
    model of ribosome-flux-dependent decay.

6.2. **Empirical coupling of translation and decay, 2018 to 2025, is
absent.**

  - Chan LY, et al. Non-invasive measurement of mRNA decay reveals
    translation initiation as the major determinant of mRNA stability.
    eLife. 2018;7:e32536.
  - Bicknell AA, et al. Attenuating ribosome load improves protein
    output from mRNA by limiting translation-dependent mRNA decay. Cell
    Rep. 2024;43(4):114098. This is the applied point the Discussion
    gestures at.
  - Koegel A, et al. Structural basis of mRNA decay by the human
    exosome-ribosome supercomplex. Nature. 2024;635:237-242. For the 3'
    pathway the model excludes.
  - Collart MA, Audebert L, Bushell M. Roles of the CCR4-Not complex in
    translation and dynamics of co-translation events. WIREs RNA.
    2023;15(1):e1827.
  - Hanson G, Coller J. Codon optimality, bias and usage in translation
    and mRNA decay. Nat Rev Mol Cell Biol. 2018;19:20-30.

6.3. **Plant co-translational decay literature is thin** for a paper
that uses Arabidopsis parameters.

  - Yu X, Willmann MR, Anderson SJ, Gregory BD. Genome-wide mapping of
    uncapped and cleaved transcripts reveals a role for the nuclear
    mRNA cap-binding complex in cotranslational RNA decay in
    Arabidopsis. Plant Cell. 2016;28(10):2385-2397.
  - Carpentier MC, et al. 5' to 3' cotranslational mRNA decay is
    dynamically regulated during Arabidopsis seedling development and
    fine-tunes translation efficiency. Plant Physiol.
    2020;184(3):1251-1262.
  - Deragon JM, Merret R. Co-translational mRNA decay in plants: recent
    advances and future directions. J Exp Bot. 2025.
    doi:10.1093/jxb/eraf146.
  - Wu HL, et al. Improved super-resolution ribosome profiling reveals
    prevalent translation of upstream ORFs and small ORFs in
    Arabidopsis. Plant Cell. 2024;36(3):510-539.

6.4. **Parameter sources for kappa and tau.** Beside Dao Duc and Song
(2018), cite Riba A, et al. Protein synthesis rates and ribosome
occupancies reveal determinants of translation elongation rates. PNAS.
2019;116(30):15023-15032, and Erdmann-Pham DD, Dao Duc K, Song YS. The
key parameters that govern translation efficiency. Cell Syst.
2020;10(2):183-192.

6.5. **Odd choices.** RN18 and RN19 (Wu and Tian, generic multistep
stochastic models) and RN20 (selenoprotein knockdown) are cited as
"models of mRNA degradation". Replace with the Valleriani line above and
Hanson and Coller (2018). Ensembl release 109 and R 3.6 should be
updated; data.table is cited at a 2021 version.

# 7. Figures

7.1. Fig 11 prints rho to seven digits. Use two significant figures and
state n and the test in the caption.

7.2. Fig 10 legend gives half-lives in minutes while the caption gives
mu in 1/s. Use one. The legend text "(quantile))" has an extra
parenthesis.

7.3. Fig 2 (parameter histograms) could move to the supplement. Figs 3
and 4 share one layout and would read better as a single figure.

# 8. Writing

8.1. A bad find-and-replace has damaged words: "Th such", "Thwever",
"Th m*_0", "Cois", "Coft" (for "As such", "However", "If", "This",
"Left"). Search the source for these.

8.2. The abstract's first sentence has two grammatical errors. "Novel"
appears four times. "Impressive" (Results) is not a result.

8.3. The second Introduction paragraph on NMD and no-go decay is flagged
in a margin note as out of place. Agreed; fold it into Limitations and
pair it with the drop-off discussion in Section 5.

# 9. Code verification (added in v1.1)

The R code lives in the repository `rurquidi/Ribosome` (public; 11
commits, 2020-05-10 to 2023-07-25). The solver is
`R/Polysome_functions_RAUC.R`; the paper's figures are produced by
`Text/Paper final figures.Rmd` (1700 lines), not by the scripts under
`R/`. `R/Manuscript_Figures.R` (3954 lines) is an earlier, partly
superseded version of the same material. I re-implemented the coded
model in base R (40 by 40 linear solve, no limSolve) and reproduced
the Fig 10 curves to three figures.

9.1. **The code matches the ODEs in the manuscript.** Initiation is
kappa (1 - i/imax), termination from class i is tau_c i/(9 imax) with
tau_c in codons per second, and the decapped classes are
m*_i = mu S_i / tau_i, with m*_0 = mu S_0 / delta. This confirms 3.2.

9.2. **The decapped share of protein output is exactly
9 mu imax/tau_c over 1 + 9 mu imax/tau_c**, independent of kappa. The
re-implementation gives the same share at every kappa from 1e-4 to 0.5.
This confirms 3.3. The paper should state the closed form and drop the
"never more than half" claim, which fails at the paper's own
high-decapping value (see 9.4).

9.3. **Fig 10 panel assignment.** The parameter grid is built with
`expand.grid`, so kappa varies fastest and each block of 5000 rows is
one mu value. The figure code plots these blocks:

| Panel as drawn | Rows plotted | Actual mu (1/s) | Half-life | Caption | Peak, code | Peak, figure |
|---|---|---|---|---|---|---|
| A top curve, panel C | 45001-50000 | 2.2e-4 | 52 min | 2.2e-4, 52 min | 1.00 | 1.0 |
| A middle curve, panel B | 40001-45000 | 3.3e-4 | 35 min | 1.7e-3, 7 min | 0.69 | 0.69 |
| A bottom curve, panel D | 20001-25000 | 2.3e-3 | 5 min | 5.7e-3, 2 min | 0.14 | 0.14 |

The low-decapping plot (peak 1.0) is panel C in the composite and the
mid-decapping plot (peak 0.69) is panel B, so the panels were swapped
when the figure was assembled. The decapped share at mu = 2.3e-3 is 45
percent, which is the "41 percent" quoted in the text.

9.4. **Every model figure assumes an elongation rate of one codon per
second.** All figure runs set tau = 1 and treat kappa as kappa prime,
but mu is entered in absolute units and is never scaled by tau. So
mu/tau_c is a second dimensionless group that the paper's
"collapse to one parameter" argument (Methods, Data Sources) does not
acknowledge, and the figures implicitly make yeast decapping 5 to 10
times faster relative to elongation than it is. The effect on the
headline:

| tau_c (codons/s) | Decapped share, mu = 5.7e-3, imax = 39 |
|---|---|
| 1 (as coded) | 67 percent |
| 5 | 29 percent |
| 10 | 17 percent |

At the captioned mu the share exceeds one half, contradicting the
Discussion. At realistic yeast elongation rates (Dao Duc and Song 2018,
Riba et al. 2019) it is a sixth to a third. The validation in Fig 11,
by contrast, uses gene-specific tau_c and kappa from Dao Duc and Song,
so the model figures and the validation are not on the same footing.
Fix: present results in terms of the two dimensionless groups
kappa/tau_c and 9 mu imax/tau_c, or rerun the figures with a
realistic tau_c and mu scaled by it.

9.5. **Fig 11 is coded as described.** Empirical load is RPF RPKM over
mRNA RPKM times length over 200, with length taken as 27 imax
nucleotides; model load uses per-gene kappa and tau_c from Dao Duc and
Song and mu from Presnyak; the test is Spearman on log10 values. The
reporting issues in 4.1 and 7.1 stand.

9.6. **Protein production in the code is total bound ribosomes**
(sum of i m_i), which is proportional to the production rate only at
fixed imax. That is correct for Fig 10 but would be wrong if reused
across the imax range.

9.7. **Reproducibility.** The input data directory was not pushed (see
1.5), limSolve is required, the solver exists in three near-duplicate
variants (`CalcMarkedClass`, `_newtau`, `NP`) of which the figure code
calls the first, `CalcMarkedClass` recomputes the capped solution with
`CalcUnmarkedClass` rather than the `_newtau` variant when none is
passed, and the figure code selects parameter sets by hard-coded row
ranges. None of this blocks publication, but the repository needs a
README that maps each figure to its chunk and parameter block, a
license, and the data files before it is cited in a submission.

# Sources consulted

PLOS Genetics submission guidelines and journal information pages
(journals.plos.org/plosgenetics), PubMed, bioRxiv, the publisher
pages for each reference listed above, and the `rurquidi/Ribosome`
repository at commit f70c00a (cloned to `~/Repositories/Ribosome` on
obed, 2026-10-09).
