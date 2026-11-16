#!/usr/bin/env python3
"""Weighted-burden discovery report: gene-level discovery and LoF-like attribution from biobank-wide results.

Two questions per gene x trait:
  1. Discovery    - is the gene associated under the combined pLoF;missenseLC burden?
  2. Attribution  - when weighted missense variants drive an association, does pLoF in the same gene act in
                    the same direction ("LoF-like")?
Attribution is asked only of genes selected on the missense-only test. pLoF and missense carriers are, in
practice, different people, so the pLoF test is independent of the selection and its p is read directly,
with no selection adjustment.

Minimum pLoF allele count (--plof-min-mac, default 5): a pLoF test on fewer than 5 minor alleles is
treated as not run and the gene is classed pLoF_MAC_lt5. Such a test rests on one to four carriers; on
height, two reached p < 2.5e-6 with a single allele each. The pipeline does not remove them (saige-slim
tests a group at any allele count, and the production v8 SKAT-O output also contains them), so the rule
lives here. It applies to real and null phenotypes and to both runs. It does not touch the missense or
combined p-values, so it changes the evidence and attribution labels of a hit, never the hit list. The
novelty label still uses the unfiltered production pLoF p.

Minimum effective carrier count (--min-neff, default 10): a FlexRV weight set whose burden test has
nEff_Burden below 10 rests on a handful of carriers, and on quantitative traits such sets produce
combined-mask p-values of 1e-8 and below that neither the pLoF nor the missense test supports. They are
dropped and each group's CCT, direction summaries and top test are recomputed, with equal weights, over
the weight sets that remain; a group left with no set counts as not tested. Tests without nEff_Burden
(the Beta(1,25) pLoF test) are kept. The per-group table holds only the CCT over all sets, so the rule
is applied from the per-weight-set exports next to it (`gene_results_<max_MAF>.txt.bgz`, written for
flexRV_native runs; about 100-200 MB per phenotype, which is what makes the report slower than with
--min-neff 0). A run without those exports (other gene modes, flat copies) is used as is, and the log says
so. The meta-analysis table is never refiltered: --create-flexrv-meta applied the same rule
(--meta-min-neff) in each ancestry before combining them. It also blanks pLoF tests on fewer than 5 pooled
alleles itself (--meta-plof-min-mac; the value stays in Pvalue_meta_cct_unfiltered, which this report
does not read), so on a meta table --plof-min-mac removes nothing further at the default.

Run
---
    python utils/discovery.py \\
        --weighted-root   gs://.../ht_results/<map_source>/<mode>/<WEIGHT>/EUR \\
        --production-root <ht_results ancestry root of the unweighted run>    # optional: novelty, comparison
        --out discovery_EUR
    python utils/discovery.py \\
        --weighted-root   gs://.../ht_results/<map_source>/flexRV_native/<WEIGHT>/META \\
        --out discovery_META                                                  # cross-ancestry meta-analysis
    python utils/discovery.py --out discovery_EUR --figures-only              # redraw the figures


Outputs (--out)
---------------
  report_log.txt                         every number of the report
  h_hits_weighted.tsv                    one row per weighted hit: evidence, novelty, attribution call
  d_attribution_{weighted,production}.tsv  step 3, one row per gene x trait in the attribution set
  g_joint_only_weighted.tsv              joint-only hits (descriptive)
  b_calibration_null.tsv, c_hits_by_evidence.tsv, c_matched_fpr.tsv, f_concordance_by_missense_bin.tsv,
  z_data_checks.tsv                      calibration, discovery, concordance, pLoF-vs-production check
  phenotypes.tsv, gene_table_null_*.parquet, candidates_*.parquet
  figures/fig1_discovery                 hits by evidence pattern, weighted vs production; matched-FPR
  figures/fig4_attribution               attribution set: missense CCT vs directional pLoF p, by call
  figures/fig6_concordance               pLoF-in-missense-direction share by missense bin, real vs null
"""
import argparse
import os
import re
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import numpy as np
import pandas as pd
from scipy import stats

ALPHA = 2.5e-6
CCT_FILE = "gene_results_cauchy_by_group.tsv.gz"
# Cross-ancestry FlexRV meta-analysis table (--create-flexrv-meta); read when a phenotype directory
# has no per-ancestry table, with its columns renamed to the per-ancestry names. One "test" of a
# group is then one weight set meta-analysed across ancestries (signed Stouffer), MAC is summed over
# the ancestries, and n_cases / n_controls are totals. It has no BETA_top / SE_top (Z_top instead).
META_FILE = "gene_results_flexrv_meta.tsv.gz"
# Cross-ancestry meta-analysis of the unweighted production run (--create-meta-analysis-ht export):
# one row per gene x mask x max MAF tier with the signed Stouffer burden statistic META_Stats_Burden and
# its two-sided Pvalue_Burden. Read last, and turned into the per-group layout by production_meta_cct().
PROD_META_FILE = "gene_results.txt.bgz"
META_RENAME = {"Pvalue_meta_cct": "Pvalue_cct", "Pvalue_meta_cct_pos": "Pvalue_cct_pos",
               "Pvalue_meta_cct_neg": "Pvalue_cct_neg", "n_sets": "n_tests"}
GROUPS = {"pLoF": "plof", "missenseLC": "mis", "pLoF;missenseLC": "comb"}
GROUP_VALUES = ["MAC", "n_tests", "Pvalue_cct", "Pvalue_cct_pos", "Pvalue_cct_neg", "n_pos", "n_neg",
                "n_pos_nom", "n_neg_nom", "BETA_top", "SE_top", "Z_top", "max_MAF_top"]
# without these a table cannot be used (a meta table written before the direction summaries were added)
REQUIRED = ["Pvalue_cct", "Pvalue_cct_pos", "Pvalue_cct_neg", "n_pos", "n_neg", "n_pos_nom", "n_neg_nom",
            "MAC", "n_tests"]
META_COLS = ["trait_type", "n_cases", "n_controls"]
# binary null phenotypes with at least this case fraction form the empirical null of the directional pLoF
# p (nominal prevalence >= 1%; the 0.5 margin absorbs sampling, e.g. 2265 / 227273 = 0.00997)
NULL_REF_MIN_PREV = 0.005
# binary null phenotypes used for the matched-FPR thresholds (drops prevalence 1e-4 nulls)
NULL_FPR_MIN_PREV = 5e-4
FPRS = (1e-4, 1e-3)
MB = [0, None, 1e-4, 1e-2, 0.1, 1.0000001]          # missense CCT bins; None is replaced by alpha
MBL = ["<alpha", "alpha-1e-4", "1e-4-0.01", "0.01-0.1", ">0.1"]
# A pLoF test on fewer than this many minor alleles is treated as not run (class pLoF_MAC_lt5). SAIGE
# itself may emit such tests (saige-slim tests a group at any allele count; production v8 SKAT-O output
# also has them), so the rule is applied here, to real and null phenotypes and to both runs alike.
PLOF_MIN_MAC = 5
# FlexRV weight sets with nEff_Burden below this are dropped before the per-group CCT (see the module
# docstring); same value as FLEXRV_META_MIN_NEFF of the pipeline's cross-ancestry meta-analysis.
MIN_NEFF = 10.0
TIER_FILE = re.compile(r"gene_results_[0-9][0-9.eE+-]*\.txt\.bgz")     # per-weight-set export of one max MAF tier
SIGN_P = 0.05                                                         # as in the pipeline's per-group table
ATTR_SET_CALLS = ["LoF_like_strict", "LoF_like_suggestive", "not_LoF_like"]


# ------------------------------------------------------------------------------------------------------
# Inputs
# ------------------------------------------------------------------------------------------------------
def _fs():
    import hailtop.fs as hfs
    return hfs


def list_cct_files(root: str) -> Dict[str, str]:
    """{phenotype: path} under a results root.

    Accepts an ht_results ancestry root (`phenotype_<X>/gene_results_cauchy_by_group.tsv.gz`), the META
    root of the cross-ancestry FlexRV meta-analysis (`phenotype_<X>/gene_results_flexrv_meta.tsv.gz`), the
    META root of the production run (`phenotype_<X>/gene_results.txt.bgz`), or a flat directory of
    `<X>.tsv.gz` copies. Existence of the per-phenotype file is checked when it is read; under
    `phenotype_<X>/` the per-ancestry table is tried first, then the two meta-analysis tables.
    """
    out = {}
    if root.startswith("gs://"):
        names = [e.path.rstrip("/") for e in _fs().ls(root)]
    else:
        names = [os.path.join(root, n) for n in os.listdir(root)]
    for p in sorted(names):
        b = os.path.basename(p)
        if b.startswith("phenotype_"):
            out[b[len("phenotype_"):]] = f"{p}/{CCT_FILE}|{p}/{META_FILE}|{p}/{PROD_META_FILE}"
        elif b.endswith(".tsv.gz"):
            out[b[:-len(".tsv.gz")]] = p
    return out


def read_cct(path: str, nrows: int = None) -> Optional[pd.DataFrame]:
    """One phenotype's per-group CCT table, restricted to the columns used; None if absent."""
    want = set(["phenoname", "gene_id", "gene_symbol", "Group"] + GROUP_VALUES + META_COLS) | set(META_RENAME)
    kw = dict(sep="\t", usecols=lambda c: c in want, dtype={"phenoname": str, "trait_type": str}, nrows=nrows)
    for one in path.split("|"):
        if os.path.basename(one) == PROD_META_FILE:
            c = production_meta_cct(one, nrows)
            if c is None:
                continue
            c.attrs["is_meta"], c.attrs["directory"] = "production", os.path.dirname(one)
            return c
        try:
            if one.startswith("gs://"):
                with _fs().open(one, "rb") as f:
                    c = pd.read_csv(f, compression="gzip", **kw)
            else:
                c = pd.read_csv(one, **kw)
        except FileNotFoundError:
            continue
        c = c.rename(columns=META_RENAME)
        c.attrs["is_meta"] = os.path.basename(one) == META_FILE
        c.attrs["directory"] = os.path.dirname(one) if os.path.basename(one) in (META_FILE, CCT_FILE) else None
        miss = [k for k in REQUIRED if k not in c.columns]
        if miss:
            raise ValueError(f"{one} lacks {', '.join(miss)}: " + (
                "this meta-analysis table predates the direction summaries, rerun --create-flexrv-meta "
                "--overwrite-meta" if os.path.basename(one) == META_FILE else
                "rerun the load (or the group-Cauchy backfill) of this phenotype"))
        return c
    return None


def _cauchy_t(p):
    """Per-test Cauchy contribution tan((0.5-p)pi), with the 1/(p pi) small-p limit (as in the pipeline)."""
    p = np.clip(np.asarray(p, dtype=float), 1e-300, 1 - 1e-15)
    small = p < 1e-16
    return np.where(small, 1 / (np.where(small, p, 1) * np.pi), np.tan((0.5 - p) * np.pi))


def _cauchy_p(T):
    """Equal-weight CCT p from the mean Cauchy statistic."""
    T = np.asarray(T, dtype=float)
    return np.where(T > 1e15, 1 / (np.where(T > 1e15, T, 1) * np.pi), stats.cauchy.sf(T))


def read_weight_sets(directory: str) -> Optional[pd.DataFrame]:
    """Per-weight-set burden tests of one phenotype (flexRV_native `gene_results_<max_MAF>.txt.bgz`
    exports in its ht_results directory); None if the directory has none."""
    try:
        if directory.startswith("gs://"):
            names = [os.path.basename(e.path.rstrip("/")) for e in _fs().ls(directory)]
        else:
            names = os.listdir(directory)
    except FileNotFoundError:
        return None
    want = {"gene_id", "gene_symbol", "Run_group", "Weight", "max_MAF", "Pvalue_Burden", "Pvalue_Burden_log10",
            "BETA_Burden", "SE_Burden", "nEff_Burden", "MAC"}
    kw = dict(sep="\t", usecols=lambda c: c in want, compression="gzip", na_values=["NA"],
              dtype={"Run_group": "category", "Weight": "category"})
    parts = []
    for n in sorted(n for n in names if TIER_FILE.fullmatch(n)):
        if directory.startswith("gs://"):
            with _fs().open(f"{directory}/{n}", "rb") as f:
                parts.append(pd.read_csv(f, **kw))
        else:
            parts.append(pd.read_csv(os.path.join(directory, n), **kw))
    if not parts or not {"Run_group", "nEff_Burden"} <= set(parts[0].columns):
        return None
    d = pd.concat(parts, ignore_index=True)
    for c in ("Run_group", "Weight"):
        d[c] = d[c].astype(str)
    return d


def production_meta_cct(path: str, nrows: int = None) -> Optional[pd.DataFrame]:
    """Per-group table from the production run's cross-ancestry meta-analysis export; None if the file
    is absent or is not that export (no META_Stats_Burden column).

    One burden test per mask and max MAF tier: p = Pvalue_Burden, direction = sign of META_Stats_Burden
    (BETA_top holds that Z, SE_top is 1), MAC = alleles summed over the ancestries. As in the per-ancestry
    production tables, a tier with the same MAC as a lower tier of the same mask is the same variant set
    and enters the CCT once.
    """
    want = {"gene_id", "gene_symbol", "annotation", "max_MAF", "MAC", "META_Stats_Burden", "Pvalue_Burden",
            "Pvalue_Burden_log10"}
    kw = dict(sep="\t", usecols=lambda c: c in want, compression="gzip", na_values=["NA"], nrows=nrows)
    try:
        if path.startswith("gs://"):
            with _fs().open(path, "rb") as f:
                d = pd.read_csv(f, **kw)
        else:
            d = pd.read_csv(path, **kw)
    except FileNotFoundError:
        return None
    if "META_Stats_Burden" not in d.columns:
        return None
    d = d[(d.annotation != "Cauchy") & (d.max_MAF > 0) & d.META_Stats_Burden.notna()]
    d = d.sort_values("max_MAF", kind="stable").drop_duplicates(["gene_id", "annotation", "MAC"])
    d = d.rename(columns={"annotation": "Run_group", "META_Stats_Burden": "BETA_Burden"}).assign(
        Weight="Beta_1_25", SE_Burden=1.0, nEff_Burden=np.nan)
    return neff_floor_cct(d, 0)[0]


def neff_floor_cct(d: pd.DataFrame, min_neff: float) -> (pd.DataFrame, dict):
    """Per-group table (columns of the pipeline's gene_results_cauchy_by_group.tsv.gz that the report
    reads) from read_weight_sets() rows, over the weight sets with nEff_Burden >= min_neff.

    Tests without nEff_Burden are kept. A gene tested in two overlapping interval chunks has two copies
    of each test; the copy with the most alleles is used (the pipeline keeps the chunk with the most).
    Returns the table and the counts of tests read / dropped.
    """
    p = d.Pvalue_Burden.to_numpy(dtype=float)
    if "Pvalue_Burden_log10" in d.columns:                     # p below the float range of the export
        under = (p == 0) & d.Pvalue_Burden_log10.notna().to_numpy()
        p = np.where(under, np.power(10.0, -d.Pvalue_Burden_log10.to_numpy(dtype=float)), p)
    d = d.assign(p=p)
    d = d[d.p.notna()]
    key = ["gene_id", "Run_group", "Weight", "max_MAF"]
    if d.duplicated(key).any():
        d = d.sort_values("MAC", ascending=False, kind="stable").drop_duplicates(key)
    sparse = (d.nEff_Burden < min_neff).to_numpy()
    counts = dict(tests=len(d), dropped=int(sparse.sum()),
                  groups=d.groupby(["gene_id", "Run_group"], observed=True).ngroups)
    d = d[~sparse]
    p, sg = d.p.to_numpy(), np.sign(d.BETA_Burden.to_numpy(dtype=float))
    nom = p < SIGN_P
    t = pd.DataFrame({
        "gene_id": d.gene_id.to_numpy(), "Group": d.Run_group.to_numpy(), "p": p, "MAC": d.MAC.to_numpy(),
        "_t": _cauchy_t(p), "_pos": sg > 0, "_neg": sg < 0, "_pos_nom": (sg > 0) & nom, "_neg_nom": (sg < 0) & nom,
        "_t_pos": _cauchy_t(np.where(sg > 0, p / 2, 1 - p / 2)),
        "_t_neg": _cauchy_t(np.where(sg < 0, p / 2, 1 - p / 2)),
    }, index=d.index)
    g = t.groupby(["gene_id", "Group"])
    out = g.agg(n_tests=("_t", "size"), MAC=("MAC", "max"), cct=("_t", "mean"),
                n_pos=("_pos", "sum"), n_neg=("_neg", "sum"), n_pos_nom=("_pos_nom", "sum"),
                n_neg_nom=("_neg_nom", "sum"), t_pos=("_t_pos", "mean"), t_neg=("_t_neg", "mean"))
    out["Pvalue_cct"] = _cauchy_p(out.cct)
    out["Pvalue_cct_pos"] = _cauchy_p(out.t_pos)
    out["Pvalue_cct_neg"] = _cauchy_p(out.t_neg)
    top = d.loc[g.p.idxmin()]
    top = top.set_index([top.gene_id, top.Run_group.rename("Group")])
    out["BETA_top"], out["SE_top"], out["max_MAF_top"] = top.BETA_Burden, top.SE_Burden, top.max_MAF
    out = out.drop(columns=["cct", "t_pos", "t_neg"]).reset_index()
    counts["groups_lost"] = counts["groups"] - len(out)
    sym = d.drop_duplicates("gene_id").set_index("gene_id").gene_symbol
    out.insert(1, "gene_symbol", out.gene_id.map(sym))
    return out, counts


def cct_meta(c: pd.DataFrame) -> dict:
    """trait_type / n_cases / n_controls written by the pipeline (empty for older tables)."""
    return {k: c[k].iloc[0] for k in META_COLS if k in c.columns and len(c) and pd.notna(c[k].iloc[0])}


def gene_table(c: pd.DataFrame, pheno: str, plof_min_mac: float = PLOF_MIN_MAC) -> pd.DataFrame:
    """Per-gene table of one phenotype from its per-group CCT table.

      p_{plof,mis,comb}_cct  Pvalue_cct of pLoF (Beta_1_25 only), missenseLC, pLoF;missenseLC
      s_mis        shared sign of the missense tests with p < 0.05 (all tests if none is nominal); 0 = ambiguous
      p_plof_dir   pLoF Pvalue_cct_pos if s_mis = +1, Pvalue_cct_neg if s_mis = -1
      plof_mac_lt5 pLoF carriers exist but there is no usable pLoF test: SAIGE emitted none (combined MAC >
                   missense MAC), or its MAC is below `plof_min_mac`
      p_plof_cct_all  pLoF Pvalue_cct before the `plof_min_mac` rule (used only for the novelty label)

    pLoF tests with MAC < `plof_min_mac` are set to missing in p_plof_cct and p_plof_dir (0 = keep all).
    The missense and combined p-values are never filtered.
    """
    c = c[c.Group.isin(GROUPS)]
    w = c[["gene_id", "gene_symbol"]].drop_duplicates("gene_id")
    for grp, tag in GROUPS.items():
        x = c[c.Group == grp].drop_duplicates("gene_id").set_index("gene_id")
        x = x[[v for v in GROUP_VALUES if v in x.columns]].add_suffix(f"_{tag}")
        w = w.merge(x, left_on="gene_id", right_index=True, how="left")
    any_nom = (w.n_pos_nom_mis + w.n_neg_nom_mis) > 0
    pos = np.where(any_nom, w.n_pos_nom_mis, w.n_pos_mis)
    neg = np.where(any_nom, w.n_neg_nom_mis, w.n_neg_mis)
    s = np.where((pos > 0) & (neg == 0), 1, np.where((neg > 0) & (pos == 0), -1, 0))
    w = w.reset_index(drop=True)

    def v(name):
        return w[name].to_numpy() if name in w.columns else np.full(len(w), np.nan)

    out = pd.DataFrame({
        "pheno": pheno, "gene_id": w.gene_id.values, "symbol": w.gene_symbol.values,
        "gene": (w.gene_id + "_" + w.gene_symbol.astype(str)).values,
        "p_plof_cct": w.Pvalue_cct_plof.values, "p_mis_cct": w.Pvalue_cct_mis.values,
        "p_comb_cct": w.Pvalue_cct_comb.values,
        "s_mis": np.where(w.n_tests_mis.notna(), s, np.nan),
        "n_sign_tests": np.where(any_nom, w.n_pos_nom_mis + w.n_neg_nom_mis, w.n_pos_mis + w.n_neg_mis),
        "MAC_plof": w.MAC_plof.values, "MAC_mis": w.MAC_mis.values, "n_plof_tests": w.n_tests_plof.values,
        "BETA_top_comb": v("BETA_top_comb"), "SE_top_comb": v("SE_top_comb"),
        "max_MAF_top_comb": v("max_MAF_top_comb"),
    })
    # Z of the most significant combined test: BETA / SE per ancestry, the Stouffer Z in a meta-analysis
    out["Z_top_comb"] = v("Z_top_comb") if "Z_top_comb" in w.columns else out.BETA_top_comb / out.SE_top_comb
    out["p_plof_dir"] = np.select([out.s_mis == 1, out.s_mis == -1],
                                  [w.Pvalue_cct_pos_plof.values, w.Pvalue_cct_neg_plof.values], np.nan)
    out["p_plof_cct_all"] = out.p_plof_cct
    low = (out.MAC_plof < plof_min_mac).to_numpy()
    out.loc[low, ["p_plof_cct", "p_plof_dir"]] = np.nan
    out["plof_mac_lt5"] = out.p_plof_cct.isna() & (low | ((w.MAC_comb - w.MAC_mis.fillna(0)).values > 0))
    return out


def cct_loader(root: str, plof_min_mac: float = PLOF_MIN_MAC,
               min_neff: float = MIN_NEFF) -> Callable[[str], Optional[pd.DataFrame]]:
    """pheno -> gene_table, reading `root`. load.meta(pheno) reads only the pipeline-written metadata.

    min_neff > 0: per-ancestry tables are rebuilt from the per-weight-set exports with the nEff_Burden
    floor (neff_floor_cct). load.neff counts what the floor did; load.neff["unfiltered"] lists the
    phenotypes without per-weight-set exports, used as written. Meta-analysis tables are used as written.
    """
    files = list_cct_files(root)
    neff = dict(min_neff=min_neff, phenotypes=0, tests=0, dropped=0, groups=0, groups_lost=0, unfiltered=[],
                meta=0, production_meta=0)

    def load(pheno):
        c = read_cct(files[pheno]) if pheno in files else None
        if c is None:
            return None
        if c.attrs.get("is_meta"):
            neff["production_meta" if c.attrs["is_meta"] == "production" else "meta"] += 1
        elif min_neff > 0:
            sets = read_weight_sets(c.attrs["directory"]) if c.attrs.get("directory") else None
            if sets is None:
                neff["unfiltered"].append(pheno)
            else:
                f, n = neff_floor_cct(sets, min_neff)
                for k in META_COLS + ["phenoname"]:
                    if k in c.columns and len(c):
                        f[k] = c[k].iloc[0]
                c = f
                neff["phenotypes"] += 1
                for k in ("tests", "dropped", "groups", "groups_lost"):
                    neff[k] += n[k]
        return gene_table(c, pheno, plof_min_mac)

    def meta(pheno):
        c = read_cct(files[pheno], nrows=1) if pheno in files else None
        return {} if c is None else cct_meta(c)

    load.phenotypes = list(files)
    load.meta = meta
    load.neff = neff
    return load


# ------------------------------------------------------------------------------------------------------
# Phenotype metadata
# ------------------------------------------------------------------------------------------------------
def is_null_name(pheno: str) -> bool:
    return pheno.startswith("random_")


def read_manifest(path: str, ancestry: Optional[str] = None) -> pd.DataFrame:
    """phenoname, trait_type (binary / continuous / quantitative), n_cases, n_controls[, is_null, ancestry]."""
    m = pd.read_csv(path, dtype={"phenoname": str})
    if ancestry and "ancestry" in m.columns:
        m = m[m.ancestry.str.upper() == ancestry.upper()]
    m = m.drop_duplicates("phenoname").set_index("phenoname")
    if "is_null" not in m.columns:
        m["is_null"] = [is_null_name(p) or c == "random_phenotype"
                        for p, c in zip(m.index, m.get("category", pd.Series("", index=m.index)))]
    return m[["trait_type", "n_cases", "n_controls", "is_null"]]


def classify(meta: pd.DataFrame) -> pd.DataFrame:
    """Adds kind (quantitative / binary / null_continuous / null_binary), binary, case_frac,
    and the role of each binary null in the held-out empirical-p check (fit / test)."""
    m = meta.copy()
    binary = m.trait_type.astype(str).str.lower().eq("binary")
    m["binary"] = binary
    m["kind"] = np.where(m.is_null, np.where(binary, "null_binary", "null_continuous"),
                         np.where(binary, "binary", "quantitative"))
    nc, nn = pd.to_numeric(m.n_cases, errors="coerce"), pd.to_numeric(m.n_controls, errors="coerce")
    m["case_frac"] = np.where(binary, nc / (nc + nn), np.nan)
    pool = m[(m.kind == "null_binary") & (m.case_frac >= NULL_REF_MIN_PREV)].sort_values("case_frac")
    m["null_ref_pool"] = m.index.isin(pool.index)
    m["null_ref_role"] = ""
    m.loc[pool.index[0::2], "null_ref_role"] = "fit"      # alternate by prevalence: fit 1/5/20%, test 2/10/50%
    m.loc[pool.index[1::2], "null_ref_role"] = "test"
    m["null_fpr_pool"] = (m.kind == "null_continuous") | ((m.kind == "null_binary") &
                                                          (m.case_frac >= NULL_FPR_MIN_PREV))
    return m


# ------------------------------------------------------------------------------------------------------
# Building blocks (steps 2-3)
# ------------------------------------------------------------------------------------------------------
def null_ref(nulls: pd.DataFrame, phenos) -> dict:
    """Binary empirical reference: sorted directional pLoF p of binary-null genes, one array per direction."""
    n = nulls[nulls.pheno.isin(phenos) & nulls.p_plof_dir.notna()]
    return {s: np.sort(n.loc[n.s_mis == s, "p_plof_dir"].values) for s in (1, -1)}


def emp_p(p, s, ref):
    """(1 + #null <= p) / (1 + n_null), using the null genes with the same missense direction."""
    p = np.asarray(p, float)
    s = np.asarray(s)
    out = np.full(len(p), np.nan)
    for sv in (1, -1):
        r = ref[sv]
        if len(r) == 0:
            continue
        k = (s == sv) & ~np.isnan(p)
        out[k] = (np.searchsorted(r, p[k], side="right") + 1) / (len(r) + 1)
    return out


def novelty(d: pd.DataFrame, alpha: float, has_ref: bool):
    if not has_ref:
        return np.full(len(d), "NA", dtype=object)
    return np.select([d.PROD_plof < alpha, d.PROD_plof_mis < alpha, d.PROD_mis < alpha],
                     ["prod_pLoF", "prod_pLoF;missenseLC", "prod_missenseLC"], "new")


def evidence(d: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Steps 1-2: combined-mask hits on real phenotypes, by which single-mask CCTs are < alpha."""
    H = d[~d.kind.str.startswith("null") & (d.p_comb_cct < alpha)].copy()
    ps, ms = H.p_plof_cct < alpha, H.p_mis_cct < alpha
    H["evidence"] = np.select([ps & ms, ps, ms], ["both", "pLoF_only", "missense_only"], "joint_only")
    return H


def attribution(d: pd.DataFrame, ref: Optional[dict], alpha: float, has_ref: bool):
    """Step 3: real-phenotype genes with missense CCT < alpha; one directional pLoF test each.

    Quantitative: p_plof_dir as is. Binary: empirical p against the binary nulls with the same missense
    direction (raw p, flagged in p_attr_source, when no binary null reference is available)."""
    A = d[~d.kind.str.startswith("null") & (d.p_mis_cct < alpha)].copy()
    use_emp = A.binary.values & (ref is not None)
    A["p_attr"] = A.p_plof_dir.values
    if ref is not None:
        A.loc[use_emp, "p_attr"] = emp_p(A.p_plof_dir.values[use_emp], A.s_mis.values[use_emp], ref)
    A["p_attr_source"] = np.where(use_emp, "binary_empirical", np.where(A.binary, "binary_raw", "raw"))
    has_plof = A.p_plof_cct.notna()
    ambiguous = has_plof & (A.s_mis == 0)
    tested = has_plof & ~ambiguous & A.p_attr.notna()
    n_test = int(tested.sum())
    A["tested"] = tested
    A["attr"] = np.select(
        [~has_plof & A.plof_mac_lt5, ~has_plof, ambiguous,
         tested & (A.p_attr < 0.05 / max(n_test, 1)), tested & (A.p_attr < 0.05)],
        ["pLoF_MAC_lt5", "no_pLoF_test", "direction_ambiguous", "LoF_like_strict", "LoF_like_suggestive"],
        "not_LoF_like")
    A["novelty"] = novelty(A, alpha, has_ref)
    return A, n_test


# ------------------------------------------------------------------------------------------------------
# Streaming over phenotypes
# ------------------------------------------------------------------------------------------------------
@dataclass
class Collected:
    nulls: Dict[str, list] = field(default_factory=dict)       # method -> [null gene tables]
    cand: Dict[str, list] = field(default_factory=dict)        # method -> [candidate rows of real phenos]
    conc: Dict[str, list] = field(default_factory=dict)        # method -> [concordance counts]
    fpr: list = field(default_factory=list)                    # matched-FPR hit counts
    check: list = field(default_factory=list)                  # weighted vs production pLoF sample
    missing: list = field(default_factory=list)


def _annotate(d, pheno, mrow):
    d = d.assign(kind=mrow.kind, binary=bool(mrow.binary), case_frac=mrow.case_frac)
    d["E_plof"] = d.MAC_plof * d.case_frac
    return d


def _cross(tables: dict) -> dict:
    """Production p-values (novelty) and the other method's combined p onto each method's table."""
    prod = tables.get("production")
    if prod is not None:
        # novelty uses the production pLoF p before the minimum-MAC rule: a gene production already
        # reports is not "new", however few alleles its test rests on
        ref = prod[["gene_id", "p_plof_cct_all", "p_mis_cct", "p_comb_cct"]].drop_duplicates("gene_id").rename(
            columns={"p_plof_cct_all": "PROD_plof", "p_mis_cct": "PROD_mis", "p_comb_cct": "PROD_plof_mis"})
    out = {}
    for m, d in tables.items():
        if prod is not None:
            d = d.merge(ref, on="gene_id", how="left")
            other = [o for o in tables if o != m]
            if other:
                o = tables[other[0]][["gene_id", "p_comb_cct"]].drop_duplicates("gene_id")
                d = d.merge(o.rename(columns={"p_comb_cct": "other_p_comb_cct"}), on="gene_id", how="left")
        else:
            d = d.assign(PROD_plof=np.nan, PROD_mis=np.nan, PROD_plof_mis=np.nan)
        out[m] = d
    return out


def _conc_counts(d, alpha):
    x = d[d.p_plof_dir.notna() & (d.s_mis != 0)]
    mb = list(MB)
    mb[1] = alpha
    b = pd.cut(x.p_mis_cct, mb, labels=MBL, right=False)
    g = (x.p_plof_dir < 0.5).groupby(b, observed=True).agg(["size", "sum"]).reset_index()
    return g.rename(columns={"p_mis_cct": "mbin"}).assign(
        tt="binary" if bool(d.binary.iloc[0]) else "quantitative", real=not d.kind.iloc[0].startswith("null"))


def collect(loaders: Dict[str, Callable], meta: pd.DataFrame, alpha: float = ALPHA, keep_p: float = None,
            check_n: int = 5000, log=print) -> (Collected, dict):
    """Two passes over the phenotypes in `meta`: nulls (kept in full), then real phenotypes (reduced).

    Returns the collected pieces and the matched-FPR thresholds {(method, trait_type, fpr): threshold}."""
    keep_p = alpha if keep_p is None else max(keep_p, alpha)
    col = Collected()
    for m in loaders:
        col.nulls[m], col.cand[m], col.conc[m] = [], [], []

    def load(pheno):
        tables = {}
        for m, f in loaders.items():
            d = f(pheno)
            if d is not None and len(d):
                tables[m] = d
        if "weighted" not in tables:
            col.missing.append(pheno)
            return None
        return {m: _annotate(d, pheno, meta.loc[pheno]) for m, d in _cross(tables).items()}

    nulls = [p for p in meta.index if meta.loc[p, "kind"].startswith("null")]
    reals = [p for p in meta.index if p not in set(nulls)]
    log(f"phenotypes: {len(nulls)} null, {len(reals)} real")
    for i, p in enumerate(nulls):
        t = load(p)
        if t is None:
            continue
        for m, d in t.items():
            col.nulls[m].append(d)
            col.conc[m].append(_conc_counts(d, alpha).assign(method=m))
    thr = {}
    for m in loaders:
        n = pd.concat(col.nulls[m]) if col.nulls[m] else None
        for tt, kind in [("quantitative", "null_continuous"), ("binary", "null_binary")]:
            if n is None:
                continue
            ok = meta.index[meta.null_fpr_pool & (meta.kind == kind)]
            nl = n[n.pheno.isin(ok)].p_comb_cct.dropna()
            for f in FPRS:
                if len(nl):
                    thr[(m, tt, f)] = (float(np.quantile(nl, f)), len(nl))

    rng = np.random.default_rng(0)
    for i, p in enumerate(reals):
        if i % 100 == 0:
            log(f"  real phenotypes read: {i}/{len(reals)}")
        t = load(p)
        if t is None:
            continue
        for m, d in t.items():
            col.cand[m].append(d[(d.p_comb_cct < keep_p) | (d.p_mis_cct < keep_p)])
            col.conc[m].append(_conc_counts(d, alpha).assign(method=m))
            tt = "binary" if bool(d.binary.iloc[0]) else "quantitative"
            for f in FPRS:
                if (m, tt, f) in thr:
                    col.fpr.append(dict(trait_type=tt, fpr=f, method=m, pheno=p,
                                        real_hits=int((d.p_comb_cct < thr[(m, tt, f)][0]).sum())))
        w = t["weighted"]
        x = w.dropna(subset=["p_plof_cct", "PROD_plof"])
        if len(x):
            idx = rng.choice(len(x), size=min(check_n, len(x)), replace=False) if check_n < len(x) else slice(None)
            col.check.append(x.iloc[idx][["pheno", "gene_id", "symbol", "kind", "p_plof_cct", "PROD_plof"]])
    # nulls go into the U2 check too
    for d in col.nulls.get("weighted", []):
        x = d.dropna(subset=["p_plof_cct", "PROD_plof"])
        if len(x):
            col.check.append(x[["pheno", "gene_id", "symbol", "kind", "p_plof_cct", "PROD_plof"]])
    if col.missing:
        log(f"no weighted CCT table for {len(col.missing)} phenotypes: {', '.join(col.missing[:20])}"
            + (" ..." if len(col.missing) > 20 else ""))
    return col, thr
class Log:
    def __init__(self, path):
        self.f = open(path, "w")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(s + "\n")
        self.f.flush()


def build_meta(args, loaders):
    """Phenotype metadata: the manifest if given, else the columns written by the pipeline."""
    phenos = loaders["weighted"].phenotypes
    if args.phenotypes:
        keep = set(args.phenotypes.split(","))
        phenos = [p for p in phenos if p in keep]
    if args.manifest:
        m = read_manifest(args.manifest, args.ancestry)
        miss = [p for p in phenos if p not in m.index]
        if miss:
            print(f"WARNING: {len(miss)} phenotypes not in the manifest are skipped: {', '.join(miss[:20])}")
        return classify(m.loc[[p for p in phenos if p in m.index]])
    # no manifest: first row of every table for its trait_type / n_cases / n_controls
    rows = {}
    for p in phenos:
        md = loaders["weighted"].meta(p)
        if "trait_type" not in md:
            print(f"WARNING: {p}: no trait_type column (table predates the pipeline metadata); pass --manifest")
            continue
        rows[p] = dict(trait_type=md["trait_type"], n_cases=md.get("n_cases"), n_controls=md.get("n_controls"),
                       is_null=is_null_name(p))
    return classify(pd.DataFrame.from_dict(rows, orient="index"))


def report(col: Collected, thr: dict, meta: pd.DataFrame, out: str, alpha: float, say):
    methods = [m for m in ("weighted", "production") if col.nulls.get(m) or col.cand.get(m)]
    has_prod = "production" in methods
    nulls = {m: pd.concat(col.nulls[m], ignore_index=True) if col.nulls[m] else None for m in methods}
    cand = {m: pd.concat(col.cand[m], ignore_index=True) for m in methods}
    pool = list(meta.index[meta.null_ref_pool])
    fit, test = list(meta.index[meta.null_ref_role == "fit"]), list(meta.index[meta.null_ref_role == "test"])
    meta.to_csv(os.path.join(out, "phenotypes.tsv"), sep="\t", index_label="pheno")
    say(f"alpha = {alpha}")
    say(meta.groupby("kind").size().to_string())
    say(f"binary null reference pool (case fraction >= cutoff): {', '.join(pool) or 'none'}")
    for m in methods:
        if nulls[m] is not None:
            nulls[m].to_parquet(os.path.join(out, f"gene_table_null_{m}.parquet"), index=False)
        cand[m].to_parquet(os.path.join(out, f"candidates_{m}.parquet"), index=False)

    # 0. U2 check: weighted-run pLoF (Beta_1_25) vs production pLoF burden
    say("\n## 0. Weighted-run pLoF CCT vs production pLoF burden CCT (sampled gene x pheno, real + null)")
    rows = []
    if col.check:
        chk = pd.concat(col.check, ignore_index=True)
        la, lb = -np.log10(chk.p_plof_cct), -np.log10(chk.PROD_plof)
        r = dict(comparison="pLoF CCT vs production pLoF burden", n=len(chk),
                 spearman=stats.spearmanr(chk.p_plof_cct, chk.PROD_plof)[0],
                 pearson_log10=np.corrcoef(la, lb)[0, 1], med_abs_dlog10=np.median(np.abs(la - lb)))
        rows.append(r)
        say(f"n={r['n']:,} spearman={r['spearman']:.4f} pearson(-log10)={r['pearson_log10']:.4f} "
            f"median|dlog10|={r['med_abs_dlog10']:.3f}")
    else:
        say("no production run given / no overlapping pLoF tests")
    pd.DataFrame(rows).to_csv(os.path.join(out, "z_data_checks.tsv"), sep="\t", index=False)

    # A. Null phenotypes: hits at alpha
    say("\n## A. Null phenotypes")
    for m in methods:
        n = nulls[m]
        if n is None:
            say(f"{m}: no null phenotypes")
            continue
        say(f"{m}: null hits at alpha  combined {(n.p_comb_cct < alpha).sum()}  missense {(n.p_mis_cct < alpha).sum()}"
            f"  pLoF {(n.p_plof_cct < alpha).sum()}   (null gene x pheno: {len(n):,})")

    # B. Calibration of the directional pLoF p on nulls
    say("\n## B. Directional pLoF p on null phenotypes (expected 0.05 / 0.01 / 0.001)")
    say("method | set | n | <0.05 | <0.01 | <0.001")
    rows, refs = [], {}
    for m in methods:
        refs[m] = None
        if nulls[m] is None:
            continue
        n = nulls[m][nulls[m].p_plof_dir.notna()]
        nc, nb = n.kind == "null_continuous", n.pheno.isin(pool)
        risk, prot = n.s_mis == 1, n.s_mis == -1
        sets = [("continuous, all", nc, "raw"), ("continuous, missense +", nc & risk, "raw"),
                ("continuous, missense -", nc & prot, "raw"),
                ("continuous, missense p<0.01", nc & (n.p_mis_cct < 0.01), "raw"),
                ("binary, risk (missense +)", nb & risk, "raw"), ("binary, protective (missense -)", nb & prot, "raw")]
        held = None
        if fit and test:
            held = emp_p(n.p_plof_dir.values, n.s_mis.values, null_ref(nulls[m], fit))
            tst = n.pheno.isin(test)
            sets += [("binary, risk, empirical p, held-out", tst & risk, "held"),
                     ("binary, protective, empirical p, held-out", tst & prot, "held"),
                     ("binary, risk, empirical p, held-out, missense p<0.01", tst & risk & (n.p_mis_cct < 0.01), "held")]
        for name, msk, src in sets:
            x = n.loc[msk, "p_plof_dir"] if src == "raw" else pd.Series(held[msk.values])
            if not len(x):
                continue
            r = dict(method=m, set=name, n=len(x), lt05=(x < .05).mean(), lt01=(x < .01).mean(),
                     lt001=(x < .001).mean())
            rows.append(r)
            say(f"{m} | {name} | {r['n']:,} | {r['lt05']:.4f} | {r['lt01']:.4f} | {r['lt001']:.5f}")
        say(f"{m}: null genes with ambiguous missense direction: {(nulls[m].s_mis == 0).mean():.3f}")
        if pool:
            refs[m] = null_ref(nulls[m], pool)
            say(f"{m}: binary empirical reference size  risk {len(refs[m][1]):,}  protective {len(refs[m][-1]):,}  "
                f"-> p floor {1 / (len(refs[m][1]) + 1):.1e}")
        else:
            say(f"{m}: WARNING no binary null phenotypes in the reference pool; binary attribution uses the raw "
                "directional p (risk tail inflated, see README)")
    pd.DataFrame(rows).to_csv(os.path.join(out, "b_calibration_null.tsv"), sep="\t", index=False)

    # C. Discovery and evidence pattern
    say("\n## C. Discovery hits (combined CCT < alpha, real phenotypes) by evidence pattern")
    ev = {m: evidence(cand[m], alpha) for m in methods}
    t = pd.concat([h.assign(method=m) for m, h in ev.items()])
    if len(t):
        say(pd.crosstab([t.method, t.kind], t.evidence, margins=True).to_string())
    t.groupby(["method", "kind", "evidence"]).size().rename("n").reset_index() \
     .to_csv(os.path.join(out, "c_hits_by_evidence.tsv"), sep="\t", index=False)
    H = ev["weighted"].copy()
    H["novelty"] = novelty(H, alpha, has_prod)
    if has_prod:
        say("\nweighted hits: novelty vs production")
        say(pd.crosstab(H.novelty, H.kind, margins=True).to_string())
        hp = ev["production"].set_index(["pheno", "gene_id"]).index
        hw = H.set_index(["pheno", "gene_id"]).index
        say(f"overlap with production combined hits: weighted {len(hw)}, production {len(hp)}, "
            f"both {len(hw.intersection(hp))}, production-only {len(hp.difference(hw))}")
        po = ev["production"][~hp.isin(hw)]
        say("production-only hits: " + ", ".join(
            f"{r.symbol}/{r.pheno} (w={r.other_p_comb_cct:.1e})" for r in po.sort_values("p_comb_cct").itertuples()))

    say("\nmatched-FPR sensitivity: threshold = FPR quantile of the null combined CCTs, per trait type")
    F = pd.DataFrame(col.fpr, columns=["trait_type", "fpr", "method", "pheno", "real_hits"])
    rows = []
    for (m, tt, f), (th, n_null) in thr.items():
        k = int(F[(F.method == m) & (F.trait_type == tt) & (F.fpr == f)].real_hits.sum())
        rows.append(dict(trait_type=tt, fpr=f, method=m, n_null=n_null, threshold=th, real_hits=k))
    R = pd.DataFrame(rows, columns=["trait_type", "fpr", "method", "n_null", "threshold", "real_hits"])
    say(R.to_string(index=False))
    R.to_csv(os.path.join(out, "c_matched_fpr.tsv"), sep="\t", index=False)

    # D. Attribution
    say("\n## D. Attribution set (missense CCT < alpha, real phenotypes)")
    att = {}
    for m in methods:
        A, n_test = attribution(cand[m], refs[m], alpha, has_prod)
        att[m] = A
        say(f"\n### {m}: {len(A)} gene x pheno; {n_test} tested; strict = 0.05/{n_test} = {0.05 / max(n_test, 1):.2e}")
        if not len(A):
            continue
        say(f"in combined hits: {(A.p_comb_cct < alpha).sum()} / {len(A)}")
        say(pd.crosstab(A.kind, A.attr, margins=True).to_string())
        if has_prod:
            say("\nnovelty x attr:")
            say(pd.crosstab(A.novelty, A.attr, margins=True).to_string())
        for cls in ["direction_ambiguous", "pLoF_MAC_lt5", "no_pLoF_test"]:
            z = A[A.attr == cls]
            if len(z):
                say(f"{cls}: {', '.join(z.symbol.astype(str) + '/' + z.pheno)}")
        cols = ["pheno", "kind", "gene", "gene_id", "symbol", "attr", "novelty", "p_mis_cct", "p_comb_cct",
                "p_plof_cct", "PROD_plof", "s_mis", "n_sign_tests", "p_plof_dir", "p_attr", "p_attr_source",
                "n_plof_tests", "MAC_plof", "E_plof"]
        A[cols].sort_values(["attr", "p_attr"]).to_csv(os.path.join(out, f"d_attribution_{m}.tsv"), sep="\t",
                                                       index=False)
    Aw = att["weighted"]
    if has_prod and len(Aw):
        for cls in ["LoF_like_strict", "LoF_like_suggestive"]:
            z = Aw[(Aw.attr == cls) & (Aw.novelty == "new")].sort_values("p_attr")
            say(f"\nweighted, new (not production pLoF/combined/missense) {cls}: "
                + ", ".join(f"{r.symbol}/{r.pheno} (p={r.p_attr:.1e})" for r in z.itertuples()))
        z = Aw[Aw.attr.str.startswith("LoF_like") & (Aw.novelty != "prod_pLoF")].sort_values("p_attr")
        say(f"\nweighted LoF-like that production pLoF burden misses ({len(z)}): "
            + ", ".join(f"{r.symbol}/{r.pheno}[{r.attr.replace('LoF_like_', '')}]" for r in z.itertuples()))
    if len(Aw):
        pb = Aw[Aw.binary & (Aw.s_mis == -1) & Aw.attr.str.startswith("LoF_like")]
        say("binary protective LoF-like calls (unconfirmed, O2): "
            + (", ".join(pb.symbol.astype(str) + "/" + pb.pheno) if len(pb) else "none"))
    if has_prod:
        say("\n## D2. Weighted vs production attribution (gene x pheno)")
        kw = Aw.set_index(["pheno", "gene_id"]).attr
        kp = att["production"].set_index(["pheno", "gene_id"]).attr
        both = pd.concat([kw.rename("weighted"), kp.rename("production")], axis=1).fillna("not_in_set")
        if len(both):
            say(pd.crosstab(both.weighted, both.production, margins=True).to_string())

    # E. Concordance
    say("\n## E. Concordance by missense bin: share of genes with p_plof_dir < 0.5 (pLoF leans the missense way)")
    C = pd.concat([c for m in methods for c in col.conc[m]], ignore_index=True)
    E = C.groupby(["method", "tt", "real", "mbin"], observed=True)[["size", "sum"]].sum().reset_index()
    E["mean"] = E["sum"] / E["size"]
    E = E[["tt", "real", "mbin", "size", "mean", "method"]]
    say(E.pivot_table(index=["method", "tt", "mbin"], columns="real", values=["size", "mean"], observed=True)
        .round(3).to_string())
    E.to_csv(os.path.join(out, "f_concordance_by_missense_bin.tsv"), sep="\t", index=False)

    # G. Joint-only hits
    say("\n## G. Joint-only weighted hits (descriptive, NOT selection-adjusted)")
    J = H[H.evidence == "joint_only"].copy()
    J["weighting_needed"] = ~(J.PROD_plof_mis < alpha) if has_prod else np.nan
    J["parts_concordant_nominal"] = (J.s_mis != 0) & (J.p_plof_dir < 0.05) & (J.p_mis_cct < 0.05)
    say(f"n={len(J)}  parts concordant (nominal): {J.parts_concordant_nominal.sum()}"
        + (f"  weighting needed: {J.weighting_needed.sum()}" if has_prod else ""))
    cols = ["pheno", "kind", "gene", "gene_id", "symbol", "p_comb_cct", "p_plof_cct", "p_mis_cct", "s_mis",
            "p_plof_dir", "parts_concordant_nominal", "PROD_plof_mis", "weighting_needed", "novelty"]
    J = J[cols].sort_values("p_comb_cct")
    if len(J):
        say(J.to_string(index=False, float_format=lambda v: f"{v:.2e}"))
    J.to_csv(os.path.join(out, "g_joint_only_weighted.tsv"), sep="\t", index=False)

    Hx = H.merge(Aw[["pheno", "gene_id", "attr", "p_attr"]], on=["pheno", "gene_id"], how="left") if len(Aw) \
        else H.assign(attr=np.nan, p_attr=np.nan)
    Hx["attr"] = Hx.attr.fillna("not_in_attribution_set")
    Hx[["pheno", "kind", "gene", "gene_id", "symbol", "evidence", "novelty", "attr", "p_attr", "p_comb_cct",
        "p_plof_cct", "p_mis_cct", "s_mis", "p_plof_dir", "BETA_top_comb", "SE_top_comb", "Z_top_comb",
        "max_MAF_top_comb",
        "PROD_plof", "PROD_mis", "PROD_plof_mis"]].sort_values("p_comb_cct") \
        .to_csv(os.path.join(out, "h_hits_weighted.tsv"), sep="\t", index=False)




# ------------------------------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------------------------------
SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
C_W, C_P = "#2a78d6", "#eb6834"                       # weighted, production
EVID = {"both": ("Both parts < α", "#2a78d6"), "pLoF_only": ("pLoF only", "#eb6834"),
        "missense_only": ("Missense only", "#1baf7a"), "joint_only": ("Joint only", "#eda100")}
CALL = {"LoF_like_strict": ("LoF-like, strict", "#1c5cab"),
        "LoF_like_suggestive": ("LoF-like, suggestive", "#86b6ef"),
        "not_LoF_like": ("Not LoF-like", AXIS)}
RC = {
    "figure.dpi": 150, "savefig.dpi": 200, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "font.family": "sans-serif", "font.size": 9, "text.color": INK,
    "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
    "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "legend.frameon": False, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlelocation": "left",
}


def _nl10(p, cap=None):
    v = -np.log10(np.clip(np.asarray(p, float), 1e-300, 1))
    return v if cap is None else np.minimum(v, cap)


def _place_labels(ax, xs, ys, texts, fontsize=7, step=1.15, max_tries=40):
    """Non-overlapping direct labels: each starts right of its point and moves up/down in steps of its own
    height until it clears the labels already placed; a thin leader joins a moved label to its point."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    placed = []
    for i in np.argsort(ys):
        x, y, t = xs[i], ys[i], texts[i]
        off = 0
        for k in range(max_tries):
            dy = (k + 1) // 2 * (1 if k % 2 else -1) if k else 0
            txt = ax.annotate(t, (x, y), xytext=(5, 2), textcoords="offset points", fontsize=fontsize, color=INK)
            bb = txt.get_window_extent(r)
            txt.remove()
            off = dy * bb.height * step
            if not any(bb.translated(0, off).overlaps(o) for o in placed):
                break
        txt = ax.annotate(t, (x, y), xytext=(5, 2 + off * 72 / fig.dpi), textcoords="offset points",
                          fontsize=fontsize, color=INK,
                          arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5, shrinkA=0, shrinkB=2)
                          if abs(off) > 1 else None)
        placed.append(txt.get_window_extent(r))


def _stacked_hbar(ax, rows, order, style, title):
    """rows: {label: {category: n}}; one horizontal stacked bar per row, segment counts labelled."""
    ax_max = max([sum(c.values()) for c in rows.values()] + [1])
    for i, (lab, cnt) in enumerate(rows.items()):
        left = 0
        for c in order:
            n = cnt.get(c, 0)
            if n == 0:
                continue
            ax.barh(i, n, left=left, color=style[c][1], height=0.6, edgecolor=SURFACE, linewidth=2)
            if n >= 0.04 * ax_max:
                ax.text(left + n / 2, i, str(n), ha="center", va="center", fontsize=8, color=INK)
            left += n
        ax.text(left + 0.01 * ax_max, i, f"{left}", va="center", fontsize=8.5, color=INK, fontweight="bold")
    ax.set_xlim(0, ax_max * 1.12)
    ax.set_yticks(range(len(rows)), list(rows))
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_title(title)


def make_figures(out: str, alpha: float = ALPHA):
    """fig1_discovery, fig4_attribution, fig6_concordance from the tables discovery report wrote in `out`."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(RC)
    fig_dir = os.path.join(out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    def save(fig, name):
        for ext in ("png", "svg"):
            fig.savefig(os.path.join(fig_dir, f"{name}.{ext}"), bbox_inches="tight")
        plt.close(fig)
        print("wrote", os.path.join(fig_dir, name))

    def rd(name):
        p = os.path.join(out, name)
        return pd.read_csv(p, sep="\t", dtype={"pheno": str}) if os.path.exists(p) else None

    n_kind = rd("phenotypes.tsv").kind.value_counts()
    ev, fpr, conc, Aw = (rd("c_hits_by_evidence.tsv"), rd("c_matched_fpr.tsv"),
                         rd("f_concordance_by_missense_bin.tsv"), rd("d_attribution_weighted.tsv"))
    methods = [m for m in ("weighted", "production") if (ev.method == m).any() or (fpr.method == m).any()]

    # fig1: hits by evidence pattern + matched-FPR sensitivity
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.4), gridspec_kw=dict(width_ratios=[1.3, 1, 1.1], wspace=0.35))
    order = list(EVID)
    for ax, kind, ttl in [(axs[0], "quantitative", f"a  Quantitative traits ({n_kind.get('quantitative', 0)})"),
                          (axs[1], "binary", f"b  Binary traits ({n_kind.get('binary', 0)})")]:
        rows = {m.capitalize(): ev[(ev.method == m) & (ev.kind == kind)].set_index("evidence").n.to_dict()
                for m in methods}
        _stacked_hbar(ax, rows, order, EVID, ttl)
        ax.set_xlabel("Hits (combined CCT < α)")
    axs[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=EVID[c][1]) for c in order],
                  labels=[EVID[c][0] for c in order], ncol=4, loc="upper left", bbox_to_anchor=(0, -0.32),
                  fontsize=8)
    ax = axs[2]
    combos = [(tt, f) for tt in ["quantitative", "binary"] for f in FPRS]
    x = np.arange(len(combos))
    width = 0.76 / max(len(methods), 1)
    for j, m in enumerate(methods):
        v = []
        for tt, f in combos:
            s = fpr[(fpr.method == m) & (fpr.trait_type == tt) & np.isclose(fpr.fpr, f)].real_hits
            v.append(int(s.iloc[0]) if len(s) else 0)
        xs = x + (j - (len(methods) - 1) / 2) * width
        ax.bar(xs, v, width=width * 0.95, color=[C_W, C_P][j], label=m.capitalize())
        for xi, vi in zip(xs, v):
            ax.text(xi, vi, str(vi), ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(x, [f"{'Quant.' if tt == 'quantitative' else 'Binary'}\nFPR {f:.0e}".replace("e-0", "e-")
                      for tt, f in combos])
    ax.set_ylabel("Real-trait hits")
    ax.set_title("c  Sensitivity: thresholds matched on null traits")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="x", visible=False)
    fig.suptitle("Gene-level discoveries by evidence pattern, at α and at matched FPR",
                 x=0.01, ha="left", fontsize=11, fontweight="bold", y=1.04)
    save(fig, "fig1_discovery")

    # fig4: attribution set, missense CCT vs directional pLoF p, coloured by call
    if Aw is not None and len(Aw):
        A = Aw[Aw.attr.isin(list(CALL))].copy()
        n_test = max(len(A), 1)
        CAPX, CAPY = 30, 12
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        for c in ["not_LoF_like", "LoF_like_suggestive", "LoF_like_strict"]:
            a = A[A.attr == c]
            for kind, mk in [("quantitative", "o"), ("binary", "s")]:
                b = a[a.kind == kind]
                ax.scatter(_nl10(b.p_mis_cct, CAPX), _nl10(b.p_attr, CAPY), s=30, c=CALL[c][1], marker=mk,
                           edgecolors=SURFACE, linewidths=0.8, zorder=3,
                           label=CALL[c][0] if kind == "quantitative" else None)
        ax.axhline(-np.log10(0.05 / n_test), color=INK2, lw=1, ls="--")
        ax.axhline(-np.log10(0.05), color=MUTED, lw=1, ls=":")
        ax.text(CAPX - 0.3, -np.log10(0.05 / n_test) + 0.15, f"strict: 0.05/{n_test}", ha="right", fontsize=7.5,
                color=INK2)
        ax.text(CAPX - 0.3, -np.log10(0.05) + 0.15, "suggestive: 0.05", ha="right", fontsize=7.5, color=INK2)
        ax.set_xlim(-np.log10(alpha) - 0.6, CAPX)
        ax.set_ylim(0, CAPY + 0.5)
        ax.set_xlabel("Missense CCT, −log10 p (selection; capped at 30)")
        ax.set_ylabel("Directional pLoF p, −log10 (capped at 12)")
        has_prod = (Aw.novelty.notna() & (Aw.novelty != "NA")).any()      # "NA" is read back as missing
        sub = (f"{n_test} tested gene × trait pairs with missense CCT < α; squares binary (empirical p)\n"
               + ("labels = LoF-like calls not significant in any production test" if has_prod else ""))
        ax.set_title("Attribution: does pLoF agree in the missense direction?", pad=8 + 11 * (sub.count("\n") + 1))
        ax.text(0, 1.012, sub, transform=ax.transAxes, fontsize=7.5, color=INK2, va="bottom")
        new = A[A.attr.str.startswith("LoF_like") & (A.novelty == "new")]
        if len(new):
            sym = new.symbol.astype(str)
            _place_labels(ax, np.minimum(_nl10(new.p_mis_cct.values), CAPX),
                          np.minimum(_nl10(new.p_attr.values), CAPY),
                          list(sym + np.where(sym.duplicated(keep=False), " / " + new.pheno, "")))
        ax.legend(loc="upper right", bbox_to_anchor=(1, 0.93), fontsize=8)
        save(fig, "fig4_attribution")

    # fig6: concordance by missense bin, real vs null
    XL = ["< α", "α–1e-4", "1e-4–0.01", "0.01–0.1", "> 0.1"]
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.8), sharey=True)
    cw = conc[conc.method == "weighted"]
    for ax, tt, ttl in [(axs[0], "quantitative", "a  Quantitative traits"), (axs[1], "binary", "b  Binary traits")]:
        ns = XL
        for real, col, lab, ls in [(True, C_W, "Real traits", "-"), (False, MUTED, "Null traits", "--")]:
            g = cw[(cw.tt == tt) & (cw.real == real)].set_index("mbin").reindex(MBL)
            xs = np.arange(len(MBL))
            ok = g["size"].notna().values
            ax.plot(xs[ok], g["mean"].values[ok], color=col, lw=2, ls=ls, marker="o", ms=6, label=lab,
                    markeredgecolor=SURFACE)
            if real:
                ns = [f"{xl}\nn={int(sz):,}" if pd.notna(sz) else xl for xl, sz in zip(XL, g["size"])]
        ax.axhline(0.5, color=AXIS, lw=1)
        ax.set_xticks(range(len(MBL)), ns, fontsize=7.5)
        ax.set_xlabel("Missense CCT bin (n = real-trait genes)")
        ax.set_title(ttl)
        ax.set_ylim(0.25, 1)
        ax.legend(loc="upper right", fontsize=8)
    axs[0].set_ylabel("Share with pLoF leaning the\nmissense direction (p_dir < 0.5)")
    fig.suptitle("pLoF agrees with missense more often the stronger the missense signal (read binary against its "
                 "null)", x=0.01, ha="left", fontsize=11, fontweight="bold", y=1.03)
    save(fig, "fig6_concordance")


# ------------------------------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weighted-root",
                    help="ht_results ancestry root of the weighted run (phenotype_<X>/ subdirectories), "
                         "or a local directory of <X>.tsv.gz copies")
    ap.add_argument("--production-root", help="same, for the unweighted production run (novelty labels and "
                                              "side-by-side comparison); omitted = weighted-only report")
    ap.add_argument("--manifest", help="phenotype CSV: phenoname, trait_type, n_cases, n_controls "
                                       "[, ancestry, category, is_null]. Default: metadata columns of the tables")
    ap.add_argument("--ancestry", help="filter the manifest's ancestry column")
    ap.add_argument("--phenotypes", help="comma-separated subset of phenotypes")
    ap.add_argument("--out", required=True, help="local output directory")
    ap.add_argument("--alpha", type=float, default=ALPHA)
    ap.add_argument("--keep-p", type=float, default=None,
                    help="also keep real-phenotype genes with combined or missense CCT below this in "
                         "candidates_*.parquet (default: alpha)")
    ap.add_argument("--check-n", type=int, default=5000,
                    help="gene x pheno pairs sampled per real phenotype for the pLoF-vs-production check")
    ap.add_argument("--plof-min-mac", type=float, default=PLOF_MIN_MAC,
                    help="pLoF tests with fewer minor alleles are treated as not run (class pLoF_MAC_lt5); "
                         "applied to both runs, real and null phenotypes. 0 keeps every test")
    ap.add_argument("--min-neff", type=float, default=MIN_NEFF,
                    help="FlexRV weight sets with nEff_Burden below this are dropped and the per-group CCTs "
                         "recomputed from the per-weight-set exports (gene_results_<max_MAF>.txt.bgz); applied "
                         "to both runs where those exports exist. 0 uses the per-group tables as written")
    ap.add_argument("--skip-figures", action="store_true", help="tables only")
    ap.add_argument("--figures-only", action="store_true", help="redraw the figures from an existing --out")
    args = ap.parse_args(argv)

    if not args.figures_only:
        if not args.weighted_root:
            ap.error("--weighted-root is required (unless --figures-only)")
        os.makedirs(args.out, exist_ok=True)
        say = Log(os.path.join(args.out, "report_log.txt"))
        loaders = {"weighted": cct_loader(args.weighted_root, args.plof_min_mac, args.min_neff)}
        if args.production_root:
            loaders["production"] = cct_loader(args.production_root, args.plof_min_mac, args.min_neff)
        meta = build_meta(args, loaders)
        say(f"weighted: {args.weighted_root}\nproduction: {args.production_root or '-'}")
        say(f"pLoF tests with MAC < {args.plof_min_mac:g} are treated as not run (--plof-min-mac)")
        say(f"FlexRV weight sets with nEff_Burden < {args.min_neff:g} are dropped (--min-neff)" if args.min_neff > 0
            else "no nEff_Burden floor (--min-neff 0): per-group tables used as written")
        col, thr = collect(loaders, meta, args.alpha, args.keep_p, args.check_n, log=say)
        for m, f in loaders.items():
            n = f.neff
            if n["phenotypes"]:
                say(f"{m}: nEff floor applied to {n['phenotypes']} phenotypes: {n['dropped']:,} of {n['tests']:,} "
                    f"weight-set tests dropped ({n['dropped'] / max(n['tests'], 1):.1%}); {n['groups_lost']:,} of "
                    f"{n['groups']:,} gene x group left with no test")
            if n["meta"]:
                say(f"{m}: {n['meta']} meta-analysis tables used as written (nEff floor applied per ancestry "
                    f"by --create-flexrv-meta; see their min_nEff column)")
            if n["production_meta"]:
                say(f"{m}: {n['production_meta']} production meta-analysis tables (Beta(1,25) burden per mask and "
                    f"max MAF tier, no weight sets: the nEff floor does not apply)")
            if args.min_neff > 0 and n["unfiltered"]:
                u = sorted(set(n["unfiltered"]))
                say(f"{m}: NO nEff floor for {len(u)} phenotypes (no per-weight-set exports next to the "
                    f"per-group table): {', '.join(u[:20])}" + (" ..." if len(u) > 20 else ""))
        report(col, thr, meta, args.out, args.alpha, say)
    if not args.skip_figures:
        make_figures(args.out, args.alpha)


if __name__ == "__main__":
    main()
