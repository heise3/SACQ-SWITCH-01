"""Original fixed module definitions and pure public-input parsers; no analysis runner."""
import gzip,hashlib,json,math,os,re,sys
from collections import defaultdict
import numpy as np,pandas as pd
MODULES = {
    "M_IFN": [
        "IFI27", "IFI44", "IFI44L", "ISG15", "IFIT1", "IFIT3",
        "OAS1", "OAS2", "OAS3", "OASL", "MX1", "RSAD2",
        "EPSTI1", "HERC5", "LY6E", "STAT1", "SIGLEC1", "USP18"
    ],
    "M_LDG": [
        "DEFA4", "DEFA3", "CAMP", "LTF", "LCN2", "MMP8",
        "MPO", "ELANE", "CTSG", "CEACAM8", "OLFM4", "BPI", "AZU1"
    ],
    "M_PB": [
        "MZB1", "JCHAIN", "TXNDC5", "XBP1", "PRDM1", "SDC1",
        "TNFRSF17", "POU2AF1", "CD38", "IGJ", "RRM2", "TYMS"
    ],
    "NC_HK": [
        "ACTB", "B2M", "GAPDH", "PGK1", "PPIA", "RPL13A",
        "RPLP0", "SDHA", "TBP", "YWHAZ", "HPRT1"
    ]
}

# Canonical HGNC synonym harmonization (IGJ is official alias of JCHAIN; DEFA1B/DEFA3)
SYNONYMS = {
    "IGJ": "JCHAIN",
    "JCHAIN": "JCHAIN",
    "TNFRSF17": "TNFRSF17",
    "BCMA": "TNFRSF17",
    "SDC1": "SDC1",
    "CD138": "SDC1",
    "CEACAM8": "CEACAM8",
    "CD66B": "CEACAM8",
    "DEFA3": "DEFA3",
    "DEFA1 /// DEFA1B /// DEFA3": "DEFA3",
}

CANONICAL_MODULES = {}
ALL_TARGET_GENES = set()
for m, genes in MODULES.items():
    canon = []
    for g in genes:
        cg = SYNONYMS.get(g, g)
        if cg not in canon:
            canon.append(cg)
        ALL_TARGET_GENES.add(cg)
        ALL_TARGET_GENES.add(g)
    CANONICAL_MODULES[m] = canon


def normalize_gene_symbol(raw_sym: str):
    if not raw_sym or raw_sym in ("---", "NA", "nan", "NULL"):
        return None
    parts = [p.strip() for p in re.split(r"///|;|,", raw_sym) if p.strip()]
    mapped = []
    for p in parts:
        p_up = p.upper()
        cg = SYNONYMS.get(p_up, p_up)
        if cg in ALL_TARGET_GENES:
            mapped.append(cg)
    uniq = sorted(set(mapped))
    if len(uniq) == 1:
        return SYNONYMS.get(uniq[0], uniq[0])
    return None


def cohens_d(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled = math.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2))
    return float((np.mean(x) - np.mean(y)) / pooled) if pooled > 0 else np.nan


def compute_module_scores(gene_df: pd.DataFrame, hc_mask: np.ndarray):
    """
    gene_df: samples x canonical_genes (already log2 scale).
    Standardize each gene to Healthy Control mean & SD (if >=5 HCs), else whole-cohort SD.
    """
    hc_idx = np.where(hc_mask)[0]
    if len(hc_idx) >= 5:
        mu = gene_df.iloc[hc_idx].mean(axis=0)
        sd = gene_df.iloc[hc_idx].std(axis=0, ddof=1).replace(0, np.nan)
        # Fallback if a gene has ~0 variance in HC
        sd_all = gene_df.std(axis=0, ddof=1).replace(0, np.nan)
        sd = sd.fillna(sd_all).where(sd > 1e-6, sd_all)
    else:
        mu = gene_df.mean(axis=0)
        sd = gene_df.std(axis=0, ddof=1).replace(0, np.nan)
    z_df = (gene_df - mu) / sd
    scores = pd.DataFrame(index=gene_df.index)
    coverage = {}
    for m, g_list in CANONICAL_MODULES.items():
        avail = [g for g in g_list if g in z_df.columns and z_df[g].notna().sum() > 0]
        coverage[m] = avail
        scores[m] = z_df[avail].mean(axis=1)
    scores["ESI"] = 0.5 * scores["M_LDG"] + 0.5 * scores["M_PB"]
    scores["Uncoupling_Delta"] = scores["ESI"] - 0.5 * scores["M_IFN"]
    return scores, coverage


# ==============================================================================
# Platform Parsers
# ==============================================================================
def load_gpl10558_map():
    path = os.path.join(CACHE_DIR, "GPL10558.annot.gz")
    probe2gene = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        in_table = False
        headers = None
        sym_idx = None
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                cols = line.split("\t")
                if headers is None:
                    headers = cols
                    sym_idx = headers.index("Gene symbol")
                    continue
                if len(cols) > sym_idx:
                    g = normalize_gene_symbol(cols[sym_idx])
                    if g:
                        probe2gene[cols[0].strip()] = g
    print(f"[GPL10558] Mapped {len(probe2gene)} target probes across {len(set(probe2gene.values()))} target genes")
    return probe2gene


def load_gpl14663_map():
    raw_map = json.load(open(os.path.join(CACHE_DIR, "GPL14663_map.json")))
    probe2gene = {}
    for pid, sym in raw_map.items():
        g = normalize_gene_symbol(str(sym))
        if g:
            probe2gene[pid.strip()] = g
            probe2gene[pid.strip() + "_at"] = g
    print(f"[GPL14663] Mapped {len(probe2gene)} target probes across {len(set(probe2gene.values()))} target genes")
    return probe2gene


def load_gpl17586_map():
    raw_map = json.load(open(os.path.join(CACHE_DIR, "GPL17586_map.json")))
    probe2gene = {}
    for pid, sym in raw_map.items():
        g = normalize_gene_symbol(str(sym))
        if g:
            probe2gene[pid.strip()] = g
    print(f"[GPL17586] Mapped {len(probe2gene)} target probesets across {len(set(probe2gene.values()))} target genes")
    return probe2gene


def load_gpl13158_map():
    raw_map = json.load(open(os.path.join(CACHE_DIR, "GPL13158_map.json")))
    probe2gene = {}
    for pid, sym in raw_map.items():
        g = normalize_gene_symbol(str(sym))
        if g:
            probe2gene[pid.strip()] = g
    print(f"[GPL13158] Mapped {len(probe2gene)} target probes across {len(set(probe2gene.values()))} target genes")
    return probe2gene


def parse_series_matrix(gz_path: str, probe2gene: dict):
    """
    Streams a GEO series_matrix.txt.gz file, extracting all sample metadata
    and ONLY expression rows belonging to probe2gene. Automatically converts linear -> log2 if needed.
    """
    meta_rows = defaultdict(list)
    sample_ids = []
    probe_vals = {}
    in_table = False
    raw_sample_vals = []

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not in_table:
                if line.startswith("!Sample_geo_accession"):
                    sample_ids = [x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_title"):
                    meta_rows["_title"] = [x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_source_name_ch1"):
                    meta_rows["_source"] = [x.strip('"') for x in line.split("\t")[1:]]
                elif line.startswith("!Sample_description"):
                    if "_desc" not in meta_rows:
                        meta_rows["_desc"] = [x.strip('"') for x in line.split("\t")[1:]]
                    else:
                        vals = [x.strip('"') for x in line.split("\t")[1:]]
                        meta_rows["_desc"] = [a + " | " + b for a, b in zip(meta_rows["_desc"], vals)]
                elif line.startswith("!Sample_supplementary_file"):
                    vals = [x.strip('"') for x in line.split("\t")[1:]]
                    if "_supp" not in meta_rows:
                        meta_rows["_supp"] = vals
                elif line.startswith("!Sample_characteristics_ch1"):
                    vals = [x.strip('"') for x in line.split("\t")[1:]]
                    meta_rows["_chars_list"].append(vals)
                elif line.startswith("!series_matrix_table_begin"):
                    in_table = True
            else:
                if line.startswith("!series_matrix_table_end"):
                    break
                if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                    hdr = [x.strip('"') for x in line.split("\t")[1:]]
                    if hdr:
                        sample_ids = hdr
                    continue
                if probe2gene is None:
                    continue
                tab_pos = line.find("\t")
                if tab_pos == -1:
                    continue
                pid = line[:tab_pos].strip('"')
                if pid in probe2gene:
                    arr = np.array(
                        [float(x) if x not in ("", "null", "NULL", "NA", "NaN") else np.nan
                         for x in line[tab_pos + 1:].split("\t")],
                        dtype=float
                    )
                    probe_vals[pid] = arr
                    if len(raw_sample_vals) < 5:
                        raw_sample_vals.append(arr[np.isfinite(arr)][:50])

    # Build metadata DataFrame
    meta = pd.DataFrame(index=sample_ids)
    for k in ("_title", "_source", "_desc", "_supp"):
        if k in meta_rows and len(meta_rows[k]) == len(sample_ids):
            meta[k] = meta_rows[k]
    char_dicts = [dict() for _ in range(len(sample_ids))]
    for row_vals in meta_rows.get("_chars_list", []):
        for idx, v in enumerate(row_vals):
            if idx < len(sample_ids) and ":" in v:
                ck, cv = v.split(":", 1)
                char_dicts[idx][ck.strip()] = cv.strip()
    char_df = pd.DataFrame(char_dicts, index=sample_ids)
    meta = pd.concat([meta, char_df], axis=1)

    if probe2gene is None:
        return meta, None, {}

    expr_probe = pd.DataFrame(probe_vals, index=sample_ids)
    finite_vals = expr_probe.values[np.isfinite(expr_probe.values)]
    q99 = float(np.percentile(finite_vals, 99)) if len(finite_vals) else np.nan
    q01 = float(np.percentile(finite_vals, 1)) if len(finite_vals) else np.nan
    scale_note = "already_log2"
    if q99 > 50.0:
        # Linear intensity -> shift if min < 0 then log2(x + 1)
        min_v = float(np.nanmin(finite_vals))
        shift = abs(min_v) + 1.0 if min_v < 0 else 1.0
        expr_probe = np.log2(expr_probe + shift)
        scale_note = f"log2_transformed(shift={shift:.2f}, raw_q01={q01:.2f}, raw_q99={q99:.2f})"
    else:
        scale_note = f"already_log2(q01={q01:.2f}, q99={q99:.2f})"

    # Collapse probes to canonical gene by mean
    gene_cols = defaultdict(list)
    for pid in expr_probe.columns:
        gene_cols[probe2gene[pid]].append(pid)
    gene_df = pd.DataFrame(
        {g: expr_probe[pids].mean(axis=1) for g, pids in gene_cols.items()},
        index=sample_ids
    )
    return meta, gene_df, {"scale": scale_note, "n_probes": len(probe_vals), "n_genes": gene_df.shape[1]}


def safe_float(s):
    try:
        return float(str(s).strip())
    except Exception:
        return np.nan


