#!/usr/bin/env python3
"""
09_reannotation_flags.py — Flag aquaporins needing reannotation + generate IGV regions

Identifies aquaporins that are partial or mis-annotated based on:
- fuente_seq == GFF3_FALLBACK (both GFF3 and EXONERATE annotations are bad)
- fuente_seq == MAKER_GFF3 (MAKER-annotated, not found by tblastn+exonerate)
- TMHs != 6 (missing transmembrane domains = partial)
- veredicto == AMBAS_MAL (both annotations bad)

Cross-references with RNA-seq expression data to prioritize candidates.
Generates BED file for IGV inspection.

FIXES applied:
  - Uses absolute paths derived from PROJECT_DIR (not relative)
  - Robust path handling
"""

import pandas as pd
import os
import sys

# ---- Configuration (absolute paths) -----------------------------------------
PROJECT_DIR   = os.environ.get("PROJECT_DIR", os.getcwd())
AQP_TABLE     = os.path.join(PROJECT_DIR, "tabla_Aquaporinas_traduccion.tabular")
GFF3_FILE     = os.path.join(PROJECT_DIR, "refs", "genome.gff3")
TPM_FILE      = os.path.join(PROJECT_DIR, "results", "basal_aquaporins", "basal_aquaporins_tpm.csv")
COUNTS_FILE   = os.path.join(PROJECT_DIR, "counts", "counts_all.csv")
RESULTS_DIR   = os.path.join(PROJECT_DIR, "results", "basal_aquaporins")
FLANKING_BP   = 5000  # Flanking region for IGV BED file

os.makedirs(RESULTS_DIR, exist_ok=True)

print("# =============================================================")
print("# AQUAPORIN REANNOTATION FLAG ANALYSIS")
print("# =============================================================")

# ---- Load aquaporin metadata ------------------------------------------------
aqp = pd.read_csv(AQP_TABLE, sep="\t")
print(f"# Total aquaporins: {len(aqp)}")

# ---- Flag candidates for reannotation ---------------------------------------
def flag_reason(row):
    reasons = []
    if row["fuente_seq"] == "GFF3_FALLBACK":
        reasons.append("GFF3_FALLBACK: both GFF3 and EXONERATE annotations unreliable")
    if row["fuente_seq"] == "MAKER_GFF3":
        reasons.append("MAKER_GFF3: only MAKER annotation, not found by tblastn+exonerate")
    if row.get("veredicto") == "AMBAS_MAL":
        reasons.append("AMBAS_MAL: both annotation sources confirmed bad")
    if pd.notna(row.get("TMHs")) and row["TMHs"] != 6:
        reasons.append(f"TMHs={row['TMHs']}: expected 6 transmembrane domains for complete aquaporin")
    if pd.notna(row.get("longitud_aa")) and row["longitud_aa"] < 200:
        reasons.append(f"Short protein ({row['longitud_aa']} aa): likely partial/truncated")
    return "; ".join(reasons) if reasons else ""

aqp["flag_reasons"] = aqp.apply(flag_reason, axis=1)
aqp["needs_reannotation"] = aqp["flag_reasons"] != ""

candidates = aqp[aqp["needs_reannotation"]].copy()
ok_genes   = aqp[~aqp["needs_reannotation"]].copy()

print(f"# Aquaporins needing reannotation: {len(candidates)}/{len(aqp)}")
print(f"# Aquaporins OK: {len(ok_genes)}/{len(aqp)}")

# ---- Cross-reference with expression data -----------------------------------
if os.path.exists(TPM_FILE):
    tpm = pd.read_csv(TPM_FILE)
    # Get sample columns (not metadata columns)
    meta_cols = ["gene_id", "aqp_family_subfamily", "subfamilia_phylo",
                 "fuente_seq", "needs_reannotation"]
    sample_cols = [c for c in tpm.columns if c not in meta_cols]

    if sample_cols:
        tpm_vals = tpm.set_index("gene_id")[sample_cols]
        # Max TPM across any tissue for each aquaporin
        max_tpm = tpm_vals.max(axis=1)
        mean_tpm = tpm_vals.mean(axis=1)

        candidates = candidates.merge(
            pd.DataFrame({"max_tpm": max_tpm, "mean_tpm": mean_tpm}),
            left_on="gene_id", right_index=True, how="left"
        )
        candidates["expressed_in_rnaseq"] = candidates["max_tpm"].fillna(0) > 0.5
    else:
        candidates["max_tpm"] = 0
        candidates["mean_tpm"] = 0
        candidates["expressed_in_rnaseq"] = False
else:
    print(f"# WARNING: TPM file not found ({TPM_FILE}). Skipping expression cross-reference.")
    candidates["max_tpm"] = None
    candidates["mean_tpm"] = None
    candidates["expressed_in_rnaseq"] = None

# ---- Extract genomic coordinates from GFF3 -----------------------------------
print("# Extracting genomic coordinates from GFF3...")
gene_coords = {}

if os.path.exists(GFF3_FILE):
    with open(GFF3_FILE) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "gene":
                continue
            attrs = parts[8]
            gene_id = None
            for attr in attrs.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr.split("=")[1]
                    break
            if gene_id:
                gene_coords[gene_id] = {
                    "chr": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6]
                }
else:
    print(f"# WARNING: GFF3 file not found ({GFF3_FILE}). Cannot extract coordinates.")

# Add coordinates to candidates
for col in ["chr", "start", "end", "strand"]:
    candidates[col] = candidates["gene_id"].map(
        lambda x: gene_coords.get(x, {}).get(col, "")
    )

# ---- Sort by priority -------------------------------------------------------
# Priority: expressed in RNA-seq > shorter protein > more reasons
candidates = candidates.sort_values(
    by=["expressed_in_rnaseq", "longitud_aa"],
    ascending=[False, True]
)

# ---- Save candidates table ---------------------------------------------------
out_cols = ["gene_id", "aqp_family_subfamily", "subfamilia_phylo",
            "fuente_seq", "veredicto", "TMHs", "longitud_aa",
            "chr", "start", "end", "strand",
            "max_tpm", "mean_tpm", "expressed_in_rnaseq",
            "flag_reasons"]
out_cols = [c for c in out_cols if c in candidates.columns]

candidates_out = candidates[out_cols]
candidates_out.to_csv(
    os.path.join(RESULTS_DIR, "reannotation_candidates.tsv"),
    sep="\t", index=False
)
print(f"# Saved: {RESULTS_DIR}/reannotation_candidates.tsv")

# ---- Generate IGV BED file ---------------------------------------------------
bed_lines = []
for _, row in candidates.iterrows():
    if row["chr"] and row["start"]:
        chrom = row["chr"]
        start = max(0, int(row["start"]) - FLANKING_BP)
        end = int(row["end"]) + FLANKING_BP
        name = f"{row['aqp_family_subfamily']}|{row['gene_id']}|{row['fuente_seq']}"
        score = 1000 if row.get("expressed_in_rnaseq") else 500
        strand = row["strand"] if row["strand"] in ["+", "-"] else "."
        bed_lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")

# Also add all aquaporin genes (not just candidates) for complete visualization
for _, row in ok_genes.iterrows():
    gene_id = row["gene_id"]
    if gene_id in gene_coords:
        gc = gene_coords[gene_id]
        start = max(0, gc["start"] - FLANKING_BP)
        end = gc["end"] + FLANKING_BP
        name = f"{row['aqp_family_subfamily']}|{gene_id}|OK"
        bed_lines.append(f"{gc['chr']}\t{start}\t{end}\t{name}\t100\t{gc['strand']}")

bed_file = os.path.join(RESULTS_DIR, "igv_aquaporin_regions.bed")
with open(bed_file, "w") as f:
    f.write("# BED file for IGV: aquaporin loci with flanking regions\n")
    f.write("# Score 1000 = candidate needing reannotation + expressed\n")
    f.write("# Score 500  = candidate needing reannotation + not expressed\n")
    f.write("# Score 100  = OK aquaporin (for reference)\n")
    for line in sorted(bed_lines):
        f.write(line + "\n")

print(f"# Saved: {bed_file} ({len(bed_lines)} regions)")

# ---- Generate summary report ------------------------------------------------
report_file = os.path.join(RESULTS_DIR, "reannotation_report.txt")
with open(report_file, "w") as f:
    f.write("=" * 70 + "\n")
    f.write("AQUAPORIN REANNOTATION CANDIDATES REPORT\n")
    f.write("=" * 70 + "\n\n")

    f.write(f"Total aquaporins analyzed: {len(aqp)}\n")
    f.write(f"Aquaporins OK (6 TMHs, good annotation): {len(ok_genes)}\n")
    f.write(f"Aquaporins needing reannotation: {len(candidates)}\n\n")

    # Breakdown by category
    f.write("--- Breakdown by category ---\n")
    for src in ["GFF3_FALLBACK", "MAKER_GFF3"]:
        n = len(candidates[candidates["fuente_seq"] == src])
        f.write(f"  {src}: {n}\n")
    n_partial = len(candidates[candidates["TMHs"] != 6])
    f.write(f"  Partial (TMHs != 6): {n_partial}\n")

    if "expressed_in_rnaseq" in candidates.columns:
        n_expressed = candidates["expressed_in_rnaseq"].sum()
        f.write(f"\nExpressed in RNA-seq (priority for reannotation): {n_expressed}\n")

    f.write("\n--- Candidates (sorted by priority) ---\n\n")
    for _, row in candidates.iterrows():
        expressed_tag = " [EXPRESSED]" if row.get("expressed_in_rnaseq") else ""
        f.write(f"{row['gene_id']} — {row['aqp_family_subfamily']}{expressed_tag}\n")
        f.write(f"  Source: {row['fuente_seq']}, Verdict: {row.get('veredicto', 'N/A')}\n")
        f.write(f"  TMHs: {row['TMHs']}, Length: {row['longitud_aa']} aa\n")
        if row["chr"]:
            f.write(f"  Location: {row['chr']}:{row['start']}-{row['end']} ({row['strand']})\n")
        if pd.notna(row.get("max_tpm")):
            f.write(f"  Max TPM: {row['max_tpm']:.2f}, Mean TPM: {row['mean_tpm']:.2f}\n")
        f.write(f"  Reasons: {row['flag_reasons']}\n\n")

    f.write("\n--- Action items ---\n")
    f.write("1. Load BAM files and this BED file in IGV\n")
    f.write("2. For each [EXPRESSED] candidate:\n")
    f.write("   a. Check read coverage (pileups) — does it extend beyond GFF3 boundaries?\n")
    f.write("   b. Look for spliced reads — they define intron/exon boundaries precisely\n")
    f.write("   c. If evidence supports a longer/different gene model, edit GFF3 coordinates\n")
    f.write("   d. Translate the new sequence — verify 6 TMHs and NPA motifs\n")
    f.write("3. For non-expressed candidates: lower priority, may be pseudogenes\n")

print(f"# Saved: {report_file}")
print("# === 09_reannotation_flags.py complete ===")
