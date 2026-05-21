#!/usr/bin/env python3
"""
12_substitute_gff3.py -- Create a corrected GFF3 for Fragaria x ananassa
                         with improved aquaporin gene models from exonerate.

SCIENTIFIC RATIONALE:
    The Benihoppe v1 genome annotation (MAKER-based) contains 109,320 gene models.
    Among the 144 aquaporin genes identified in our analysis, 20 had defective gene
    structures in the original GFF3: truncated CDS, missing exons, or merged loci.
    These 20 genes were re-annotated using exonerate protein2genome alignment against
    curated aquaporin protein templates, producing structurally correct models with
    proper exon-intron boundaries validated by multiple sequence alignment and
    transmembrane helix prediction (TMHMM).

    This script performs a surgical substitution: the 20 exonerate-corrected gene
    blocks replace their defective original counterparts, while all other genes
    (109,300 non-aquaporin + 124 correctly-annotated aquaporin) are passed through
    with metadata attributes indicating aquaporin identity and annotation provenance.

STRATEGY:
    1. Read tabla_Aquaporinas_traduccion.tabular to identify the 20 EXONERATE genes.
    2. Read homeolog_groups.tsv for annotation_source and is_partial metadata.
    3. Parse consensus GFF3 (consenso_aqp_fixed.gff3) to extract exonerate blocks.
    4. Stream the original GFF3 (gzipped, ~1.3M lines):
       - Non-aquaporin genes: pass through unchanged.
       - Non-exonerate aquaporin genes: pass through with added attributes.
       - Exonerate aquaporin genes: SKIP entirely.
    5. Insert normalized exonerate blocks at their correct genomic positions.
    6. Write the output sorted by chromosome and position as gzipped GFF3.

ID NORMALIZATION (exonerate blocks):
    Consensus format:           Normalized format (original-compatible):
    Gene  mRNA_XXXXX-FxaID_gene  ->  FxaID
    mRNA  mRNA_XXXXX-FxaID       ->  FxaID-mRNA-1
    Exon  mRNA_XXXXX-FxaID.exonN ->  FxaID-mRNA-1.exonN
    CDS   mRNA_XXXXX-FxaID.cds.N ->  FxaID-mRNA-1.cds.N

INPUT FILES:
    - Fragaria_ananassa_Benihoppe_Genome.gff3.gz   (original genome annotation)
    - consenso_aqp_fixed.gff3                       (consensus aquaporin GFF3)
    - tabla_Aquaporinas_traduccion.tabular           (aquaporin classification table)
    - homeolog_groups.tsv                            (homeolog group metadata)

OUTPUT:
    - Fragaria_ananassa_Benihoppe_Genome_AQP_corrected.gff3.gz

USAGE:
    python 12_substitute_gff3.py
"""

import gzip
import csv
import re
import os
import sys
from collections import OrderedDict

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
WORK_DIR = os.path.dirname(os.path.abspath(__file__))

ORIGINAL_GFF3 = os.path.join(WORK_DIR,
    "Fragaria_ananassa_Benihoppe_Genome.gff3.gz")
CONSENSUS_GFF3 = os.path.join(WORK_DIR, "consenso_aqp_fixed.gff3")
TABULAR_FILE = os.path.join(WORK_DIR, "tabla_Aquaporinas_traduccion.tabular")
HOMEOLOG_FILE = os.path.join(WORK_DIR, "homeolog_groups.tsv")
OUTPUT_GFF3 = os.path.join(WORK_DIR,
    "Fragaria_ananassa_Benihoppe_Genome_AQP_corrected.gff3.gz")

# Chromosome sort order (1-7, each with A/B/C/D subgenomes)
CHR_ORDER = []
for num in range(1, 8):
    for sub in ["A", "B", "C", "D"]:
        CHR_ORDER.append(f"chr_{num}{sub}")
CHR_RANK = {c: i for i, c in enumerate(CHR_ORDER)}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def parse_gff3_id(attributes_col):
    """Extract the ID value from a GFF3 attributes column (column 9)."""
    m = re.search(r'ID=([^;\s]+)', attributes_col)
    return m.group(1) if m else None


def extract_base_gene_id_from_exonerate(consensus_gene_id):
    """
    Extract the base gene ID from an exonerate consensus gene ID.

    Example: 'mRNA_54367-Fxa2Ag00184_gene' -> 'Fxa2Ag00184'
    """
    # Strip _gene suffix
    s = consensus_gene_id
    if s.endswith("_gene"):
        s = s[:-5]
    # Strip mRNA_XXXXX- prefix
    m = re.match(r'mRNA_\d+-(.+)', s)
    if m:
        return m.group(1)
    return s


def normalize_exonerate_block(lines, base_gene_id):
    """
    Normalize IDs in an exonerate gene block to match original GFF3 format.

    Transforms:
      Gene ID:  mRNA_XXXXX-FxaID_gene  ->  FxaID
      mRNA ID:  mRNA_XXXXX-FxaID       ->  FxaID-mRNA-1
      Exon ID:  mRNA_XXXXX-FxaID.exonN ->  FxaID-mRNA-1.exonN
      CDS ID:   mRNA_XXXXX-FxaID.cds.N ->  FxaID-mRNA-1.cds.N

    Also changes the source column from exonerate:protein2genome:local to maker,
    and adds custom annotation attributes.
    """
    normalized = []
    # Build the prefix that appears in consensus IDs (e.g., "mRNA_54367-Fxa2Ag00184")
    # We detect it from the gene line
    consensus_prefix = None  # e.g., "mRNA_54367-Fxa2Ag00184"

    for line in lines:
        if line.startswith("#") or not line.strip():
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        feature_type = parts[2]
        attrs = parts[8]

        # Change source column to 'maker' for GTF compatibility
        parts[1] = "maker"

        if feature_type == "gene":
            # Extract the consensus prefix from gene ID
            gene_id = parse_gff3_id(attrs)
            # gene_id is like: mRNA_54367-Fxa2Ag00184_gene
            if gene_id and gene_id.endswith("_gene"):
                consensus_prefix = gene_id[:-5]  # mRNA_54367-Fxa2Ag00184

            # Normalize: ID=base_gene_id;Name=base_gene_id + custom attrs
            new_attrs = (
                f"ID={base_gene_id};Name={base_gene_id}"
                f";annotation_source=exonerate;is_aquaporin=true;is_partial=false"
            )
            parts[8] = new_attrs

        elif feature_type == "mRNA":
            mRNA_id = f"{base_gene_id}-mRNA-1"
            new_attrs = (
                f"ID={mRNA_id};Parent={base_gene_id};Name={base_gene_id}"
                f";annotation_source=exonerate;is_aquaporin=true"
            )
            parts[8] = new_attrs

        elif feature_type == "exon":
            # Extract exon number from original ID
            feat_id = parse_gff3_id(attrs)
            exon_num = ""
            if feat_id:
                m = re.search(r'\.exon(\d+)$', feat_id)
                if m:
                    exon_num = m.group(1)
            mRNA_id = f"{base_gene_id}-mRNA-1"
            exon_id = f"{mRNA_id}.exon{exon_num}"
            new_attrs = f"ID={exon_id};Parent={mRNA_id}"
            parts[8] = new_attrs

        elif feature_type == "CDS":
            # Extract CDS number from original ID
            feat_id = parse_gff3_id(attrs)
            cds_num = ""
            if feat_id:
                m = re.search(r'\.cds\.(\d+)$', feat_id)
                if m:
                    cds_num = m.group(1)
            mRNA_id = f"{base_gene_id}-mRNA-1"
            cds_id = f"{mRNA_id}.cds.{cds_num}"
            new_attrs = f"ID={cds_id};Parent={mRNA_id}"
            parts[8] = new_attrs

        else:
            # Other feature types (five_prime_UTR, three_prime_UTR, etc.)
            # Normalize their IDs too
            mRNA_id = f"{base_gene_id}-mRNA-1"
            feat_id = parse_gff3_id(attrs)
            if feat_id and consensus_prefix:
                # Replace consensus prefix with normalized one
                new_feat_id = feat_id.replace(consensus_prefix, mRNA_id)
                new_attrs = f"ID={new_feat_id};Parent={mRNA_id}"
            else:
                new_attrs = f"Parent={mRNA_id}"
            parts[8] = new_attrs

        normalized.append("\t".join(parts))

    return normalized


def chr_sort_key(chrom, pos):
    """Return a sort key for chromosome + position ordering."""
    rank = CHR_RANK.get(chrom, 999)
    return (rank, int(pos))


# ---------------------------------------------------------------------------
# STEP 1: Read the tabular file to identify EXONERATE genes
# ---------------------------------------------------------------------------

def load_exonerate_gene_ids():
    """
    Read tabla_Aquaporinas_traduccion.tabular and return the set of gene_id
    values where fuente_seq == 'EXONERATE'.
    """
    exonerate_ids = set()
    all_aquaporin_ids = set()
    with open(TABULAR_FILE, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"].strip()
            all_aquaporin_ids.add(gene_id)
            if row["fuente_seq"].strip() == "EXONERATE":
                exonerate_ids.add(gene_id)

    print(f"[Step 1] Loaded {len(all_aquaporin_ids)} aquaporin gene IDs "
          f"from tabular file.")
    print(f"         {len(exonerate_ids)} genes with fuente_seq=EXONERATE.")
    return exonerate_ids, all_aquaporin_ids


# ---------------------------------------------------------------------------
# STEP 2: Read homeolog_groups.tsv for metadata
# ---------------------------------------------------------------------------

def load_homeolog_metadata():
    """
    Read homeolog_groups.tsv and return a dict:
      gene_id -> {annotation_source, is_partial, family, subfamily, ...}
    """
    metadata = {}
    with open(HOMEOLOG_FILE, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"].strip()
            metadata[gene_id] = {
                "annotation_source": row["annotation_source"].strip(),
                "is_partial": row["is_partial"].strip(),
                "family": row["family"].strip(),
                "subfamily": row["subfamily"].strip(),
                "homeolog_group": row["homeolog_group"].strip(),
                "is_tandem_duplicate": row["is_tandem_duplicate"].strip(),
            }
    print(f"[Step 2] Loaded metadata for {len(metadata)} genes from "
          f"homeolog_groups.tsv.")
    return metadata


# ---------------------------------------------------------------------------
# STEP 3: Parse consensus GFF3 for exonerate gene blocks
# ---------------------------------------------------------------------------

def load_exonerate_blocks(exonerate_ids):
    """
    Parse consenso_aqp_fixed.gff3 and extract all exonerate-source gene blocks.

    Returns a dict: base_gene_id -> list of normalized GFF3 lines.
    Each block is identified by source column containing 'exonerate'.
    Gene blocks are separated by '###' lines in the consensus.
    """
    blocks = {}
    current_block = []
    current_is_exonerate = False
    current_gene_id_consensus = None

    with open(CONSENSUS_GFF3, "r") as fh:
        for line in fh:
            stripped = line.strip()

            # Block separator
            if stripped == "###" or stripped == "":
                if current_block and current_is_exonerate and current_gene_id_consensus:
                    base_id = extract_base_gene_id_from_exonerate(
                        current_gene_id_consensus)
                    if base_id in exonerate_ids:
                        normalized = normalize_exonerate_block(
                            current_block, base_id)
                        blocks[base_id] = normalized
                current_block = []
                current_is_exonerate = False
                current_gene_id_consensus = None
                continue

            # Skip header lines
            if stripped.startswith("#"):
                continue

            parts = stripped.split("\t")
            if len(parts) < 9:
                continue

            current_block.append(stripped)

            # Detect if this block is exonerate-sourced
            if "exonerate" in parts[1]:
                current_is_exonerate = True

            # Capture the gene ID
            if parts[2] == "gene" and current_gene_id_consensus is None:
                current_gene_id_consensus = parse_gff3_id(parts[8])

    # Handle the last block if file does not end with ###
    if current_block and current_is_exonerate and current_gene_id_consensus:
        base_id = extract_base_gene_id_from_exonerate(
            current_gene_id_consensus)
        if base_id in exonerate_ids:
            normalized = normalize_exonerate_block(current_block, base_id)
            blocks[base_id] = normalized

    print(f"[Step 3] Extracted and normalized {len(blocks)} exonerate gene "
          f"blocks from consensus GFF3.")

    # Verify all expected exonerate genes are found
    missing = exonerate_ids - set(blocks.keys())
    if missing:
        print(f"  WARNING: Missing exonerate blocks for: {sorted(missing)}")
    else:
        print(f"  All {len(exonerate_ids)} exonerate genes found in consensus.")

    return blocks


# ---------------------------------------------------------------------------
# STEP 4: Get chromosome and position for each exonerate block
# ---------------------------------------------------------------------------

def get_exonerate_positions(blocks):
    """
    For each exonerate block, extract (chromosome, start_position) from the
    gene line for sorted insertion.

    Returns: dict of base_gene_id -> (chromosome, start_position)
    """
    positions = {}
    for gene_id, lines in blocks.items():
        for line in lines:
            parts = line.split("\t")
            if len(parts) >= 5 and parts[2] == "gene":
                positions[gene_id] = (parts[0], int(parts[3]))
                break
    return positions


# ---------------------------------------------------------------------------
# STEP 5: Process the original GFF3 and produce the corrected output
# ---------------------------------------------------------------------------

def process_original_gff3(exonerate_ids, all_aquaporin_ids, homeolog_meta,
                          exonerate_blocks, exonerate_positions):
    """
    Stream through the original GFF3 (gzipped), and produce a list of
    (chromosome, position, line_group) tuples for all gene blocks.

    - Non-aquaporin genes: passed through verbatim.
    - Non-exonerate aquaporin genes: passed through with added attributes
      on the gene line.
    - Exonerate aquaporin genes: SKIPPED (will be replaced by consensus blocks).
    """
    # We will collect gene blocks as: (chrom, start_pos, [lines])
    output_blocks = []

    # State variables for tracking gene blocks
    current_gene_id = None
    current_chrom = None
    current_start = None
    current_block_lines = []
    current_is_exonerate_skip = False
    current_is_aquaporin = False

    def flush_block():
        """Save the current block to output (unless it should be skipped)."""
        nonlocal current_gene_id, current_chrom, current_start
        nonlocal current_block_lines, current_is_exonerate_skip
        nonlocal current_is_aquaporin

        if current_block_lines and not current_is_exonerate_skip:
            output_blocks.append(
                (current_chrom, current_start, current_block_lines))

        # Reset
        current_gene_id = None
        current_chrom = None
        current_start = None
        current_block_lines = []
        current_is_exonerate_skip = False
        current_is_aquaporin = False

    lines_read = 0
    genes_total = 0
    genes_skipped = 0
    genes_annotated = 0

    print(f"[Step 5] Processing original GFF3: {ORIGINAL_GFF3}")

    with gzip.open(ORIGINAL_GFF3, "rt") as fh:
        for line in fh:
            lines_read += 1
            stripped = line.rstrip("\n")

            if lines_read % 200000 == 0:
                print(f"  ... processed {lines_read:,} lines, "
                      f"{genes_total:,} genes ...")

            # Handle comment and separator lines
            if stripped.startswith("#") or stripped == "":
                # If inside a block, just add; otherwise ignore
                if current_gene_id is not None:
                    current_block_lines.append(stripped)
                continue

            parts = stripped.split("\t")
            if len(parts) < 9:
                # Malformed line, pass through if inside block
                if current_gene_id is not None:
                    current_block_lines.append(stripped)
                continue

            feature_type = parts[2]
            chrom = parts[0]
            start_pos = int(parts[3])

            # Detect start of a new gene block
            if feature_type == "gene":
                # Flush previous block
                if current_gene_id is not None:
                    flush_block()

                genes_total += 1
                gene_id = parse_gff3_id(parts[8])
                current_gene_id = gene_id
                current_chrom = chrom
                current_start = start_pos

                if gene_id in exonerate_ids:
                    # This gene will be replaced by exonerate block
                    current_is_exonerate_skip = True
                    genes_skipped += 1

                elif gene_id in all_aquaporin_ids:
                    # Non-exonerate aquaporin: add metadata attributes
                    current_is_aquaporin = True
                    genes_annotated += 1
                    meta = homeolog_meta.get(gene_id, {})
                    ann_src = meta.get("annotation_source", "unknown")
                    is_partial = meta.get("is_partial", "no")
                    is_partial_bool = "true" if is_partial == "yes" else "false"

                    # Append custom attributes to the gene line
                    existing_attrs = parts[8].rstrip()
                    new_attrs = (
                        f"{existing_attrs}"
                        f";annotation_source={ann_src}"
                        f";is_aquaporin=true"
                        f";is_partial={is_partial_bool}"
                    )
                    parts[8] = new_attrs
                    stripped = "\t".join(parts)

                # else: non-aquaporin gene, pass through unchanged

                current_block_lines.append(stripped)
            else:
                # Non-gene line: part of current block
                current_block_lines.append(stripped)

    # Flush the last block
    if current_gene_id is not None:
        flush_block()

    print(f"  Read {lines_read:,} lines, {genes_total:,} gene blocks total.")
    print(f"  Skipped {genes_skipped} exonerate gene blocks (to be replaced).")
    print(f"  Annotated {genes_annotated} non-exonerate aquaporin gene blocks.")
    print(f"  Passed through {genes_total - genes_skipped - genes_annotated:,} "
          f"non-aquaporin gene blocks unchanged.")

    # Now add the exonerate replacement blocks
    print(f"[Step 6] Inserting {len(exonerate_blocks)} exonerate replacement "
          f"blocks ...")
    for gene_id, block_lines in exonerate_blocks.items():
        chrom, start = exonerate_positions[gene_id]
        output_blocks.append((chrom, start, block_lines))

    return output_blocks


# ---------------------------------------------------------------------------
# STEP 6: Sort and write the output GFF3
# ---------------------------------------------------------------------------

def write_output_gff3(output_blocks):
    """
    Sort all gene blocks by chromosome and position, then write as gzipped GFF3.
    """
    print(f"[Step 7] Sorting {len(output_blocks):,} gene blocks by "
          f"chromosome and position ...")

    # Sort by chromosome rank, then start position
    output_blocks.sort(key=lambda b: (CHR_RANK.get(b[0], 999), b[1]))

    print(f"[Step 8] Writing output: {OUTPUT_GFF3}")

    total_lines = 0
    total_genes = 0

    with gzip.open(OUTPUT_GFF3, "wt") as out:
        # Write GFF3 header
        out.write("##gff-version 3\n")

        for chrom, start, block_lines in output_blocks:
            for bline in block_lines:
                out.write(bline + "\n")
                total_lines += 1
                # Count gene lines for verification
                parts = bline.split("\t")
                if len(parts) >= 3 and parts[2] == "gene":
                    total_genes += 1

    print(f"  Wrote {total_lines:,} lines ({total_genes:,} genes) "
          f"to {os.path.basename(OUTPUT_GFF3)}.")
    return total_genes


# ---------------------------------------------------------------------------
# STEP 7: Verification
# ---------------------------------------------------------------------------

def verify_output(total_genes, exonerate_ids, all_aquaporin_ids):
    """
    Quick verification of the output file.
    """
    print(f"\n[Verification]")
    # Expected: 109,320 original genes
    #         - 20 exonerate genes removed
    #         + 20 exonerate genes added back (corrected)
    #         = 109,320 total genes
    expected_genes = 109320
    print(f"  Expected genes: {expected_genes}")
    print(f"  Actual genes:   {total_genes}")

    if total_genes == expected_genes:
        print(f"  PASS: Gene count matches.")
    else:
        print(f"  WARNING: Gene count mismatch! "
              f"Difference: {total_genes - expected_genes}")

    # Verify exonerate genes are present with correct attributes
    exo_found = set()
    aqp_found = set()
    with gzip.open(OUTPUT_GFF3, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            gene_id = parse_gff3_id(parts[8])
            if gene_id in exonerate_ids:
                exo_found.add(gene_id)
                if "annotation_source=exonerate" in parts[8]:
                    pass  # correct
                else:
                    print(f"  WARNING: Exonerate gene {gene_id} missing "
                          f"annotation_source=exonerate attribute.")
            if gene_id in all_aquaporin_ids:
                aqp_found.add(gene_id)

    print(f"  Exonerate genes in output: {len(exo_found)}/20")
    print(f"  Total aquaporin genes in output: {len(aqp_found)}/144")

    missing_exo = exonerate_ids - exo_found
    if missing_exo:
        print(f"  WARNING: Missing exonerate genes: {sorted(missing_exo)}")

    missing_aqp = all_aquaporin_ids - aqp_found
    if missing_aqp:
        print(f"  WARNING: Missing aquaporin genes: {sorted(missing_aqp)}")
    else:
        print(f"  PASS: All 144 aquaporin genes present in output.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("12_substitute_gff3.py")
    print("Creating corrected GFF3 with exonerate aquaporin gene models")
    print("=" * 72)
    print()

    # Verify input files exist
    for path, desc in [
        (ORIGINAL_GFF3, "Original genome GFF3"),
        (CONSENSUS_GFF3, "Consensus aquaporin GFF3"),
        (TABULAR_FILE, "Aquaporin translation table"),
        (HOMEOLOG_FILE, "Homeolog groups metadata"),
    ]:
        if not os.path.isfile(path):
            print(f"ERROR: {desc} not found: {path}")
            sys.exit(1)
        size_mb = os.path.getsize(path) / (1024 * 1024)
        print(f"  Found {desc}: {os.path.basename(path)} ({size_mb:.1f} MB)")
    print()

    # Step 1: Identify exonerate genes
    exonerate_ids, all_aquaporin_ids = load_exonerate_gene_ids()
    print()

    # Step 2: Load homeolog metadata
    homeolog_meta = load_homeolog_metadata()
    print()

    # Step 3: Extract and normalize exonerate blocks from consensus
    exonerate_blocks = load_exonerate_blocks(exonerate_ids)
    print()

    # Step 4: Get positions for exonerate blocks
    exonerate_positions = get_exonerate_positions(exonerate_blocks)
    print()

    # Step 5-6: Process original GFF3
    output_blocks = process_original_gff3(
        exonerate_ids, all_aquaporin_ids, homeolog_meta,
        exonerate_blocks, exonerate_positions)
    print()

    # Step 7-8: Sort and write output
    total_genes = write_output_gff3(output_blocks)
    print()

    # Step 9: Verification
    verify_output(total_genes, exonerate_ids, all_aquaporin_ids)

    print()
    print("=" * 72)
    print("Done.")
    print("=" * 72)


if __name__ == "__main__":
    main()
