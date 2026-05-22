#!/usr/bin/env python3
"""
11_homeolog_grouping.py — Algorithmic Homeolog Group Identification
for Fragaria × ananassa aquaporins from phylogenetic tree topology.

SCIENTIFIC BASIS:
    Fragaria × ananassa is an allo-octoploid (2n=8x=56) derived from four
    diploid progenitors contributing subgenomes A, B, C, D. Each ancestral
    gene can have up to four homeologous copies, one per subgenome, located
    on homeologous chromosomes (same chromosome number, different subgenome
    letter). Homeologs diverged at (or shortly after) the polyploidization
    events and form monophyletic clades in gene trees, separated from
    paralogs (ancient duplications) by outgroup orthologs.

ALGORITHM (3 phases):
    Phase 1 — Maximal Fxa-only clade detection
        Parse the IQ-TREE3 maximum-likelihood phylogeny (Q.PLANT+R7 model,
        1000 ultrafast bootstraps). Root on the monocot outgroup (Oryza
        sativa, Os). Traverse the tree bottom-up (postorder) and mark each
        node as "Fxa-only" if ALL descendant leaves are Fragaria × ananassa
        genes. A node is a "maximal Fxa-only clade" if it is Fxa-only but
        its parent is not. These clades correspond to groups of Fxa genes
        not interrupted by any outgroup ortholog.

    Phase 2 — Chromosome-based splitting
        Some maximal Fxa-only clades contain genes from different chromosome
        numbers (i.e., paralogs on different chromosomes that lack an
        outgroup gene between them). These are split into per-chromosome
        groups, first attempting a tree-topology-based split (if the binary
        tree structure cleanly separates chromosomes), then falling back to
        direct chromosome assignment.

    Phase 3 — Distance-aware merging
        Due to incomplete lineage sorting, gene tree discord, or stochastic
        outgroup placement, some true homeolog groups may be split into
        multiple fragments (separate Fxa-only clades or singletons) by a
        small number of outgroup insertions. These fragments are merged if:
          (a) they share the same dominant chromosome number,
          (b) their MRCA in the tree contains ≤ MAX_OUTGROUP_INSERTIONS
              non-Fxa genes between the fragments, AND
          (c) the minimum inter-fragment phylogenetic distance is ≤ the
              empirically determined homeolog distance threshold (95th
              percentile of within-group distances × safety margin).
        The threshold is estimated from Phase 2 groups that already contain
        multiple subgenomes on a single chromosome (clear homeolog signal).

SUBFAMILY ASSIGNMENT:
    For each homeolog group, the aquaporin subfamily (e.g., PIP2, TIP1) is
    determined by finding the closest non-Fxa gene in the tree for a
    representative member and parsing the subfamily from the outgroup
    gene's name (e.g., MdPIP2_6 → PIP2). This is cross-validated against
    the subfamily annotations in the metadata table for GFF3-confirmed
    members. Discrepancies between the tree-derived and metadata subfamilies
    for MAKER_GFF3 genes are flagged as reclassifications.

TANDEM DUPLICATE DETECTION:
    Within each homeolog group, if multiple genes originate from the same
    subgenome, they are flagged as tandem duplicates. This is validated
    by checking for adjacent gene numbers on the same chromosome arm.

PCA VALIDATION:
    When PCA_Coordenadas_Finales.csv is available, protein-level PCA
    coordinates are used as an independent validation of homeolog groups.
    Homeologs encode nearly identical proteins and should cluster tightly
    in PCA space. Groups with high PCA spread may contain misassigned
    genes. PCA distances are also used as an additional criterion in the
    merge step to prevent merging of phylogenetically adjacent but
    functionally divergent genes.

INPUT:
    fxa_final.treefile                      — IQ-TREE3 ML tree (304 seqs)
    tabla_Aquaporinas_traduccion.tabular    — Gene metadata table
    PCA_Coordenadas_Finales.csv             — PCA coordinates (optional)

OUTPUT:
    homeolog_groups.tsv         — Master table (144 genes × 16 columns)
    homeolog_groups_summary.tsv — Per-group summary

USAGE:
    python 11_homeolog_grouping.py

REFERENCE:
    Tree inferred with IQ-TREE3 v3.0 (Minh et al., 2020) using the
    Q.PLANT+R7 substitution model. Outgroups: Malus domestica (Md),
    Hevea brasiliensis (Hb), Arabidopsis thaliana (At), Oryza sativa (Os).
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

try:
    from Bio import Phylo
except ImportError:
    sys.exit("ERROR: BioPython is required. Install with: pip install biopython")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CONFIGURATION
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Outgroup species prefix for rooting (monocot, most distant from Fragaria)
ROOT_OUTGROUP_PREFIX = "Os"

# Maximum number of non-Fxa genes tolerated between two fragments
# when evaluating merge candidates. Biologically motivated: allows for
# 1-3 outgroup insertions due to ILS or tree estimation error.
MAX_OUTGROUP_INSERTIONS = 3

# Safety factor applied to the estimated homeolog distance threshold.
# The threshold is: 95th percentile of within-group distances × this factor.
DISTANCE_THRESHOLD_FACTOR = 1.5

# Minimum number of genes with known subfamily (from metadata) needed
# to establish a group's subfamily by majority vote (used as fallback
# when outgroup neighbor method is ambiguous).
MIN_GENES_FOR_MAJORITY_VOTE = 2

# Regex patterns
FXA_GENE_PATTERN = re.compile(r"Fxa(\d)([ABCD])g(\d+)")
FXA_BASE_ID_PATTERN = re.compile(r"(Fxa\d[ABCD]g\d+)")
OUTGROUP_SUBFAMILY_PATTERN = re.compile(
    r"^[A-Za-z]{2,3}((?:PIP|TIP|NIP|SIP|XIP)\d+)"
)
# Pattern to extract sub-subfamily: e.g. MdPIP2_6 → ("PIP2", "6") → PIP2;6
OUTGROUP_SUB_SUBFAMILY_PATTERN = re.compile(
    r"^([A-Za-z]{2,3})((?:PIP|TIP|NIP|SIP|XIP)(\d+))_(\d+)"
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GENE NAME PARSING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def is_fxa_gene(name):
    """Test whether a tree leaf name corresponds to a Fragaria × ananassa gene."""
    return bool(FXA_GENE_PATTERN.search(name or ""))


def parse_fxa_info(tree_name):
    """Extract structured metadata from an Fxa tree leaf name.

    Handles all naming conventions in the tree:
        Fxa1Ag01329                   (standard GFF3)
        mRNA_54367-Fxa2Ag00184        (exonerate prefix)
        Fxa1Bg02622-partial-maker     (partial + maker suffix)
        Fxa4Ag01484-partial           (partial suffix)
        Fxa6Bg00715-maker             (maker suffix)

    Returns dict with: base_id, tree_name, chr, subgenome, gene_num
    or None if name does not match the Fxa pattern.
    """
    m = FXA_GENE_PATTERN.search(tree_name)
    if not m:
        return None
    base_match = FXA_BASE_ID_PATTERN.search(tree_name)
    return {
        "base_id": base_match.group(1),
        "tree_name": tree_name,
        "chr": int(m.group(1)),
        "subgenome": m.group(2),
        "gene_num": int(m.group(3)),
    }


def parse_outgroup_subfamily(name):
    """Parse subfamily from an outgroup gene name.

    Examples:
        MdPIP2_6  → ("PIP2", "PIP")
        AtTIP3_1  → ("TIP3", "TIP")
        OsNIP5_3  → ("NIP5", "NIP")
        FaPIP1_1  → ("PIP1", "PIP")

    Returns (subfamily, family) or (None, None) if unparseable.
    """
    m = OUTGROUP_SUBFAMILY_PATTERN.match(name or "")
    if m:
        subfamily = m.group(1)
        family = re.match(r"(PIP|TIP|NIP|SIP|XIP)", subfamily).group(1)
        return subfamily, family
    return None, None


def parse_outgroup_sub_subfamily(name):
    """Parse full sub-subfamily from an outgroup gene name.

    Examples:
        MdPIP2_6  → ("PIP2;6", "PIP2", "PIP", "Md")
        AtTIP3_1  → ("TIP3;1", "TIP3", "TIP", "At")
        OsNIP5_3  → ("NIP5;3", "NIP5", "NIP", "Os")
        FaPIP1_1  → ("PIP1;1", "PIP1", "PIP", "Fa")

    Returns (sub_subfamily, subfamily, family, species) or (None,None,None,None).
    """
    m = OUTGROUP_SUB_SUBFAMILY_PATTERN.match(name or "")
    if m:
        species = m.group(1)
        subfamily = m.group(2)  # e.g. PIP2
        family = re.match(r"(PIP|TIP|NIP|SIP|XIP)", subfamily).group(1)
        isoform = m.group(4)    # e.g. 6
        sub_subfamily = f"{subfamily};{isoform}"  # e.g. PIP2;6
        return sub_subfamily, subfamily, family, species
    return None, None, None, None


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# TREE OPERATIONS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def build_parent_map(tree):
    """Build a child_id → parent_clade mapping for the entire tree.

    Uses id() of clade objects as keys because Bio.Phylo Clade objects
    are not hashable by default.
    """
    parent_map = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parent_map[id(child)] = clade
    return parent_map


def root_tree(tree, outgroup_prefix):
    """Root the tree on the specified outgroup species.

    Outgroup rooting with the monocot O. sativa places the root between
    monocots and dicots, which is the accepted phylogenetic position for
    plant aquaporin trees. Falls back to midpoint rooting if outgroup
    terminals are absent.
    """
    outgroup_terminals = [
        t for t in tree.get_terminals() if t.name.startswith(outgroup_prefix)
    ]
    if not outgroup_terminals:
        print(f"  WARNING: No '{outgroup_prefix}' outgroup found; using midpoint rooting")
        tree.root_at_midpoint()
        return "midpoint"

    try:
        tree.root_with_outgroup(*outgroup_terminals)
        return f"outgroup ({outgroup_prefix}, n={len(outgroup_terminals)})"
    except Exception as e:
        print(f"  WARNING: Outgroup rooting failed ({e}); using midpoint rooting")
        tree.root_at_midpoint()
        return "midpoint (fallback)"


def get_fxa_leaf_names(clade):
    """Return a list of Fxa tree leaf names under a clade."""
    if clade.is_terminal():
        return [clade.name] if is_fxa_gene(clade.name) else []
    return [t.name for t in clade.get_terminals() if is_fxa_gene(t.name)]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PHASE 1: MAXIMAL FXA-ONLY CLADE DETECTION
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def find_maximal_fxa_clades(tree):
    """Identify maximal monophyletic clades containing only Fxa genes.

    A clade is "Fxa-only" if every leaf descendant is an Fxa gene.
    A clade is "maximal" if it is Fxa-only but its parent is not
    (i.e., adding any more nodes would include a non-Fxa gene).

    These clades are the primary signal for homeolog groups: genes within
    the same Fxa-only clade diverged more recently than from any outgroup
    ortholog, consistent with post-polyploidization divergence.

    Also collects singleton Fxa terminals (those not inside any multi-gene
    Fxa-only clade), which arise when the gene's immediate sister in the
    tree is a non-Fxa gene.

    Returns:
        List of Bio.Phylo Clade objects (internal nodes or terminal singletons).
    """
    parent_map = build_parent_map(tree)

    # Bottom-up labeling
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            clade._fxa_only = is_fxa_gene(clade.name)
        else:
            clade._fxa_only = all(c._fxa_only for c in clade.clades)

    # Collect maximal Fxa-only clades (multi-gene)
    maximal_clades = []
    for clade in tree.find_clades(order="postorder"):
        if clade._fxa_only and not clade.is_terminal():
            parent = parent_map.get(id(clade))
            if parent is None or not parent._fxa_only:
                maximal_clades.append(clade)

    # Identify Fxa terminals not covered by any multi-gene clade
    covered = set()
    for clade in maximal_clades:
        for t in clade.get_terminals():
            covered.add(t.name)

    singletons = []
    for t in tree.get_terminals():
        if is_fxa_gene(t.name) and t.name not in covered:
            singletons.append(t)

    return maximal_clades, singletons


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PHASE 2: CHROMOSOME-BASED SPLITTING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def split_clade_by_chromosome(clade, gene_info):
    """Split a multi-chromosome Fxa clade into per-chromosome groups.

    In allopolyploids, homeologs reside on the same chromosome number
    across subgenomes (syntenic position). A maximal Fxa-only clade
    containing genes from different chromosome numbers contains paralogs
    (genes duplicated before polyploidization) that lack an outgroup
    separator in the tree.

    Strategy:
        1. If the clade is binary and each child subtree has a distinct
           set of chromosome numbers → split at the root (tree-guided).
        2. Otherwise → assign genes to groups by chromosome number
           (chromosome-guided fallback).

    Returns list of lists of tree_names (each list = one chromosome group).
    """
    leaves = get_fxa_leaf_names(clade)
    if not leaves:
        return []

    chr_map = {}
    for name in leaves:
        info = gene_info.get(name)
        if info:
            chr_map[name] = info["chr"]

    unique_chrs = set(chr_map.values())
    if len(unique_chrs) <= 1:
        return [leaves]

    # Attempt tree-topology-based splitting (recursive)
    if not clade.is_terminal() and len(clade.clades) >= 2:
        child_data = []
        for child in clade.clades:
            child_leaves = get_fxa_leaf_names(child)
            if child_leaves:
                child_chrs = set(chr_map.get(n) for n in child_leaves if n in chr_map)
                child_data.append((child, child_leaves, child_chrs))

        if len(child_data) >= 2:
            # Check for non-overlapping chromosome sets across children
            all_separate = True
            for i in range(len(child_data)):
                for j in range(i + 1, len(child_data)):
                    if child_data[i][2] & child_data[j][2]:
                        all_separate = False
                        break
                if not all_separate:
                    break

            if all_separate:
                result = []
                for child, _, _ in child_data:
                    result.extend(split_clade_by_chromosome(child, gene_info))
                return result

    # Fallback: group directly by chromosome number
    by_chr = defaultdict(list)
    for name in leaves:
        if name in chr_map:
            by_chr[chr_map[name]].append(name)
    return list(by_chr.values())


def phase2_split(clades, singletons, gene_info):
    """Apply chromosome-based splitting to all Phase 1 clades.

    Returns a flat list of gene name lists (each = one candidate group).
    """
    groups = []
    split_count = 0

    for clade in clades:
        sub_groups = split_clade_by_chromosome(clade, gene_info)
        if len(sub_groups) > 1:
            split_count += 1
        groups.extend(sub_groups)

    # Add singletons as individual groups
    for s in singletons:
        groups.append([s.name])

    return groups, split_count


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PHASE 3: DISTANCE-AWARE MERGING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def estimate_homeolog_distance_threshold(tree, groups, gene_info, terminal_map):
    """Estimate the maximum phylogenetic distance expected within a homeolog group.

    Calibration strategy:
        Collect pairwise distances within groups that (a) span ≥2 subgenomes
        and (b) are restricted to a single chromosome number. These represent
        confident homeolog pairs. The threshold is set at the 95th percentile
        of these distances, multiplied by DISTANCE_THRESHOLD_FACTOR to provide
        a safety margin for groups with missing subgenomes or longer branches
        due to relaxed selection.

    Returns:
        threshold (float): maximum distance for merge eligibility.
        calibration_distances (list): the distances used for calibration.
    """
    calibration_distances = []

    for group in groups:
        if len(group) < 2:
            continue

        chrs = set()
        subgenomes = set()
        for name in group:
            info = gene_info.get(name)
            if info:
                chrs.add(info["chr"])
                subgenomes.add(info["subgenome"])

        if len(chrs) == 1 and len(subgenomes) >= 2:
            # This is a confident homeolog group — collect within-group distances
            for i, n1 in enumerate(group):
                for n2 in group[i + 1 :]:
                    t1 = terminal_map.get(n1)
                    t2 = terminal_map.get(n2)
                    if t1 and t2:
                        d = tree.distance(t1, t2)
                        calibration_distances.append(d)

    if not calibration_distances:
        print("  WARNING: No calibration distances found; using default threshold 0.5")
        return 0.5, []

    calibration_distances.sort()
    idx_95 = min(int(0.95 * len(calibration_distances)), len(calibration_distances) - 1)
    p95 = calibration_distances[idx_95]
    threshold = p95 * DISTANCE_THRESHOLD_FACTOR

    return threshold, calibration_distances


def count_outgroup_between(tree, group1, group2, gene_info, terminal_map):
    """Count non-Fxa genes under the MRCA of two groups.

    This measures how many outgroup genes are 'between' the two groups
    in the tree. A low count (≤ MAX_OUTGROUP_INSERTIONS) suggests the
    groups were artificially split by a small number of outgroup insertions.

    Other Fxa genes under the MRCA are handled as follows:
      - Same chromosome as the merge candidate: block the merge (these may
        belong to a different homeolog group on the same chromosome, and
        merging could incorrectly combine distinct homeolog sets).
      - Different chromosome: do not block (these are paralogs that happen
        to be phylogenetically adjacent; their presence is expected in
        gene families with tandem duplications across chromosomes).

    Returns:
        n_outgroup (int): number of non-Fxa leaves under the MRCA.
        Returns 999 to block merge when same-chromosome Fxa interference
        is detected.
    """
    all_terminals = []
    for name in group1 + group2:
        t = terminal_map.get(name)
        if t:
            all_terminals.append(t)

    if len(all_terminals) < 2:
        return 999

    try:
        mrca = tree.common_ancestor(all_terminals)
    except Exception:
        return 999

    # Determine the chromosomes of the merge candidate
    merge_chrs = set()
    for name in group1 + group2:
        info = gene_info.get(name)
        if info:
            merge_chrs.add(info["chr"])

    all_leaves = mrca.get_terminals()
    fxa_in_groups = set(group1 + group2)
    n_outgroup = 0

    for leaf in all_leaves:
        if leaf.name in fxa_in_groups:
            continue
        if not is_fxa_gene(leaf.name):
            n_outgroup += 1
        else:
            # Other Fxa gene: only block if same chromosome
            other_info = gene_info.get(leaf.name)
            if other_info and other_info["chr"] in merge_chrs:
                return 999  # Same-chromosome Fxa interference
            # Different chromosome → tolerate (paralog, not a competing homeolog)

    return n_outgroup


def min_inter_group_distance(tree, group1, group2, terminal_map):
    """Minimum phylogenetic distance between any gene in group1 and group2."""
    min_d = float("inf")
    for n1 in group1:
        for n2 in group2:
            t1 = terminal_map.get(n1)
            t2 = terminal_map.get(n2)
            if t1 and t2:
                d = tree.distance(t1, t2)
                if d < min_d:
                    min_d = d
    return min_d


def phase3_merge(tree, groups, gene_info, terminal_map, pca_coords=None):
    """Merge fragments of the same homeolog group separated by outgroup insertions.

    Iteratively evaluates all pairs of groups on the same chromosome.
    Pairs are ranked by inter-group distance (closest first) and merged
    greedily if they satisfy the outgroup count, distance, and (optionally)
    PCA proximity criteria.

    When PCA coordinates are available, an additional check ensures that
    candidate merge groups are close in PCA space (protein sequence
    similarity), preventing merges of phylogenetically adjacent but
    functionally divergent genes.

    Returns merged groups and a log of merge operations.
    """
    # Estimate distance threshold from confident homeolog groups
    threshold, cal_dists = estimate_homeolog_distance_threshold(
        tree, groups, gene_info, terminal_map
    )
    merge_log = []
    merge_log.append(
        f"Phylogenetic distance threshold: {threshold:.4f} "
        f"(95th pctl of {len(cal_dists)} calibration distances "
        f"x {DISTANCE_THRESHOLD_FACTOR})"
    )

    # PCA threshold (if coordinates available)
    pca_threshold = None
    if pca_coords:
        pca_threshold = estimate_pca_threshold(groups, gene_info, pca_coords)
        merge_log.append(f"PCA distance threshold: {pca_threshold:.4f}")

    def dominant_chr(group):
        """Most frequent chromosome number in a group."""
        chr_counts = defaultdict(int)
        for name in group:
            info = gene_info.get(name)
            if info:
                chr_counts[info["chr"]] += 1
        if not chr_counts:
            return None
        return max(chr_counts, key=chr_counts.get)

    changed = True
    n_merges = 0
    while changed:
        changed = False

        # Build list of all candidate merge pairs (same chromosome)
        candidates = []
        for i in range(len(groups)):
            chr_i = dominant_chr(groups[i])
            for j in range(i + 1, len(groups)):
                chr_j = dominant_chr(groups[j])
                if chr_i is not None and chr_i == chr_j:
                    d = min_inter_group_distance(
                        tree, groups[i], groups[j], terminal_map
                    )
                    candidates.append((d, i, j))

        # Sort by distance (merge closest pairs first)
        candidates.sort()

        for d, i, j in candidates:
            # Criterion 1: distance within homeolog range
            if d > threshold:
                continue

            # Criterion 2: few outgroup insertions between groups
            n_out = count_outgroup_between(
                tree, groups[i], groups[j], gene_info, terminal_map
            )
            if n_out > MAX_OUTGROUP_INSERTIONS:
                continue

            # Criterion 3: Subgenome compatibility
            # If both groups have unique (non-tandem) subgenome representation,
            # block the merge when their subgenomes overlap. Overlapping
            # unique subgenomes indicate the groups are distinct homeolog
            # sets (paralogs), not fragments of the same homeolog group.
            # Groups with existing tandem duplicates (>1 gene per subgenome)
            # are exempt, since tandem arrays naturally have overlapping
            # subgenomes when split by outgroup insertions.
            sg1_counts = defaultdict(int)
            sg2_counts = defaultdict(int)
            for n in groups[i]:
                _info = gene_info.get(n)
                if _info:
                    sg1_counts[_info["subgenome"]] += 1
            for n in groups[j]:
                _info = gene_info.get(n)
                if _info:
                    sg2_counts[_info["subgenome"]] += 1
            has_tandems_1 = any(c > 1 for c in sg1_counts.values())
            has_tandems_2 = any(c > 1 for c in sg2_counts.values())
            if not has_tandems_1 and not has_tandems_2:
                common_sgs = set(sg1_counts.keys()) & set(sg2_counts.keys())
                if common_sgs:
                    continue  # Would create subgenome duplicates → distinct groups

            # Criterion 4 (optional): PCA proximity validation
            if pca_coords:
                pca_d = pca_inter_group_distance(
                    groups[i], groups[j], gene_info, pca_coords
                )
                if pca_d != float("inf") and pca_d > pca_threshold:
                    continue

            # Merge j into i
            base_ids_i = sorted(
                gene_info[n]["base_id"] for n in groups[i] if n in gene_info
            )
            base_ids_j = sorted(
                gene_info[n]["base_id"] for n in groups[j] if n in gene_info
            )
            merge_log.append(
                f"MERGE: [{','.join(base_ids_i)}] + [{','.join(base_ids_j)}] "
                f"(dist={d:.4f}, outgroup_between={n_out})"
            )

            groups[i] = groups[i] + groups[j]
            groups.pop(j)
            changed = True
            n_merges += 1
            break

    merge_log.append(f"Total merges: {n_merges}")
    return groups, merge_log


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# PCA COORDINATES (PROTEIN SEQUENCE SPACE VALIDATION)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def load_pca_coordinates(filepath):
    """Load PCA coordinates from PCA_Coordenadas_Finales.csv.

    The PCA was computed on protein sequence features (motifs, diagnostic
    scores) and captures functional similarity between aquaporins. Genes
    in the same homeolog group should cluster tightly in PCA space because
    they encode nearly identical proteins (recent divergence).

    Returns dict: gene_id → (PC1, PC2).
    """
    coords = {}
    if not os.path.exists(filepath):
        return coords
    with open(filepath, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gid = row.get("ID", "").strip()
            try:
                pc1 = float(row.get("PC1", 0))
                pc2 = float(row.get("PC2", 0))
                coords[gid] = (pc1, pc2)
            except (ValueError, TypeError):
                continue
    return coords


def pca_euclidean_distance(coord1, coord2):
    """Euclidean distance between two PCA coordinate pairs."""
    return ((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2) ** 0.5


def pca_group_coherence(group, gene_info, pca_coords):
    """Compute PCA coherence metrics for a homeolog group.

    Returns:
        centroid (PC1, PC2): mean PCA position of group members.
        max_dist: maximum distance from centroid to any member.
        mean_dist: mean distance from centroid to members.
        n_with_coords: number of genes with PCA data.
    """
    points = []
    for name in group:
        info = gene_info.get(name)
        if info:
            base_id = info["base_id"]
            if base_id in pca_coords:
                points.append(pca_coords[base_id])

    if not points:
        return None, None, None, 0

    centroid = (
        sum(p[0] for p in points) / len(points),
        sum(p[1] for p in points) / len(points),
    )
    distances = [pca_euclidean_distance(p, centroid) for p in points]
    max_dist = max(distances) if distances else 0
    mean_dist = sum(distances) / len(distances) if distances else 0

    return centroid, max_dist, mean_dist, len(points)


def pca_inter_group_distance(group1, group2, gene_info, pca_coords):
    """Minimum PCA distance between any gene in group1 and group2."""
    min_d = float("inf")
    for n1 in group1:
        info1 = gene_info.get(n1)
        if not info1 or info1["base_id"] not in pca_coords:
            continue
        c1 = pca_coords[info1["base_id"]]
        for n2 in group2:
            info2 = gene_info.get(n2)
            if not info2 or info2["base_id"] not in pca_coords:
                continue
            c2 = pca_coords[info2["base_id"]]
            d = pca_euclidean_distance(c1, c2)
            if d < min_d:
                min_d = d
    return min_d


def estimate_pca_threshold(groups, gene_info, pca_coords):
    """Estimate PCA distance threshold from confident homeolog groups.

    Uses within-group PCA distances from groups with ≥2 subgenomes on
    a single chromosome. Returns 95th percentile × safety factor.
    """
    within_distances = []
    for group in groups:
        if len(group) < 2:
            continue
        chrs = set()
        sgs = set()
        for name in group:
            info = gene_info.get(name)
            if info:
                chrs.add(info["chr"])
                sgs.add(info["subgenome"])
        if len(chrs) == 1 and len(sgs) >= 2:
            for i, n1 in enumerate(group):
                info1 = gene_info.get(n1)
                if not info1 or info1["base_id"] not in pca_coords:
                    continue
                c1 = pca_coords[info1["base_id"]]
                for n2 in group[i + 1 :]:
                    info2 = gene_info.get(n2)
                    if not info2 or info2["base_id"] not in pca_coords:
                        continue
                    c2 = pca_coords[info2["base_id"]]
                    within_distances.append(pca_euclidean_distance(c1, c2))

    if not within_distances:
        return 5.0  # Default generous threshold
    within_distances.sort()
    idx_95 = min(int(0.95 * len(within_distances)), len(within_distances) - 1)
    return within_distances[idx_95] * DISTANCE_THRESHOLD_FACTOR


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# SUBFAMILY ASSIGNMENT
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def assign_subfamily_from_tree(tree, group, terminal_map):
    """Determine subfamily by finding the closest outgroup gene in the tree.

    For each gene in the group, compute the distance to every non-Fxa
    terminal. The closest non-Fxa gene's name is parsed for its subfamily
    annotation (e.g., MdPIP2_6 → PIP2).

    Uses up to 3 representative genes from the group and takes the
    consensus closest outgroup to be robust against long-branch artifacts.

    Returns (subfamily, family, closest_outgroup_name, distance) or Nones.
    """
    # Cache non-Fxa terminals for reuse
    outgroup_terminals = [t for t in tree.get_terminals() if not is_fxa_gene(t.name)]
    if not outgroup_terminals:
        return None, None, None, None

    # Use up to 3 representatives for robustness
    reps = group[:min(3, len(group))]
    subfamily_votes = defaultdict(int)
    best_overall = None
    best_dist = float("inf")

    for rep_name in reps:
        rep_terminal = terminal_map.get(rep_name)
        if not rep_terminal:
            continue

        min_dist = float("inf")
        closest = None
        for t in outgroup_terminals:
            d = tree.distance(rep_terminal, t)
            if d < min_dist:
                min_dist = d
                closest = t.name

        if closest:
            sf, fam = parse_outgroup_subfamily(closest)
            if sf:
                subfamily_votes[(sf, fam, closest)] += 1
            if min_dist < best_dist:
                best_dist = min_dist
                best_overall = closest

    # Pick the consensus subfamily (most votes)
    if subfamily_votes:
        (sf, fam, outg), _ = max(subfamily_votes.items(), key=lambda x: x[1])
        return sf, fam, outg, best_dist

    # Fallback: use closest outgroup even if subfamily parsing failed
    if best_overall:
        sf, fam = parse_outgroup_subfamily(best_overall)
        return sf, fam, best_overall, best_dist

    return None, None, None, None


def assign_all_sub_subfamilies(tree, named_groups, gene_info, terminal_map, group_subfamilies, gff3_coords=None):
    """Global sub-subfamily assignment using optimal 1:1 matching per subfamily.

    ALGORITHM:
        1. Collect all outgroup (reference) genes and parse their sub-subfamily
           annotations (e.g., MdPIP2_6 → PIP2;6).
        2. Group HGs by their assigned subfamily (e.g., all PIP2 HGs together).
        3. Detect inter-HG tandem duplicates using GFF3 coordinates: if genes
           from two HGs of the same subfamily are within 50kb on the same
           chromosome in ≥2 subgenomes, they are tandem duplicates and MUST
           share the same sub-subfamily (e.g., both are NIP1;1).
        4. For each subfamily, compute the average phylogenetic distance from
           each HG (using ALL members) to each reference sub-subfamily.
        5. Use greedy optimal assignment: assign the closest (HG, sub-subfamily)
           pair first, then remove both from the pool and repeat.
           Tandem clusters count as a single unit for assignment.
        6. When there are more HG groups than unique reference sub-subfamilies
           (i.e., Fragaria has more paralogs than the reference), create new
           isoform numbers (e.g., NIP1;3), ordered by phylogenetic distance.
        7. Species weighting: Md (Rosaceae, closest) and Fa weight 3x,
           At 2x, Hb 1x, Os 0.5x for distance calculations.

    Returns dict: group_id → (sub_subfamily, confidence, detail_string).
    """
    SPECIES_WEIGHT = {"Md": 3.0, "Fa": 3.0, "At": 2.0, "Hb": 1.0, "Os": 0.5}
    TANDEM_MAX_DISTANCE_BP = 50000  # 50kb threshold for tandem detection
    TANDEM_MIN_SUBGENOMES = 2  # Need adjacency in at least 2 subgenomes

    # Collect all reference genes by subfamily
    outgroup_terminals = [t for t in tree.get_terminals() if not is_fxa_gene(t.name)]
    ref_by_subfamily = defaultdict(list)  # subfamily → [(ref_name, sub_subfamily, terminal)]

    for t in outgroup_terminals:
        ssf, sf, fam, species = parse_outgroup_sub_subfamily(t.name)
        if ssf:
            ref_by_subfamily[sf].append((t.name, ssf, t, species))

    # Build lookup: base_id → terminal
    fxa_terminals_by_base = {}
    for t in tree.get_terminals():
        info = parse_fxa_info(t.name)
        if info:
            fxa_terminals_by_base[info["base_id"]] = t

    # Group HGs by subfamily
    hgs_by_subfamily = defaultdict(list)  # subfamily → [(group_id, group)]
    for group_id, group in named_groups:
        sf = group_subfamilies.get(group_id, "unknown")
        hgs_by_subfamily[sf].append((group_id, group))

    results = {}

    for sf, hg_list in hgs_by_subfamily.items():
        refs = ref_by_subfamily.get(sf, [])

        if not refs:
            # No references at all — assign sequentially
            for i, (gid, group) in enumerate(hg_list, 1):
                results[gid] = (f"{sf};{i}", "none", "no_reference")
            continue

        # ── Detect inter-HG tandem duplicates using GFF3 coordinates ──
        # Two HGs are tandem if their genes are within 50kb on the same
        # chromosome in at least 2 subgenomes.
        tandem_clusters = []  # list of sets of group_ids that are tandems
        tandem_assigned = set()

        if gff3_coords and len(hg_list) > 1:
            # Build per-subgenome coordinate map for each HG
            hg_coords = {}  # gid → {subgenome → (chr, min_start, max_end)}
            for gid, group in hg_list:
                sg_coords = {}
                for member_name in group:
                    info = gene_info.get(member_name)
                    if not info:
                        continue
                    base_id = info["base_id"]
                    sg = info["subgenome"]
                    if base_id in gff3_coords:
                        c = gff3_coords[base_id]
                        if sg not in sg_coords:
                            sg_coords[sg] = (c["chr"], c["start"], c["end"])
                        else:
                            old = sg_coords[sg]
                            sg_coords[sg] = (
                                old[0],
                                min(old[1], c["start"]),
                                max(old[2], c["end"]),
                            )
                hg_coords[gid] = sg_coords

            # Pairwise comparison of HGs within this subfamily
            hg_ids = [gid for gid, _ in hg_list]
            for i in range(len(hg_ids)):
                for j in range(i + 1, len(hg_ids)):
                    gid_a = hg_ids[i]
                    gid_b = hg_ids[j]
                    coords_a = hg_coords.get(gid_a, {})
                    coords_b = hg_coords.get(gid_b, {})
                    shared_sgs = set(coords_a.keys()) & set(coords_b.keys())
                    adjacent_count = 0
                    for sg in shared_sgs:
                        chr_a, start_a, end_a = coords_a[sg]
                        chr_b, start_b, end_b = coords_b[sg]
                        if chr_a != chr_b:
                            continue
                        # Distance between the two genes
                        gap = max(0, max(start_a, start_b) - min(end_a, end_b))
                        if gap <= TANDEM_MAX_DISTANCE_BP:
                            adjacent_count += 1
                    if adjacent_count >= TANDEM_MIN_SUBGENOMES:
                        print(f"    TANDEM: {gid_a} <-> {gid_b} ({sf}, adjacent in {adjacent_count} subgenomes)")
                        # Merge into clusters
                        found_cluster = None
                        for cluster in tandem_clusters:
                            if gid_a in cluster or gid_b in cluster:
                                cluster.add(gid_a)
                                cluster.add(gid_b)
                                found_cluster = cluster
                                break
                        if not found_cluster:
                            tandem_clusters.append({gid_a, gid_b})
                        tandem_assigned.add(gid_a)
                        tandem_assigned.add(gid_b)

        # Collect unique sub-subfamilies from references, averaging distances
        # per sub-subfamily (not per individual ref gene).
        unique_ssfs = sorted(set(r[1] for r in refs))

        # Compute weighted average distance: each HG → each unique sub-subfamily
        # Weight by species relevance (Md closest to Fragaria).
        dist_matrix = {}  # (group_id, ssf) → weighted_avg_distance

        for gid, group in hg_list:
            for ssf in unique_ssfs:
                # Find all reference genes with this ssf
                ssf_refs = [r for r in refs if r[1] == ssf]
                weighted_sum = 0.0
                weight_sum = 0.0
                for ref_name, _, ref_term, species in ssf_refs:
                    sp_weight = SPECIES_WEIGHT.get(species, 1.0)
                    # Average distance from ALL group members to this ref gene
                    member_dists = []
                    for member_name in group:
                        mt = terminal_map.get(member_name)
                        if mt:
                            d = tree.distance(mt, ref_term)
                            member_dists.append(d)
                    if member_dists:
                        avg_d = sum(member_dists) / len(member_dists)
                        weighted_sum += avg_d * sp_weight
                        weight_sum += sp_weight
                if weight_sum > 0:
                    dist_matrix[(gid, ssf)] = weighted_sum / weight_sum
                else:
                    dist_matrix[(gid, ssf)] = float("inf")

        # GREEDY OPTIMAL 1:1 ASSIGNMENT
        # Sort all (group_id, ssf) pairs by distance.
        # Tandem cluster members share one slot: use the cluster representative
        # (the one closest to any ref) for assignment, then propagate.

        # Build effective HG list: replace tandem cluster members with cluster reps
        cluster_rep = {}  # gid → representative gid
        for cluster in tandem_clusters:
            # The representative is the member closest to any reference
            best_gid = None
            best_dist = float("inf")
            for gid in cluster:
                for ssf in unique_ssfs:
                    d = dist_matrix.get((gid, ssf), float("inf"))
                    if d < best_dist:
                        best_dist = d
                        best_gid = gid
            for gid in cluster:
                cluster_rep[gid] = best_gid

        # Build the set of "effective" HGs for assignment (one per cluster)
        effective_gids = set()
        for gid, _ in hg_list:
            rep = cluster_rep.get(gid, gid)
            effective_gids.add(rep)

        all_pairs = sorted(dist_matrix.items(), key=lambda x: x[1])

        assigned_groups = set()
        assigned_ssfs = set()
        assignments = {}  # group_id → (ssf, distance)

        for (gid, ssf), dist in all_pairs:
            rep = cluster_rep.get(gid, gid)
            if rep not in effective_gids:
                continue  # not a representative
            if rep in assigned_groups:
                continue
            if ssf in assigned_ssfs:
                continue
            assignments[rep] = (ssf, dist)
            assigned_groups.add(rep)
            assigned_ssfs.add(ssf)

        # Handle unassigned effective groups (more HGs than unique sub-subfamilies)
        unassigned = [gid for gid in effective_gids if gid not in assigned_groups]
        if unassigned:
            # Find the max existing isoform number in this subfamily
            existing_nums = set()
            for ssf in unique_ssfs:
                # Parse the number after the semicolon
                parts = ssf.split(";")
                if len(parts) == 2 and parts[1].isdigit():
                    existing_nums.add(int(parts[1]))
            for ssf in assigned_ssfs:
                parts = ssf.split(";")
                if len(parts) == 2 and parts[1].isdigit():
                    existing_nums.add(int(parts[1]))

            next_num = max(existing_nums) + 1 if existing_nums else 1

            # Sort unassigned by their minimum distance to any reference
            # (to assign them in phylogenetic order)
            unassigned_with_dist = []
            for gid in unassigned:
                min_d = min(dist_matrix.get((gid, ssf), float("inf")) for ssf in unique_ssfs)
                unassigned_with_dist.append((min_d, gid))
            unassigned_with_dist.sort()

            for _, gid in unassigned_with_dist:
                new_ssf = f"{sf};{next_num}"
                # Find closest already-assigned ssf to determine the new isoform
                closest_ssf_dist = float("inf")
                closest_ssf = None
                for ssf in unique_ssfs:
                    d = dist_matrix.get((gid, ssf), float("inf"))
                    if d < closest_ssf_dist:
                        closest_ssf_dist = d
                        closest_ssf = ssf
                assignments[gid] = (new_ssf, closest_ssf_dist)
                next_num += 1

        # Propagate assignments to all tandem cluster members
        for cluster in tandem_clusters:
            rep = None
            for gid in cluster:
                if gid in assignments:
                    rep = gid
                    break
            if rep:
                for gid in cluster:
                    if gid != rep:
                        assignments[gid] = assignments[rep]

        # Build results with confidence
        for gid, group in hg_list:
            if gid in assignments:
                ssf, dist = assignments[gid]

                # Confidence based on how clearly it separates from alternatives
                dists_to_all = sorted(
                    (dist_matrix.get((gid, s), float("inf")), s)
                    for s in unique_ssfs
                )
                if len(dists_to_all) >= 2:
                    best_d = dists_to_all[0][0]
                    second_d = dists_to_all[1][0]
                    if second_d > 0 and (second_d - best_d) / second_d > 0.15:
                        confidence = "high"
                    elif second_d > 0 and (second_d - best_d) / second_d > 0.05:
                        confidence = "medium"
                    else:
                        confidence = "low"
                else:
                    confidence = "high"

                # Build detail string
                detail_parts = []
                for d, s in dists_to_all[:3]:
                    detail_parts.append(f"{s}={d:.4f}")
                detail_str = "; ".join(detail_parts)

                results[gid] = (ssf, confidence, detail_str)
            else:
                results[gid] = (f"{sf};?", "none", "unassigned")

    return results


def assign_subfamily_from_metadata(group, gene_info, metadata):
    """Determine subfamily by majority vote from metadata table.

    Only considers genes with reliable annotations (GFF3 or EXONERATE
    sources, not MAKER_GFF3 which may be misclassified).

    Returns (subfamily, family) or (None, None) if no consensus.
    """
    votes = defaultdict(int)
    for name in group:
        info = gene_info.get(name)
        if not info:
            continue
        base_id = info["base_id"]
        meta = metadata.get(base_id)
        if not meta:
            continue
        # Only use reliable sources for voting
        source = meta.get("annotation_source", "")
        if source in ("MAKER_GFF3",):
            continue
        subfam = meta.get("original_subfamily", "").replace("Fa", "")
        fam = meta.get("original_family", "")
        if subfam and fam:
            votes[(subfam, fam)] += 1

    if not votes:
        return None, None

    (best_subfam, best_fam), count = max(votes.items(), key=lambda x: x[1])
    if count >= MIN_GENES_FOR_MAJORITY_VOTE or count == len(votes):
        return best_subfam, best_fam
    return None, None


def assign_subfamily(tree, group, gene_info, metadata, terminal_map):
    """Hybrid subfamily assignment: tree-based primary, metadata validation.

    Returns (subfamily, family, method, details).
    """
    # Method 1: closest outgroup neighbor in tree
    tree_subfam, tree_fam, outgroup_name, dist = assign_subfamily_from_tree(
        tree, group, terminal_map
    )

    # Method 2: metadata majority vote
    meta_subfam, meta_fam = assign_subfamily_from_metadata(
        group, gene_info, metadata
    )

    # Decision logic
    if tree_subfam and meta_subfam:
        if tree_subfam == meta_subfam:
            return tree_subfam, tree_fam, "tree+metadata", outgroup_name
        else:
            # Prefer metadata when they disagree (metadata has expert curation)
            # But log the discrepancy
            return (
                meta_subfam, meta_fam, "metadata (tree_disagrees)",
                f"tree={tree_subfam} via {outgroup_name}"
            )
    elif tree_subfam:
        return tree_subfam, tree_fam, "tree_only", outgroup_name
    elif meta_subfam:
        return meta_subfam, meta_fam, "metadata_only", ""
    else:
        return "unknown", "unknown", "none", ""


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# METADATA & ANNOTATION
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def load_metadata(filepath):
    """Load gene metadata from tabla_Aquaporinas_traduccion.tabular."""
    metadata = {}
    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = row["gene_id"].strip()
            if gid:
                metadata[gid] = {
                    "annotation_source": row.get("fuente_seq", "").strip(),
                    "veredicto": row.get("veredicto", "").strip(),
                    "original_subfamily": row.get("aqp_family_subfamily", "").strip(),
                    "original_family": row.get("subfamilia_phylo", "").strip(),
                    "confianza": row.get("confianza", "").strip(),
                    "longitud_aa": row.get("longitud_aa", "").strip(),
                    "TMHs": row.get("TMHs", "").strip(),
                }
    return metadata


def determine_annotation_type(base_id, metadata):
    """Classify annotation source for a gene."""
    if base_id not in metadata:
        return "unknown"
    m = metadata[base_id]
    source = m["annotation_source"]
    veredicto = m["veredicto"]
    if source == "EXONERATE" or veredicto == "EXONERATE":
        return "exonerate"
    elif source == "MAKER_GFF3":
        return "maker_partial"
    elif source == "GFF3_FALLBACK":
        return "gff3_fallback"
    elif source == "GFF3":
        return "gff3_original"
    return source.lower() if source else "unknown"


def determine_is_partial(base_id, metadata):
    """Determine if a gene model is partial/incomplete."""
    if base_id not in metadata:
        return False
    m = metadata[base_id]
    source = m["annotation_source"]
    veredicto = m["veredicto"]
    tmhs = int(m["TMHs"]) if m["TMHs"].isdigit() else 0
    length = int(m["longitud_aa"]) if m["longitud_aa"].isdigit() else 0
    if source == "MAKER_GFF3":
        return True
    if source == "GFF3_FALLBACK" and veredicto == "AMBAS_MAL":
        return True
    if tmhs < 6 and length < 230:
        return True
    return False


def load_gff3_coordinates(gff3_path):
    """Load gene start/end coordinates from the consensus GFF3.

    Returns dict: base_gene_id → {chr, start, end, strand}.
    """
    coords = {}
    if not os.path.exists(gff3_path):
        return coords
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            m = FXA_GENE_PATTERN.search(parts[8])
            if m:
                base_id = f"Fxa{m.group(1)}{m.group(2)}g{m.group(3)}"
                coords[base_id] = {
                    "chr": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                }
    return coords


def detect_overlapping_fragments(group, gene_info, metadata, gff3_coords):
    """Detect two types of problematic annotations within a homeolog group:

    1. overlapping_fragment: a gene whose genomic coordinates overlap with
       another gene on the same chromosome+subgenome. This indicates the
       gene finder annotated the same locus twice.

    2. putative_pseudogene: a gene with fewer than 6 TMHs. Functional
       aquaporins require all 6 transmembrane helices to form the water
       channel pore. Genes with <6 TMHs are truncated and non-functional.

    Returns dict: tree_name → (flag, overlaps_with_base_id_or_empty)
    """
    flags = {}  # tree_name → (flag_type, overlaps_with)

    # Group genes by subgenome for overlap checks
    by_subgenome = defaultdict(list)
    for name in group:
        info = gene_info.get(name)
        if info:
            by_subgenome[info["subgenome"]].append(name)

    # Check for genomic overlaps within each subgenome
    for sg, names in by_subgenome.items():
        if len(names) < 2:
            continue

        # Build coordinate list
        genes_with_coords = []
        for name in names:
            info = gene_info.get(name)
            if not info:
                continue
            base_id = info["base_id"]
            if base_id in gff3_coords:
                c = gff3_coords[base_id]
                md = metadata.get(base_id, {})
                tmhs = int(md.get("TMHs", "0") or "0")
                aa_len = int(md.get("longitud_aa", "0") or "0")
                genes_with_coords.append({
                    "name": name,
                    "base_id": base_id,
                    "start": c["start"],
                    "end": c["end"],
                    "tmhs": tmhs,
                    "aa_len": aa_len,
                })

        # Sort by start position
        genes_with_coords.sort(key=lambda x: x["start"])

        # Check pairwise overlaps
        for i in range(len(genes_with_coords)):
            for j in range(i + 1, len(genes_with_coords)):
                g1 = genes_with_coords[i]
                g2 = genes_with_coords[j]
                # Do they overlap? (g2.start <= g1.end)
                if g2["start"] <= g1["end"]:
                    # Mark the smaller/worse one as the fragment
                    if g1["tmhs"] >= g2["tmhs"] and g1["aa_len"] >= g2["aa_len"]:
                        flags[g2["name"]] = ("overlapping_fragment", g1["base_id"])
                    else:
                        flags[g1["name"]] = ("overlapping_fragment", g2["base_id"])

    # Check for putative pseudogenes (not already flagged as overlapping)
    for name in group:
        if name in flags:
            continue
        info = gene_info.get(name)
        if not info:
            continue
        base_id = info["base_id"]
        md = metadata.get(base_id, {})
        tmhs = int(md.get("TMHs", "0") or "0")
        if tmhs < 6:
            flags[name] = ("putative_pseudogene", "")

    return flags


def detect_tandem_duplicates(group, gene_info):
    """Detect tandem duplicates: multiple genes from the same subgenome.

    Within a homeolog group, each subgenome should contribute at most one
    gene. Additional genes from the same subgenome indicate tandem (or
    local) duplications that occurred independently of the polyploidization.

    Returns a set of tree_names that are tandem duplicates (all copies
    from a duplicated subgenome are flagged, not just the 'extra' ones).
    """
    sg_counts = defaultdict(list)
    for name in group:
        info = gene_info.get(name)
        if info:
            sg_counts[info["subgenome"]].append(name)

    tandem_names = set()
    for sg, names in sg_counts.items():
        if len(names) > 1:
            for name in names:
                tandem_names.add(name)

    return tandem_names


def detect_reclassifications(group, group_subfamily, group_family, gene_info, metadata):
    """Detect genes whose metadata classification disagrees with tree-derived assignment.

    Only evaluates MAKER_GFF3 genes, which were classified before the
    phylogenetic tree was available and may have incorrect assignments.

    Two categories:
      - Family-level reclassification: the aquaporin family changed
        (e.g., TIP → NIP). These represent genuine classification errors.
      - Subfamily clarification: the family is correct but the subfamily
        is now specified (e.g., NIP → NIP1). These represent new
        information from the tree, not errors.

    Returns list of (base_id, old_fam, old_subfam, new_subfam, reclass_type).
    """
    reclassifications = []
    for name in group:
        info = gene_info.get(name)
        if not info:
            continue
        base_id = info["base_id"]
        meta = metadata.get(base_id)
        if not meta:
            continue
        if meta["annotation_source"] != "MAKER_GFF3":
            continue

        old_subfam = meta["original_subfamily"].replace("Fa", "")
        old_fam = meta["original_family"]

        # Check family-level disagreement
        if old_fam and old_fam != group_family:
            reclassifications.append(
                (base_id, old_fam, old_subfam, group_subfamily, "family_change")
            )
        # Check subfamily clarification (same family, more specific now)
        elif old_subfam and old_subfam != group_subfamily:
            if old_subfam == old_fam:
                # Original had only family-level annotation (e.g., "NIP")
                reclassifications.append(
                    (base_id, old_fam, old_subfam, group_subfamily, "subfamily_clarification")
                )

    return reclassifications


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GROUP NAMING & SORTING
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def assign_group_ids(groups, gene_info):
    """Assign systematic group IDs sorted by chromosome number and position.

    Groups are sorted by:
        1. Chromosome number (1-7)
        2. Minimum gene number within the group (genomic position proxy)

    Regular groups: HG-01, HG-02, ...
    Singleton/doublet groups with incomplete subgenome representation
    that appear to be orphans: HG-S1, HG-S2, ...
    """
    def sort_key(group):
        chrs = []
        gene_nums = []
        for name in group:
            info = gene_info.get(name)
            if info:
                chrs.append(info["chr"])
                gene_nums.append(info["gene_num"])
        if not chrs:
            return (99, 99999)
        # Use most common chromosome, then minimum gene number
        chr_counts = defaultdict(int)
        for c in chrs:
            chr_counts[c] += 1
        dom_chr = max(chr_counts, key=chr_counts.get)
        min_num = min(gene_nums) if gene_nums else 99999
        return (dom_chr, min_num)

    sorted_groups = sorted(groups, key=sort_key)

    # Separate regular groups (≥3 genes or ≥2 subgenomes) from orphans
    regular = []
    orphans = []
    for group in sorted_groups:
        subgenomes = set()
        for name in group:
            info = gene_info.get(name)
            if info:
                subgenomes.add(info["subgenome"])
        if len(group) >= 3 or len(subgenomes) >= 2:
            regular.append(group)
        else:
            orphans.append(group)

    # Assign IDs
    named_groups = []
    for i, group in enumerate(regular, 1):
        named_groups.append((f"HG-{i:02d}", group))
    for i, group in enumerate(orphans, 1):
        named_groups.append((f"HG-S{i}", group))

    return named_groups


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# OUTPUT
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def determine_completeness(subgenomes):
    """Classify group completeness based on subgenome representation."""
    n = len(subgenomes)
    if n == 4:
        return "complete"
    elif n == 3:
        return "triplet"
    elif n == 2:
        return "doublet"
    else:
        return "singleton"


def generate_output_tables(
    named_groups, gene_info, metadata, tree, terminal_map, pca_coords, output_dir,
    gff3_coords=None, output_suffix=""
):
    """Generate homeolog_groups.tsv and homeolog_groups_summary.tsv."""

    fieldnames = [
        "homeolog_group", "gene_id", "chromosome", "chr_number", "subgenome",
        "family", "subfamily", "sub_subfamily", "sub_subfamily_confidence",
        "annotation_source", "is_partial",
        "is_tandem_duplicate", "fragment_flag", "overlaps_with",
        "group_completeness", "group_size",
        "n_subgenomes", "subgenomes_present", "was_reclassified",
        "old_classification",
    ]

    # ── PASS 1: Assign subfamilies to all groups ──
    group_subfamilies = {}   # group_id → subfamily
    group_families = {}      # group_id → family
    group_methods = {}       # group_id → (method, detail)
    for group_id, group in named_groups:
        subfamily, family, method, detail = assign_subfamily(
            tree, group, gene_info, metadata, terminal_map
        )
        group_subfamilies[group_id] = subfamily
        group_families[group_id] = family
        group_methods[group_id] = (method, detail)

    # ── PASS 2: Global sub-subfamily assignment (detects inter-HG tandems) ──
    print("  [Sub-subfamily] Running global optimal assignment...")
    ssf_results = assign_all_sub_subfamilies(
        tree, named_groups, gene_info, terminal_map, group_subfamilies,
        gff3_coords=gff3_coords
    )
    # Print assignments
    for gid in sorted(ssf_results.keys()):
        ssf, conf, detail = ssf_results[gid]
        print(f"    {gid}: {group_subfamilies[gid]} → {ssf} ({conf})")

    # ── PASS 3: Generate output rows ──
    all_rows = []
    summary_rows = []
    all_reclassifications = []

    for group_id, group in named_groups:
        subfamily = group_subfamilies[group_id]
        family = group_families[group_id]
        method, detail = group_methods[group_id]

        # Sub-subfamily from global assignment
        sub_subfamily, ssf_confidence, ssf_detail = ssf_results.get(
            group_id, (f"{subfamily};?", "none", "")
        )

        # Tandem duplicates
        tandem_set = detect_tandem_duplicates(group, gene_info)

        # Fragment / pseudogene detection
        if gff3_coords:
            fragment_flags = detect_overlapping_fragments(
                group, gene_info, metadata, gff3_coords
            )
        else:
            fragment_flags = {}

        # Reclassifications
        reclass = detect_reclassifications(group, subfamily, family, gene_info, metadata)
        all_reclassifications.extend(reclass)
        reclass_ids = {r[0] for r in reclass}

        # Subgenome analysis
        subgenomes = sorted(set(
            gene_info[n]["subgenome"] for n in group if n in gene_info
        ))
        completeness = determine_completeness(subgenomes)

        # Per-gene rows
        group_rows = []
        for name in sorted(group, key=lambda n: (
            gene_info.get(n, {}).get("subgenome", "Z"),
            gene_info.get(n, {}).get("gene_num", 99999),
        )):
            info = gene_info.get(name)
            if not info:
                continue
            base_id = info["base_id"]

            # Check reclassification (only flag family-level changes in TSV)
            was_reclass = False
            old_class = ""
            for r in reclass:
                if r[0] == base_id:
                    if r[4] == "family_change":
                        was_reclass = True
                        old_class = f"{r[1]}({r[2]})"
                    break

            # Fragment flag
            frag_flag = ""
            overlaps_with = ""
            if name in fragment_flags:
                frag_flag, overlaps_with = fragment_flags[name]

            row = {
                "homeolog_group": group_id,
                "gene_id": base_id,
                "chromosome": f"chr_{info['chr']}{info['subgenome']}",
                "chr_number": info["chr"],
                "subgenome": info["subgenome"],
                "family": family,
                "subfamily": subfamily,
                "sub_subfamily": sub_subfamily,
                "sub_subfamily_confidence": ssf_confidence,
                "annotation_source": determine_annotation_type(base_id, metadata),
                "is_partial": "yes" if determine_is_partial(base_id, metadata) else "no",
                "is_tandem_duplicate": "yes" if name in tandem_set else "no",
                "fragment_flag": frag_flag,
                "overlaps_with": overlaps_with,
                "group_completeness": completeness,
                "group_size": len(group),
                "n_subgenomes": len(subgenomes),
                "subgenomes_present": ",".join(subgenomes),
                "was_reclassified": "yes" if was_reclass else "no",
                "old_classification": old_class,
            }
            group_rows.append(row)
            all_rows.append(row)

        # Summary row
        tandem_count = sum(1 for r in group_rows if r["is_tandem_duplicate"] == "yes")
        partial_count = sum(1 for r in group_rows if r["is_partial"] == "yes")
        exonerate_count = sum(
            1 for r in group_rows if r["annotation_source"] == "exonerate"
        )
        # PCA coherence
        pca_centroid, pca_max_dist, pca_mean_dist, pca_n = pca_group_coherence(
            group, gene_info, pca_coords
        )

        summary_rows.append({
            "homeolog_group": group_id,
            "family": family,
            "subfamily": subfamily,
            "sub_subfamily": sub_subfamily,
            "sub_subfamily_confidence": ssf_confidence,
            "chr_number": group_rows[0]["chr_number"] if group_rows else "",
            "n_genes": len(group_rows),
            "n_subgenomes": len(subgenomes),
            "subgenomes": ",".join(subgenomes),
            "completeness": completeness,
            "n_tandem_duplicates": tandem_count,
            "n_partial": partial_count,
            "n_exonerate": exonerate_count,
            "gene_ids": ",".join(r["gene_id"] for r in group_rows),
            "subfamily_method": method,
            "sub_subfamily_detail": ssf_detail,
            "pca_max_dist": f"{pca_max_dist:.4f}" if pca_max_dist is not None else "NA",
            "pca_mean_dist": f"{pca_mean_dist:.4f}" if pca_mean_dist is not None else "NA",
        })

    # Write main table
    out_path = os.path.join(output_dir, f"homeolog_groups{output_suffix}.tsv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    # Write summary table
    summary_path = os.path.join(output_dir, f"homeolog_groups_summary{output_suffix}.tsv")
    summary_fields = list(summary_rows[0].keys()) if summary_rows else []
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    return all_rows, summary_rows, all_reclassifications


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# MAIN
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def main():
    # ── Argument parsing ──
    parser = argparse.ArgumentParser(
        description="Homeolog grouping for Fragaria × ananassa aquaporins."
    )
    parser.add_argument(
        "--mode", choices=["clean", "full"], default="clean",
        help=(
            "clean (default): Use pruned tree (arbol_acuaporinas.treefile) "
            "and cleaned tabular (121 complete AQPs). "
            "full: Use original tree with partials (old_tree_with_partials/fxa_final.treefile) "
            "and backup tabular (144 AQPs including partials/MAKER). "
            "Output files are suffixed with '_full' in full mode."
        )
    )
    args = parser.parse_args()
    mode = args.mode

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # ── Select input files based on mode ──
    if mode == "clean":
        tree_path = os.path.join(script_dir, "arbol_acuaporinas.treefile")
        tabular_path = os.path.join(script_dir, "tabla_Aquaporinas_traduccion.tabular")
        output_suffix = ""
        mode_label = "CLEAN (no partials, no MAKER)"
    else:  # full
        tree_path = os.path.join(script_dir, "old_tree_with_partials", "fxa_final.treefile")
        tabular_path = os.path.join(script_dir, "tabla_Aquaporinas_traduccion.tabular.backup")
        output_suffix = "_full"
        mode_label = "FULL (includes partials + MAKER sequences)"

    pca_path = os.path.join(script_dir, "PCA_Coordenadas_Finales.csv")
    gff3_path = os.path.join(script_dir, "consenso_aqp_fixed.gff3")

    for path, label in [(tree_path, "Tree"), (tabular_path, "Metadata")]:
        if not os.path.exists(path):
            sys.exit(f"ERROR: {label} file not found: {path}")

    print("=" * 70)
    print(f"Homeolog Grouping — Algorithmic Identification [{mode.upper()} mode]")
    print(f"  Mode: {mode_label}")
    print(f"  Tree: {os.path.basename(tree_path)}")
    print(f"  Table: {os.path.basename(tabular_path)}")
    print("Fragaria × ananassa Aquaporins (IQ-TREE3, Q.PLANT+R7)")
    print("=" * 70)

    # ── Step 1: Parse tree ────────────────────────────────────────────
    print("\n[Step 1] Parsing phylogenetic tree...")
    tree = Phylo.read(tree_path, "newick")
    n_terminals = tree.count_terminals()
    fxa_terminals = [t for t in tree.get_terminals() if is_fxa_gene(t.name)]
    print(f"  Total terminals: {n_terminals}")
    print(f"  Fxa terminals:   {len(fxa_terminals)}")

    # Build gene info map (tree_name → metadata)
    gene_info = {}
    tree_to_base = {}
    for t in fxa_terminals:
        info = parse_fxa_info(t.name)
        if info:
            gene_info[t.name] = info
            tree_to_base[t.name] = info["base_id"]

    # ── Step 2: Root tree ─────────────────────────────────────────────
    print("\n[Step 2] Rooting tree on outgroup...")
    rooting_method = root_tree(tree, ROOT_OUTGROUP_PREFIX)
    print(f"  Method: {rooting_method}")

    # Terminal lookup map (needed for distance calculations)
    terminal_map = {t.name: t for t in tree.get_terminals()}

    # ── Step 3: Load metadata ─────────────────────────────────────────
    print("\n[Step 3] Loading gene metadata...")
    metadata = load_metadata(tabular_path)
    print(f"  Genes in metadata: {len(metadata)}")

    # ── Step 3b: Load PCA coordinates ─────────────────────────────────
    print("\n[Step 3b] Loading PCA coordinates...")
    pca_coords = load_pca_coordinates(pca_path)
    if pca_coords:
        print(f"  Genes with PCA coords: {len(pca_coords)}")
    else:
        print("  WARNING: PCA file not found; skipping PCA validation")

    # ── Step 3c: Load GFF3 coordinates ────────────────────────────────
    print("\n[Step 3c] Loading GFF3 coordinates...")
    gff3_coords = load_gff3_coordinates(gff3_path)
    if gff3_coords:
        print(f"  Genes with GFF3 coords: {len(gff3_coords)}")
    else:
        print("  WARNING: GFF3 file not found; skipping overlap detection")

    # ── Phase 1: Find maximal Fxa-only clades ─────────────────────────
    print("\n[Phase 1] Finding maximal Fxa-only clades...")
    clades, singletons = find_maximal_fxa_clades(tree)
    n_in_clades = sum(len(get_fxa_leaf_names(c)) for c in clades)
    print(f"  Multi-gene clades:  {len(clades)} ({n_in_clades} genes)")
    print(f"  Singletons:         {len(singletons)} genes")

    # Check for multi-chromosome clades
    multi_chr = 0
    for c in clades:
        names = get_fxa_leaf_names(c)
        chrs = set(gene_info[n]["chr"] for n in names if n in gene_info)
        if len(chrs) > 1:
            multi_chr += 1
    print(f"  Multi-chromosome:   {multi_chr} (will be split in Phase 2)")

    # ── Phase 2: Chromosome-based splitting ───────────────────────────
    print("\n[Phase 2] Splitting multi-chromosome clades...")
    groups, n_splits = phase2_split(clades, singletons, gene_info)
    print(f"  Clades split:     {n_splits}")
    print(f"  Groups after P2:  {len(groups)}")

    # Verify all genes assigned
    all_assigned = set()
    for g in groups:
        for n in g:
            base = gene_info.get(n, {}).get("base_id", n)
            all_assigned.add(base)
    print(f"  Genes assigned:   {len(all_assigned)} / {len(fxa_terminals)}")

    # ── Phase 3: Distance-aware merging ───────────────────────────────
    print("\n[Phase 3] Merging fragments separated by outgroup insertions...")
    groups, merge_log = phase3_merge(tree, groups, gene_info, terminal_map, pca_coords)
    for entry in merge_log:
        print(f"  {entry}")
    print(f"  Groups after P3:  {len(groups)}")

    # ── Assign group IDs ──────────────────────────────────────────────
    print("\n[Assigning group IDs]...")
    named_groups = assign_group_ids(groups, gene_info)
    regular = sum(1 for gid, _ in named_groups if not gid.startswith("HG-S"))
    orphans = sum(1 for gid, _ in named_groups if gid.startswith("HG-S"))
    print(f"  Regular groups:   {regular}")
    print(f"  Orphan groups:    {orphans}")

    # ── Generate output ───────────────────────────────────────────────
    print("\n[Generating output tables]...")
    all_rows, summary_rows, reclassifications = generate_output_tables(
        named_groups, gene_info, metadata, tree, terminal_map, pca_coords, script_dir,
        gff3_coords=gff3_coords, output_suffix=output_suffix
    )

    # ── Summary statistics ────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("RESULTS SUMMARY")
    print(f"{'=' * 70}")
    print(f"  Total genes:      {len(all_rows)}")
    print(f"  Total groups:     {len(named_groups)}")

    # Fragment/pseudogene summary
    n_overlapping = sum(1 for r in all_rows if r.get("fragment_flag") == "overlapping_fragment")
    n_pseudogene = sum(1 for r in all_rows if r.get("fragment_flag") == "putative_pseudogene")
    if n_overlapping or n_pseudogene:
        print(f"\n  Flagged annotations:")
        if n_overlapping:
            print(f"    Overlapping fragments:  {n_overlapping}")
            for r in all_rows:
                if r.get("fragment_flag") == "overlapping_fragment":
                    print(f"      {r['gene_id']} overlaps with {r['overlaps_with']}")
        if n_pseudogene:
            print(f"    Putative pseudogenes:   {n_pseudogene}")
            for r in all_rows:
                if r.get("fragment_flag") == "putative_pseudogene":
                    print(f"      {r['gene_id']}")

    completeness_counts = defaultdict(int)
    for s in summary_rows:
        completeness_counts[s["completeness"]] += 1
    for cat in ["complete", "triplet", "doublet", "singleton"]:
        print(f"  {cat.capitalize():16s} {completeness_counts.get(cat, 0)}")

    if reclassifications:
        family_changes = [r for r in reclassifications if r[4] == "family_change"]
        clarifications = [r for r in reclassifications if r[4] == "subfamily_clarification"]
        if family_changes:
            print(f"\n  Family-level reclassifications ({len(family_changes)} genes):")
            for base_id, old_fam, old_subfam, new_subfam, _ in family_changes:
                print(f"    {base_id}: {old_fam}({old_subfam}) -> {new_subfam}")
        if clarifications:
            print(f"\n  Subfamily clarifications ({len(clarifications)} genes):")
            for base_id, old_fam, old_subfam, new_subfam, _ in clarifications:
                print(f"    {base_id}: {old_fam} -> {new_subfam}")

    print(f"\n  Output: homeolog_groups{output_suffix}.tsv ({len(all_rows)} rows)")
    print(f"  Output: homeolog_groups_summary{output_suffix}.tsv ({len(summary_rows)} rows)")

    # ── Algorithm report ──────────────────────────────────────────────
    report_path = os.path.join(script_dir, f"homeolog_algorithm_report{output_suffix}.txt")
    with open(report_path, "w") as f:
        f.write("HOMEOLOG GROUPING — ALGORITHM REPORT\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Tree file:       {os.path.basename(tree_path)}\n")
        f.write(f"Metadata file:   {os.path.basename(tabular_path)}\n")
        f.write(f"Rooting method:  {rooting_method}\n")
        f.write(f"Total terminals: {n_terminals}\n")
        f.write(f"Fxa terminals:   {len(fxa_terminals)}\n\n")

        f.write("CONFIGURATION\n")
        f.write(f"  ROOT_OUTGROUP_PREFIX:       {ROOT_OUTGROUP_PREFIX}\n")
        f.write(f"  MAX_OUTGROUP_INSERTIONS:    {MAX_OUTGROUP_INSERTIONS}\n")
        f.write(f"  DISTANCE_THRESHOLD_FACTOR:  {DISTANCE_THRESHOLD_FACTOR}\n")
        f.write(f"  MIN_GENES_FOR_MAJORITY_VOTE: {MIN_GENES_FOR_MAJORITY_VOTE}\n\n")

        f.write("PHASE 1: MAXIMAL FXA-ONLY CLADES\n")
        f.write(f"  Multi-gene clades: {len(clades)}\n")
        f.write(f"  Singletons:        {len(singletons)}\n")
        if singletons:
            f.write("  Singleton genes:\n")
            for s in singletons:
                info = gene_info.get(s.name, {})
                f.write(f"    {info.get('base_id', s.name)} (chr{info.get('chr','?')}{info.get('subgenome','?')})\n")
        f.write(f"  Multi-chromosome clades: {multi_chr}\n\n")

        f.write("PHASE 2: CHROMOSOME SPLITTING\n")
        f.write(f"  Clades split:    {n_splits}\n")
        f.write(f"  Groups after P2: {len(groups) + sum(1 for ml in merge_log if 'MERGE' in ml)}\n\n")

        f.write("PHASE 3: DISTANCE-AWARE MERGING\n")
        for entry in merge_log:
            f.write(f"  {entry}\n")
        f.write(f"  Groups after P3: {len(groups)}\n\n")

        f.write("FINAL GROUPS\n")
        for gid, group in named_groups:
            info_list = [gene_info.get(n) for n in group if n in gene_info]
            base_ids = [i["base_id"] for i in info_list if i]
            chrs = set(i["chr"] for i in info_list if i)
            sgs = sorted(set(i["subgenome"] for i in info_list if i))
            f.write(f"  {gid}: chr={chrs} sg={sgs} genes={','.join(base_ids)}\n")

        if reclassifications:
            family_changes = [r for r in reclassifications if r[4] == "family_change"]
            clarifications = [r for r in reclassifications if r[4] == "subfamily_clarification"]
            if family_changes:
                f.write(f"\nFAMILY-LEVEL RECLASSIFICATIONS ({len(family_changes)} genes)\n")
                for base_id, old_fam, old_subfam, new_subfam, _ in family_changes:
                    f.write(f"  {base_id}: {old_fam}({old_subfam}) -> {new_subfam}\n")
            if clarifications:
                f.write(f"\nSUBFAMILY CLARIFICATIONS ({len(clarifications)} genes)\n")
                for base_id, old_fam, old_subfam, new_subfam, _ in clarifications:
                    f.write(f"  {base_id}: {old_fam} -> {new_subfam}\n")

        # Gene coverage check
        meta_ids = set(metadata.keys())
        tree_base_ids = set(tree_to_base.values())
        missing_from_tree = meta_ids - tree_base_ids
        if missing_from_tree:
            f.write(f"\nGENES IN METADATA BUT MISSING FROM TREE ({len(missing_from_tree)})\n")
            for gid in sorted(missing_from_tree):
                f.write(f"  {gid}\n")

    print(f"  Output: homeolog_algorithm_report.txt")

    # Print gene coverage warning
    meta_ids = set(metadata.keys())
    tree_base_ids = set(tree_to_base.values())
    missing = meta_ids - tree_base_ids
    if missing:
        print(f"\n  WARNING: {len(missing)} genes in metadata but absent from tree:")
        for gid in sorted(missing):
            print(f"    {gid}")

    print("\nDone.")


if __name__ == "__main__":
    main()
