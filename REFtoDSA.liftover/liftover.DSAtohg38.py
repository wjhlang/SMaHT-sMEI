#!/usr/bin/env python3
import pysam
from collections import defaultdict

# ======================
# Parameters
# ======================
hap1_bam = "Hap1.verkko.sorted.bam"
hap2_bam = "Hap2.verkko.sorted.bam"
window = 100
OUTPUT_FILE = "output_dsa_to_hg38.txt"
INPUT_FILE = "input_dsa.txt"
allow_partial = False


# ======================
# Map query/contig window to reference (strand-aware)
# ======================
def map_window_query_to_ref(aln, win_query_start, win_query_end):
    ref_pos = aln.reference_start
    query_pos = 0
    cigar = aln.cigartuples
    contig_length = aln.query_length
    ref_positions = []
    n_M, n_eq, n_X = 0, 0, 0

    if aln.is_reverse:
        consumed_query_bases = 0
        for op, length in cigar:
            if op in (0, 7, 8):
                for _ in range(length):
                    query_1b = contig_length - consumed_query_bases
                    if win_query_start <= query_1b <= win_query_end:
                        ref_positions.append(ref_pos + 1)
                        if op == 0: n_M += 1
                        elif op == 7: n_eq += 1
                        elif op == 8: n_X += 1
                    ref_pos += 1
                    consumed_query_bases += 1
            elif op in (1, 4):
                consumed_query_bases += length
                query_pos += length
            elif op in (2, 3):
                ref_pos += length
            elif op == 5:
                continue
    else:
        for op, length in cigar:
            if op in (0, 7, 8):
                for _ in range(length):
                    query_1b = query_pos + 1
                    if win_query_start <= query_1b <= win_query_end:
                        ref_positions.append(ref_pos + 1)
                        if op == 0: n_M += 1
                        elif op == 7: n_eq += 1
                        elif op == 8: n_X += 1
                    ref_pos += 1
                    query_pos += 1
            elif op == 1:
                query_pos += length
            elif op == 2:
                ref_pos += length
            elif op == 3:
                ref_pos += length
            elif op == 4:
                query_pos += length
            elif op == 5:
                continue

    if not ref_positions:
        return (None, None, 0, 0, 0)

    return (min(ref_positions), max(ref_positions), n_M, n_eq, n_X)


# ======================
# Process one BAM in a single pass
# ======================
def process_bam_single_pass(bam_path, site_dict, fout, hap):
    """
    site_dict: {contig_name: [(pos, label), ...]}
    Returns set of (contig, pos) that were hit.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    n_aln = 0
    n_hit = 0
    hit_set = set()  # track (contig, pos) pairs that got at least one hit

    for aln in bam.fetch(until_eof=True):
        n_aln += 1
        if n_aln % 100000 == 0:
            print(f"  [{hap}] {n_aln} alignments scanned, {n_hit} hits so far...")

        if aln.mapping_quality != 60:
            continue

        contig = aln.query_name
        if contig not in site_dict:
            continue

        if aln.is_reverse:
            q_start_1b = aln.query_length - aln.query_alignment_end + 1
            q_end_1b   = aln.query_length - aln.query_alignment_start
        else:
            q_start_1b = aln.query_alignment_start + 1
            q_end_1b   = aln.query_alignment_end

        for pos, label in site_dict[contig]:
            win_start = max(1, pos - window)
            win_end   = pos + window

            if not allow_partial:
                if not (q_start_1b <= win_start and q_end_1b >= win_end):
                    continue
            else:
                if q_end_1b < win_start or q_start_1b > win_end:
                    continue

            ref_start, ref_end, n_M, n_eq, n_X = map_window_query_to_ref(aln, win_start, win_end)
            if ref_start is None or ref_end is None:
                continue

            chrom    = aln.reference_name
            strand   = '-' if aln.is_reverse else '+'
            flag     = aln.flag
            mapq     = aln.mapping_quality
            ref_mean = (ref_start + ref_end) // 2
            n_hit   += 1
            hit_set.add((contig, pos))

            fout.write(f"{label}\t{contig}\t{pos}\t{chrom}\t{ref_start}\t{ref_end}\t{ref_mean}\t{strand}\t{flag}\t{mapq}\t{n_M}\t{n_eq}\t{n_X}\t{hap}\n")

    bam.close()
    print(f"  [{hap}] Done. {n_aln} alignments scanned, {n_hit} hits written.")

    # Write NA for any site that never got a hit
    n_miss = 0
    for contig, sites in site_dict.items():
        for pos, label in sites:
            if (contig, pos) not in hit_set:
                fout.write(f"{label}\t{contig}\t{pos}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{hap}\n")
                n_miss += 1
    print(f"  [{hap}] {n_miss} sites with no hits written as NA.")

    return hit_set


# ======================
# Main
# ======================

# Load input sites grouped by contig and haplotype
hap1_sites = defaultdict(list)
hap2_sites = defaultdict(list)
na_lines   = []

with open(INPUT_FILE) as fin:
    for line in fin:
        if not line.strip():
            continue
        parts  = line.strip().split()
        contig = parts[0]
        pos    = int(parts[1])
        label  = parts[2] if len(parts) > 2 else "."

        if contig.startswith("haplotype1-"):
            hap1_sites[contig].append((pos, label))
        elif contig.startswith("haplotype2-"):
            hap2_sites[contig].append((pos, label))
        else:
            na_lines.append((contig, pos, label))

print(f"Loaded {sum(len(v) for v in hap1_sites.values())} hap1 sites across {len(hap1_sites)} contigs")
print(f"Loaded {sum(len(v) for v in hap2_sites.values())} hap2 sites across {len(hap2_sites)} contigs")

with open(OUTPUT_FILE, "w") as fout:
    fout.write(f"# Window size: ±{window}bp (center±{window}bp)\n")
    fout.write("label\tcontig\tpos\tgrch38_chr\tgrch38_start\tgrch38_end\tgrch38_mean\tstrand\tflag\tmapq\tn_M\tn_eq\tn_X\thap\n")

    # Write unknown haplotype lines as NA
    for contig, pos, label in na_lines:
        fout.write(f"{label}\t{contig}\t{pos}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tunknown\n")

    print("\nProcessing Hap1 BAM...")
    process_bam_single_pass(hap1_bam, hap1_sites, fout, "hap1")

    print("\nProcessing Hap2 BAM...")
    process_bam_single_pass(hap2_bam, hap2_sites, fout, "hap2")

print(f"\nDone. Output: {OUTPUT_FILE}")
