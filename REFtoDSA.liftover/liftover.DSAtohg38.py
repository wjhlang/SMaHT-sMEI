#!/usr/bin/env python3
import pysam

# ======================
# Parameters
# ======================
hap1_bam = "Hap1.chr10.bam"
hap2_bam = "Hap2.chr10.bam"
window = 100
OUTPUT_FILE = "output_dsa_to_hg38.txt"
INPUT_FILE = "input_dsa.txt"
allow_partial = False  # Set True to allow partial coverage of ±window


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
# Get hits from BAM for DSA contig + pos
# ======================
def dsa_region_to_grch38_hits(bam_path, contig, pos_1based, debug=False):
    hits = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    win_start = max(1, pos_1based - window)
    win_end   = pos_1based + window

    for aln in bam.fetch(until_eof=True):
        if aln.query_name != contig:
            continue
        # Include primary, secondary, and supplementary
        if aln.mapping_quality != 60:
            continue

        # Strand-aware contig alignment start/end
        if aln.is_reverse:
            q_start_1b = aln.query_length - aln.query_alignment_end + 1
            q_end_1b   = aln.query_length - aln.query_alignment_start
        else:
            q_start_1b = aln.query_alignment_start + 1
            q_end_1b   = aln.query_alignment_end

        # Check full coverage (or partial)
        if not allow_partial:
            if not (q_start_1b <= win_start and q_end_1b >= win_end):
                continue
        else:
            # accept if any overlap with ±window
            if q_end_1b < win_start or q_start_1b > win_end:
                continue

        ref_start, ref_end, n_M, n_eq, n_X = map_window_query_to_ref(aln, win_start, win_end)
        if ref_start is None or ref_end is None:
            continue

        chrom  = aln.reference_name
        strand = '-' if aln.is_reverse else '+'
        flag   = aln.flag
        mapq   = aln.mapping_quality
        ref_mean = (ref_start + ref_end) // 2

        if debug:
            print(f"DEBUG: {contig}:{win_start}-{win_end} -> {chrom}:{ref_start}-{ref_end} ({strand})")

        hits.append((chrom, ref_start, ref_end, strand, flag, mapq, ref_mean, n_M, n_eq, n_X))
    bam.close()
    return hits


# ======================
# Main
# ======================
with open(OUTPUT_FILE, "w") as fout, open(INPUT_FILE) as fin:
    fout.write(f"# Window size: ±{window}bp (center±{window}bp)\n")
    fout.write("contig\tpos\tgrch38_chr\tgrch38_start\tgrch38_end\tgrch38_mean\tstrand\tflag\tmapq\tn_M\tn_eq\tn_X\thap\n")

    for line in fin:
        if not line.strip():
            continue
        contig, pos = line.strip().split()[:2]
        pos = int(pos)

        if contig.startswith("haplotype1-"):
            hits = dsa_region_to_grch38_hits(hap1_bam, contig, pos)
            hap = "hap1"
        elif contig.startswith("haplotype2-"):
            hits = dsa_region_to_grch38_hits(hap2_bam, contig, pos)
            hap = "hap2"
        else:
            hits = []
            hap = "unknown"

        if hits:
            for chrom, rS, rE, strand, flag, mapq, rM, nM, nEq, nX in hits:
                fout.write(f"{contig}\t{pos}\t{chrom}\t{rS}\t{rE}\t{rM}\t{strand}\t{flag}\t{mapq}\t{nM}\t{nEq}\t{nX}\t{hap}\n")
        else:
            fout.write(f"{contig}\t{pos}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{hap}\n")

