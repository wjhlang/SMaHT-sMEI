import pysam

hap1_bam = "Hap1.verkko.sorted.bam"
hap2_bam = "Hap2.verkko.sorted.bam"
window = 100
OUTPUT_FILE="output_hg38_to_dsa.txt"
INPUT_FILE="input_hg38.txt"


def map_window_ref_to_query(aln, win_ref_start, win_ref_end):
    """
    For a given pysam AlignedSegment, and 0-based inclusive window on the reference,
    return: [contig_start_1based, contig_end_1based] = minimum and maximum contig/query positions
    that map to the window.

    Returns (contig_start, contig_end) or (None, None) if not fully covered.
    """
    ref_pos = aln.reference_start  # On reference
    query_pos = 0  # On query/contig

    cigar = aln.cigartuples
    n_M = 0   # op code 0 (M)
    n_eq = 0  # op code 7 (=)
    n_X = 0   # op code 8 (X)
    contig_positions = []
    contig_length = aln.query_length  # Total length of the contig

    # For reverse: need to know how many query bases are used within aligned region
    if aln.is_reverse:
        consumed_query_bases = 0

        for op, length in cigar:
            if op in (0, 7, 8):  # M, =, X
                for i in range(length):
                    if win_ref_start <= ref_pos <= win_ref_end:
                        # Compute position from RIGHT (reverse)
                        query_1b = contig_length - consumed_query_bases
                        contig_positions.append(query_1b)
                        if op == 0:
                            n_M += 1
                        elif op == 7:
                            n_eq += 1
                        elif op == 8:
                            n_X += 1
                    ref_pos += 1
                    consumed_query_bases += 1
            elif op in (1, 4):  # I, S (insertion, soft clip): consume query only
                consumed_query_bases += length
                query_pos += length
            elif op in (2, 3):  # D, N : consume reference only
                ref_pos += length
            elif op == 5:  # H: hard clip, do nothing; doesn't consume anything
                continue

    else:
        for op, length in cigar:
            if op in (0, 7, 8):  # M, =, X
                for i in range(length):
                    if win_ref_start <= ref_pos <= win_ref_end:
                        query_1b = query_pos + 1
                        contig_positions.append(query_1b)
                        if op == 0:
                            n_M += 1
                        elif op == 7:
                            n_eq += 1
                        elif op == 8:
                            n_X += 1
                    ref_pos += 1
                    query_pos += 1
            elif op == 1:   # I
                query_pos += length
            elif op == 2:   # D
                ref_pos += length
            elif op == 3:   # N
                ref_pos += length
            elif op == 4:   # S
                query_pos += length
            elif op == 5:   # H
                continue

    if not contig_positions:
        return (None, None, 0, 0, 0)

    return (min(contig_positions), max(contig_positions), n_M, n_eq, n_X)



def grch38_region_to_contig_hits(bam_path, chrom, pos_start_1based, pos_end_1based, debug=False):
    """
    For a given window, return list of tuples:
    (contig, contig_start, contig_end, strand, flag, mapq, contig_mean, [optional DEBUG])
    """
    hits = []
    bam = pysam.AlignmentFile(bam_path, "rb")
    pos_start_0based = pos_start_1based - 1
    pos_end_0based = pos_end_1based - 1

    for aln in bam.fetch(chrom, pos_start_0based, pos_end_0based+1):
        print(f"\n--- Alignment being considered ---")
        print(f"Query name: {aln.query_name}")
        print(f"Reference start: {aln.reference_start}")
        print(f"Reference end: {aln.reference_end}")
        print(f"Mapping quality: {aln.mapping_quality}")
        #print(f"CIGAR: {aln.cigarstring}")
        print(f"Is reverse: {aln.is_reverse}")
        
        
        
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        strand = '-' if aln.is_reverse else '+'

        # Only accept alignments that FULLY cover the region
        if not (ref_start <= pos_start_0based and ref_end > pos_end_0based):
            print(">> Alignment does not fully span the window, skipping.")
            continue
        if aln.mapping_quality != 60:
            print(">> Alignment mapping quality != 60, skipping.")
            continue

        contig = aln.query_name
        flag = aln.flag
        mapq = aln.mapping_quality
        cigar = aln.cigartuples

        # Get left soft-clip
        left_softclip = 0
        if cigar and cigar[0][0] == 4:
            left_softclip = cigar[0][1]

        # Get right soft-clip
        right_softclip = 0
        if cigar and cigar[-1][0] == 4:
            right_softclip = cigar[-1][1]

        print(">> Alignment passes filters. Mapping reference window to query...")
        contig_start, contig_end, n_M, n_eq, n_X = map_window_ref_to_query(aln, pos_start_0based, pos_end_0based)
        print(f"Mapped contig_start: {contig_start}, contig_end: {contig_end}")


        debug_info = {}
        if debug:
            debug_info = {
                'qname': contig,
                'ref_start': ref_start,
                'ref_end': ref_end,
                'query_start': aln.query_alignment_start + 1, # 1-based
                'query_end': aln.query_alignment_end,
                
                'left_softclip': left_softclip,
                'right_softclip': right_softclip,

                'window_ref_start': pos_start_0based,
                'window_ref_end': pos_end_0based,
                'contig_start': contig_start,
                'contig_end': contig_end,
                'strand': strand,
                'n_M': n_M,
                'n_eq': n_eq,
                'n_X': n_X,
            }
           # print(f"DEBUG: {debug_info}")

        if contig_start is not None and contig_end is not None:
            contig_mean = (contig_start + contig_end) / 2
            hits.append((contig, contig_start, contig_end, strand, flag, mapq, contig_mean, n_M, n_eq, n_X, debug_info))
            print(">> Alignment added to HIT list.")
        else:
            print(">> Alignment does not map the window fully on query, skipping.")
    bam.close()
    return hits

# ======= Main processing =======


with open(OUTPUT_FILE, "w") as fout, open(INPUT_FILE) as fin:
    fout.write(f"# Window size: ±{window}bp (region: center±{window}bp)\n")
    fout.write("chrom\tpos\tlabel\tgrch38_start\tgrch38_end\thap\tn_contig_hits\tcontig\tcontig_start\tcontig_end\tcontig_mean\tstrand\tflag\tmapq\tn_M\tn_eq\tn_X\n")
    for line in fin:
        if not line.strip():
            continue
        chrom, pos, label = line.strip().split()[:3]
        pos = int(pos)
        grch38_start = pos - window
        grch38_end = pos + window

        hits1 = grch38_region_to_contig_hits(hap1_bam, chrom, grch38_start, grch38_end, debug=True)
        hits2 = grch38_region_to_contig_hits(hap2_bam, chrom, grch38_start, grch38_end, debug=True)

        # Hap1 output
        n1 = len(hits1)
        if n1:
            for h in hits1:
                fout.write(f"{chrom}\t{pos}\t{label}\t{grch38_start}\t{grch38_end}\thap1\t{n1}\t{h[0]}\t{h[1]}\t{h[2]}\t{h[6]}\t{h[3]}\t{h[4]}\t{h[5]}\t{h[7]}\t{h[8]}\t{h[9]}\n")

        else:
            fout.write(f"{chrom}\t{pos}\t{label}\t{grch38_start}\t{grch38_end}\thap1\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")

        # Hap2 output
        n2 = len(hits2)
        if n2:
            for h in hits2:
                fout.write(f"{chrom}\t{pos}\t{label}\t{grch38_start}\t{grch38_end}\thap2\t{n1}\t{h[0]}\t{h[1]}\t{h[2]}\t{h[6]}\t{h[3]}\t{h[4]}\t{h[5]}\t{h[7]}\t{h[8]}\t{h[9]}\n")
        else:
            fout.write(f"{chrom}\t{pos}\t{label}\t{grch38_start}\t{grch38_end}\thap2\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
