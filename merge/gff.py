import csv


TEST_AUGUSTUS = "aug.gff3"
TEST_STRINGTIE = "stringtie.gff3"


def parse_row(row):
    # chrom, start, end, strand, attr, feat
    return row[0], row[3], row[4], row[6], row[8], row[2]


def parse_attrs(attr_str, attr_name):
    # Find start and ending index of the attribute we want
    if not attr_str:
        return ""
    start_idx = attr_str.find(attr_name)
    end_idx = attr_str[start_idx:].find(";")

    # Slice string to the attribute we want
    value = attr_str[start_idx:end_idx]
    value = value.replace(f"{attr_name}=", "")

    return value


def parse_file(gff, target="gene"):
    ann_store = dict()
    meta_data = dict()
    stream = csv.reader(open(gff), delimiter="\t")
    key = None
    for row in stream:
        if row.startswith("#"):
            continue

        chrom, start, end, strand, attr, feat = parse_row(row)
        if feat == target:
            key = parse_attrs(attr, "ID")
            # Fill meta data one for each key
            meta_data[key] = dict(chrom=chrom, start=start, end=end, attr=attr, strand=strand)
        if key:
            ann_store.setdefault(key, []).append(row)

    return ann_store, meta_data


def calculate_coverage(tmeta, gmeta):

    ch1, ch2 = tmeta['chrom'], gmeta['chrom']
    en1, en2 = tmeta['end'], gmeta['end']
    st1, st2 = tmeta['start'], gmeta['start']

    cov = 0
    len1 = abs(en1 - st1) + 1
    len2 = abs(en2 - st2) + 1

    within = st1 >= st2 and en1 <= en2
    within2 = st1 <= st2 and en1 >= en2
    # Parameters to evaluate left and right transcript coverage.
    left = st1 < st2 and en1 >= st2 and en1 < en2
    right = st1 > st2 and st2 < en2 and en1 > en2

    if within:
        cov = ((en1 - st1) / len1) * 100
    if within2:
        cov = ((en2 - st2) / len2) * 100
    if left:
        cov = ((en1 - st2) / len1) * 100
    if right:
        cov = ((en2 - st1) / len1) * 100

    return cov, ch1, ch2


def is_overlap(tmeta, gmeta):
    # Parse relevant coordinates from meta data
    coverage, ch1, ch2 = calculate_coverage(tmeta=tmeta, gmeta=gmeta)

    if ch1 != ch2:
        return False
        # if strand1 != strand2:
        #    return False

    if ch1 == ch2 and coverage >= 60:
        return True

    return False


def check_overlap(gmeta, tmeta):
    overlap_store = []
    for gene, gene_meta in gmeta.items():
        # Append gene to list of overlapping genes
        if is_overlap(tmeta=tmeta, gmeta=gmeta):
            overlap_store.append(gene)
    return overlap_store


def get_count(rows, target_str, target_idx):

    # Filter for rows that have the target_str in the index ( target_idx )
    target_lst = list(filter(lambda l: l[target_idx] == target_str, rows))
    count = len(target_lst)
    return count


def handle_gene_overlap(trans_meta, gmeta_data, gstore, count, overlap_store, prev_gene=""):
    # Get the transcript counts for this gene
    gene_overlap = len(overlap_store) == 1
    g_key = overlap_store[0] if gene_overlap else ""
    gene_meta = gmeta_data.get(g_key, {})
    t_best = parse_attrs(trans_meta['attr'], "best_match")
    g_best = parse_attrs(gene_meta.get('attr', ""), "best_match")

    g_best = parse_attrs(gene_meta.get('attr', ""), "best_match")

    if t_best != "stringtie" and g_best != "augustus":
        return


    return


def create_merged_gff(gene_file, transcript_file):
    gstore, gmeta_data = parse_file(gene_file, "gene")
    tstore, tmeta_data = parse_file(transcript_file, "transcript")
    prev_gene, count = "", 0

    for tkey, tvals in tstore.items():

        # Get the transcript meta data
        trans_meta = tmeta_data.get(tkey)
        trans_strand = trans_meta['strand']

        # Keep track of genes overlapping with this transcript.
        overlap_store = check_overlap(gmeta_data, trans_meta)
        is_overlapped = len(overlap_store) > 1
        is_gene_overlap = len(overlap_store) == 1
        g_key = overlap_store[0] if is_gene_overlap else ""
        gene_meta = gmeta_data.get(g_key, {})
        gene_strand = gene_meta.get('strand', '')
        gene = gstore[g_key]
        # CASE 1 : # transcript overlaps with more than one gene - remove transcript.
        if is_overlapped:
            continue

        # Discard if the overlapped gene and transcript do not have the same orientation
        if is_gene_overlap and gene_strand != trans_strand:
            continue

        if is_gene_overlap:

            tcount = get_count(rows=gene, target_str='transcript', target_idx=2)
            count = count + 1 if g_key == prev_gene else 1
            tcount += count
            prev_gene = g_key

            handle_gene_overlap(prev_gene=prev_gene)

        1 / 0
        pass

    return


def main():

    create_merged_gff(gene_file=TEST_AUGUSTUS, transcript_file=TEST_STRINGTIE)

    return


if __name__ == "__main__":
    main()