import csv


def parse_row(row):
    # chrom, start, end, strand, attr, feat
    return row[0], row[3], row[4], row[6], row[8], row[2]


def check_overlap(gmeta_data, tmeta_data):
    return


def extract_attrs(attr_str, attr_name):
    # Find star and ending index
    start_idx = attr_str.find(attr_name)

    end_idx = attr_str[start_idx:].find(";")

    # Slice string to the attribute we want
    value = attr_str[start_idx:end_idx]

    value = value.replace(f"{attr_name}=", "")

    return value


def is_overlap(chrom1="", chrom2="",end1=0, start1=0,end2=0, start2=0):
    if chrom1 != chrom2:
        return False
        # if strand1 != strand2:
        #    return False

    cov = 0
    len1 = abs(end1 - start1) + 1
    len2 = abs(end2 - start2) + 1

    overlap_within = start1 >= start2 and end1 <= end2
    overlap_within2 = start1 <= start2 and end1 >= end2

    # check overlap - within a transcript
    if overlap_within:
        cov = ((end1 - start1) / len1) * 100

    if overlap_within2:
        cov = ((end2 - start2) / len2) * 100

    if start1 < start2 and (end1 >= start2 and end1 < end2):
        # overlap left
        cov = ((end1 - start2) / len1) * 100
    if start1 > start2 and (start1 < end2 and end1 > end2):
        # overlap right
        cov = ((end2 - start1) / len1) * 100
    if chrom1 == chrom2 and cov >= 60:
        return True

    return False


def parse_cords(mdata1, mdata2):
    chrom1 = mdata1['chrom']
    chrom2 = mdata2['chrom']
    end1 = mdata1['end']
    start1 = mdata1['start']
    end2 = mdata2['end']
    start2 = mdata2['start']

    return chrom1, chrom2, end1, end2, start1, start2


def create_merged_gff(gstore, gmeta_data, tstore, tmeta_data):

    for tkey, tvals in tstore.items():

        # Get the transcript meta data
        trans_meta = tmeta_data.get(tkey)

        # Keep track of genes overlapping with this transcript.
        overlap_store = []

        for gene, gene_meta in gmeta_data.items():
            # Parse relevant coordinates from meta data
            ch1, ch2, en1, en2, st1, st2 = parse_cords(mdata1=trans_meta, mdata2=gene_meta)
            # Check for overlap
            if is_overlap(chrom1=ch1, chrom2=ch2, end1=en1, start1=st1, end2=en2, start2=st2):
                # Append gene to list of overlapping genes
                overlap_store.append(gene)

            # CASE 1 : # transcript overlaps with more than one gene - remove transcript.
            if len(overlap_store) > 1:
                continue

        print(chrom, attr)
        1 / 0
        pass

    return


def parse_file(gff, target_feat="gene"):
    ann_store = dict()
    meta_data = dict()
    stream = csv.reader(open(gff), delimiter="\t")
    feat_id = None

    for row in stream:

        chrom, start, end, strand, attr, feat = parse_row(row)

        if chrom.startswith("#"):
            continue

        if feat == target_feat:
            feat_id = extract_attrs(attr, "ID")
            mdata = dict(chrom=chrom, start=start, end=end, attr=attr, strand=strand)
            meta_data[feat_id] = mdata

        ann_store.setdefault(feat_id, []).append(row)

    return ann_store, meta_data


def main():
    aug_gff = "test_aug.gff3"
    aug_gff_store = parse_file(aug_gff, "gene")

    #
    # print(len(aug_gff_store))
    # for k, val in aug_gff_store.items():
    #     for v in val:
    #         v = list(map(str, v))
    #         print("\t".join(v))

    return


if __name__ == "__main__":
    main()
