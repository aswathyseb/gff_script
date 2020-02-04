import csv


def parse_row(row):
    # chrom, start, end, strand, attr, feat
    return row[0], row[3], row[4], row[6], row[8], row[2]


def extract_attrs(attr_str, attr_name):
    # Find star and ending index
    start_idx = attr_str.find(attr_name)

    end_idx = attr_str[start_idx:].find(";")

    # Slice string to the attribute we want
    value = attr_str[start_idx:end_idx]

    value = value.replace(attr_name, "")

    return value


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
            feat_id = extract_attrs(attr, "ID=")
            meta_data[feat_id] = dict(chrom=chrom, start=start, end=end, attr=attr, strand=strand)
            ann_store.setdefault(feat_id, []).append(row)
            continue

        if feat_id:
            ann_store[feat_id].append(row)
        #print(key)
        1 / 0


def main():
    aug_gff = "test_aug.gff3"

    return


if __name__ == "__main__":
    main()


