# This script merges stringtie annotated with augustus annoatated gff3.
# A transcript that overlaps with more than one gene is removed.
#

import sys, csv
from collections import defaultdict


def extract_attr(attr_str, attr_name):
    attrs = attr_str.split(";")
    feat_id = [val for val in attrs if val.startswith(attr_name)][0]
    feat_id = feat_id.replace(attr_name, "")
    return feat_id


def extract_tid(attr):
    tid = attr.split(";")[1].split(" ")
    tid = list(filter(lambda x: x != '', tid))[1]
    tid = tid.replace("\"", "")
    return tid

def read_gff(gff, feat="gene"):
    """
    Reads augustus annotated gff3 file.
    """

    ann_store = defaultdict(list)
    stream = csv.reader(open(gff), delimiter="\t")

    for row in stream:
        if row[0].startswith("#"):
            continue
        if row[2] == feat:
            chrom, start, end, strand, attr = row[0], row[3], row[4], row[6], row[8]
            feat_id = extract_attr(attr, "ID=")
            key = ":".join([chrom, start, end, strand, feat_id])
            ann_store[key].append(row)
            continue
        ann_store[key].append(row)

    return ann_store


def read_ann(fname):
    """
    Stores augustus annotated genes
    """
    ann_store = defaultdict(list)
    stream = csv.DictReader(open(fname), delimiter="\t")
    for row in stream:
        key = ":".join([row['chrom'],row['start'],row['end'], row['strand']])
        ann_store[key] = row
    return ann_store


def store_recp_map(fname, species):
    """
    Reads a tab-delimited idmap file and returns a dictionary of list.
    Eg: giraffe_recp_protein_map_stringtie.txt
    """
    trans_store = defaultdict(dict)
    stream = csv.DictReader(open(fname), delimiter="\t")
    for row in stream:
        protein_id = species + "_pid"
        pid = row.get(protein_id)
        trans_store[pid] = row
    return trans_store


def modify_gene_attr(gid, gname, source, qcovs, pident):
    if gname is None or gname == "":
        new_attr = f'ID={gid};gene_id={gid};gene_name={gid};Name={gid};best_match=stringtie'
    else:
        new_attr = f'ID={gid};gene_id={gid};Name={gname};gene_name={gname};best_match={source};query_coverage={qcovs};percentage_identity={pident}'

    return new_attr


def modify_transcript_attr(tid, gid, gname, source, qcovs, pident, evidence):
    transcript_no = ".".join(tid.split(".")[-2:])
    parent = gid

    if evidence == "stringtie":
        new_attr = f'ID={tid};Parent={parent};transcript_id={tid};gene_name={parent};Name={tid};best_match={evidence}'
    elif evidence == "locus":
        trans_id = ".".join([gname, transcript_no])
        # trans_id = ".".join([gname, tid])
        new_attr = f'ID={tid};Parent={parent};transcript_id={trans_id};gene_name={gname};Name={trans_id};best_match={evidence}'
    else:
        trans_id = ".".join([gname, transcript_no])
        # trans_id = ".".join([gname, tid])
        new_attr = f'ID={tid};Parent={parent};transcript_id={trans_id};gene_name={gname};Name={trans_id};best_match={source};query_coverage={qcovs};percentage_identity={pident}'
    return new_attr


def modify_attr(attr, tid, gname):
    parent = tid
    transcript_no = ".".join(tid.split(".")[-2:])

    if gname is None or gname == "":
        new_attr = f'Parent={parent};transcript_id={parent}'
    else:
        trans_id = ".".join([gname, transcript_no])
        # trans_id = ".".join([gname, tid])
        new_attr = f'Parent={parent};transcript_id={trans_id};gene_name={gname}'

    if attr.startswith("ID"):
        feat = attr.split(";")[0].split(".")[-1]
        id_ = ".".join([parent, feat])
        new_attr = ";".join([id_, new_attr])

    return new_attr


def is_overlap(chrom1,start1,end1, strand1, chrom2, start2,end2,strand2):
    """
    Checks coordinate overlap.
    """

    if chrom1 !=chrom2:
       return False
    cov  = 0
    #if strand1 == "-":
    #   start1,end1 = end1,start1
    #if strand2 == "-":
    #   start2, end2 = end2, start2

    len1 = abs(end1-start1) +1
    len2 = abs(end2-start2) + 1
    # check overlap - within a transcript
    if start1 >=start2 and end1 <=end2:
       cov = ((end1-start1)/len1) * 100
    if start1 <= start2 and end1 >= end2:
       cov = ((end2-start2)/len2) * 100
    if  start1 < start2 and (end1 >=start2 and end1 <end2):
       # overlap left
       cov = ((end1-start2)/len1)*100
    if start1 > start2 and (start1 <end2 and end1 >end2):
       # overlap right
       cov = ((end2-start1)/len1)*100
    if chrom1 == chrom2 and cov >=60:
        return True
    return False


def check_gene_overlap(store, transcript):
    """
    Check if a stringtie transcript overlap with an augustus predicted gene.
    """
    gene_overlap = list()

    tr_chrom, tr_start, tr_end, tr_strand, tr_id = transcript.split(":")
    tr_start, tr_end = int(tr_start), int(tr_end)

    for gene in store.keys():
        gene_chrom, gene_start, gene_end, gene_strand, gene_id = gene.split(":")
        gene_start, gene_end = int(gene_start), int(gene_end)

        if is_overlap(tr_chrom, tr_start, tr_end, tr_strand, gene_chrom, gene_start, gene_end, gene_strand):
            gene_overlap.append(gene)
    return gene_overlap


def modify_current_ann(vals ,tid, gid, gname, source, qcovs, pident, evidence):

    for item in vals:
        if item[2] == "gene":
            gene_attr = modify_gene_attr(gid, gname, source, qcovs, pident)
            item[8] = gene_attr
        elif item[2] == "transcript":
            tran_attr = modify_transcript_attr(tid, gid, gname, source, qcovs, pident, evidence)
            item[8] = tran_attr
        else:
            new_attr = modify_attr(item[8], tid, gname)
            item[8] = new_attr
    return vals


def annotate_existing(attr, best_match, existing_vals):
    tid = extract_attr(attr, "ID=")
    gid = ".".join(tid.split(".")[:-1])
    gname = extract_attr(attr, "gene_name=")
    source = best_match
    evidence = get_evidence(source)
    pident, qcovs = get_evidence_vals(evidence, attr)
    # annotate augustus gene
    modified = modify_current_ann(existing_vals, tid, gid, gname, source, qcovs, pident,
                                       evidence)
    return modified


def get_evidence(source):
    if source != "stringtie" and source != "locus":
        evidence = "species"
    else:
        evidence = source
    return evidence


def get_evidence_vals(evidence, attrs):
    """
    returns pident, qcovs
    """
    pident, qcovs = 0, 0
    if evidence == "species":
        pident = extract_attr(attrs, "percentage_identity=")
        qcovs = extract_attr(attrs, "query_coverage=")
    return pident, qcovs


def get_merged_keys(merged):

    store = list()
    for item in merged:
        if item[2] != "gene":
            continue
        gid = extract_attr(item[8], "ID=")
        chrom, start, end, strand = item[0], item[3], item[4], item[6]
        key = ":".join([chrom, start,end,strand,gid])
        store.append(key)
    return store


def get_unchanged_augustus():
    pass



def created_merged_gff(aug_store, transcript_store):

    #transcript_store_corrected = transcript_gff_store
    new_trans_store = list()
    new_gene_store = defaultdict(list)
    merged_store = list()
    seen = set()

    i =0
    for transcript_key, transcript_vals in transcript_store.items():

        i+=1
        attrs = transcript_vals[0][8]
        tr_best_match = extract_attr(attrs, "best_match=")
        print("iter",transcript_key)

        #  check transcript overlap with augustus annotated gene
        overlap_store = check_gene_overlap(aug_store, transcript_key)

        # CASE 1 : # transcript overlaps with more than one gene - remove transcript.
        if len(overlap_store) > 1:
            # del transcript_store_corrected[transcript_key]
            continue

        print("AUG store")
        print(aug_store)
        print("merged")
        print(merged_store)


        # Overlap with an augustus gene
        if len(overlap_store) == 1:
            print("OVERLAP")
            #print(overlap_store)
            gene_key = overlap_store[0]
            gene_attr = aug_store[gene_key][0][8]
            gene_best_match = extract_attr(gene_attr, "best_match=")

            # CASE2 : transcript is annotated and gene is annotated - nothing to do
            if tr_best_match != "stringtie" and gene_best_match != "augustus":
                print("OVERLAP with annotated", overlap_store)
                if gene_key not in seen:
                    merged_store.extend(aug_store[gene_key])
                    seen.add(gene_key)
                merged_store.extend(transcript_store[transcript_key])
                continue

            # CASE3 : transcript is unannotated and gene is unannotated - nothong to do
            if tr_best_match == "stringtie" and gene_best_match == "augustus":
                print("OVERLAP with unannotated", overlap_store)
                print(f"iteration{i}")
                #print(transcript_key, gene_key)
                #print(seen)
                if gene_key not in seen:
                    merged_store.extend(aug_store[gene_key])
                    seen.add(gene_key)
                merged_store.extend(transcript_store[transcript_key])
                print("case 3: aug", aug_store)
                continue

            # CASE 4 : transcript is annotated gene is unannotated - annotate gene with transcript ann
            if tr_best_match != "stringtie"  and gene_best_match == "augustus":
                print("OVERLAP with unannotated gene", overlap_store)
                print(f"iteration{i}")
                #print(transcript_key, gene_key)
                #print(seen)

                print(" case 4: aug", aug_store)
                unann_vals = aug_store[gene_key]
                print(type(aug_store))
                print(type(unann_vals))

                modified_gene  =  annotate_existing(attrs, tr_best_match, unann_vals)
                print(type(modified_gene))
                #1/0
                print("modified gene", modified_gene)
                print("aug", aug_store)
                #1/0

                if gene_key not in seen:
                    merged_store.extend(modified_gene)
                    seen.add(gene_key)
                merged_store.extend(transcript_vals)
                print(f"iteration{i}")
                print("aug", aug_store)
                #print("mmmm",merged_store)
                #print("----")
                continue

            # CASE 5 : transcript is unannotated and gene is annotated - annotate transcript with gene ann
            if tr_best_match == "stringtie" and gene_best_match != "augustus":
                print("OVERLAP with unannotated transcript", overlap_store)
                #print(transcript_key, gene_key)
                #print(seen)
                modified_transcript =  annotate_existing(gene_attr, gene_best_match, transcript_vals)
                if gene_key not in seen:
                    merged_store.extend(aug_store[gene_key])
                    seen.add(gene_key)
                merged_store.extend(modified_transcript)
                continue

        # No overlap - add gene row to the largest transcript
        if len(overlap_store) == 0:
            print("NO overlap")
            print(transcript_key)
            # # check if it is overlaps with any other transcript seen so far.
            # # specify gene length as the length of the largest transcript.
            tr_chrom, tr_start, tr_end, tr_strand, tr_id = transcript_key.split(":")
            tr_start, tr_end = int(tr_start), int(tr_end)
            gchrom, gstart, gend, gstrand = tr_chrom, tr_start, tr_end, tr_strand
            if new_trans_store:
                for item in new_trans_store:
                    chrom, start, end, strand, tid = item.split(":")
                    start, end = int(start), int(end)

                    if is_overlap(chrom, start, end, strand, tr_chrom, tr_start, tr_end, tr_strand):
                        seen_length = end- start + 1
                        curr_length = tr_end - tr_start + 1

                        gstart, gend, gstrand = (start,end, strand ) if seen_length > curr_length else (tr_start,tr_end, strand)
                    break


            best_match =extract_attr(attrs, "best_match=")
            source = best_match
            tid = extract_attr(attrs, "ID=")
            gid = ".".join(tid.split(".")[:-1])
            gname = extract_attr(attrs, "gene_name=")
            evidence = get_evidence(best_match)
            pident, qcovs = get_evidence_vals(evidence, attrs)
            gene_attr = modify_gene_attr(gid, gname, source, qcovs, pident)
            gene_row = [gchrom, "Stringtie", "gene" , gstart, gend, ".", gstrand, gene_attr]

            #print("###",gene_row)

            if gname not in new_gene_store:
                new_gene_store[gname].append(gene_row)
            new_gene_store[gname].extend(transcript_vals)
            #merged_store.extend(gene_row)

            new_trans_store.append(transcript_key)
            continue

    return merged_store, new_gene_store



if __name__ == "__main__":
    aug_gff = "aug_ann_test.gff3" #"genes_ann_any.gff3"
    stringtie_gff = "stringtie_ann_test.gff3" #"stringtie_uniq_ann.gff3"
    #stringtie_ann = "giraffe_recp_protein_map_stringtie.txt"
    #augustus_ann = "augustus_annotated_any.txt"
    species = "giraffe"

    aug_gff_store = read_gff(aug_gff, "gene")
    transcript_gff_store = read_gff(stringtie_gff, "transcript")

    #aug_genes = read_ann(augustus_ann)
    #recp_map = store_recp_map(stringtie_ann, species)

    merged, added = created_merged_gff(aug_gff_store, transcript_gff_store)
    # Get augustus genes that did not have any intersection with stringtie, ie, augustus uniq

    #print(len(merged))
    #print(len(added))

    current_keys = list()
    # get merged
    if merged:
        merged_keys = get_merged_keys(merged)
        #print("***",merged_keys)
        current_keys.extend(merged_keys)
    if added:
        new_keys =added.keys()
        #print(new_keys)
        current_keys.extend(new_keys)
    #print(merged_keys)
    #print("current",current_keys)


    if not current_keys:
        current_keys = aug_gff_store.keys()

    unchanged = [key for key in aug_gff_store if key not in current_keys]

    #print(unchanged)

    # print unchanged
    for key in unchanged:
        for item in aug_gff_store[key]:
            item = list(map(str, item))
            print("\t".join(item))

    # print merged:
    for item in merged:
        item = list(map(str, item))
        print("\t".join(item))

    # print additional from stringtie
    for gene, childeren in added.items():
        for child in childeren:
            child= list(map(str, child))
            print("\t".join(child))
