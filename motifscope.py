import pandas as pd
import numpy as np
import os
import argparse
import sys
import pandas as pd
from Bio import SeqIO, pairwise2
#from Bio.Seq import Seq
from multiprocess import Pool
from itertools import repeat
import math
from itertools import combinations_with_replacement
import umap
import copy
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib
from random import random
import scipy.cluster.hierarchy as shc
import Levenshtein
import binascii
import logging
from suffix_tree import Tree
#import matplotlib.gridspec as gridspec
from itertools import groupby


class Kmer:
    def __init__(self, kmer, count, nseq, start_pos_list):
        self.kmer = kmer
        self.all_counts = count
        self.start_pos_list = start_pos_list
        indel_size, indel_start, indel_end = get_indel_size(self.kmer, self.start_pos_list)
        self.max_kmer_mask = {(nseq, kmer): max(indel_size or [0])}
        indel_pos_list = [(indel_start[i], indel_end[i], self.kmer) for i in range(len(indel_start))]
        self.consecutive_region = {(nseq, kmer): indel_pos_list}
        self.sum_kmer_consecutive_mask_per_seq = {(nseq, kmer): sum(y for y in indel_size if y > len(self.kmer))}
        self.sum_kmer_mask_per_seq = {(nseq, kmer): sum(indel_size)}
        self.sum_kmer_consecutive_mask = {kmer: sum(y for y in indel_size if y > len(self.kmer))}
                    
    def merge_seq(self, obj):
        self.max_kmer_mask = {**self.max_kmer_mask, **obj.max_kmer_mask}
        self.consecutive_region = {**self.consecutive_region, **obj.consecutive_region}
        self.sum_kmer_consecutive_mask_per_seq = {**self.sum_kmer_consecutive_mask_per_seq, **obj.sum_kmer_consecutive_mask_per_seq}
        self.sum_kmer_mask_per_seq = {**self.sum_kmer_mask_per_seq, **obj.sum_kmer_mask_per_seq}
        for shift_kmer in self.sum_kmer_consecutive_mask:
            self.sum_kmer_consecutive_mask[shift_kmer] += obj.sum_kmer_consecutive_mask[shift_kmer]
    
    def perm_kmer(self, another_kmer, nseq, another_start_pos_list):
            if another_kmer != self.kmer:
                another_kmer_indel_size, another_kmer_indel_start, another_kmer_indel_end = get_indel_size(another_kmer, another_start_pos_list)
                self.max_kmer_mask[(nseq, another_kmer)] = max(another_kmer_indel_size or [0])
                indel_pos_list = [(another_kmer_indel_start[i], another_kmer_indel_end[i], another_kmer) for i in range(len(another_kmer_indel_start))]
                self.consecutive_region[(nseq, another_kmer)] = indel_pos_list
                self.sum_kmer_consecutive_mask_per_seq[(nseq, another_kmer)] = sum(y for y in another_kmer_indel_size if y > len(another_kmer))
                self.sum_kmer_mask_per_seq[(nseq, another_kmer)] = sum(another_kmer_indel_size)
                self.sum_kmer_consecutive_mask[another_kmer] = sum(y for y in another_kmer_indel_size if y > len(another_kmer))

def get_indel_size(kmer, start_pos_list):
    start_pos_no_overlap = []
    start_pos_list = list(start_pos_list)
    start_pos_list.sort()
    for x in start_pos_list:
        if not start_pos_no_overlap:
            start_pos_no_overlap.append(x)
        elif x - start_pos_no_overlap[-1] >= len(kmer):
            start_pos_no_overlap.append(x)

    if len(start_pos_no_overlap) == 1:
        return [len(kmer)], [start_pos_no_overlap[0]], [start_pos_no_overlap[0] + len(kmer)]
    merged_start = []
    merged_end = []
    size = []
    for i in range(len(start_pos_no_overlap)):
        if i == 0:
            start_pos = start_pos_no_overlap[i]
            masked_size = len(kmer)
        elif start_pos_no_overlap[i] - start_pos_no_overlap[i - 1] == len(kmer):
            masked_size += len(kmer)
            if i == len(start_pos_no_overlap) - 1:
                merged_start += [start_pos]
                merged_end += [start_pos + masked_size]
                size += [masked_size]
        else:
            merged_start += [start_pos]
            merged_end += [start_pos + masked_size]
            size += [masked_size]
            start_pos = start_pos_no_overlap[i]
            masked_size = len(kmer)
            if i == len(start_pos_no_overlap) - 1:
                merged_start += [start_pos]
                merged_end += [start_pos + masked_size]
                size += [masked_size]
    return size, merged_start, merged_end

def parse_fasta(filename):
    """Parse a FASTA file and return a dictionary of sequences"""
    sequences = {}
    for record in SeqIO.parse(filename, 'fasta'):
        if record.seq != "":
            sequences[record.description] = str(record.seq)
    return sequences

def count_kmer(all_seq, fasta_file, k_min, k_max):
    shorted_len = min(len(value) for value in all_seq.values())
    output_cnt_file = fasta_file.replace(".fa", "." + str(k_min) + "." + str(k_max) + "mer.txt")
    os.system("/project/holstegelab/Share/yaran/bin/aardvark/build/aardvark kc -k %d -m %d %s > %s " %(min(shorted_len - 1, k_min), min(shorted_len, k_max), fasta_file, output_cnt_file))
    df = pd.read_csv(output_cnt_file, sep='\t', header=None)
    os.system("rm %s" % (output_cnt_file))
    count_dict = dict(zip(df[0], df[1]))
    count_dict = remove_redundant_kmer(count_dict)
    return count_dict



def count_kmer_per_seq(record, seq_identifier, fasta_file, k_min, k_max):
    #all_rounds_kmer = pd.DataFrame()
    #seq_used_motif = pd.DataFrame()
    #all_kmer_df = pd.DataFrame()
    if str(record) != "":
        seq_length = max([len(x) for x in record.split("x")])
        k_max = min([k_max, seq_length])
        tmp_file = fasta_file.replace(".fa", "_") + "tmp_seq_" + str(seq_identifier).replace(" ", "_") + ".fa"
        with open (tmp_file, 'w') as tmp_fasta:
            tmp_fasta.write(">" + seq_identifier + "\n")
            tmp_fasta.write(record + "\n")
           
        output_cnt_file = tmp_file.replace(".fa", "." + str(k_min) + "." + str(k_max) + "mer.txt")
        os.system("/project/holstegelab/Share/yaran/bin/aardvark/build/aardvark kc -k %d -m %d %s > %s " %(min(len(record) - 1, k_min), min(len(record), k_max), tmp_file, output_cnt_file))
        df = pd.read_csv(output_cnt_file, header = None, sep = "\t")
        df.columns = ["kmer", "count"]
        df["ckmer"] = df.apply(lambda x: min([x["kmer"][i:] + x["kmer"][:i] for i in range(len( x["kmer"]))]), axis = 1)
        df["count"] = df.apply(lambda x: math.ceil(x["count"]/len(x["ckmer"])), axis = 1)
        os.system("rm %s" % (output_cnt_file))
        os.system("rm %s" % (tmp_file))
        count_dict = dict(zip(df["ckmer"], df["count"]))
        count_dict = remove_redundant_kmer(count_dict)
        return count_dict
    else:
        return {}

def count_all_kmers(all_seq, fasta_file, k_min, k_max):
    all_count = {}
    for seq in all_seq:
        all_count[seq] = count_kmer_per_seq(all_seq[seq], seq, fasta_file, k_min, k_max)
    return all_count

def count_all_kmers_multi(all_seq, fasta_file, k_min, k_max):
    
    count_list = pool.starmap(count_kmer_per_seq, zip(all_seq.values(), all_seq.keys(), repeat(fasta_file), repeat(k_min), repeat(k_max)))
    seq_list = list(all_seq.keys())
    all_count = {}
    for i in range(len(count_list)):
        all_count[seq_list[i]] = count_list[i]
    return all_count

def is_smallest_kmer(kmer):
    k = 1
    while k <= int((len(kmer) + 1) / 2):
        if len(kmer) % k == 0:
            sub = [kmer[i:i+k] for i in range(0, len(kmer), k)]
            if sub.count(sub[0]) == len(sub):
                unit = sub[0]
                break
            else:
                k += 1
        else:
            k += 1
    if k == int((len(kmer) + 1) / 2) + 1:
        unit = kmer
    return kmer == unit

def remove_redundant_kmer(kmer_dict):
    keys_to_keep = [key for key in kmer_dict if is_smallest_kmer(key)]
    keys_to_keep = [key for key in keys_to_keep if kmer_dict[key] > 1]
    return {key:value for key,value in kmer_dict.items() if key in keys_to_keep}


def find_positions(input_string, substring_list):
    st = Tree({"seq": input_string})
    kmer_start = {}
    for substring in substring_list:
        start_pos = []
        for id_, path in st.find_all(substring):
            start_pos += [len(input_string)-len(str(path).strip("$").replace(" ", ""))]
        kmer_start[substring] = start_pos
    return kmer_start

def find_sporadic(seq, kmer_list):
    if seq == '':
        return []
    else:
        start_pos_dict = find_positions(seq, kmer_list)
        start_pos_dict = {key:value for key, value in start_pos_dict.items() if value != {}}
        return start_pos_dict


def select_kmer_fast(kmer_per_seq):
    kmer_perm_dict = {}

    for kmer_one_seq in kmer_per_seq:
        for kmer_set in kmer_one_seq:
            if not kmer_perm_dict:
                for a in kmer_set.perm_kmer:
                    kmer_perm_dict[a] = kmer_set
            else:
                if kmer_set.kmer in kmer_perm_dict:
                    kmer_perm_dict[kmer_set.kmer].merge_seq(kmer_set)
                else:
                    for a in kmer_set.perm_kmer:
                        kmer_perm_dict[a] = kmer_set

    kmer_all_seq = [kmer for kmer in set(kmer_perm_dict.values()) if max(kmer.max_kmer_mask.values()) / len(kmer.kmer) > 1]
    return sorted(kmer_all_seq, key = lambda x: [max(x.max_kmer_mask.values()), -max(x.max_kmer_mask.values())/len(x.kmer), x.count], reverse = True)

def pairwise_alignment(seq1, seq2):
    from Bio import pairwise2
    alignments = pairwise2.align.globalms(seq1, seq2, 0, -2, -5, -5)
    scores = []
    for alignment in alignments:
        scores.append(alignment.score)
    return max(scores)

def get_best_perm_kmer(kmer_obj, previous_kmer):
    perm_kmer = []
    regions = list(kmer_obj.consecutive_region.values())
    regions = [item for sublist in regions for item in sublist]
    for region in regions:
        if region[1] - region[0] > len(kmer_obj.kmer):
            if region[2] not in perm_kmer:
                perm_kmer += [region[2]]
    if len(perm_kmer) == 0:
        for region in regions:
            if region[2] not in perm_kmer:
                perm_kmer += [region[2]]

    perm_kmer_score = list(map(lambda i: sum([pairwise_alignment(a, i) for a in previous_kmer]), perm_kmer))

    kmer_idx = perm_kmer_score.index(max(perm_kmer_score))
    return perm_kmer[kmer_idx]

def consecutive_mask_seq(seq, seq_number, kmer, kmer_obj):
    regions = kmer_obj.consecutive_region[(seq_number, kmer)]
    masked_region = []
    for region in regions:
        if region[1] - region[0] > len(kmer):
            if "x" not in seq[region[0]:region[1]]:
                seq = seq[:region[0]] + "x" * (region[1] - region[0]) + seq[region[1]:]
                masked_region += [region]
    return seq, {(seq_number, kmer): masked_region}

def sporadic_mask_seq(seq, seq_number, kmer, kmer_obj):
    regions = kmer_obj.consecutive_region[(seq_number, kmer)]
    masked_region = []
    for region in regions:
        if region[1] - region[0] == len(kmer):
            if "x" not in seq[region[0]:region[1]]:
                seq = seq[:region[0]] + "x" * (region[1] - region[0]) + seq[region[1]:]
                masked_region += [region]
    return seq, {(seq_number, kmer): masked_region}

def sporadic_mask_seq_new(seq, seq_number, kmer_dict):
    to_be_masked_dict = {}
    for kmer in kmer_dict:
        regions = get_indel_size(kmer, list(kmer_dict[kmer]))
        #print(regions)
        masked_region = []
        for i in range(len(regions[0])):
            if regions[0][i] == len(kmer):
                if "x" not in seq[regions[1][i]:regions[2][i]]:
                    seq = seq[:regions[1][i]] + "x" * (regions[2][i] - regions[1][i]) + seq[regions[2][i]:]
                    masked_region += [(regions[1][i], regions[2][i], kmer)]
        to_be_masked_dict[(seq_number, kmer)] = masked_region

    return seq, to_be_masked_dict

def drop_kmer(selected_kmer, kmer_list_per_seq, df):
    for i in range(len(kmer_list_per_seq)):
        if selected_kmer == kmer_list_per_seq[i][0].kmer:
            df = df.drop(df[df["nseq"] == kmer_list_per_seq[i][0].nseq & df["cnt"] >= df.loc[df['ckmer'] == selected_kmer, 'cnt'].values[0]].index)
            
    return df


def kmer_selection(all_seq, kmer_counts):
    all_kmer_counts  = {}
    for seq in kmer_counts:
        for kmer in kmer_counts[seq]:
            if kmer in all_kmer_counts:
                all_kmer_counts[kmer] += kmer_counts[seq][kmer]
            else:
                all_kmer_counts[kmer] = kmer_counts[seq][kmer]
    max_kmer_counts = {}
    for seq in kmer_counts:
        for kmer in kmer_counts[seq]:
            if kmer in max_kmer_counts:
                max_kmer_counts[kmer] = max([kmer_counts[seq][kmer], max_kmer_counts[kmer]])
            else:
                max_kmer_counts[kmer] = kmer_counts[seq][kmer]
    

    kmer_to_screen_list = sorted(max_kmer_counts, key=lambda k: [max_kmer_counts[k], -len(k)], reverse = True)

    st = Tree(all_seq)
    kmer_set_dict = {}
    all_start_pos = {}
    max_mask = 0
    for kmer in kmer_to_screen_list:
        if len(kmer) * max_kmer_counts[kmer] > max_mask:
            for shift_kmer in [kmer[x:] + kmer[:x] for x in range(len(kmer))]:
                start_pos = {seq:[] for seq in all_seq}
                for id_, path in st.find_all(shift_kmer):
                    start_pos[id_] += [len(all_seq[id_])-len(str(path).strip("$").replace(" ", ""))]
                all_start_pos[shift_kmer] = start_pos
            
            for seq in all_seq:
                kmer_object = Kmer(kmer, all_kmer_counts[kmer], seq, all_start_pos[kmer][seq])
                for shift_kmer in [kmer[x:] + kmer[:x] for x in range(len(kmer))]:
                    kmer_object.perm_kmer(shift_kmer, seq, all_start_pos[shift_kmer][seq])
                if kmer in kmer_set_dict:
                    kmer_set_dict[kmer].merge_seq(kmer_object)
                else:
                    kmer_set_dict[kmer] = kmer_object
            if max(kmer_set_dict[kmer].max_kmer_mask.values()) / len(kmer) > 1:
                max_mask = max(max_mask, max(kmer_set_dict[kmer].max_kmer_mask.values()))

        else:
            continue
    st = None
    kmer_all_seq = [kmer for kmer in set(kmer_set_dict.values()) if max(kmer.max_kmer_mask.values()) / len(kmer.kmer) > 1]
    return sorted(kmer_all_seq, key = lambda x: [max(x.max_kmer_mask.values()), -len(x.kmer)], reverse = True)



def mask_all_seq(input_fasta, min_k, max_k):
    n = 0
    all_seq = parse_fasta(input_fasta)
    seq_to_be_masked = all_seq.copy()
    seq_list = [str(seq) for seq in all_seq.values()]
    selected_kmer_list = []
    all_seq_masked_dict = {}
    output_fasta = input_fasta.replace(".fa", ".tmp.fa")
    os.system("cp %s %s" %(input_fasta, output_fasta))
    while max([len(x) for x in seq_list]) != 1:
        print("counting kmer")
        all_seq_dict = parse_fasta(output_fasta)
        kmer_counts = count_all_kmers_multi(all_seq_dict, output_fasta, min_k, max_k)
        print("selecting kmer")
        kmer_set_all_seq = kmer_selection(all_seq_dict, kmer_counts)
        if len(kmer_set_all_seq) == 0:
            break
        else:
            selected_kmer = kmer_set_all_seq[0]
            if len(selected_kmer.kmer) == 1:
                break
            else: 
                if n == 0: 
                    sum_max_perm_kmer_region = {}
                    for kmer_region in selected_kmer.max_kmer_mask:
                        if kmer_region[1] in sum_max_perm_kmer_region:
                            sum_max_perm_kmer_region[kmer_region[1]] += selected_kmer.max_kmer_mask[kmer_region]
                        else:
                            sum_max_perm_kmer_region[kmer_region[1]] = selected_kmer.max_kmer_mask[kmer_region]
                    kmer = sorted(sum_max_perm_kmer_region, key = sum_max_perm_kmer_region.get, reverse = True)[0]
                else:
                    logging.debug("getting best perm kmer")
                    kmer = get_best_perm_kmer(selected_kmer, selected_kmer_list)
                    logging.debug("best kmer selected")
                logging.debug(kmer)
                print(kmer)
                selected_kmer_list += [kmer]
                one_kmer_masked = {}
                for region in selected_kmer.max_kmer_mask:
                    if region[1] == kmer:
                        if selected_kmer.max_kmer_mask[region] != 0:
                            seq_to_be_masked[region[0]], masked_region = consecutive_mask_seq(seq_to_be_masked[region[0]], region[0], kmer, selected_kmer)
                            one_kmer_masked = {**one_kmer_masked, **masked_region}
                        
                if all(value == [] for value in one_kmer_masked.values()):
                    break
                
                else:
                    all_seq_masked_dict = {**all_seq_masked_dict, **one_kmer_masked}
                seq_list = [x.split("x") for x in seq_to_be_masked.values()]
                seq_list = [item for sublist in seq_list for item in sublist]
                with open (output_fasta, 'w') as output_fasta_file:
                    for seq in seq_to_be_masked:
                        output_fasta_file.write(">" + seq + "\n")
                        output_fasta_file.write(seq_to_be_masked[seq] + "\n")
                max_k = min([max_k, max([len(x) for x in seq_list])])
                n += 1
    
    os.system("rm %s" %(output_fasta))
    all_seq_masked_sporadic_seq = {}
    selected_kmer_list = sorted(selected_kmer_list, key = len, reverse = True)
    print("finding unconsecutive motifs")
    sporadic_motifs = {}
    all_seq_masked_sporadic_seq = {}
    for seq in seq_to_be_masked:
        sporadic_motifs[seq] = find_sporadic(seq_to_be_masked[seq], selected_kmer_list)

        seq_to_be_masked[seq], masked_region = sporadic_mask_seq_new(seq_to_be_masked[seq], seq, sporadic_motifs[seq])
        all_seq_masked_sporadic_seq = {**all_seq_masked_sporadic_seq, **masked_region}

    
    return all_seq_masked_dict, all_seq_masked_sporadic_seq

def mask_all_seq_ref_motif(input_fasta, min_k, max_k, ref_motifs_list):
    ref_motifs_dict = {}
    for motif in ref_motifs_list:
        ref_motifs_dict[motif] = min([motif[x:] + motif[:x] for x in range(len(motif))])
    ref_motifs_dict_r = {value: key for key, value in ref_motifs_dict.items()}
    n = 0
    all_seq = parse_fasta(input_fasta)
    seq_to_be_masked = all_seq.copy()
    seq_list = [str(seq) for seq in all_seq.values()]
    selected_kmer_list = []
    all_seq_masked_dict = {}
    output_fasta = input_fasta.replace(".fa", ".tmp.fa")
    os.system("cp %s %s" %(input_fasta, output_fasta))
    while max([len(x) for x in seq_list]) != 1:
        print("counting kmer")
        all_seq_dict = parse_fasta(output_fasta)
        kmer_counts = count_all_kmers_multi(all_seq_dict, output_fasta, min_k, max_k)
        print("selecting kmer")
        kmer_set_all_seq = kmer_selection(all_seq_dict, kmer_counts)
        if len(kmer_set_all_seq) == 0:
            break
        else:
            selected_kmer = kmer_set_all_seq[0]
            if len(selected_kmer.kmer) == 1:
                break
            else: 
                if n == 0:
                    if selected_kmer.kmer in ref_motifs_dict_r:
                        kmer = ref_motifs_dict_r[selected_kmer.kmer]
                    else:
                        sum_max_perm_kmer_region = {}
                        for kmer_region in selected_kmer.max_kmer_mask:
                            if kmer_region[1] in sum_max_perm_kmer_region:
                                sum_max_perm_kmer_region[kmer_region[1]] += selected_kmer.max_kmer_mask[kmer_region]
                            else:
                                sum_max_perm_kmer_region[kmer_region[1]] = selected_kmer.max_kmer_mask[kmer_region]
                        kmer = sorted(sum_max_perm_kmer_region, key = sum_max_perm_kmer_region.get, reverse = True)[0]
                else:
                    if selected_kmer.kmer in ref_motifs_dict_r:
                        kmer = ref_motifs_dict_r[selected_kmer.kmer]
                    else:
                        logging.debug("getting best perm kmer")
                        kmer = get_best_perm_kmer(selected_kmer, selected_kmer_list)
                        logging.debug("best kmer selected")
                logging.debug(kmer)
                print(kmer)
                selected_kmer_list += [kmer]
                one_kmer_masked = {}
                for region in selected_kmer.max_kmer_mask:
                    if region[1] == kmer:
                        if selected_kmer.max_kmer_mask[region] != 0:
                            seq_to_be_masked[region[0]], masked_region = consecutive_mask_seq(seq_to_be_masked[region[0]], region[0], kmer, selected_kmer)
                            one_kmer_masked = {**one_kmer_masked, **masked_region}
                        
                if all(value == [] for value in one_kmer_masked.values()):
                    break
                
                else:
                    all_seq_masked_dict = {**all_seq_masked_dict, **one_kmer_masked}
                seq_list = [x.split("x") for x in seq_to_be_masked.values()]
                seq_list = [item for sublist in seq_list for item in sublist]
                with open (output_fasta, 'w') as output_fasta_file:
                    for seq in seq_to_be_masked:
                        output_fasta_file.write(">" + seq + "\n")
                        output_fasta_file.write(seq_to_be_masked[seq] + "\n")
                n += 1

    os.system("rm %s" %(output_fasta))
    all_seq_masked_sporadic_seq = {}
    selected_kmer_list = list(set(selected_kmer_list + ref_motifs_list))
    selected_kmer_list = sorted(selected_kmer_list, key = len, reverse = True)
    print("finding unconsecutive motifs")
    sporadic_motifs = {}
    all_seq_masked_sporadic_seq = {}
    for seq_number in list(range(len(seq_to_be_masked))):
        sporadic_motifs[seq_number] = find_sporadic(seq_to_be_masked[seq_number], selected_kmer_list)

        seq_to_be_masked[seq_number], masked_region = sporadic_mask_seq_new(seq_to_be_masked[seq_number], seq_number, sporadic_motifs[seq_number])
        all_seq_masked_sporadic_seq = {**all_seq_masked_sporadic_seq, **masked_region}
    
    return all_seq_masked_dict, all_seq_masked_sporadic_seq

def tag_all_seq(all_seq, all_seq_masked_dict, all_seq_masked_sporadic_seq):
    seq_df = pd.DataFrame()
    for seq in all_seq:
        seq_array = pd.DataFrame(list(str(all_seq[seq])), columns = [seq])
        masked_region_one_seq = [all_seq_masked_dict[region] for region in all_seq_masked_dict if region[0] == seq]
        masked_region_one_seq = [item for sublist in masked_region_one_seq for item in sublist]
        for region in masked_region_one_seq:
            seq_array.iloc[region[0]:region[1], 0] = region[2]
        masked_region_one_seq = [all_seq_masked_sporadic_seq[region] for region in all_seq_masked_sporadic_seq if region[0] == seq]
        masked_region_one_seq = [item for sublist in masked_region_one_seq for item in sublist]
        for region in masked_region_one_seq:
            seq_array.iloc[region[0]:region[1], 0] = region[2]
            
        seq_df = pd.concat([seq_df, seq_array], axis = 1)
    return seq_df

def generate_hex_chr_list():
    usable_hex = list(range(1,256))
    usable_hex.remove(60)
    usable_hex.remove(61)
    usable_hex.remove(62)
    usable_hex.remove(45)
    usable_hex.remove(32)
    usable_hex.remove(10)
    usable_hex.remove(13)
    usable_hex = ['{:02x}'.format(x) for x in usable_hex]
    return usable_hex

def assign_chr_to_motif(all_unique_motif, printable_asii, input_fasta_file_name):
    if len(all_unique_motif) > len(printable_asii):
        return None, "too many motifs"        
    else:
        printable_asii_to_include = printable_asii[:len(all_unique_motif)]
        motif_dict = dict(zip(all_unique_motif, printable_asii_to_include))
        pd.DataFrame(motif_dict.items(), columns=['motif', 'ascii']).to_csv(input_fasta_file_name.replace(".fa", "_motif_to_ascii.txt"), sep = "\t", index = False)
        return motif_dict


def write_seq_in_hex_chr(all_seq_motif_df, motif_dict, input_fasta_file_name, random_num):
    all_seq_motif_chr_df = all_seq_motif_df.replace(motif_dict)
    # write into fasta file
    chr_fasta_file = input_fasta_file_name.replace(".fa", "_motif_in_hex_" + str(random_num) + ".hex")
    seq_chr_dict = {}    
    with open (chr_fasta_file, 'w') as chr_fasta:
        for column in all_seq_motif_df:
            seq_description = column
            seq_motif_in_chr = all_seq_motif_chr_df[column].dropna()
            seq_motif = all_seq_motif_df[column].dropna()
            seq_chr = []
            i = 0
            while i < len(seq_motif):
                seq_chr += [seq_motif_in_chr[i]]
                i += len(seq_motif[i])
            seq_chr_string = " ".join(seq_chr)
            chr_fasta.write(">" + seq_description + "\n")
            chr_fasta.write(seq_chr_string + "\n")
            seq_chr_dict[column] = seq_chr_string
            
    motif_fasta_file = input_fasta_file_name.replace(".fa", "_motif_to_hex_" + str(random_num) + ".fa")
    with open (motif_fasta_file, 'w') as dict_fasta:
        for item in motif_dict:
            dict_fasta.write(">" + motif_dict[item] + "\n")
            dict_fasta.write(item + "\n")        
    return seq_chr_dict
            


def compress_string(string, motif_dict):
    parts = string.split(" ")
    transformed_string = ""
    for key, group in groupby(parts):
        count = len(list(group))
        transformed_string += f"{motif_dict[key]}{count} "
    return transformed_string.strip()


def write_compressed_seq(seq_chr_dict, motif_dict, input_fasta_file_name, random_num):
    compressed_fasta_file = input_fasta_file_name.replace(".fa", "_compressed_representation_" + str(random_num) + ".fa")
    motif_dict_r = {value: key for key, value in motif_dict.items()}
    with open (compressed_fasta_file, 'w') as compressed_fasta_file:
        for seq in seq_chr_dict:
            compressed_fasta_file.write(">" + seq + "\n")
            compressed_fasta_file.write(compress_string(seq_chr_dict[seq], motif_dict_r) + "\n")

def hex_to_bytes(hex_string):
    return binascii.unhexlify(hex_string)


def normalize_edit_distance(distance, seq_1_length, seq_2_length):
    if max(seq_1_length, seq_2_length) != 0:
        return distance / max(seq_1_length, seq_2_length)
    else:
        return 0


def edit_distance_between_seq_byte(input_seq_chr_dict):
    seq_chr_dict = {}
    for seq in input_seq_chr_dict:
        if input_seq_chr_dict[seq] == "":
            seq_chr_dict[seq] = str.encode(" ")
        else:
            seq_chr_dict[seq] = hex_to_bytes(input_seq_chr_dict[seq].replace(" ", ""))
    
    seq_distance_df = pd.DataFrame(list(combinations_with_replacement(list(seq_chr_dict.keys()), 2)), columns=['nseq_1', 'nseq_2'])
    seq_distance_df["seq_1"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_1"]], axis = 1)
    seq_distance_df["seq_2"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_2"]], axis = 1)
    ls = pool.starmap(Levenshtein.distance, zip(seq_distance_df.seq_1, seq_distance_df.seq_2))
    seq_distance_df['distance'] = ls
    seq_distance_df['distance'] = seq_distance_df.apply(lambda x: normalize_edit_distance(x["distance"], len(x["seq_1"]), len(x["seq_2"])), axis = 1)

    seq_distance_df = seq_distance_df.drop(columns = ["seq_1", "seq_2"])

    seq_distance_df_2 = seq_distance_df.copy()
    seq_distance_df_2["nseq_2"] = seq_distance_df["nseq_1"]
    seq_distance_df_2["nseq_1"] = seq_distance_df["nseq_2"]
    seq_distance_df = pd.concat([seq_distance_df, seq_distance_df_2], ignore_index = True)
    seq_distance_df = seq_distance_df.drop_duplicates().pivot(index = 'nseq_1', columns = 'nseq_2', values = 'distance').astype(float)

    return seq_distance_df

def get_motif_pairwise_distance(motif_list):
    score_df = pd.DataFrame(list(combinations_with_replacement(motif_list, 2)), columns=['motif_1', 'motif_2'])
    
    score_df['tuple'] = list(zip(score_df.motif_1, score_df.motif_2))
    ls = pool.starmap(pairwise_alignment, score_df['tuple'].tolist())
    score_df['score'] = ls
    score_df = score_df.drop(columns = ['tuple'])

    score_df_2 = score_df.copy()
    score_df_2["motif_2"] = score_df["motif_1"]
    score_df_2["motif_1"] = score_df["motif_2"]
    score_df = pd.concat([score_df, score_df_2], ignore_index = True)
    score_df = score_df.drop_duplicates().pivot(index = 'motif_1', columns = 'motif_2', values = 'score').astype(float)

    return score_df

def edit_distance_between_seq_byte(input_seq_chr_dict):
    seq_chr_dict = {}
    for seq in input_seq_chr_dict:
        if input_seq_chr_dict[seq] == "":
            seq_chr_dict[seq] = str.encode(" ")
        else:
            seq_chr_dict[seq] = hex_to_bytes(input_seq_chr_dict[seq].replace(" ", ""))
    
    seq_distance_df = pd.DataFrame(list(combinations_with_replacement(list(seq_chr_dict.keys()), 2)), columns=['nseq_1', 'nseq_2'])
    seq_distance_df["seq_1"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_1"]], axis = 1)
    seq_distance_df["seq_2"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_2"]], axis = 1)
    
    ls = pool.starmap(Levenshtein.distance, zip(seq_distance_df.seq_1, seq_distance_df.seq_2))
    seq_distance_df['distance'] = ls
    seq_distance_df['distance'] = seq_distance_df.apply(lambda x: normalize_edit_distance(x["distance"], len(x["seq_1"]), len(x["seq_2"])), axis = 1)

    seq_distance_df = seq_distance_df.drop(columns = ["seq_1", "seq_2"])

    seq_distance_df_2 = seq_distance_df.copy()
    seq_distance_df_2["nseq_2"] = seq_distance_df["nseq_1"]
    seq_distance_df_2["nseq_1"] = seq_distance_df["nseq_2"]
    seq_distance_df = pd.concat([seq_distance_df, seq_distance_df_2], ignore_index = True)
    seq_distance_df = seq_distance_df.drop_duplicates().pivot(index = 'nseq_1', columns = 'nseq_2', values = 'distance').astype(float)

    return seq_distance_df

def get_motif_pairwise_distance_new(motif_list):
    score_df = pd.DataFrame(list(combinations_with_replacement(motif_list, 2)), columns=['motif_1', 'motif_2'])
    
    score_df['tuple'] = list(zip(score_df.motif_1, score_df.motif_2))
    ls = pool.starmap(pairwise_alignment, score_df['tuple'].tolist())
    score_df['score'] = ls
    score_df = score_df.drop(columns = ['tuple'])
    
    score_df_2 = score_df.copy()
    score_df_2["motif_2"] = score_df["motif_1"]
    score_df_2["motif_1"] = score_df["motif_2"]
    score_df = pd.concat([score_df, score_df_2], ignore_index = True)
    score_df = score_df.drop_duplicates().pivot(index = 'motif_1', columns = 'motif_2', values = 'score').astype(float)

    return score_df

def run_umap(dm):
    umap_f = umap.UMAP(n_components = 1, metric = "manhattan", n_neighbors = 10, min_dist = 0.5, random_state = 0)
    X_transform_L2 = pd.DataFrame(umap_f.fit_transform(dm))
    X_transform_L2.columns = ["dimension_reduction"]
    X_transform_L2["motif"] = dm.index
    X_transform_L2 = X_transform_L2.astype({"dimension_reduction": float})

    return X_transform_L2

def map_score_to_alignment(all_seq_motifs, X_transform_L2):
    motif_dict = dict(zip(X_transform_L2.motif, map(float, X_transform_L2.dimension_reduction)))
    motif_dict['nan'] = np.nan
    df = all_seq_motifs.applymap(lambda x: motif_dict[str(x)])
    return df.transpose()

def summarize_motif(all_seq_motifs, reads_name, input_fasta, random_num, seq_distance_df, figtitle):
    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary["index"] = summary.apply(lambda x: x["index"].split("_")[-1], axis = 1)
    summary = summary.merge(reads_name, left_on = "index", right_on = "nseq", how = "left").drop(columns = ["index", "nseq"]).set_index("read_name").fillna(0)
    summary.to_csv(input_fasta.replace(".fa", "_umap_" + str(random_num) + "_only_consecutive_same_string_motif_count" + ".txt"), sep = "\t")
    seq_distance_df = seq_distance_df.reset_index().merge(reads_name, left_on = "nseq_1", right_on = "nseq", how = "left").drop(columns = ["nseq_1", "nseq"]).set_index("read_name")
    clusters = shc.linkage(seq_distance_df, method='average', metric="euclidean")
    
    plt.figure(figsize=(25, 50), dpi = 200)
    d = shc.dendrogram(Z = clusters, labels = seq_distance_df.index, orientation = "right")
    plt.title(figtitle)
    plt.savefig(input_fasta.replace(".fa", "_umap_" + str(random_num) + "_only_consecutive_same_string_cluster" + ".png"), bbox_inches = "tight")
    plt.clf()
    return list(reversed(d["ivl"]))


def plot_df(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, population_df):
    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[4,1,40,10], height_ratios=[1], wspace=0.02)

    for col in range(3):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary = summary.set_index("index")
    clusters = shc.linkage(seq_distance_df, method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.transpose()
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")
    cbaxes = fig.add_axes([0.925, 0.7, .02, 0.1], title="population")
    

    sns.heatmap(df, cmap=cmap2, ax=all_axes[2], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[2].set(ylabel="")
    all_axes[2].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[2].tick_params(axis ='x', which ='major', rotation=90)
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    
    all_seq_motifs = all_seq_motifs.transpose()
    
    for i in range(all_seq_motifs.shape[1]):
        j = 0
        while j < all_seq_motifs.shape[0]:
            if str(all_seq_motifs.iloc[j, i]) != "nan":
                all_axes[2].add_patch(Rectangle((j, i), len(all_seq_motifs.iloc[j, i]), 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
                j += len(all_seq_motifs.iloc[j, i])
            else:
                break


    seq_population = pd.DataFrame({'seq': seq_order_list})
    seq_population["Sample"] = seq_population.apply(lambda x: x["seq"].split("#")[0], axis = 1)
        
    seq_population = seq_population.merge(population_df[["Sample", "Group"]], on = "Sample", how = "left")
    g = len(seq_population.Group.unique())
    color = ["skyblue", "silver", "lightgreen", "orange", "violet"]
    color_df = pd.DataFrame({'Group': list(seq_population.Group.unique()), 'color': color[:g], 'n': list(range(0, 2*g, 2))})
    seq_population = seq_population.merge(color_df, on = "Group", how = "left")       
    cmap4 = (matplotlib.colors.ListedColormap(color[:g]))
    sns.heatmap(seq_population.set_index("seq")[["n"]], ax=all_axes[1], cmap=cmap4, cbar=False, xticklabels = False, yticklabels = False, linewidths=0.1, linecolor='black')
    all_axes[1].set(ylabel="")

    bounds = list(range(-1, 2*g, 2))
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap4.N)
    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap4, norm=norm), cax=cbaxes, ticks=list(range(0, 2*g, 2)), spacing='uniform', orientation='vertical')
    cbaxes.set_yticklabels(list(seq_population.Group.unique()))
    
    
    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')

    all_axes[2].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")


def plot_df_reads(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle):
    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=3, nrows=1, width_ratios=[4,40,5], height_ratios=[1], wspace=0.02)

    for col in range(2):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary = summary.set_index("index")
    clusters = shc.linkage(seq_distance_df, method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.transpose()
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")    

    sns.heatmap(df, cmap=cmap2, ax=all_axes[1], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[1].set(ylabel="")
    all_axes[1].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[1].tick_params(axis ='x', which ='major', rotation=90)
    cbar_ax.set_yticklabels(dm.motif.tolist()) 
    
    all_seq_motifs = all_seq_motifs.transpose()
    
    for i in range(all_seq_motifs.shape[1]):
        j = 0
        while j < all_seq_motifs.shape[0]:
            if str(all_seq_motifs.iloc[j, i]) != "nan":
                all_axes[1].add_patch(Rectangle((j, i), len(all_seq_motifs.iloc[j, i]), 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
                j += len(all_seq_motifs.iloc[j, i])
                
            else:
                break

    all_axes[1].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")



def msa_with_characters(input_fasta_file_name, random_num):
    chr_fasta_file = input_fasta_file_name.replace(".fa", "_motif_in_hex_" + str(random_num) + ".hex")
    os.system("/project/holstegelab/Share/yaran/bin/libexec/mafft/hex2maffttext %s > %s" %(chr_fasta_file, chr_fasta_file.replace(".hex", ".ASCII")))
    os.system("mafft --text --op 2.0 --ep 0.1 %s > %s" %(chr_fasta_file.replace(".hex", ".ASCII"), chr_fasta_file.replace(".hex", "_mafft_output.ASCII")))
    os.system("/project/holstegelab/Share/yaran/bin/libexec/mafft/maffttext2hex %s > %s" % (chr_fasta_file.replace(".hex", "_mafft_output.ASCII"), chr_fasta_file.replace(".hex", "_mafft_output.hex")))
    msa_result = list(SeqIO.parse(chr_fasta_file.replace(".hex", "_mafft_output.hex"), "fasta"))
    
    return msa_result

def prepare_for_plotting(msa_result, motif_dict, dm):
    motif_dict_r = {y: x for x, y in motif_dict.items()}
    motif_dict_r["--"] = ''
    
    msa_df = pd.DataFrame()
    for i in range(len(msa_result)):
        single_seq = str(msa_result[i].seq)
        msa_df = pd.concat([msa_df, pd.DataFrame([single_seq[i:i+2] for i in range(0, len(single_seq), 2)], columns = [msa_result[i].description])], axis = 1)
    
    msa_df["motifs"] = msa_df.apply(lambda x: x.unique().tolist(), axis = 1)
    msa_df["max_length"] = msa_df.apply(lambda x: max(len(motif_dict_r[i]) for i in x["motifs"]), axis = 1)
    msa_df = msa_df.drop(columns = ["motifs"])
    max_motif_length = msa_df["max_length"].tolist()
    msa_adjusted_df = msa_df.loc[msa_df.index.repeat(msa_df.max_length)]
    msa_adjusted_df = msa_adjusted_df.drop(columns = ["max_length"])
    msa_adjusted_df = msa_adjusted_df.transpose()
    msa_adjusted_df = msa_adjusted_df.applymap(lambda x: motif_dict_r[x])
    msa_adjusted_df.columns = range(msa_adjusted_df.columns.size)

    score_dict = dict(zip(dm.motif, dm.dimension_reduction))
    score_dict[''] = np.nan
    df = msa_adjusted_df.applymap(lambda x: score_dict[x])
    
    return df, max_motif_length
    

def plot_msa_df(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, population_df, max_motif_length):

    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[4,1,40,10], height_ratios=[1], wspace=0.01)

    for col in range(3):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    clusters = shc.linkage(seq_distance_df, method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))


    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    cmap2.set_over('none')

    
    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")
    cbaxes = fig.add_axes([0.925, 0.7, .02, 0.1], title="population")

    sns.heatmap(df, cmap=cmap2, ax=all_axes[2], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[2].set(ylabel="")
    all_axes[2].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[2].tick_params(axis ='x', which ='major', rotation=90)
    cbar_ax.set_yticklabels(dm.motif.tolist()) 

    for i in range(df.shape[0]):
        j = 0
        m = 0
        while j < df.shape[1]:
            if str(df.iloc[i, j]) != "nan":
                all_axes[2].add_patch(Rectangle((j, i), max_motif_length[m], 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
            j += max_motif_length[m]
            m += 1

    seq_population = pd.DataFrame({'seq': seq_order_list})
    seq_population["Sample"] = seq_population.apply(lambda x: x["seq"].split("#")[0], axis = 1)
    
    seq_population = seq_population.merge(population_df[["Sample", "Group"]], on = "Sample", how = "left")
    g = len(seq_population.Group.unique())
    color = ["skyblue", "silver", "lightgreen", "orange", "violet"]
    color_df = pd.DataFrame({'Group': list(seq_population.Group.unique()), 'color': color[:g], 'n': list(range(0, 2*g, 2))})
    seq_population = seq_population.merge(color_df, on = "Group", how = "left")       
    cmap4 = (matplotlib.colors.ListedColormap(color[:g]))
    sns.heatmap(seq_population.set_index("seq")[["n"]], ax=all_axes[1], cmap=cmap4, cbar=False, xticklabels = False, yticklabels = False, linewidths=0.1, linecolor='black')
    all_axes[1].set(ylabel="")

    bounds = list(range(-1, 2*g, 2))
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap4.N)
    cb = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap4, norm=norm), cax=cbaxes, ticks=list(range(0, 2*g, 2)), spacing='uniform', orientation='vertical')
    cbaxes.set_yticklabels(list(seq_population.Group.unique()))
    
    cb.outline.set_linewidth(0.1)
    cb.outline.set_edgecolor('black')

    all_axes[2].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")


def plot_msa_df_reads(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, max_motif_length):

    fig = plt.figure(figsize=(min(max(50, 0.015 * df.shape[1]), 120), min(120,0.5 * df.shape[0])), dpi = 300)
    spec = fig.add_gridspec(ncols=3, nrows=1, width_ratios=[4,40,5], height_ratios=[1], wspace=0.01)

    for col in range(2):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    clusters = shc.linkage(seq_distance_df, method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.reindex(seq_order_list)

    cmap2 = copy.copy(plt.get_cmap('YlGnBu_r'))
    cmap2.set_over('none')

    cbar_ax = fig.add_axes([0.975, .2, .02, .6], title="motif")

    sns.heatmap(df, cmap=cmap2, ax=all_axes[1], cbar_ax = cbar_ax, cbar_kws={"ticks": list(map(float, dm.dimension_reduction))})
    all_axes[1].set(ylabel="")
    all_axes[1].tick_params(right=True, left = False, top=False, labelright=True, labelleft=False,labeltop=False,rotation=0)
    all_axes[1].tick_params(axis ='x', which ='major', rotation=90)
    cbar_ax.set_yticklabels(dm.motif.tolist()) 


    for i in range(df.shape[0]):
        j = 0
        m = 0
        while j < df.shape[1]:
            if str(df.iloc[i, j]) != "nan":
                all_axes[1].add_patch(Rectangle((j, i), max_motif_length[m], 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
            j += max_motif_length[m]
            m += 1

    all_axes[1].title.set_text(figtitle)
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")

#######################


parser = argparse.ArgumentParser(description='MotifScope')
parser.add_argument('--analysis_type', dest = 'analysis_type',
                    help='type of input sequences [assembly / reads].', type = str,
                    required=True)

parser.add_argument('-i', '--input', default = None, dest='input_fasta_to_count',
                    metavar="aligned_None_None_None_rfc1.centenarian.blood.hpc.one.line_1014_4129_one_line.fa", type=str,
                    help='input fasta file to count')

parser.add_argument('-mink', '--min_kmer', default = None, dest='min_kmer_size',
                    metavar=5, type=int,
                    help='minimum k to count')

parser.add_argument('-maxk', '--max_kmer', default = None, dest='max_kmer_size',
                    metavar=5, type=int,
                    help='maximum k to count')

parser.add_argument('-t', '--title', default = None, dest='title',
                    metavar="RFC11", type=str,
                    help='title of the plot')

parser.add_argument('-msa', '--msa', dest = 'run_msa', type = str, 
                    help = 'Boolean (True/False).', default = 'False')

parser.add_argument('-m', '--m', dest = 'motif_guided', type = str, 
                    help = 'Boolean (True/False).', default = 'False')

parser.add_argument('-motif', '--motif', dest = 'ref_motifs', type = str, required = False,
                    help = 'file with ref motifs separated by tabs.', metavar = 'motifs.txt')

parser.add_argument('-p', '--population', default = None, dest='population',
                    metavar="metadata.txt", type=str,
                    help='population metadata file')

sys.getrecursionlimit()
args = parser.parse_args()

title = args.title
analysis_type = args.analysis_type
input_fasta_to_count = args.input_fasta_to_count
max_kmer_size = args.max_kmer_size
min_kmer_size = args.min_kmer_size
run_msa = args.run_msa
motif_guided = args.motif_guided
ref_motifs = args.ref_motifs
original_limit = sys.getrecursionlimit()

new_limit = 3000
sys.setrecursionlimit(new_limit)

random_num = str(random()).replace('.', '')
all_seq = list(SeqIO.parse(input_fasta_to_count, "fasta"))
all_seq_dict = parse_fasta(input_fasta_to_count)

pool = Pool()

if motif_guided == "False":
    all_seq_masked, all_seq_masked_sporadic = mask_all_seq(input_fasta_to_count, min_kmer_size, max_kmer_size)
elif motif_guided == "True":
    with open(ref_motifs, 'r') as ref:
        ref_motifs = ref.readlines()
    ref_motifs_list = [i.strip().split("\t") for i in ref_motifs]
    ref_motifs_list = [i for l in ref_motifs_list for i in l]
    all_seq_masked, all_seq_masked_sporadic = mask_all_seq_ref_motif(input_fasta_to_count, min_kmer_size, max_kmer_size, ref_motifs_list)

sys.setrecursionlimit(original_limit)

all_seq_df = tag_all_seq(all_seq_dict, all_seq_masked, all_seq_masked_sporadic)

all_seq_masked.clear()
all_seq_masked_sporadic.clear()

unique_motifs = pd.concat([all_seq_df[col] for col in all_seq_df.columns]).unique().tolist()
if np.nan in unique_motifs:
    unique_motifs.remove(np.nan)

printable_asii_list = generate_hex_chr_list()
motif_chr_dict = assign_chr_to_motif(unique_motifs, printable_asii_list, input_fasta_to_count)
all_seq_chr_dict = write_seq_in_hex_chr(all_seq_df, motif_chr_dict, input_fasta_to_count, random_num)
write_compressed_seq(all_seq_chr_dict, motif_chr_dict, input_fasta_to_count, random_num)

all_seq_distance_df = edit_distance_between_seq_byte(all_seq_chr_dict)
all_seq_chr_dict.clear()
alignment_score_matrix = get_motif_pairwise_distance(unique_motifs)
dimension_reduction_result = run_umap(alignment_score_matrix)


if analysis_type == "assembly":
    population = args.population
    population_metadata = pd.read_csv(population, sep = "\t")
    if run_msa == "True":
        msa = msa_with_characters(input_fasta_to_count, random_num)
        msa_df_to_plot, motif_length = prepare_for_plotting(msa, motif_chr_dict, dimension_reduction_result)
        plot_msa_df(msa_df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, input_fasta_to_count.replace(".fa", "_" + str(random_num) + "_msa" + ".png"), title, population_metadata, motif_length)
    elif run_msa == "False":
        df_to_plot = map_score_to_alignment(all_seq_df, dimension_reduction_result)
        plot_df(df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, input_fasta_to_count.replace(".fa", "_" + str(random_num) + ".png"), title, population_metadata)

elif analysis_type == "reads":
    if run_msa == "True":
        msa = msa_with_characters(input_fasta_to_count, random_num)
        msa_df_to_plot, motif_length = prepare_for_plotting(msa, motif_chr_dict, dimension_reduction_result)
        plot_msa_df_reads(msa_df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, input_fasta_to_count.replace(".fa", "_reads_" + str(random_num) + "_msa" + ".png"), title, motif_length)
    elif run_msa == "False":
        df_to_plot = map_score_to_alignment(all_seq_df, dimension_reduction_result)
        plot_df_reads(df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, input_fasta_to_count.replace(".fa", "_reads_" + str(random_num) + ".png"), title)
