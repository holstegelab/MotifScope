import pandas as pd
import numpy as np
import os
import argparse
import sys
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
from suffix_trees import STree
import matplotlib.gridspec as gridspec

class Kmer:
    def __init__(self, kmer, count, seq, nseq, start_pos_list):
        self.kmer = kmer
        self.count = count
        self.seq = seq
        self.nseq = nseq
        self.perm_kmer = [self.kmer[x:] + self.kmer[:x] for x in range(len(self.kmer))]
        self.kmer_count = {(nseq, kmer): count}
        self.sum_count = count
        self.start_pos_list = start_pos_list
        indel_size, indel_start, indel_end = get_indel_size(self.kmer, self.start_pos_list)
        self.max_kmer_mask = {(nseq, kmer): max(indel_size)}
        indel_pos_list = [(indel_start[i], indel_end[i], self.kmer) for i in range(len(indel_start))]
        self.consecutive_region = {(nseq, kmer): indel_pos_list}
        self.sum_kmer_consecutive_mask = {(nseq, kmer): sum(y for y in indel_size if y > len(self.kmer))}
        self.sum_kmer_mask = {(nseq, kmer): sum(indel_size)}
                
        
    def is_equal(self, another_kmer):
        if another_kmer in self.perm_kmer:
            return True
        else:
            return False
    
    def append_kmer_set(self, another_kmer, another_count, another_start_pos_list):
        if self.is_equal(another_kmer):
            self.kmer_count[(self.nseq, another_kmer)] = another_count
            another_kmer_indel_size, another_kmer_indel_start, another_kmer_indel_end = get_indel_size(another_kmer, another_start_pos_list)
            self.max_kmer_mask[(self.nseq, another_kmer)] = max(another_kmer_indel_size)
            indel_pos_list = [(another_kmer_indel_start[i], another_kmer_indel_end[i], another_kmer) for i in range(len(another_kmer_indel_start))]
            self.consecutive_region[(self.nseq, another_kmer)] = indel_pos_list
            self.sum_kmer_consecutive_mask[(self.nseq, another_kmer)] = sum(y for y in another_kmer_indel_size if y > len(another_kmer))
            self.sum_kmer_mask[(self.nseq, another_kmer)] = sum(another_kmer_indel_size)
            self.sum_count += another_count
    
    def merge_seq(self, obj):
        self.kmer_count = {**self.kmer_count, **obj.kmer_count}
        self.max_kmer_mask = {**self.max_kmer_mask, **obj.max_kmer_mask}
        self.consecutive_region = {**self.consecutive_region, **obj.consecutive_region}
        self.sum_kmer_consecutive_mask = {**self.sum_kmer_consecutive_mask, **obj.sum_kmer_consecutive_mask}
        self.sum_kmer_mask = {**self.sum_kmer_mask, **obj.sum_kmer_mask}
        self.sum_count += obj.sum_count



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
            sequences[record.id] = str(record.seq)
    return sequences


def get_all_seq_info(all_seq):
    read_name_df = pd.DataFrame()
    nseq = 0
    for record in all_seq:
        read_name_df = pd.concat([read_name_df, pd.DataFrame({"nseq": "n" + str(nseq), "read_name": record}, index = [nseq])])
        nseq += 1
    return read_name_df

def process_record(fasta_file, k_min, k_max, record, seq_number):
    #all_rounds_kmer = pd.DataFrame()
    #seq_used_motif = pd.DataFrame()
    #all_kmer_df = pd.DataFrame()
    if str(record.seq) != "":
        #logging.debug(record.id)
        tmp_file = fasta_file.replace(".fa", "_") + "tmp_seq_" + str(seq_number) + ".fa"
        with open (tmp_file, 'w') as tmp_fasta:
            SeqIO.write(record, tmp_fasta, 'fasta')
            
        seq_list = list(SeqIO.parse(tmp_file, "fasta"))
        
        output_cnt_file = tmp_file.replace(".fa", "." + str(k_min) + "." + str(k_max) + "mer.txt")
        os.system("~/bin/aardvark/build/aardvark kc -k %d -m %d %s > %s " %(min(len(record.seq) - 1, k_min), min(len(record.seq), k_max), tmp_file, output_cnt_file))
        kmer_df = pd.read_csv(output_cnt_file, header = None, sep = "\t")
        kmer_df.columns = ["ckmer", "cnt"]
        kmer_df["seq"] = record.id
        kmer_df["nseq"] = seq_number
    
        return kmer_df
    else:
        return pd.DataFrame()
    
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


def find_positions(input_string, substring_list):
    st = STree.STree(input_string)
    kmer_start = {}
    for substring in substring_list:
        kmer_start[substring] = st.find_all(substring)
    return kmer_start


def calculation_per_kmer_pool(seq, seq_number, kmer_df):

    kmer_df = kmer_df[kmer_df["nseq"] == seq_number]
    kmer_counts = dict(zip(kmer_df.ckmer, kmer_df.cnt))
    kmer_list = kmer_df.loc[:,"ckmer"].tolist()
    all_kmer = {}
    for kmer in kmer_list:
        for x in range(len(kmer)):
            all_kmer[kmer[x:] + kmer[:x]] = kmer
   
    if seq == '':
        return []
    else:
        start_pos_dict = find_positions(seq, list(all_kmer.keys()))
        kmer_set_dict = {}
        for kmer in all_kmer:
            if kmer in seq:
                if all_kmer[kmer] in kmer_set_dict:
                    kmer_set_dict[all_kmer[kmer]].append_kmer_set(kmer, kmer_counts[all_kmer[kmer]], start_pos_dict[kmer])
                else:
                    kmer_set_dict[all_kmer[kmer]] = Kmer(kmer, kmer_counts[all_kmer[kmer]], seq, seq_number, start_pos_dict[kmer])
                
        kmer_set_list = [kmer for kmer in set(kmer_set_dict.values())]
        return sorted(kmer_set_list, key = lambda x: max(x.max_kmer_mask.values()), reverse = True)
        
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
    return sorted(kmer_all_seq, key = lambda x: max(x.max_kmer_mask.values()), reverse = True)

def pairwise_alignment(seq1, seq2):
    from Bio import pairwise2
    alignments = pairwise2.align.globalms(seq1, seq2, 0, -2, -5, -5)
    #alignments = pairwise2.align.globalxx(seq1, seq2)
    scores = []
    for alignment in alignments:
        scores.append(alignment.score)
    #return max(scores) / len(alignments[scores.index(max(scores))].seqA)
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
            seq = seq[:region[0]] + "x" * (region[1] - region[0]) + seq[region[1]:]
            masked_region += [region]
    return seq, {(seq_number, kmer): masked_region}

def sporadic_mask_seq(seq, seq_number, kmer, kmer_obj):
    regions = kmer_obj.consecutive_region[(seq_number, kmer)]
    masked_region = []
    for region in regions:
        if region[1] - region[0] == len(kmer):
            seq = seq[:region[0]] + "x" * (region[1] - region[0]) + seq[region[1]:]
            masked_region += [region]
    return seq, {(seq_number, kmer): masked_region}



def mask_all_seq(all_seq, df):
    n = 0
    seq_to_be_masked = [str(seq) for seq in all_seq.values()]
    seq_list = seq_to_be_masked
    selected_kmer_list = []
    all_seq_masked_dict = {}
    while max([len(x) for x in seq_list]) != 1:
        #logging.debug("start calculation")
        print("start calculation")
        kmer_set_per_seq = pool.starmap(calculation_per_kmer_pool, zip(seq_to_be_masked, list(range(len(seq_to_be_masked))), repeat(df)))
        #logging.debug("calculation done")
        print("calculation done")
        #logging.debug("merging kmer object")
        #kmer_set_all_seq = select_kmer(kmer_set_per_seq) #can be faster
        kmer_set_all_seq = select_kmer_fast(kmer_set_per_seq) #can be faster
        #logging.debug("merging finished")
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
                    #logging.debug("getting best perm kmer")
                    kmer = get_best_perm_kmer(selected_kmer, selected_kmer_list)
                    #logging.debug("best kmer selected")
                #logging.debug(kmer)
                print(kmer)
                #print(sum_max_perm_kmer_region[kmer])
                selected_kmer_list += [kmer]
                one_kmer_masked = {}
                for region in selected_kmer.max_kmer_mask:
                    if region[1] == kmer:
                        #print(region)
                        seq_to_be_masked[region[0]], masked_region = consecutive_mask_seq(seq_to_be_masked[region[0]], region[0], kmer, selected_kmer)
                        one_kmer_masked = {**one_kmer_masked, **masked_region}
                        
                if all(value == [] for value in one_kmer_masked.values()):
                    break
                
                else:
                    all_seq_masked_dict = {**all_seq_masked_dict, **one_kmer_masked}
                seq_list = [x.split("x") for x in seq_to_be_masked]
                seq_list = [item for sublist in seq_list for item in sublist]
                n += 1
        
    all_seq_masked_sporadic_seq = {}
    while max([len(x) for x in seq_list]) != 1:
        print("finding unconsecutive motifs")
        kmer_set_per_seq = []
        kmer_round = []
        kmer_set_per_seq = pool.starmap(calculation_per_kmer_pool, zip(seq_to_be_masked, list(range(len(seq_to_be_masked))), repeat(df)))
        for used_one_kmer in selected_kmer_list:
            one_kmer_masked = {}
            for kmer_set in kmer_set_per_seq:
                for one_kmer in kmer_set:
                    if one_kmer.is_equal(used_one_kmer):
                        selected_kmer = one_kmer
                        
                        #logging.debug(used_one_kmer)
                        print(kmer)
                        for region in selected_kmer.max_kmer_mask:
                            if region[1] == used_one_kmer:
                                seq_to_be_masked[region[0]], masked_region = sporadic_mask_seq(seq_to_be_masked[region[0]], region[0], used_one_kmer, selected_kmer)
                                if masked_region.values() != [[]]:
                                    kmer_round += [used_one_kmer]
                                    one_kmer_masked = {**one_kmer_masked, **masked_region}
                                
        
                        if all(value == [] for value in one_kmer_masked.values()):
                            break

           
            seq_list = [x.split("x") for x in seq_to_be_masked]
            seq_list = [item for sublist in seq_list for item in sublist]
            if one_kmer_masked != {}:
                all_seq_masked_sporadic_seq = {**all_seq_masked_sporadic_seq, **one_kmer_masked}
                break
                        
                        
        if kmer_round == []:
            break
        
    
    return all_seq_masked_dict, all_seq_masked_sporadic_seq

def tag_all_seq(all_seq, all_seq_masked_dict, all_seq_masked_sporadic_seq):
    seq_df = pd.DataFrame()
    i = 0
    for seq in all_seq:
        seq_array = pd.DataFrame(list(str(all_seq[seq])), columns = ["n" + str(i)])
        masked_region_one_seq = [all_seq_masked_dict[region] for region in all_seq_masked_dict if region[0] == i]
        #masked_region_one_seq = [new_all_seq_masked[region] for region in new_all_seq_masked if region[0] == i]
        masked_region_one_seq = [item for sublist in masked_region_one_seq for item in sublist]
        for region in masked_region_one_seq:
            seq_array.iloc[region[0]:region[1], 0] = region[2]
        masked_region_one_seq = [all_seq_masked_sporadic_seq[region] for region in all_seq_masked_sporadic_seq if region[0] == i]
        #masked_region_one_seq = [new_all_seq_masked[region] for region in new_all_seq_masked if region[0] == i]
        masked_region_one_seq = [item for sublist in masked_region_one_seq for item in sublist]
        for region in masked_region_one_seq:
            seq_array.iloc[region[0]:region[1], 0] = region[2]
            
        seq_df = pd.concat([seq_df, seq_array], axis = 1)
        i += 1
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


def write_seq_in_hex_chr(all_seq_motif_df, motif_dict, read_name_df, input_fasta_file_name, random_num):
    all_seq_motif_chr_df = all_seq_motif_df.replace(motif_dict)
    # write into fasta file
    chr_fasta_file = input_fasta_file_name.replace(".fa", "_motif_in_hex_" + str(random_num) + ".hex")
    seq_chr_dict = {}    
    with open (chr_fasta_file, 'w') as chr_fasta:
        for index, row in read_name_df.iterrows():
            seq_description = row["read_name"]
            seq_motif_in_chr = all_seq_motif_chr_df.loc[:, row["nseq"]].dropna()
            seq_motif = all_seq_motif_df.loc[:, row["nseq"]].dropna()
            seq_chr = []
            i = 0
            while i < len(seq_motif):
                seq_chr += [seq_motif_in_chr[i]]
                i += len(seq_motif[i])
            

            seq_chr_string = " ".join(seq_chr)
            chr_fasta.write(">" + seq_description + "\n")
            chr_fasta.write(seq_chr_string + "\n")
            seq_chr_dict[row["nseq"]] = seq_chr_string
            
    motif_fasta_file = input_fasta_file_name.replace(".fa", "_motif_to_hex_" + str(random_num) + ".fa")
    with open (motif_fasta_file, 'w') as dict_fasta:
        for item in motif_dict:
            dict_fasta.write(">" + motif_dict[item] + "\n")
            dict_fasta.write(item + "\n")        
    return seq_chr_dict


def hex_to_bytes(hex_string):
    """Convert a hexadecimal string to bytes"""
    return binascii.unhexlify(hex_string)


def normalize_edit_distance(distance, seq_1_length, seq_2_length):
    if max(seq_1_length, seq_2_length) != 0:
        return distance / max(seq_1_length, seq_2_length)
    else:
        return 0


def edit_distance_between_seq_byte(input_seq_chr_dict):

    #seq_distance_df = pd.DataFrame(list(combinations_with_replacement(list(all_seq_chr_dict.keys()), 2)), columns=['nseq_1', 'nseq_2'])
    seq_chr_dict = {}
    for seq in input_seq_chr_dict:
        if input_seq_chr_dict[seq] == "":
            seq_chr_dict[seq] = str.encode(" ")
        else:
            seq_chr_dict[seq] = hex_to_bytes(input_seq_chr_dict[seq].replace(" ", ""))
    
    seq_distance_df = pd.DataFrame(list(combinations_with_replacement(list(seq_chr_dict.keys()), 2)), columns=['nseq_1', 'nseq_2'])
    seq_distance_df["seq_1"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_1"]], axis = 1)
    seq_distance_df["seq_2"] = seq_distance_df.apply(lambda x: seq_chr_dict[x["nseq_2"]], axis = 1)
    #seq_distance_df['tuple'] = list(zip(seq_chr_dict[seq_distance_df.seq_1], seq_chr_dict[seq_distance_df.seq_2]))
    pool = Pool()
    
    ls = pool.starmap(Levenshtein.distance, zip(seq_distance_df.seq_1, seq_distance_df.seq_2))
    seq_distance_df['distance'] = ls
    seq_distance_df['distance'] = seq_distance_df.apply(lambda x: normalize_edit_distance(x["distance"], len(x["seq_1"]), len(x["seq_2"])), axis = 1)

    #seq_distance_df['distance'] = seq_distance_df.apply(lambda x: x["distance"] / len(x["seq_1"]) if len(x["seq_1"]) >= len(x["seq_2"]) else x["distance"] / len(x["seq_2"]), axis = 1)
    seq_distance_df = seq_distance_df.drop(columns = ["seq_1", "seq_2"])
    #seq_distance_df = seq_distance_df.drop(columns = ['tuple'])
    
    #score_df["score"] = score_df.apply(lambda x: pairwise_alignment(x["motif_1"], x["motif_2"]), axis = 1)
    seq_distance_df_2 = seq_distance_df.copy()
    seq_distance_df_2["nseq_2"] = seq_distance_df["nseq_1"]
    seq_distance_df_2["nseq_1"] = seq_distance_df["nseq_2"]
    seq_distance_df = pd.concat([seq_distance_df, seq_distance_df_2], ignore_index = True)
    seq_distance_df = seq_distance_df.drop_duplicates().pivot(index = 'nseq_1', columns = 'nseq_2', values = 'distance').astype(float)

    return seq_distance_df

def get_motif_pairwise_distance(motif_list):
    score_df = pd.DataFrame(list(combinations_with_replacement(motif_list, 2)), columns=['motif_1', 'motif_2'])
    
    score_df['tuple'] = list(zip(score_df.motif_1, score_df.motif_2))
    pool = Pool()
    ls = pool.starmap(pairwise_alignment, score_df['tuple'].tolist())
    score_df['score'] = ls
    score_df = score_df.drop(columns = ['tuple'])
    
    #score_df["score"] = score_df.apply(lambda x: pairwise_alignment(x["motif_1"], x["motif_2"]), axis = 1)
    score_df_2 = score_df.copy()
    score_df_2["motif_2"] = score_df["motif_1"]
    score_df_2["motif_1"] = score_df["motif_2"]
    score_df = pd.concat([score_df, score_df_2], ignore_index = True)
    score_df = score_df.drop_duplicates().pivot(index = 'motif_1', columns = 'motif_2', values = 'score').astype(float)

    return score_df

def run_umap(dm):
        
    #umap_f = umap.UMAP(n_components = 1, metric = "manhattan", n_neighbors = 15, random_state = 0)
    umap_f = umap.UMAP(n_components = 1, metric = "manhattan", n_neighbors = 10, min_dist = 0.5, random_state = 0)
    #tsne = TSNE(random_state = 0, n_components = 1)
    X_transform_L2 = pd.DataFrame(umap_f.fit_transform(dm))
    X_transform_L2.columns = ["dimension_reduction"]
    X_transform_L2["motif"] = dm.index
    X_transform_L2 = X_transform_L2.astype({"dimension_reduction": float})

    return X_transform_L2

def map_score_to_alignment(all_seq_motifs, X_transform_L2, reads_name_df):
    motif_dict = dict(zip(X_transform_L2.motif, map(float, X_transform_L2.dimension_reduction)))
    motif_dict['nan'] = np.nan
    df = all_seq_motifs.applymap(lambda x: motif_dict[str(x)])
    df = df.rename(columns = dict(zip(reads_name_df.nseq, reads_name_df.read_name)))
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

def plot_df(df, dm, all_seq_motifs, seq_distance_df, reads_name_df, figname, figtitle, population_df):
    fig = plt.figure(figsize=(max(50, 0.015 * df.shape[1]), 0.5 * df.shape[0]), dpi = 300)
    spec = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[4,1,40,10], height_ratios=[1], wspace=0.02)

    for col in range(3):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    summary = all_seq_motifs.apply(lambda x: x.value_counts(), axis = 0).transpose().reset_index()
    summary["index"] = summary.apply(lambda x: x["index"].split("_")[-1], axis = 1)
    summary = summary.merge(reads_name, left_on = "index", right_on = "nseq", how = "left").drop(columns = ["index", "nseq"]).set_index("read_name").fillna(0)
    seq_distance_df = seq_distance_df.reset_index().merge(reads_name, left_on = "nseq_1", right_on = "nseq", how = "left").drop(columns = ["nseq_1", "nseq"]).set_index("read_name")
    clusters = shc.linkage(seq_distance_df, method='average', metric="euclidean")

    d = shc.dendrogram(Z = clusters, ax = all_axes[0], labels = seq_distance_df.index, orientation = "left")
    all_axes[0].tick_params(right=False, left = False, top=False, labelright=False, labelleft=False,labeltop=False)
    seq_order_list = list(reversed(d["ivl"]))

    '''df = df.sort_index().reset_index()
    df = df.set_index("read_name")
    df = df.sort_index()'''
    df = df.reindex(seq_order_list)
    all_seq_motifs = all_seq_motifs.transpose()
    all_seq_motifs["nseq"] = all_seq_motifs.index.to_series().apply(lambda x: x.split("_")[-1])
    all_seq_motifs = all_seq_motifs.merge(reads_name_df, on = "nseq", how = "left").drop(columns = ["nseq"]).set_index("read_name")
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
        #print(i)
        j = 0
        while j < all_seq_motifs.shape[0]:
            if str(all_seq_motifs.iloc[j, i]) != "nan":
                all_axes[2].add_patch(Rectangle((j, i), len(all_seq_motifs.iloc[j, i]), 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
                j += len(all_seq_motifs.iloc[j, i])
                
                #print(j)
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

    
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")


def msa_with_characters(input_fasta_file_name, random_num):
    chr_fasta_file = input_fasta_file_name.replace(".fa", "_motif_in_hex_" + str(random_num) + ".hex")
    os.system("~/miniconda3/envs/py38/libexec/mafft/hex2maffttext %s > %s" %(chr_fasta_file, chr_fasta_file.replace(".hex", ".ASCII")))
    os.system("mafft --text --op 2.0 --ep 0.1 %s > %s" %(chr_fasta_file.replace(".hex", ".ASCII"), chr_fasta_file.replace(".hex", "_mafft_output.ASCII")))
    os.system("~/miniconda3/envs/py38/libexec/mafft/maffttext2hex %s > %s" % (chr_fasta_file.replace(".hex", "_mafft_output.ASCII"), chr_fasta_file.replace(".hex", "_mafft_output.hex")))
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
    #msa_adjusted_df = msa_adjusted_df.reindex(seq_order_list)
    msa_adjusted_df = msa_adjusted_df.applymap(lambda x: motif_dict_r[x])
    msa_adjusted_df.columns = range(msa_adjusted_df.columns.size)

    score_dict = dict(zip(dm.motif, dm.dimension_reduction))
    score_dict[''] = np.nan
    df = msa_adjusted_df.applymap(lambda x: score_dict[x])
    
    return df, max_motif_length
    

def plot_msa_df(df, dm, all_seq_motifs, seq_distance_df, figname, figtitle, population_df, msa_result, max_motif_length):

    fig = plt.figure(figsize=(max(50, 0.015 * sum(max_motif_length)), 0.5 * df.shape[0]), dpi = 300)
    spec = fig.add_gridspec(ncols=4, nrows=1, width_ratios=[4,1,40,10], height_ratios=[1], wspace=0.01)

    for col in range(3):
        axs = fig.add_subplot(spec[0, col])
    all_axes = fig.get_axes()    

    seq_distance_df = seq_distance_df.reset_index().merge(reads_name, left_on = "nseq_1", right_on = "nseq", how = "left").drop(columns = ["nseq_1", "nseq"]).set_index("read_name")
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

    all_seq_motifs = all_seq_motifs.transpose()

    for i in range(all_seq_motifs.shape[1]):
        #print(i)
        j = 0
        m = 0
        while j < all_seq_motifs.shape[0]:
            if str(all_seq_motifs.iloc[j, i]) not in ["nan", ""]:
                all_axes[2].add_patch(Rectangle((j, i), max_motif_length[m], 1, fill=False, edgecolor='black', lw=0.05, clip_on=False))
            j += max_motif_length[m]
            m += 1
                

            #logging.debug(j)


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

    
    plt.yticks(rotation=0)
    plt.savefig(figname, bbox_inches = "tight")



#######################


parser = argparse.ArgumentParser(description='MotifScope')

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

parser.add_argument('-p', '--population', default = None, dest='population',
                    metavar="metadata.txt", type=str,
                    help='population metadata file')

parser.add_argument('-e', '--editdistance', default = None, dest='edit_distance',
                    metavar=2, type=int,
                    help='edit distance to allow in motif screening')

args = parser.parse_args()
title = args.title
input_fasta_to_count = args.input_fasta_to_count
max_kmer_size = args.max_kmer_size
min_kmer_size = args.min_kmer_size

population = args.population
max_edit_distance = args.edit_distance

random_num = str(random()).replace('.', '')
all_seq = list(SeqIO.parse(input_fasta_to_count, "fasta"))
all_seq_dict = parse_fasta(input_fasta_to_count)
reads_name = get_all_seq_info(all_seq_dict)
population_metadata = pd.read_csv(population, sep = "\t")
#logging.basicConfig(filename='%s.test.msa.e%s.kmer.2.%s.out' % (title, str(max_edit_distance), str(max_kmer_size)), level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

pool = Pool()
#logging.debug("start counting")
kmer_counts_per_seq = pool.starmap(process_record, zip(repeat(input_fasta_to_count), repeat(min_kmer_size), repeat(max_kmer_size), all_seq, list(range(len(all_seq)))))
#dfs = pool.starmap(process_record, zip(repeat(input_fasta_to_count), repeat(kmer), all_seq, list(range(len(all_seq)))))
all_kmer_counts = pd.concat([kmer_counts for kmer_counts in kmer_counts_per_seq], axis = 0, ignore_index = True)
#logging.debug("done")
all_kmer_counts["smallest_kmer"] = all_kmer_counts.apply(lambda x: is_smallest_kmer(x[0]), axis = 1)
all_kmer_counts = all_kmer_counts[all_kmer_counts["smallest_kmer"] == True].drop(columns = ["smallest_kmer"])
all_kmer_counts["sum_cnt"] = all_kmer_counts["cnt"].groupby(all_kmer_counts["ckmer"]).transform('sum')

all_seq_masked, all_seq_masked_sporadic = mask_all_seq(all_seq_dict, all_kmer_counts)
all_seq_df = tag_all_seq(all_seq_dict, all_seq_masked, all_seq_masked_sporadic)

unique_motifs = pd.concat([all_seq_df[col] for col in all_seq_df.columns]).unique().tolist()
if np.nan in unique_motifs:
    unique_motifs.remove(np.nan)

printable_asii_list = generate_hex_chr_list()
motif_chr_dict = assign_chr_to_motif(unique_motifs, printable_asii_list, input_fasta_to_count)
all_seq_chr_dict = write_seq_in_hex_chr(all_seq_df, motif_chr_dict, reads_name, input_fasta_to_count, random_num)
all_seq_distance_df = edit_distance_between_seq_byte(all_seq_chr_dict)
alignment_score_matrix = get_motif_pairwise_distance(unique_motifs)
dimension_reduction_result = run_umap(alignment_score_matrix)

### no msa
df_to_plot = map_score_to_alignment(all_seq_df, dimension_reduction_result, reads_name)
seq_order = summarize_motif(all_seq_df, reads_name, input_fasta_to_count, random_num, all_seq_distance_df, title)
plot_df(df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, reads_name, input_fasta_to_count.replace(".fa", "_umap_" + str(random_num) + "_" + str(max_edit_distance) + ".png"), title, population_metadata)


### msa
msa = msa_with_characters(input_fasta_to_count, random_num)
msa_df_to_plot, motif_length = prepare_for_plotting(msa, motif_chr_dict, dimension_reduction_result)
plot_msa_df(msa_df_to_plot, dimension_reduction_result, all_seq_df, all_seq_distance_df, input_fasta_to_count.replace(".fa", "_umap_" + str(random_num) + "_msa_" + str(max_edit_distance) + ".png"), title, population_metadata, msa, motif_length)
