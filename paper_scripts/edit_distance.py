import pandas as pd
import os
from Bio import SeqIO
import Levenshtein
import numpy as np
import argparse
import re


def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as f:
        current_header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_header = line[1:]
                sequences[current_header] = []
            else:
                sequences[current_header].append(line)
    return sequences


def process_ms(ms_seq):
    motifs = set()
    seq = ""
    for motif in ms_seq:
        parts = motif.split()
        for part in parts:
            # Extract the motif and its count
            motif_part = ''.join(filter(str.isalpha, part))
            count_part = ''.join(filter(str.isdigit, part))
            motifs.add(motif_part)
            seq += motif_part*int(count_part)
    return list(motifs), seq

def translate_vamos(seq_in_charcter, ms_motif_dict_all):
    translated_seq = ""
    for hex in seq_in_charcter:
        if hex != '':
            if hex in ms_motif_dict_all:
                translated_seq += ms_motif_dict_all[hex]
        else:
            return ''
    return translated_seq

def parse_utr(filename):
    """Parse a FASTA file and return a dictionary of sequences"""
    sequences = {}
    for record in SeqIO.parse(filename, 'fasta'):
        sequences[record.description] = str(record.seq)
    return sequences

def extract_string_from_header(header):
    pattern = r"#Pat([^#]+)#"
    match = re.search(pattern, header)
    if match:
        return match.group(1)
    return None

def extract_substring_with_integers(header):
    pattern = r"<([^<>]+)>(\d+)"
    matches = re.findall(pattern, header)
    return matches

def modify_fasta_dict(original_dict):
    modified_dict = {}
    motif_dict = {}
    for header, sequence in original_dict.items():
        
        extracted_string = extract_string_from_header(header).strip(" ")
        matches = extract_substring_with_integers(extracted_string)
        extracted_string_new = ""
        for substring, repeat_count in matches:
            repeat_count = int(repeat_count)
            extracted_string_new += repeat_count * substring
        if extracted_string_new:
            modified_dict[header.split("#Annotation")[1].strip(" ")] = extracted_string_new
            motif_dict[header.split("#Annotation")[1].strip(" ")] = matches
    return modified_dict, motif_dict



parser = argparse.ArgumentParser(description='motifs and edit distance')
parser.add_argument('-b', '--bed', default = None, dest='bed',
                    type=str,
                    help='bed file')
parser.add_argument('-ms', '--motifscope_seq', default = None, dest='motifscope_seq',
                    type=str,
                    help='compressed representation from MotifScope MotifScope seq')
parser.add_argument('-t', '--trf_file', default = None, dest='trf_file',
                    type=str,
                    help='trf output')
parser.add_argument('-u', '--utr', default = None, dest='utr',
                    metavar="utr output", type=str,
                    help='uTR')
parser.add_argument('-s', '--seq_file', default = None, dest='seq_file',
                    type=str,
                    help='fasta file contains the two haplotype')
parser.add_argument('-vo', '--vamos_original', default = None, dest='vamos_original',
                    type=str,
                    help='vamos original output')
parser.add_argument('-ve', '--vamos_efficient', default = None, dest='vamos_efficient',
                    type=str,
                    help='vamos efficient output')



args = parser.parse_args()
bed = args.bed
motifscope = args.motifscope
motifscope_seq = args.motifscope_seq
trf_file = args.trf_file
utr = args.utr
seq_file = args.seq_file
vamos_original = args.vamos_original
vamos_efficient = args.vamos_efficient


coor = pd.read_csv(bed, sep = "\t", header = None)
original = pd.read_csv(vamos_original, comment = "#", sep = "\t")
chr = coor.iloc[0,0]
r_start = coor.iloc[0,1]
r_end = coor.iloc[0,2]
original = original[original["CHROM"] == chr]
original = original[original["POS"] == int(r_start)].transpose()

d_0_2 = pd.read_csv(vamos_efficient, comment = "#", sep = "\t")
d_0_2 = d_0_2[d_0_2["CHROM"] == chr]
d_0_2 = d_0_2[d_0_2["POS"] == int(r_start)].transpose()

all_seq = parse_fasta(seq_file)
h1 = {list(all_seq.keys())[0]: all_seq[list(all_seq.keys())[0]]}
h2 = {list(all_seq.keys())[1]: all_seq[list(all_seq.keys())[1]]}

#################
###MotifScope####
#################
sorted_keys = list(h1.keys()) + list(h2.keys())
sorted_all_seq = {key: all_seq[key] for key in sorted_keys}
all_seq_concat = ''.join(sorted_all_seq.values())

motif_scope_seq = parse_fasta(motifscope_seq)
motif_scope_final = []
translated_seq_dict = {}
ms_edit_distance_dict = {}
for header, motifs in motif_scope_seq.items():
    ms_motifs_seq, ms_seq = process_ms(motifs)
    motif_scope_final += ms_motifs_seq
    translated_seq_dict[header] = ms_seq
    ms_edit_distance_dict[header] = Levenshtein.distance(all_seq[header], translated_seq_dict[header])

ms_edit_distance = sum(value for value in ms_edit_distance_dict.values()) / (len(all_seq_concat)-2)
ms_motif_number = len(list(set(motif_scope_final)))

#################
#vamos original##
#################
if original.shape[1] != 0:
    alleles = original.iloc[7,0].split(";")[-1].split("=")[-1].split(",")
    all_motifs = original.iloc[7,0].split(";")[1].split("=")[-1].split(",")
    all_motifs_dict = {}
    for i in range(len(all_motifs)):
        all_motifs_dict[str(i)] = all_motifs[i]

    #vamos_original = [all_motifs_dict[i] for i in all_motifs_dict]

    hap_dict = {}
    for i in range(1, len(original.iloc[7,0].split(";")[-1].split("=")[-1].split(","))+1):
        hap_dict[str(i)] = original.iloc[7,0].split(";")[-1].split("=")[-1].split(",")[i-1]

    original = original.iloc[9:,:]
    original = original.sort_index()
    original.columns = ["hap"]
    vamos_all_hap_dict = {}
    vamos_all_seq_dict = {}
    all_actual_seq_dict = {}

    for i, r in original.iterrows():
        if "/" in r["hap"]:
            if i.split(".")[-1] == "1":
                if len(h1) == 1: 
                    vamos_all_hap_dict[list(h1.keys())[0]] = hap_dict[str(r["hap"].split("/")[0])]
            elif i.split(".")[-1] == "2":
                if len(h2) == 1: 
                    vamos_all_hap_dict[list(h2.keys())[0]] = hap_dict[str(r["hap"].split("/")[0])]
    

    for seq in vamos_all_hap_dict:
        vamos_all_seq_dict[seq] = []
        for m in vamos_all_hap_dict[seq].split("-"):
            print(m)
            vamos_all_seq_dict[seq] += [m]
    
    vamos_original_final = {}
    for seq in vamos_all_seq_dict:
        for m in vamos_all_seq_dict[seq]:
            if m in vamos_original_final:
                continue
            else:
                vamos_original_final[m] = all_motifs_dict[m]
    
    vamos_original_edit_distance_dict = {}
    translated_seq_dict = {}
    for seq in vamos_all_seq_dict:
        translated_seq_dict[seq] = translate_vamos(vamos_all_seq_dict[seq], vamos_original_final)
        vamos_original_edit_distance_dict[seq] = Levenshtein.distance(all_seq[seq], translated_seq_dict[seq])
    
    vamos_original_edit_distance = sum(value for value in vamos_original_edit_distance_dict.values()) / (len(all_seq_concat))
    vamos_original_motif_number = len(vamos_original_final)


#################
#vamos efficient#
#################
if d_0_2.shape[1] != 0:
    alleles = d_0_2.iloc[7,0].split(";")[-1].split("=")[-1].split(",")
    all_efficient_motifs = d_0_2.iloc[7,0].split(";")[1].split("=")[-1].split(",")
    all_efficient_motifs_dict = {}
    for i in range(len(all_efficient_motifs)):
        all_efficient_motifs_dict[str(i)] = all_efficient_motifs[i]
    hap_dict = {}
    for i in range(1, len(d_0_2.iloc[7,0].split(";")[-1].split("=")[-1].split(","))+1):
        hap_dict[str(i)] = d_0_2.iloc[7,0].split(";")[-1].split("=")[-1].split(",")[i-1]

    d_0_2 = d_0_2.iloc[9:,:]
    d_0_2 = d_0_2.sort_index()
    d_0_2.columns = ["hap"]
    vamos_all_hap_dict = {}
    vamos_all_seq_dict = {}
    all_actual_seq_dict = {}

    for i, r in d_0_2.iterrows():
        if "/" in r["hap"]:
            if i.split(".")[-1] == "1":
                if len(h1) == 1: 
                    vamos_all_hap_dict[list(h1.keys())[0]] = hap_dict[str(r["hap"].split("/")[0])]
            elif i.split(".")[-1] == "2":
                if len(h2) == 1: 
                    vamos_all_hap_dict[list(h2.keys())[0]] = hap_dict[str(r["hap"].split("/")[0])]
    
    for seq in vamos_all_hap_dict:
        vamos_all_seq_dict[seq] = []
        for m in vamos_all_hap_dict[seq].split("-"):
            vamos_all_seq_dict[seq] += [m]
    
    vamos_efficient_final = {}
    for seq in vamos_all_seq_dict:
        for m in vamos_all_seq_dict[seq]:
            if m in vamos_efficient_final:
                continue
            else:
                vamos_efficient_final[m] = all_motifs_dict[m]
    
    vamos_efficient_edit_distance_dict = {}
    translated_seq_dict = {}
    for seq in vamos_all_seq_dict:
        translated_seq_dict[seq] = translate_vamos(vamos_all_seq_dict[seq], vamos_efficient_final)
        vamos_efficient_edit_distance_dict[seq] = Levenshtein.distance(all_seq[seq], translated_seq_dict[seq])
    
    vamos_efficient_edit_distance = sum(value for value in vamos_efficient_edit_distance_dict.values()) / (len(all_seq_concat))
    vamos_efficient_motif_number = len(vamos_efficient_final)


#################
#######utr#######
#################
utr_seq = parse_utr(utr)
utr_dict, utr_motifs = modify_fasta_dict(utr_seq)
utr_dict = {key: utr_dict[key] for key in sorted_keys}
utr_motifs = {key: utr_motifs[key] for key in sorted_keys}

utr_final = []
for seq in utr_motifs:
    for m in utr_motifs[seq]:
        utr_final += [m[0]]
utr_final = list(set(utr_final))

utr_edit_distance_dict = {}
for seq in utr_dict:
    utr_edit_distance_dict[seq] = Levenshtein.distance(all_seq[seq], utr_dict[seq])


utr_edit_distance = sum(value for value in utr_edit_distance_dict.values()) / (len(all_seq_concat))
utr_motif_number = len(utr_final)


#################
#######trf#######
#################
trf_sorted_keys = [key.split(":")[-1] for key in sorted_keys]
trf_sorted_all_seq = {}
for seq in sorted_all_seq:
    trf_sorted_all_seq[seq.split(":")[-1]] = sorted_all_seq[seq]

output_fname = trf_file.replace(".txt", ".processed.txt")
#output_fname = "/trf_new/%s_new.txt"
f = open(output_fname, "w")
f.write("read_name\ttr_id\tstart\tend\ttr_length\tmotif_size\tcopy_number\tconsensus_size\tpercent_matches\tpercent_indels\tscore\tA\tC\tG\tT\tentropy\tmotif\tsequence\n")
file = open(trf_file, 'r')
#file = open(os.path.join(path, fname), 'r')
lines = file.readlines()
file.close()
n = 0
while n < len(lines):
    if "@" in lines[n]:
        m = n
        n += 1
    elif "@" not in lines[n]:
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (lines[m].split(":")[-1].strip(), lines[m].split(" ")[0][1:].strip(), int(lines[n].split(" ")[0]), int(lines[n].split(" ")[1]), int(lines[n].split(" ")[1]) - int(lines[n].split(" ")[0]) + 1,\
                                                                                            lines[n].split(" ")[2], lines[n].split(" ")[3], lines[n].split(" ")[4], lines[n].split(" ")[5], lines[n].split(" ")[6], \
                                                                                                lines[n].split(" ")[7], lines[n].split(" ")[8], lines[n].split(" ")[9], \
                                                                                                    lines[n].split(" ")[10], lines[n].split(" ")[11], lines[n].split(" ")[12], \
                                                                                                        lines[n].split(" ")[13], lines[n].split(" ")[14].strip()))
        n += 1
f.close()
    
trf = pd.read_csv(output_fname, sep = "\t")
os.system("rm %s" %(output_fname))
trf_edit_distance = 0
if trf.shape[0] != 0:
    trf = trf[trf.read_name.isin(trf_sorted_keys)]
    trf_sum = trf.groupby(by = "motif")["copy_number"].sum().reset_index().sort_values(by = "copy_number", ascending = False)

    trf_final = trf_sum["motif"].tolist()
    trf_motif_number = len(trf_final)
    trf_edit_distance_dict = {}
    translated_seq_dict = {}
    trf_edit_distance_per_seq_dict = {}
    for index, row in trf.iterrows():
        trf_motif_translated = row["motif"] * (int(row["copy_number"]) + 1)
        actual_trf_motif_translated = trf_motif_translated[:int(row["end"])-int(row["start"]) + 1]
        translated_distance = Levenshtein.distance(trf_sorted_all_seq[row["read_name"]][int(row["start"])-1:int(row["end"])], actual_trf_motif_translated)

        if row["read_name"] in translated_seq_dict:
            translated_seq_dict[row["read_name"]] += [actual_trf_motif_translated]
            trf_edit_distance_per_seq_dict[row["read_name"]] += translated_distance
        else:
            translated_seq_dict[row["read_name"]] = [actual_trf_motif_translated]
            trf_edit_distance_per_seq_dict[row["read_name"]] = translated_distance

    all_trf_translated_seq_concat = ''.join([''.join(x) for x in translated_seq_dict.values()])
    trf_edit_distance = sum(x for x in trf_edit_distance_per_seq_dict.values()) / len(all_seq_concat)


ms = pd.DataFrame({"motif": [",".join(motif_scope_final.values())], "edit_distance": ms_edit_distance, "motif_number": ms_motif_number, "method": "MotifScope"})
utr_results = pd.DataFrame({"motif": [",".join(utr_final)], "edit_distance": utr_edit_distance, "motif_number": utr_motif_number, "method": "uTR"})

if original.shape[1] != 0:
    vamos_o = pd.DataFrame({"motif": [",".join(vamos_original_final.values())], "edit_distance": vamos_original_edit_distance, "motif_number": vamos_original_motif_number, "method": "vamos original"})
    vamos_e = pd.DataFrame({"motif": [",".join(vamos_efficient_final.values())], "edit_distance": vamos_efficient_edit_distance, "motif_number": vamos_efficient_motif_number, "method": "vamos efficient"})
else:
    vamos_o = pd.DataFrame()
    vamos_e = pd.DataFrame()
if trf.shape[0] != 0:
    trf_results = pd.DataFrame({"motif": [",".join(trf_final)], "edit_distance": trf_edit_distance, "motif_number": utr_motif_number, "method": "TRF"})
else:
    trf_results = pd.DataFrame()
df = pd.concat([ms, vamos_o, vamos_e, trf_results, utr_results], ignore_index=True, axis = 0)
df["locus"] = "%s:%s-%s" %(chr, r_start, r_end)
df.to_csv("%s.%s.%s.stats.all.txt" % (chr, r_start, r_end), sep = "\t", index = False)
