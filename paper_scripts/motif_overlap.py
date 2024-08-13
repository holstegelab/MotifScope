import pandas as pd
import numpy as np

def calculate_intersection_ratio(set1, set2):
    c_A = []
    c_B = []
    for m in set1:
        c_m = min([m[i:] + m[:i] for i in range(len(m))])
        c_A += [c_m]
    for m in set2:
        c_m = min([m[i:] + m[:i] for i in range(len(m))])
        c_B += [c_m]

    intersection = set(c_A).intersection(set(c_B))
    ratio = len(intersection) / len(set(c_B))
    
    return ratio

def method_used(number):
    if str(number) != "nan":
        return 1
    
def find_smallest_units(strings):
    def smallest_unit(s):
        n = len(s)
        for length in range(1, n + 1):
            if n % length == 0:
                substring = s[:length]
                if substring * (n // length) == s:
                    return substring
        return s

    smallest_units = [smallest_unit(s) for s in strings]
    return smallest_units


#stats_file: merged output of edit_distance.py
stats = pd.read_csv(stats_file, sep = "\t")
stats = stats[stats["motif"] != "motif"]
results = []


method_motif_dict = {}
for index, row in stats.iterrows():
    loci = row["locus"]
    method = row["method"]
    motifs = row["motif"].split(',')
    motifs = find_smallest_units(motifs)
    
    if loci not in method_motif_dict:
        method_motif_dict[loci] = {}
    method_motif_dict[loci][method] = set(motifs)


for loci, method_dict in method_motif_dict.items():
    methods = list(method_dict.keys())
    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            method1 = methods[i]
            method2 = methods[j]
            
            if method1 in method_dict and method2 in method_dict:
                ratio1_to_2 = calculate_intersection_ratio(method_dict[method1], method_dict[method2])
                ratio2_to_1 = calculate_intersection_ratio(method_dict[method2], method_dict[method1])
                
                results.append({
                    "locus": loci,
                    "method1": method1,
                    "method2": method2,
                    "ratio1_to_2": ratio1_to_2,
                    "ratio2_to_1": ratio2_to_1
                })


results_df = pd.DataFrame(results)
