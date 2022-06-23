from Bio import SeqIO
import numpy as np

from tqdm import tqdm, trange

import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import csv

from multiprocessing import Pool

with open(snakemake.input[0]) as reference_h, open(snakemake.input[1]) as quasispecies_h:
    reference: str = str(list(SeqIO.parse(reference_h, "fasta"))[0].seq)
    reference_length = len(reference)
    quasispecies: list = []
    for record in tqdm(SeqIO.parse(quasispecies_h, "fasta")):
        quasispecies.append(record.seq)
        if len(record.seq) != reference_length:
            raise Exception("Not same length")

reference_gapless_length = len(reference.replace("-",""))

def get_maximum_number_contiguous_gaps():
    cur_number = 0
    maximum_number = 0
    for nuc in reference:
        if nuc == "-":
            cur_number+=1
            continue
        if cur_number > maximum_number:
            maximum_number = cur_number
            cur_number = 0
    return maximum_number

maximum_number_insertions = get_maximum_number_contiguous_gaps() + 1

def group_by_position(cur_quasispecies):
    grouped = np.empty(reference_gapless_length, dtype=f"U{maximum_number_insertions}")
    cur_pos_object = ""
    cur_pos_done = False
    gapless_pos = 0
    for gapped_pos,(ref_n, quas_n) in enumerate(zip(reference, cur_quasispecies)):
        if ref_n == "-" and quas_n == "-":
            pass
        elif ref_n == "-":
            cur_pos_object += quas_n
            if gapped_pos == len(reference) - 1 or reference[gapped_pos + 1] != "-":
                cur_pos_done = True
        elif quas_n == "-":
            cur_pos_object = ""
            cur_pos_done = True
        else:
            if cur_pos_object != "":
                grouped[gapless_pos] = cur_pos_object
                cur_pos_object = ""
                cur_pos_done = False
                gapless_pos+=1
            
            cur_pos_object = quas_n
            if gapped_pos == len(reference) - 1 or (reference[gapped_pos + 1] != "-"):
                cur_pos_done = True
        
        if cur_pos_done:
            grouped[gapless_pos] = cur_pos_object
            cur_pos_object = ""
            cur_pos_done = False
            gapless_pos+=1
    if cur_pos_object != "":
        grouped[gapless_pos] = cur_pos_object
        cur_pos_object = ""
        cur_pos_done = False
        gapless_pos+=1
    
    return grouped

def group_by_position(cur_quasispecies):
    grouped = []
    cur_pos_object = ""
    cur_pos_done = False
    gapless_pos = 0
    for gapped_pos,(ref_n, quas_n) in enumerate(zip(reference, cur_quasispecies)):
        if ref_n == "-" and quas_n == "-":
            continue
        elif ref_n == "-":
            cur_pos_object += quas_n
            if gapped_pos == len(reference) - 1 or reference[gapped_pos + 1] != "-":
                cur_pos_done = True
        elif quas_n == "-":
            cur_pos_object = ""
            cur_pos_done = True
        else:
            if cur_pos_object != "":
                grouped.append(cur_pos_object)
                cur_pos_object = ""
                cur_pos_done = False
                gapless_pos+=1
            
            cur_pos_object = quas_n
            if gapped_pos == len(reference) - 1 or (reference[gapped_pos + 1] != "-"):
                cur_pos_done = True
        
        if cur_pos_done:
            grouped.append(cur_pos_object)
            cur_pos_object = ""
            cur_pos_done = False
            gapless_pos+=1
    if cur_pos_object != "":
        grouped.append(cur_pos_object)
        cur_pos_object = ""
        cur_pos_done = False
        gapless_pos+=1
    print((grouped))
    return grouped



def group_by_position(cur_quasispecies):
    grouped = []
    cur_pos_object = ""
    cur_pos_done = False
    gapless_pos = 0
    while gapless_pos < len(reference):
        cur_q_nuc = cur_quasispecies[gapless_pos]
        cur_r_nuc = reference[gapless_pos]
        
        if cur_q_nuc == "-" and cur_r_nuc == "-": #Doesn't matter
            if gapless_pos == len(reference) - 1 and cur_pos_object != "":
                grouped[-1]+=cur_pos_object
                cur_pos_object = ""
        elif cur_q_nuc == "-": #Deletion
            cur_pos_object += ""
            cur_pos_done = True
        elif cur_r_nuc == "-": # Insertion
            cur_pos_object += cur_q_nuc
            if gapless_pos == len(reference) - 1: # Because we are adding insertions from back to forward, if an insertion occurs in the last nucleotide, need to add it to the last element in the list
                grouped[-1] += cur_pos_object
                cur_pos_object = ""
                break
        else: #Normal
            cur_pos_object += cur_q_nuc
            cur_pos_done = True
        
        if cur_pos_done:
            cur_pos_done = False
            grouped.append(cur_pos_object)
            cur_pos_object = ""
        gapless_pos+=1
    if cur_pos_object != "":
        cur_pos_done = False
        grouped.append(cur_pos_object)
        cur_pos_object = ""
    if len(grouped) != reference_gapless_length:
        print(cur_quasispecies)
        print(f"LENGTH: {len(grouped)}")
        print(grouped)
        0/0
    return grouped














def get_entropies(group_by_position_array):
    out = np.empty(reference_gapless_length, dtype=np.float64)
    for i in trange(reference_gapless_length):
        frequencies = np.unique(group_by_position_array[:,i], return_counts=True)[1]
        normalized_freq = frequencies / frequencies.sum()
        entropy = np.sum(-normalized_freq * np.log2(normalized_freq))
        out[i] = entropy
    return out
            
def get_entropies(group_by_position_array):
    out = np.empty(reference_gapless_length, dtype=np.float64)
    for i in trange(reference_gapless_length):
        frequencies = np.unique(group_by_position_array[:,i], return_counts=True)[1]
        normalized_freq = frequencies / frequencies.sum()
        entropy = np.sum(-normalized_freq * np.log2(normalized_freq))
        out[i] = entropy
    return out

def main():
    cur_output = list(tqdm(map(group_by_position, quasispecies), total=len(quasispecies), desc="Creating Group By Position Array..."))    
    with open(snakemake.output[0], "w", newline="") as out_file:
        writer = csv.writer(out_file)
        writer.writerows(cur_output)


if __name__ == "__main__":
    main()
else:
    raise Exception("Not in '__main__' mode")

