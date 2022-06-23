import inspect
import numpy as np
import scipy
import pandas as pd
from tqdm import trange
import seaborn as sns
import matplotlib.pyplot as plt
import math

from tqdm import tqdm
import matplotlib.collections as mc

NUCS = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
from functools import reduce
def get_all_possible_sequences_length(sequence):
    return reduce(lambda x, y: x * y, [len(NUCS[nuc]) for nuc in sequence])
    
def get_all_possible_sequences(sequence):
    compressed_all_sequences_non_ambiguous = [NUCS[nuc] for nuc in sequence]
    indices = [0] * len(compressed_all_sequences_non_ambiguous)
    cur_p = 0
    while True:
        cur_p+=1
        #print(f"{cur_p}/{total}")
        # Get the current sequence
        cur_seq = ""
        for position, cur_index in enumerate(indices):
            cur_seq += compressed_all_sequences_non_ambiguous[position][cur_index]
        yield cur_seq
        # Increment the indices
        for i in range(len(indices))[::-1]:
            if indices[i] + 1 >= len(compressed_all_sequences_non_ambiguous[i]):
                indices[i] = 0
            else:
                indices[i] += 1
                break
        else:
            break
    return
    

#Fix for IUPAC Ambiguous nucleotides
# 
with open(snakemake.input[0]) as group_by_position_h:
    data = pd.read_csv(group_by_position_h, header=None, dtype=str)
    for i in trange(data.shape[0], desc="Finding ends of sequences"):
        fvi = data.iloc[i].first_valid_index()
        lvi = data.iloc[i].last_valid_index()
        data.iloc[i, fvi:lvi] = data.iloc[i, fvi:lvi].fillna("")
    data2 = data.applymap(lambda x: np.NaN if (type(x) != str and math.isnan(x)) or len(x) != 1 else x)


def get_entropies(pd_data, IUPAC_changed_frequencies):
    # Get the entropy for each position
    out = np.empty(pd_data.shape[1], dtype=np.float64)
    for i, cur_frequencies in tqdm(enumerate(IUPAC_changed_frequencies), desc="Generating entropies", total=len(IUPAC_changed_frequencies)):
        frequencies = np.array(list(cur_frequencies.values()))
        normalized_freq = frequencies / frequencies.sum()
        entropy = np.sum(-normalized_freq * np.log2(normalized_freq))
        out[i] = entropy
    return out

def get_insertion_deletion_frequencies(pd_data, IUPAC_changed_frequencies):
    # Get the insertion and deletion frequencies
    out = np.empty(pd_data.shape[1], dtype=np.float64)
    for i, cur_frequencies in tqdm(enumerate(IUPAC_changed_frequencies), desc="Generating insertion/deletion frequencies", total=len(IUPAC_changed_frequencies)):
        
        insertions_deletions_count = 0
        for index,occurences in cur_frequencies.items():
            if len(index) > 1:
                insertions_deletions_count += occurences
        normalized_freq = insertions_deletions_count / sum(cur_frequencies.values())
        out[i] = normalized_freq
    return out

def get_total_occurences(pd_data, IUPAC_changed_frequencies):
    # Get the total nucleotide occurences for each position
    out = np.empty(pd_data.shape[1], dtype=np.float64)
    for i, cur_frequencies in tqdm(enumerate(IUPAC_changed_frequencies), desc="Generating total occurences", total=len(IUPAC_changed_frequencies)):
        out[i] = sum(cur_frequencies.values())
    return out

def get_IUPAC_changed_frequencies(pd_data) -> list[pd.Series]:
    out = []
    for i in trange(pd_data.shape[1], desc="Generating IUPAC changed frequencies"):
        
        total_frequencies = pd_data.iloc[:,i].value_counts()
        dict_total_frequencies = {}
        to_drop = []
        # Remove the IUPAC nucleotides
        for cur_nuc in total_frequencies.index:
            if len(cur_nuc) >= 1:
                # Substitute the IUPAC nucleotide for non-ambiguous nucleotides
                all_possibilities = get_all_possible_sequences(cur_nuc)
                all_possibilities_length = get_all_possible_sequences_length(cur_nuc)
                if all_possibilities_length > 10_000:
                    all_possibilities = [cur_nuc]
                    all_possibilities_length = 1
                for cur_possibility in all_possibilities:
                    if cur_possibility not in dict_total_frequencies:
                        dict_total_frequencies[cur_possibility] = 0
                    dict_total_frequencies[cur_possibility] += total_frequencies[cur_nuc] / all_possibilities_length
        if "" in total_frequencies.index:
            dict_total_frequencies[""] = total_frequencies[""]
        out.append(dict_total_frequencies)
    return out
                

def main():
    frequencies_insertions_deletions = get_IUPAC_changed_frequencies(data)
    frequencies_no_insertions_deletions = get_IUPAC_changed_frequencies(data2)
    entropies = get_entropies(data, frequencies_insertions_deletions)
    entropies_no_insertions_deletions = get_entropies(data, frequencies_no_insertions_deletions)
    
    insertion_deletion_frequencies = get_insertion_deletion_frequencies(data, frequencies_insertions_deletions)
    occurences = get_total_occurences(data, frequencies_insertions_deletions)
    
    
    fig, axes = plt.subplots(nrows=5, figsize=(16,16))
    dataframe = pd.DataFrame({"Entropy": entropies, "Entropy (No Insertions/Deletions)": entropies_no_insertions_deletions, "Change in Entropy": entropies - entropies_no_insertions_deletions, "Occurences": occurences, "Insertion Deletion Frequencies": insertion_deletion_frequencies, "Position": range(1, data.shape[1] + 1)})
    #dataframe = pd.DataFrame({"Entropy": entropies[2253:2550], "Position": range(2253, 2550)})
    sns.lineplot(x="Position", y="Entropy", data=dataframe, ax=axes[0])
    sns.lineplot(x="Position", y="Entropy (No Insertions/Deletions)", data=dataframe, ax=axes[1])
    sns.lineplot(x="Position", y="Change in Entropy", data=dataframe, ax=axes[2])
    sns.lineplot(x="Position", y="Insertion Deletion Frequencies", data=dataframe, ax=axes[3])
    sns.lineplot(x="Position", y="Occurences", data=dataframe, ax=axes[4])
    #print(inspect.getsource(axes[0]._children))
    #base_color = axes[0]._children[-1]._color
    
    """lines = []
    colors = []
    for i, row in enumerate(dataframe.itertuples()):
        if i != 0:
            lines[-1].append((row.Position, row.Entropy))
        lines.append([(row.Position, row.Entropy)])
        colors.append((0,0,0,row.Occurences/dataframe["Occurences"].max()))
    
    axes[0].add_collection(mc.LineCollection(lines, colors=colors))
    axes[0].autoscale()"""
    sns.despine()
    
    
    plt.savefig(snakemake.output[0], dpi=600, bbox_inches="tight")
    dataframe.to_csv(snakemake.output[1], index=False)

main()