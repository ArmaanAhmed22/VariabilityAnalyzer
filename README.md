# VariabilityAnalyzer

![frontcover](https://raw.githubusercontent.com/ArmaanAhmed22/VariabilityAnalyzer/master/Assets/front_cover.png)
A tool to help visualize the variability of a nucleic acid sequence (genome, gene, etc).

## Analysis

The variability of a region is assessed by calculating the Shannon entropy of each position in the region. Usually, insertions and deletions are ignored and only non-ambiguous nucleotides are considered. However, this tool can include both insertions and deletions in the entropy calculation. This is simply done by determining the "object" at each position (insertion, deletion, or single nucleotide), calculating the frequency of each object (`f`), and running the frequencies through the Shannon entropy formula: `-sum(f*log(f))`.

Ambiguous nucleotides are proportionally converted to all their non-ambiguous nucleotides. For example, if at a position, there is a `N` nucleotide observation, then 1 `N` is counted as 0.25 `A`, 0.25 `C`, 0.25 `G`, and 0.25 `T`. For insertions with ambiguous nucleotides, all non-ambiguous permutations are proportionally considered. (If the number of permutations grow beyond 10000, then the entropy calculation will simply use the ambiguous insertion).

## Installation

Simply download this repository!
`git clone https://github.com/ArmaanAhmed22/VariabilityAnalyzer.git`

## Usage

An input, aligned reference sequence and quasispecies files should be placed within `Input/{file}`, where `{file}` is the name of the region-to-be-analyzed (`protease`, `SARS-CoV-2`, `HIV-1`, etc, etc).

Next, run `snakemake Output/{file}/variability_map.png -c1` to generate the output files.

## Example

Directory:

```text
Input/
    Example/
        reference.fasta
        quasispecies.fasta
Pipeline/**
Snakefile
```

Run `snakemake Output/Example/variability_map.png -c1` to generate `variability_map.png` and `variability_map.csv`

### Example Output: *Protease*

![variability_map.png](https://raw.githubusercontent.com/ArmaanAhmed22/VariabilityAnalyzer/master/Output/Protease/variability_map.png)
