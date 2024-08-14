#!/usr/bin/env python3

"""

########################## INFO ################################

This helper script converts the output of snp-sites to an
olivar (https://github.com/treangenlab/Olivar) compatible output.

The aim is a head-to-head comparison of varVAMP and olivar. However,
olivar requires a csv with all variants and can not handle the initial
input alignment used by varVAMP. Therefore, we need to calculate all variants
of the respective alignments.
To make a head-to-head comparison fair, we use the gap cleaned alignment
of the varVAMP output. Therefore, the consensus sequence comprising the majority
nucleotides that was calculated by varVAMP has the identical ref nucleotide
as the most frequent variant in the called variants. This allows to use the
 majority consensus sequence as the input for olivar.

Variants were called with snp-sites (https://github.com/sanger-pathogens/snp-sites):
snp-sites -v -o snp-sites/[alignment_name.vcf] alignments/[respective alignment]

Afterwards run the script with:
python3 snp-sites_to_olivar.py

This will create olivar compatible tab sep csv files in a new folder:
snp-sites/olivar_input

########################## COPYRIGHT ################################

Dr. Jonas Fuchs (jonas.fuchs@uniklinik-freiburg.de)
Institute for Virology, Freiburg, Germany

"""


import os

import numpy as np

from in_silico_eval import get_files, get_file_names, read_fasta


def main():
    """
    reads in vcf files and majority consensus sequences and converts the
    information to a olivar compatible output.
    """
    vcf_files = sorted(get_files("snp-sites"))
    sequence_files = sorted(get_files("consensus_files"))
    names = get_file_names(vcf_files)
    # check if folder has been deleted
    if not os.path.exists("snp-sites/olivar_input"):
        os.makedirs("snp-sites/olivar_input")

    for vcf_file, name, sequence_file in zip(vcf_files, names, sequence_files):
        seq_id, seq = read_fasta(sequence_file)
        with open(f"snp-sites/olivar_input/{name}_olivar.csv", "w") as olivar_csv:
            print("START,STOP,FREQ,REF,ALT", file=olivar_csv)
            with open(vcf_file, "r") as vcf:
                for line in vcf:
                    if line.startswith("#"):
                        continue
                    line = line[:-1]
                    values = line.split("\t")
                    start_pos = int(values[1])
                    nucs = [values[3]]+values[4].split(",")
                    numeric_nucs, counts = np.unique(values[9:], return_counts=True)
                    total_counts = sum(counts)
                    # sort counts and nucleotides
                    counts_sorted = sorted(counts, reverse=True)
                    # the vcf file from snp-sites reports the respective nucleotide
                    # for each sequence id. to calculate the FREQ the most freq nuc
                    # is the same as the nuc of the consensus sequence. all other
                    # nucleotides are therefore variants of the consensus sequence
                    # the next section is to demultiplex the different variants at the
                    # same position
                    nucs_sorted = [x for _, x in sorted(zip(counts, nucs), reverse=True)]
                    nucs_sorted = list(map(lambda x: x.replace('*', '-'), nucs_sorted))
                    for alt_nuc, count in zip(nucs_sorted[1:], counts_sorted[1:]):
                        nuc_freq = count/total_counts*100
                        print(start_pos, start_pos,nuc_freq,seq[start_pos-1],alt_nuc,sep=",", file=olivar_csv)

if __name__ == "__main__":
    main()
