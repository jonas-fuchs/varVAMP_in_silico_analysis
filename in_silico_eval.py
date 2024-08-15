#!/usr/bin/env python3

"""
########################## INFO ################################

this script reproduces all plots for the in silico analysis of:

"varVAMP: automated pan-specific primer design for tiled
full genome sequencing and qPCR of highly diverse viral pathogens."

All relevant data is given within this repo. The data was produced with:

- adapted bed primer files:
    bed annotations for the primer locations of the varVAMP output were adapted to the appropriate references
    used for mapping the raw reads.
- alignments:
    Output of varVAMP. The initial non-gap masked MAFFT alignment used as the varVAMP input can be found at
    https://github.com/jonas-fuchs/ViralPrimerSchemes
- consensus files:
    varVAMP output
- coverages per amplicon:
    bamDASH output from mapped *.bam files and bed files for the amplicons used as track
    (--dump to dump the data to tabular)
- new sequences:
    output of the mapping and consensus Galaxy pipeline. The regions of the most left and right primers were
    masked with 'N'. Additionally low covered regions (<20x) and mutations between a frequency of 0.3 and 0.7
    were also masked with 'N'.
- per base coverages:
    output of bamQC (https://github.com/s-andrews/BamQC) generated from mapped *.bam files
- primer bed files:
    initial primer bed files from the varVAMP output. Locations are relative to the consensus sequences
- primer tsv files:
    varVAMP output
- reference seq:
    reference sequences used for mapping
- sequence identity:
    pairwise identities calculated with https://github.com/BioinformaticsToolsmith/Identity using either the initial
    sequences used to generate the MAFFT alignment of the varVAMP input or the new sequences
- variant tables:
    tabular files extracted from vcf files (variant callings on *.bam files) with SnpSift Extract Fields.
    for medaka variant calls (ONT SARS-CoV-2 data), the AF field was artificially added and set to 1 for all mutations,
    to ensure compatibility with the Illumina data.

########################## INSTALLATION AND RUNNING THE ANALYSIS ################################

requires python 3.11

Run:
    pip install .
    python3 extract_as_fasta.py

will produce a new "output" dir with some tabular files and all plots shown in the publication

########################## COPYRIGHT ################################

Dr. Jonas Fuchs (jonas.fuchs@uniklinik-freiburg.de)
Institute for Virology, Freiburg, Germany


"""

# BUILT-INS
import shutil
import math
import warnings
import itertools
import statistics
from os import listdir, makedirs
from os.path import isfile, join, basename, splitext, exists, isdir
# LIBS
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import primer3 as p3
from Bio.Seq import Seq
from Bio import AlignIO

AMBIG_NUCS = {
    "r": ["a", "g"],
    "y": ["c", "t"],
    "s": ["g", "c"],
    "w": ["a", "t"],
    "k": ["g", "t"],
    "m": ["a", "c"],
    "b": ["c", "g", "t"],
    "d": ["a", "g", "t"],
    "h": ["a", "c", "t"],
    "v": ["a", "c", "g"],
    "n": ["a", "c", "g", "t"]
}


def get_files(mypath):
    """
    list alignment files
    """
    return [f"{mypath}/{f}" for f in listdir(mypath) if isfile(join(mypath, f))]


def get_file_names(files):
    """
    list all file names
    """
    name_list = []
    for name in [splitext(basename(f))[0].replace("_", " ") for f in files]:
        # allows to manually sort files by a prepending 1_ etc
        if name[0].isdigit():
            name = name[2:]
        name_list.append(name)

    return name_list


def list_folder_names(folder):
    """
    list all subfolders
    """
    return [name for name in listdir(folder) if isdir(join(folder, name))]


def read_fasta(c_file):
    """
    reads lines and return the line after the header
    """
    with open(c_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].rstrip("\n")
                seq_id = seq_id.split(" ")[0]
            else:
                sequence = Seq(line.rstrip("\n"))

    return seq_id, sequence


def read_alignment(alignment_path):
    """
    read alignment with AlignIO and
    convert to list of lists
    """
    alignment_list = []

    for sequence in AlignIO.read(alignment_path, "fasta"):
        alignment_list.append([sequence.id, str(sequence.seq)])

    return alignment_list


def calc_hairpin(seq):
    """
    calculates hairpins
    """
    return p3.calc_hairpin(
            seq.upper(),
            mv_conc=100,
            dv_conc=2,
            dntp_conc=0.8,
            dna_conc=15
        )


def calc_dimer(seq1, seq2, structure=False):
    """
    Calculate the heterodimerization thermodynamics of two DNA sequences.
    Return primer3 thermo object.
    """
    return p3.calc_heterodimer(
        seq1.upper(),
        seq2.upper(),
        mv_conc=100,
        dv_conc=2,
        dna_conc=15,
        dntp_conc=0.8,
        output_structure=structure
    )


def get_permutations(seq):
    """
    get all permutations of an ambiguous sequence. needed to
    correctly report the gc and the temperature.
    """
    groups = itertools.groupby(seq, lambda char: char not in AMBIG_NUCS)
    splits = []
    for b, group in groups:
        if b:
            splits.extend([[g] for g in group])
        else:
            for nuc in group:
                splits.append(AMBIG_NUCS[nuc])
    return[''.join(p) for p in itertools.product(*splits)]


def alignment_stats(alignment_folder):
    """
    basic alignment stats
    """
    alignment_files = get_files(alignment_folder)
    names = get_file_names(alignment_files)
    for file, name in zip(alignment_files, names):
        alignment = read_alignment(file)
        gc_content = round(np.mean([(x[1].count("g")+x[1].count("c"))/(x[1].count("g")+x[1].count("c")+x[1].count("a")+x[1].count("t"))*100 for x in alignment]), 2)
        print(f"\t{name}: n of sequences: {len(alignment)}, GC content: {gc_content} %")


def pairwise_evaluation(filepath):
    """
    analysis of the pairwise sequence identities
    produced by https://github.com/BioinformaticsToolsmith/Identity
    """
    files = get_files(filepath)
    names = get_file_names(files)
    results = []

    for file, name in zip(files, names):
        if file.endswith(".txt"):
            identities = []
            with open(file, "r") as identity_file:
                for line in identity_file:
                    identities.append(float(line.strip().split("\t")[2])*100)
            if len(identities) > 1:
                results.append([name, statistics.mean(identities), statistics.stdev(identities), identities])
            else:
                results.append([name, identities[0], 0, identities])

    return pd.DataFrame(results, columns=["virus", "mean", "std", "all_values"])


def entropy(chars, states):
    """
    input is a list of characters or numbers.
    calculate relative shannon's entropy. relative values are
    achieved by using the number of possible states as the base
    """
    ent = 0
    n_chars = len(chars)
    # only one char is in the list
    if n_chars <= 1:
        return ent
    # calculate the number of unique chars and their counts
    value, counts = np.unique(chars, return_counts=True)
    probs = counts/n_chars
    if np.count_nonzero(probs) <= 1:
        return ent

    for prob in probs:
        ent -= prob*math.log(prob, states)

    return ent


def alignment_entropy(alignment, states):
    """
    calculate the entropy for every position in an alignment.
    return pandas df.
    """
    pos_index = list()
    entropys = list()
    # iterate over alignment positions and the sequences
    for nuc_pos in range(0, len(alignment[0][1])):
        pos = []
        for seq_number in range(0, len(alignment)):
            pos.append(alignment[seq_number][1][nuc_pos])
        entropys.append(entropy(pos, states))
        pos_index.append(nuc_pos)
    # create df
    entropy_df = pd.DataFrame()
    entropy_df["position"] = pos_index
    entropy_df["normalized position"] = [i/max(pos_index)*100 for i in pos_index]
    entropy_df["entropy"] = entropys
    entropy_df["normalized Shannon's entropy (average)"] = entropy_df["entropy"].rolling(int(max(pos_index)/100), center=True).mean()

    return entropy_df


def calc_permutations(primer):
    """
    get all permutations of a primer with ambiguous nucleotides
    """

    permutations = 1

    for nuc in primer:
        if nuc in AMBIG_NUCS:
            n = len(AMBIG_NUCS[nuc])
            if permutations != 1:
                permutations = permutations*n
            else:
                permutations = n

    return permutations


def rev_complement(seq):
    """
    reverse complement a sequence
    """
    return str(Seq(seq).reverse_complement())


def calculate_mismatches(primer_starts, primer_stops, primer_seqs, primer_names, alignment):
    """
    calculates all the following statistics for one primer scheme:
      - count mismatches per primer pos and normalize to number of sequences in alignment.
      - count mismatches per primer and sequence and return mean.
    """

    mean_mismatches, mismatches, per_pos_mismatches_per_scheme = [], [], []

    for start, stop, primer, primer_name in zip(primer_starts, primer_stops, primer_seqs, primer_names):
        mismatches_per_position = len(primer) * [0]
        mismatches_per_sequence = []
        primer = primer.lower()
        if "RW" in primer_name or "RIGHT" in primer_name:
            primer = rev_complement(primer)
        for sequence in alignment:
            mismatch = 0
            seq_slice = sequence[1][start:stop]
            for idx, nuc in enumerate(seq_slice):
                primer_nuc = primer[idx]
                if nuc == primer_nuc or nuc == "-":
                    continue
                if nuc in AMBIG_NUCS:
                    # check if the kmer has an amb pos
                    if primer_nuc in AMBIG_NUCS:
                        slice_nuc_set = set(AMBIG_NUCS[nuc])
                        pri_set = set(AMBIG_NUCS[primer_nuc])
                        # check if these sets have no overlap
                        # -> mismatch
                        if len(slice_nuc_set.intersection(pri_set)) == 0:
                            mismatches_per_position[idx] += 1
                            mismatch += 1
                    # if no amb pos is in kmer then check if kmer nuc
                    # is part of the amb slice nuc
                    elif primer_nuc not in AMBIG_NUCS[nuc]:
                        mismatches_per_position[idx] += 1
                        mismatch += 1
                    # check if kmer has an amb pos but the current
                    # slice_nuc is not part of this amb nucleotide
                elif primer_nuc in AMBIG_NUCS:
                    if nuc not in AMBIG_NUCS[primer_nuc]:
                        mismatches_per_position[idx] += 1
                        mismatch += 1
                    # mismatch
                else:
                    mismatches_per_position[idx] += 1
                    mismatch += 1
            mismatches_per_sequence.append(mismatch)
        # append to lists
        if "FW" in primer_name:
            mismatches_per_position.reverse()
        # norm to % alignment
        per_pos_mismatches_per_scheme.append([x/len(alignment)*100 for x in mismatches_per_position])
        # reverse FW primers so positions can be given as nt distance from 5' end
        mean_mismatches.append(np.mean(mismatches_per_sequence))
        mismatches.append(mismatches_per_sequence)

    return mismatches, mean_mismatches, per_pos_mismatches_per_scheme


def calc_dimer_and_hairpin_temps(bed_df, tsv_df, consensus_seq):
    """
    calculate melting temps and homo-dimers for primers and all permutations
    """

    melting_temp_opt, melting_temp_per, dimer_temp_opt, dimer_temp_per = [], [], [], []

    for start, stop, p_name, p_seq in zip(bed_df[1], bed_df[2], bed_df[3], tsv_df["seq"]):
        primer = consensus_seq[start:stop]
        permut_p = get_permutations(p_seq)
        if "RW" in p_name:
            primer = primer.reverse_complement()
            permut_p = [Seq(x).reverse_complement() for x in permut_p]
        primer = str(primer)
        permut_p = [str(x) for x in permut_p]
        permut_tmp = np.mean([calc_hairpin(x).tm for x in permut_p])
        permut_dimer = np.mean([calc_dimer(x, x).tm for x in permut_p])

        melting_temp_opt.append(calc_hairpin(primer).tm), melting_temp_per.append(permut_tmp)
        dimer_temp_opt.append(calc_dimer(primer, primer).tm), dimer_temp_per.append(permut_dimer)

    return melting_temp_opt + melting_temp_per, dimer_temp_opt + dimer_temp_per


def calculate_and_plot_entropy(input_folder, output_folder):
    """
    create the entropy plots
    """
    files = get_files(input_folder)
    names = get_file_names(files)
    # set colormap
    # ini fig
    fig = plt.figure(figsize=(16, 4.5))
    # cal entropy and plot
    for file, name in zip(files, names):
        alignment = read_alignment(file)
        entropies = alignment_entropy(alignment, 4)
        entropies.to_csv(f"{output_folder}/{name}_entropy.tabular", sep="\t", header=True, index=False)
        plt.plot(
            entropies["normalized position"],
            entropies["normalized Shannon\'s entropy (average)"],
            linestyle="-",
            label=name,
            linewidth=1.5,
        )
    # set plotting stuff
    plt.legend(frameon=False)
    plt.xlabel("% multiple sequence alignment")
    plt.ylabel("moving average (normalized Shannon's entropy)")
    plt.ylim([0, 1])
    plt.xlim([0, 100])
    plt.xticks(np.arange(0, 100 + 1, 10))
    sns.despine()
    fig.savefig(f"{output_folder}/entropy.pdf")


def calculate_and_plot_degeneracy(input_folder, output_folder, colorscheme):
    """
    calc the degeneracy of the individual primers of each scheme and plot
    """

    tsv_files = get_files(input_folder)
    names = get_file_names(tsv_files)

    permutations = []

    for tsv_file in tsv_files:
        permutations_per_scheme = []
        tsv_df = pd.read_csv(tsv_file, sep="\t", header=0)
        for primer_seq in tsv_df["seq"]:
            permutations_per_scheme.append(calc_permutations(primer_seq))
        permutations.append(permutations_per_scheme)

    plt.figure(figsize=(3, 4.5))
    ax = sns.stripplot(permutations, palette=colorscheme)
    sns.boxplot(showmeans=True,
                meanline=True,
                meanprops={'color': 'k', 'ls': '-', 'lw': 2},
                medianprops={'visible': False},
                whiskerprops={'visible': False},
                zorder=10,
                data=permutations,
                showfliers=False,
                showbox=False,
                showcaps=False,
                ax=ax)

    ax.set_xticklabels(labels=range(0, len(names)), rotation=45, ha="right")
    ax.set_xticklabels(names)
    sns.despine()
    plt.ylabel("n permutations per primer")
    plt.savefig(f"{output_folder}/degeneracy.pdf", bbox_inches='tight')


def calculate_and_plot_mismatches(alignment_folder, bed_folder, tsv_folder, output_folder):
    """
    calculate mismatches at each primer position and number of mismatches of
    each primer to each sequence of an alignment
    """

    flattened_mismatches, per_pos_mismatch_list_all = [], []

    # load data
    alignment_files = get_files(alignment_folder)
    bed_files = get_files(bed_folder)
    tsv_files = get_files(tsv_folder)
    names = get_file_names(bed_files)

    # calculate all mismatches
    for bed_file, alignment_file, tsv_file, name in zip(bed_files, alignment_files, tsv_files, names):
        bed_df = pd.read_csv(bed_file, sep="\t", header=None)
        tsv_df = pd.read_csv(tsv_file, sep="\t", header=0)
        alignment = read_alignment(alignment_file)
        mismatches, mean_mismatches, per_pos_mismatches_per_scheme = calculate_mismatches(
            bed_df[1], bed_df[2], tsv_df["seq"],tsv_df["primer_name"], alignment
        )
        pd.DataFrame(
            data=[mismatches, mean_mismatches, per_pos_mismatches_per_scheme],
            columns=tsv_df["primer_name"],
            index=["mismatches", "mean mismatches", "per position mismatches (3' to 5')"]
        ).to_csv(
            f"output/{name}_mismatches.tabular",
            sep="\t",
            header=True
        )
        # flatten all mismatches to one list
        flattened_mismatches.append([item for sublist in mismatches for item in sublist])
        per_pos_mismatch_list_all.append(per_pos_mismatches_per_scheme)

    # ini figure for mismatch bubble plot
    fig, ax = plt.subplots(figsize=(4, 4.5))

    for idx, mismatches in enumerate(flattened_mismatches):
        value, counts = np.unique(mismatches, return_counts=True)
        # normalize to total number
        counts_norm = [x / sum(counts) * 100 for x in counts]
        # plot to ax
        ax.scatter(x=[idx] * len(counts_norm), y=value, s=counts_norm)
    # pseudodata for legend
    legend_dot_size = ax.scatter([0, 0, 0], [0, 0, 0], s=[1, 10, 100], edgecolor=None, alpha=0.1)
    legend = plt.legend(
        *legend_dot_size.legend_elements("sizes"),
        title="Percent",
        frameon=False,
        loc='upper center',
        ncols=3
    )
    ax.add_artist(legend)
    ax.set_xticklabels(rotation=45, ha="right", labels=range(-1, len(counts_norm)))
    ax.set_xticklabels([""] + names)
    ax.set_ylabel("nt mismatches between primers and sequences")
    sns.despine()
    fig.savefig(f"{output_folder}/mismatches.pdf", bbox_inches='tight')


    # ini figure for distance mismatch plot
    fig = plt.figure(figsize=(7, 4.5))

    for i, (per_pos_mismatches_per_scheme, name) in enumerate(zip(per_pos_mismatch_list_all, names)):
        # maximum length of list
        max_len = max([len(x) for x in per_pos_mismatches_per_scheme])
        distance_5_prime = [range(1, max_len+1)]
        mismatches = [0]*max_len
        counts = [0]*max_len

        for per_pos_mismatches in per_pos_mismatches_per_scheme:
            for j, pos in enumerate(per_pos_mismatches):
                # count mismatches
                mismatches[j] += pos
                counts[j] += 1
        normalized_mismatches = [x/y for x,y in zip(mismatches, counts)]
        # add None so they have same length
        normalized_mismatches = normalized_mismatches + [None]*(25-len(normalized_mismatches))
        sns.lineplot(x=list(range(0, 25)), y=normalized_mismatches, label=name)

    plt.legend(frameon=False)
    plt.xlabel("distance from 3' primer end")
    plt.ylabel("% mismatch with alignment sequences")
    sns.despine()
    plt.xticks(np.arange(0, 25, 2.0))
    fig.savefig(f"{output_folder}/mismatches_distance_from_3_prime.pdf")


def set_size(w, h, ax=None):
    """
    set size of the axis rather than the plot - circumvents weird plotting issues that automatically scales white space
    w, h: width, height in inches
    """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


def plot_primer_stats(output_folder, tsv_folder, bed_folder, consensus_folder):
    """
    plot remaining parameters for the primers
    """

    tsv_files = get_files(tsv_folder)
    bed_files = get_files(bed_folder)
    consensus_files = get_files(consensus_folder)
    names = get_file_names(tsv_files)

    all_dfs = []

    for tsv_file, bed_file, consensus_seq, name in zip(tsv_files, bed_files, consensus_files, names):
        df_new = pd.DataFrame()
        tsv_df = pd.read_csv(tsv_file, sep="\t", header=0)
        bed_df = pd.read_csv(bed_file, sep="\t", header=None)
        consensus_seq = read_fasta(consensus_seq)[1]
        hairpin_temps = calc_dimer_and_hairpin_temps(bed_df, tsv_df, consensus_seq)
        df_new["% GC content"] = list(tsv_df["gc_best"])+list(tsv_df["mean_gc"])
        df_new["°C melting temp"] = list(tsv_df["temp_best"]) + list(tsv_df["mean_temp"])
        df_new["°C hairpin temp"] = hairpin_temps[0]
        df_new["°C homo-dimer temp"] = hairpin_temps[1]
        df_new["type"] = ["opt (most common nt)"]*len(tsv_df["temp_best"]) + ["mean (all permutations)"]*len(tsv_df["temp_best"])
        df_new["scheme"] = name
        all_dfs.append(df_new)
    # cat to one big df
    df_for_plot = pd.concat(all_dfs, ignore_index=True)

    # plot temp
    fig, ax = plt.subplots(figsize=(3, 4.5))
    sns.stripplot(
        data=df_for_plot,
        x="scheme",
        y="°C melting temp",
        hue="type",
        dodge=True,
        palette=["sienna", "slategrey"],
        ax=ax)
    sns.despine()
    ax.set(xlabel=None)
    plt.ylim([45, 75])
    plt.axhline(y=60, linestyle="--", color="black", label="target")
    plt.axhline(y=56, linestyle=":", color="grey", label="min")
    plt.axhline(y=63, linestyle=":", color="grey", label="max")
    plt.legend(frameon=False,  loc="upper right")
    plt.xticks(rotation=45, ha="right")
    fig.savefig(f"{output_folder}/temp.pdf", bbox_inches='tight')

    # plot GC
    fig, ax = plt.subplots(figsize=(3, 4.5))
    sns.stripplot(
        data=df_for_plot,
        x="scheme",
        y="% GC content",
        hue="type",
        dodge=True,
        palette=["sienna", "slategrey"],
        ax=ax)
    sns.despine()
    ax.set(xlabel=None)
    plt.ylim([0,100])
    plt.axhline(y=50, linestyle="--", color="black", label="target")
    plt.axhline(y=35, linestyle=":", color="grey", label="min")
    plt.axhline(y=65, linestyle=":", color="grey", label="max")
    plt.legend(frameon=False,  loc="upper right")
    plt.xticks(rotation=45, ha="right")
    fig.savefig(f"{output_folder}/gc_content.pdf", bbox_inches='tight')

    # plot melting temp
    fig, ax = plt.subplots(figsize=(3, 4.5))
    sns.stripplot(
        data=df_for_plot,
        x="scheme",
        y="°C hairpin temp",
        hue="type",
        dodge=True,
        palette=["sienna", "slategrey"],
        ax=ax)
    sns.despine()
    ax.set(xlabel=None)
    plt.ylim([0, 100])
    plt.axhline(y=47, linestyle="--", color="black", label="max")
    plt.legend(frameon=False, loc="upper right")
    plt.xticks(rotation=45, ha="right")
    fig.savefig(f"{output_folder}/melting_temp.pdf", bbox_inches='tight')

    # plot homo-dimer temp
    fig, ax = plt.subplots(figsize=(3, 4.5))
    sns.stripplot(
        data=df_for_plot,
        x="scheme",
        y="°C homo-dimer temp",
        hue="type",
        dodge=True,
        palette=["sienna", "slategrey"],
        ax=ax)
    sns.despine()
    ax.set(xlabel=None)
    plt.ylim(top=100)
    plt.axhline(y=47, linestyle="--", color="black", label="max")
    plt.legend(frameon=False, loc="upper right")
    plt.xticks(rotation=45, ha="right")
    fig.savefig(f"{output_folder}/dimer_temp.pdf", bbox_inches='tight')


def plot_sequence_identity_comparison(identity_all_df, identity_folder, output_folder):
    """
    plot the pairwise sequences identities of the alignment and the newly produced sequences
    as a dumbbell plot
    """
    # get one df with sequence identities and alignment identities
    identity_comparison = pairwise_evaluation(identity_folder)

    alignment_identity, text, ttest_result = [], [], []

    for name in identity_comparison["virus"]:
        if name in list(identity_all_df["virus"]):
            alignment_identity_temp = float(identity_all_df[identity_all_df["virus"] == name]["mean"])
            identity_diff = float(identity_comparison[identity_comparison["virus"] == name]["mean"] - alignment_identity_temp)
            alignment_identity.append(alignment_identity_temp)
            text.append(f"{round(identity_diff)}%")

            # t-test for the alignment seq identities and new seq identities
            all_values_new = list(identity_comparison[identity_comparison["virus"] == name]["all_values"])[0]
            all_values_aln = list(identity_all_df[identity_all_df["virus"] == name]["all_values"])[0]
            # check if n at least 3
            if all([len(all_values_new) >= 3, len(all_values_aln) >= 3]):
                p_value = round(stats.ttest_ind(all_values_new, all_values_aln, equal_var=False).pvalue, 3)
                # associate stars
                current_significance = "n.s."
                for sign_n, stars in zip((0.05, 0.01, 0.001), ("*", "**", "***")):
                    if p_value <= sign_n:
                        current_significance = stars
                if current_significance == "***":
                    ttest_result.append(
                        f"$^{{{current_significance}}}$p < 0.001"
                    )
                else:
                    ttest_result.append(
                        f"$^{{{current_significance}}}$p = {p_value}"
                    )
            else:
                ttest_result.append("n.d.")
        else:
            alignment_identity.append("NA"), text.append("NA"), ttest_result.append("NA")

    # add columns
    identity_comparison["alignment_mean"], identity_comparison["change"], identity_comparison["ttest_pvalue"] = alignment_identity, text, ttest_result
    # dumbbell plot
    plt.figure(figsize=(9, 4.5))
    colors = np.where(identity_comparison["mean"] < identity_comparison["alignment_mean"], "#d9d9d9", "#ed98a2")
    # plot hlines with colors dependent if the mean identity is higher or lower
    plt.hlines(y=identity_comparison["virus"],
               xmin=identity_comparison["mean"],
               xmax=identity_comparison["alignment_mean"],
               color=colors,
               lw=10)
    # generate x,y tuple for annotation
    xy_values = [[(x - 0.5, y + 0.25), (x2 + 1, y)] for x, y, x2 in zip(identity_comparison[["mean", "alignment_mean"]].sum(axis=1)/2, np.arange(0, len(identity_comparison.index) + 1), identity_comparison[["mean", "alignment_mean"]].max(axis=1))]
    # annotate change
    for change, p_value, xy in zip(list(identity_comparison["change"]), list(identity_comparison["ttest_pvalue"]), xy_values):
        plt.annotate(change, xy[0], verticalalignment="center")
        plt.annotate(p_value, xy[1], verticalalignment="center")
    plt.scatter(identity_comparison["mean"],
                identity_comparison["virus"],
                color="#0096d7",
                s=200,
                label="new sequences",
                zorder=3)
    plt.scatter(identity_comparison["alignment_mean"],
                identity_comparison["virus"],
                color="grey",
                s=200,
                label="alignment",
                zorder=3)
    sns.despine()
    plt.legend(ncol=2,
               bbox_to_anchor=(1., 1.01),
               loc="lower right",
               frameon=False)
    plt.xlim(right=100)
    plt.xlabel("% mean pairwise sequence identity")
    plt.tight_layout()
    plt.savefig(f"{output_folder}/identity_comparision.pdf", bbox_inches='tight')


def plot_per_amplicon_coverages(coverages, output_folder):
    """
    plot for all viruses the per amplicon mean coverages
    """

    for virus in list_folder_names(coverages):
        files = get_files(f"{coverages}/{virus}")
        names = get_file_names(files)

        # get normalized coverages
        all_dfs = []
        for file, name in zip(files, names):
            per_amplicon_cov = pd.read_csv(file, sep="\t", header=0)
            max_coverage = max(per_amplicon_cov["mean coverage"]/(per_amplicon_cov["stop"]-per_amplicon_cov["start"]))
            per_amplicon_cov["normalized_coverages"] = list((per_amplicon_cov["mean coverage"]/(per_amplicon_cov["stop"]-per_amplicon_cov["start"]))/max_coverage*100)
            per_amplicon_cov["scheme_name"] = [name]*len(per_amplicon_cov["normalized_coverages"])
            all_dfs.append(per_amplicon_cov)
        final_df = pd.concat(all_dfs)
        # plot
        fig, (ax1, ax2) = plt.subplots(nrows=2, squeeze=True, sharex=True)
        palette = sns.color_palette("copper", len(set(per_amplicon_cov["normalized_coverages"])))  # define colours
        for index, amplicon_name in enumerate(list(set(final_df["name"]))):
            sns.lineplot(
                data=final_df[final_df["name"] == amplicon_name],
                x="scheme_name",
                y="recovery",
                marker="o",
                color=palette[index],
                legend=False,
                ax=ax1
            )
            sns.lineplot(
                data=final_df[final_df["name"] == amplicon_name],
                x="scheme_name",
                y="normalized_coverages",
                marker="o",
                color=palette[index],
                legend=False,
                ax=ax2
            )
        sns.despine()
        ax1.set_ylabel("% recovery (>= 20x)")
        ax2.set_ylabel("% normalized coverage")
        ax2.set_yscale("log")
        ax2.xaxis.set_label_text("")
        ax2.set_xlim(left=-0.5, right=len(set(final_df["scheme_name"]))-0.5)  # overwrite autospacing so it matches barplot
        ax1.set_ylim([0, 105])
        plt.tight_layout()
        plt.xticks(rotation=45, ha="right")
        set_size(len(files)*0.35, 4.5)
        plt.savefig(f"{output_folder}/{virus}_per_amplicon_coverage.pdf", bbox_inches='tight')


def analyse_and_plot_primer_binding(adapted_bed_folder, ref_folder, variant_folder, per_base_coverages_folder, tsv_folder, output_folder):
    """
    analyse how well the primers bind to each newly generated sequence and plot
    """

    primer_mismatch_dictionary = {}

    tsv_files = get_files(tsv_folder)

    for virus_name in list_folder_names("adapted_bed_primer_files"):
        # read in files
        bed_files = get_files(f"{adapted_bed_folder}/{virus_name}")
        ref_files = get_files(f"{ref_folder}/{virus_name}")
        variant_files = get_files(f"{variant_folder}/{virus_name}")
        coverage_files = get_files(f"{per_base_coverages_folder}/{virus_name}")
        sample_names = get_file_names(variant_files)
        # ini stuff
        primer_mismatch_dictionary[virus_name] = {}
        mismatch_df_per_virus = pd.DataFrame()
        # create a dictionary with all relevant information from bed and tsv file
        for bed_file, ref_file in zip(bed_files, ref_files):
            ref_id, ref_seq = read_fasta(ref_file)
            bed_df = pd.read_csv(bed_file, sep="\t", header=None)
            tsv_file = [filename for filename in tsv_files if virus_name in filename][0]  # a bit hacky way to get the correct file
            tsv_df = pd.read_csv(tsv_file, sep="\t", header=0)
            primer_mismatch_dictionary[virus_name][ref_id] = {}
            primer_mismatch_dictionary[virus_name][ref_id]["ref seq"] = ref_seq
            for ref, start, stop, name, direction in zip(bed_df[0], bed_df[1], bed_df[2], bed_df[3], bed_df[5]):
                primer_mismatch_dictionary[virus_name][ref][name] = {
                    "start": start,
                    "stop": stop,
                    "direction": direction
                }
                primer_mismatch_dictionary[virus_name][ref][name]["primer_seq"] = tsv_df[tsv_df["primer_name"] == name]["seq"].to_string(header=False, index=False)
        # read in variant files
        for variant_file, coverage_file, sample_name in zip(variant_files, coverage_files, sample_names):
            # read in variants and define seq to mutate
            variant_df = pd.read_csv(variant_file, sep="\t", header=0)
            coverage_df = pd.read_csv(coverage_file, sep="\t", header=0)
            seq_id = coverage_df["#chr"][0]
            seq_mut = list(str(primer_mismatch_dictionary[virus_name][seq_id]["ref seq"]))
            # mutate the ref seq
            for row in variant_df.iterrows():
                if row[1]["AF"] >= 0.7 and len(row[1]["REF"]) == 1 and len(row[1]["ALT"]) == 1:
                    seq_mut[row[1]["POS"]-1] = row[1]["ALT"]
            seq_mut = Seq("".join(seq_mut).lower())
            # remember for each primer the number of mismatches
            mismatch_list = []
            # check mismatches between primers and mutated seq
            for attribute in primer_mismatch_dictionary[virus_name][seq_id]:
                # get the sequences
                if attribute == "ref seq":
                    continue
                start, stop = primer_mismatch_dictionary[virus_name][seq_id][attribute]["start"], primer_mismatch_dictionary[virus_name][seq_id][attribute]["stop"]
                primer_coverage = list(coverage_df[(coverage_df["pos"] > start) & (coverage_df["pos"] <= stop)]["coverage"])  # tsv is one based and bed file is 0 based
                sample_seq = seq_mut[start:stop]
                # check if the primer region is sufficiently covered, if not the primer is not included in the analysis
                if not all(x >= 20 for x in primer_coverage) or len(primer_coverage) != len(sample_seq):
                    continue
                if primer_mismatch_dictionary[virus_name][seq_id][attribute]["direction"] == "-":
                    sample_seq = rev_complement(sample_seq)
                primer_seq = primer_mismatch_dictionary[virus_name][seq_id][attribute]["primer_seq"]
                # compare each position and count mismatches
                mismatches = 0
                for primer_nuc, seq_nuc in zip(primer_seq, sample_seq):
                    if primer_nuc in AMBIG_NUCS:
                        primer_nuc = AMBIG_NUCS[primer_nuc]
                    if seq_nuc not in primer_nuc:
                        mismatches += 1
                mismatch_list.append(mismatches)

            mismatch_df_per_virus = pd.concat([mismatch_df_per_virus, pd.Series(mismatch_list, name=sample_name).value_counts().to_frame()], axis=1)
        # plot results
        mismatch_df_per_virus = mismatch_df_per_virus.T.plot(
            kind="bar",
            stacked=True,
        )
        sns.despine()
        plt.xticks(rotation=45, ha="right")
        set_size(len(variant_files) * 0.35, 4.5)
        plt.legend(loc="lower left", title="number of mismatches", ncol=3, bbox_to_anchor=(0,1))
        plt.ylabel("primer target sequences covered >= 20x")
        plt.savefig(f"{output_folder}/{virus_name}_primer_mismatches.pdf", bbox_inches='tight')


def main(color_scheme, output_folder):
    """
    main function
    """
    warnings.filterwarnings("ignore")
    if exists(output_folder):
        shutil.rmtree(output_folder)
    makedirs(output_folder)
    print("### Starting the varVAMP in silico analysis ###\n")
    print("- Alignments to be analysed:")
    alignment_stats("alignments")
    print("-pairwise sequence identity:")
    identity_all = pairwise_evaluation("sequence_identity/all")
    print(identity_all[["virus", "mean", "std"]])
    sns.set_palette(palette=color_scheme, n_colors=7)
    print("- Plotting entropy...")
    calculate_and_plot_entropy(input_folder="alignments",
                               output_folder=output_folder)
    print("- Plotting primer degeneracy...")
    calculate_and_plot_degeneracy(input_folder="primer_tsv_files",
                                  output_folder=output_folder,
                                  colorscheme=color_scheme)
    print("- Plotting mismatches with alignment sequences...")
    calculate_and_plot_mismatches(alignment_folder="alignments",
                                  bed_folder="primer_bed_files",
                                  tsv_folder="primer_tsv_files",
                                  output_folder=output_folder)
    print("- Plotting primer stats...")
    plot_primer_stats(output_folder=output_folder,
                      tsv_folder="primer_tsv_files",
                      bed_folder="primer_bed_files",
                      consensus_folder="consensus_files")
    print("- Plotting identity comparison...")
    plot_sequence_identity_comparison(identity_all_df=identity_all,
                                      identity_folder="sequence_identity/new_seq",
                                      output_folder=output_folder)
    print("- Plotting per amplicon coverages...")
    plot_per_amplicon_coverages(coverages="coverages_per_amplicon",
                                output_folder=output_folder)
    print("- Plotting how well primers match to new sequences...")
    analyse_and_plot_primer_binding(adapted_bed_folder="adapted_bed_primer_files",
                                    ref_folder="reference_seq",
                                    variant_folder="variant_tables",
                                    per_base_coverages_folder="per_base_coverages",
                                    tsv_folder="primer_tsv_files",
                                    output_folder=output_folder)
    print("\n###         Finished the analysis          ###")


#primer_starts, primer_stops, primer_seqs, primer_names, alignment

bed_file = "primer_bed_files_olivar/7_ratHEV.bed"
alignment_file = "alignments/7_ratHEV.fasta"

tsv_file = "primer_tsv_files/7_ratHEV.tsv"
bed_file_2 = "primer_bed_files/7_ratHEV.bed"


alignment = read_alignment(alignment_file)
bed_df = pd.read_csv(bed_file, sep="\t", header=None)

mismatches, mean_mismatches, per_pos_mismatches_per_scheme = calculate_mismatches(bed_df[1]-1, bed_df[2], bed_df[6], bed_df[3], alignment)


tsv_df = pd.read_csv(tsv_file, sep="\t", header=0)
bed_df = pd.read_csv(bed_file_2, sep="\t", header=None)
mismatches_2, mean_mismatches_2, per_pos_mismatches_per_scheme_2 = calculate_mismatches(bed_df[1], bed_df[2], tsv_df["seq"],tsv_df["primer_name"], alignment)

print(f"\nvarvamp mismatches:\n {sum(mean_mismatches_2)/len(mean_mismatches_2)},\n\n olivar mismatches:\n {sum(mean_mismatches)/len(mean_mismatches)}")



alignment_folder = "alignments"
bed_folder = "primer_bed_files_olivar"

flattened_mismatches, per_pos_mismatch_list_all = [], []

# load data
alignment_files = get_files(alignment_folder)
bed_files = get_files(bed_folder)
names = get_file_names(bed_files)

# calculate all mismatches
for bed_file, alignment_file, name in zip(bed_files, alignment_files, names):
    bed_df = pd.read_csv(bed_file, sep="\t", header=None)
    alignment = read_alignment(alignment_file)
    mismatches, mean_mismatches, per_pos_mismatches_per_scheme = calculate_mismatches(
        bed_df[1]-1, bed_df[2], bed_df[6],bed_df[3], alignment
    )
    # flatten all mismatches to one list
    flattened_mismatches.append([item for sublist in mismatches for item in sublist])
    per_pos_mismatch_list_all.append(per_pos_mismatches_per_scheme)

# ini figure for mismatch bubble plot
fig, ax = plt.subplots(figsize=(4, 4.5))

for idx, mismatches in enumerate(flattened_mismatches):
    value, counts = np.unique(mismatches, return_counts=True)
    # normalize to total number
    counts_norm = [x / sum(counts) * 100 for x in counts]
    # plot to ax
    ax.scatter(x=[idx] * len(counts_norm), y=value, s=counts_norm)
# pseudodata for legend
legend_dot_size = ax.scatter([0, 0, 0], [0, 0, 0], s=[1, 10, 100], edgecolor=None, alpha=0.1)
legend = plt.legend(
    *legend_dot_size.legend_elements("sizes"),
    title="Percent",
    frameon=False,
    loc='upper center',
    ncols=3
)
ax.add_artist(legend)
ax.set_xticklabels(rotation=45, ha="right", labels=range(-1, len(counts_norm)))
ax.set_xticklabels([""] + names)
ax.set_ylabel("nt mismatches between primers and sequences")
sns.despine()
#fig.savefig(f"{output_folder}/mismatches.pdf", bbox_inches='tight')
plt.show()





# run the analysis
if __name__ == "__main__":
    main(color_scheme="coolwarm", output_folder="output")


