"""
this script contains the code for the in silico varVAMP evaluation
"""

# import libraries
import shutil
import math
import warnings
from os import listdir, makedirs
from os.path import isfile, join, basename, splitext, exists
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
import primer3 as p3
import itertools
import statistics


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
    return [splitext(basename(f))[0].replace("_", " ") for f in files]


def read_fasta(c_file):
    """
    reads lines and return the line after the header
    """
    with open(c_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            return Seq(line.rstrip("\n"))


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
                results.append([name, statistics.mean(identities), statistics.stdev(identities)])
            else:
                results.append([name, identities[0], 0])

    return pd.DataFrame(results, columns=["virus", "mean", "std"])


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
        if "RW" in primer_name:
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
            linewidth=1,
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


def calculate_and_plot_degeneracy(input_folder, output_folder):
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
    ax = sns.stripplot(permutations)
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
    plt.ylabel("degeneracy per primer")
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
        consensus_seq = read_fasta(consensus_seq)
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

    alignment_identity, alignment_std, text = [], [], []

    for name in identity_comparison["virus"]:
        if name in list(identity_all_df["virus"]):
            alignment_identity_temp = float(identity_all_df[identity_all_df["virus"] == name]["mean"])
            identity_diff = alignment_identity_temp - float(identity_comparison[identity_comparison["virus"] == name]["mean"])
            alignment_identity.append(alignment_identity_temp)
            text.append(f"{round(identity_diff)}%")
        else:
            alignment_identity.append("NA"), alignment_std.append("NA"), text.append("NA")
    # add columns
    identity_comparison["alignment_mean"], identity_comparison["change"] = alignment_identity, text

    # dumbbell plot
    plt.figure(figsize=(9, 4.5))
    colors = np.where(identity_comparison["mean"] < identity_comparison["alignment_mean"], "#d9d9d9", "#d57883")
    # plot hlines with colors dependent if the mean identity is higher or lower
    plt.hlines(y=identity_comparison["virus"], xmin=identity_comparison["mean"],
               xmax=identity_comparison["alignment_mean"],
               color=colors, lw=10)
    # generate x,y tuple for annotation
    xy_values = [(x+1, y) for x, y in zip(identity_comparison[["mean", "alignment_mean"]].max(axis=1), np.arange(0, len(identity_comparison.index) + 1))]
    # annotate change
    for annotation, xy in zip(list(identity_comparison["change"]), xy_values):
        plt.annotate(annotation, xy, verticalalignment="center")
    plt.scatter(identity_comparison["mean"], identity_comparison["virus"], color="#0096d7", s=200,
                label="new sequences", zorder=3)
    plt.scatter(identity_comparison["alignment_mean"], identity_comparison["virus"], color="grey", s=200,
                label="alignment", zorder=3)
    sns.despine()
    plt.legend(ncol=2, bbox_to_anchor=(1., 1.01), loc="lower right", frameon=False)
    plt.xlim(right=100)
    plt.xlabel("% mean pairwise identity")
    plt.tight_layout()

    plt.savefig(f"{output_folder}/identity_comparision.pdf", bbox_inches='tight')


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
    print(identity_all)
    sns.set_palette(palette=color_scheme)
    print("- Plotting entropy...")
    calculate_and_plot_entropy("alignments", output_folder)
    print("- Plotting primer degeneracy...")
    calculate_and_plot_degeneracy("tsv_files", output_folder)
    print("- Plotting mismatches...")
    calculate_and_plot_mismatches("alignments", "bed_files", "tsv_files", output_folder)
    print("- Plotting primer stats...")
    plot_primer_stats(output_folder, "tsv_files", "bed_files", "consensus_files")
    print("- Plotting identity comparison...")
    plot_sequence_identity_comparison(identity_all, "sequence_identity/new_seq", output_folder)

    print("\n###         Finished the analysis          ###")


# run the analysis
if __name__ == "__main__":
    main(color_scheme="PuOr", output_folder="output")
