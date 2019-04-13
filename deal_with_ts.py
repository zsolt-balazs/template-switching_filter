#!/usr/bin/env python3

import pandas as pd
import os
from argparse import ArgumentParser
from ast import literal_eval
from Bio import SeqIO

def get_cov(position, cov_dict):
    """
    Calculates coverage based on a cov_dict
    """
    if position in cov_dict:
        return cov_dict.get(position)
    else:
        return 0

def coverage(pos_list, args, contig):
    """
    Calculates average coverages from a given distance for a list of positions
    """
    cols = ["contig", "pos", "count"]
    cov = pd.read_csv(args.coverage_file, sep = "\t", names= cols)
    cov = cov.loc[cov.contig == contig]
    cov["pos"] = cov["pos"].astype(int)
    cov["count"] = cov["count"].astype(int)
    cov.set_index("pos", drop=True, inplace=True)
    cov_dict = cov.to_dict()["count"]
    coverages = []
    for pos in pos_list:
        to_avg = []
        if args.distance > args.cov_sample:
            for position in range(pos + args.distance - args.cov_sample,
                                  pos + args.distance):
                to_avg.append(get_cov(position, cov_dict))
        elif args.distance < args.cov_sample:
            for position in range(pos + args.distance,
                                  pos + args.distance - args.cov_sample):
                to_avg.append(get_cov(position, cov_dict))
        else:
            to_avg.append(get_cov(pos + args.distance, cov_dict))
        coverages.append(sum(to_avg)/len(to_avg))
    return coverages

def check_if_qualified(df, minimum, ratio):
    """
    Checks whether the feature position satisfies the minimum count and minimum
    ratio of coverage requirements.
    """
    qual_list = []
    for index, row in df.iterrows():
        is_qual = (row["count"] >= minimum
                   and row["ratio"] >= ratio
                   and row["is_picked"])
        qual_list.append(is_qual)
    return qual_list

def pick_from_greatests(dictionary, wobble):
    """
    Picks the left- or rightmost positions of the greatests list in a window
    determined by the wobble size. Whether the left or the rightmost positions
    are desired can be set by the user, and the list is ordered accordingly.
    """
    previous = -100
    is_picked_list = []
    for pos, is_greatest in dictionary.items():
        is_picked = False
        if is_greatest:
            if previous not in range(pos - wobble, pos + wobble + 1):
                is_picked = True
            previous = pos
        is_picked_list.append(is_picked)
    return is_picked_list

def get_As(seq_record, position, args):
    if args.strand < 0:
        c = 1
        d = 19
    else:
        c = 20
        d = 0
    seq = seq_record.seq[position - c:position + d].lower()
    if args.strand < 0:
        seq = seq.reverse_complement()
    countAs = 0
    count_list = []
    for char in seq[::-1]:
        if char == "a":
            countAs += 1
            count_list.append(countAs)
        else:
            countAs += -1
            count_list.append(countAs)
        if countAs < 0:
            break
    countAs = max(count_list)
    if countAs < 0:
        countAs = 0
    return countAs

def check_if_greatest(tuples, wobble):
    """
    Finds the feature position with the highest read support in a window
    determined by the wobble size.
    """
    is_greatest_list = []
    for tup in tuples:
        pos, count = tup
        for position, count2 in tuples:
            is_greatest = True
            if position in range(pos - wobble, pos + wobble + 1):
                is_greatest = count >= count2
                if not is_greatest:
                    break
        is_greatest_list.append(is_greatest)
    return is_greatest_list

def count_average(tuples, window):
    in_window = []
    for tup in tuples:
        pos, count = tup
        hundred = 0
        for position, count2 in tuples:
            if position in range(pos - window, pos + window + 1):
                hundred += count2
        in_window.append(hundred/(2 * window + 1))
    return in_window

def get10(df, countfile, args):
    count_dict = {}
    with open(countfile) as cfile:
        for line in cfile:
            count_dict[eval(line.split("\t")[0])] = int(line.split("\t")[1])
    dictio = {}
    sums = []
    for index, row in df.iterrows():
        a = []
        for pos in range(row["pos"] - args.check_surroundings, row["pos"] 
                         + args.check_surroundings + 1):
            if (row["contig"], pos) in count_dict:
                a.append(count_dict[row["contig"], pos])
            else:
                a.append(0)
        if args.feature in ["r5", "l3"]:
                a = a[::-1]
        dictio[row["pos"]] = a
        sums.append(sum(a))
    dictio = pd.DataFrame(dictio).T
    if countfile == args.feature_file:
        mark = "f"
    else:
        mark = "r"
    dictio = dictio.rename(columns=lambda x: mark + str(x - args.check_surroundings))
    dictio["pos"] = dictio.index
    df = df.merge(dictio, how="left", on="pos")
    df["{}sum".format(mark)] = sums
    return df

def contig_ends(df, args, contig):
    """
    Processes the dataframe for one contig looking for TSSs or TESs.
    """
    df["pos"] = pd.to_numeric(df["pos"])
    if args.feature == "l5" or args.feature == "l3":
    # This makes sure that the leftmost position is taken for each left feature
        df = df.sort_values(by="pos")
    else:
        df = df.sort_values(by="pos", ascending=False)
    df["coverage_before"] = coverage(df["pos"], args, contig)
    args.distance = args.distance * -1
    args.cov_sample = args.cov_sample * -1
    df["coverage_after"] = coverage(df["pos"], args, contig)
    pos_count = list(zip(df["pos"], df["count"]))
    df["average"] = count_average(pos_count, 50)
    df["is_greatest"] = check_if_greatest(pos_count, args.wobble)
    greatests_dict = dict(zip(df["pos"], df["is_greatest"]))
    df["is_picked"] = pick_from_greatests(greatests_dict, args.wobble)
    df["ratio"] = df["count"] / df["coverage_before"]
    df["is_qualified"] = check_if_qualified(df, args.minimum, args.ratio)
    df = get10(df, args.feature_file, args)
    df = get10(df, args.feature_file.replace("_ts", ""), args)
    args.distance = args.distance * -1
    args.cov_sample = args.cov_sample * -1
    return df

def find_features(args):
    """
    Reads in a dataframe from csv, chops it and processes it on contigs.
    """
    dictionary_df = pd.read_csv(args.dictionary, sep = "\t", )
    dictionary = dictionary_df.to_dict()
    df = pd.read_csv(args.feature_file, sep = "\t", names = ["pos", "count"])
    df["pos"] = df["pos"].apply(literal_eval)
    df[["contig", "pos"]] = df["pos"].apply(pd.Series)
    contig_set = list(set(df.contig))
    new_df = ()
    for contig in contig_set:
        current_df = df.loc[df.contig == contig].copy()
        current_df = contig_ends(current_df, args, contig)
        if len(new_df) != 0:
            new_df = pd.concat([new_df, current_df])
        else:
            new_df = current_df
    feat_list = []
    cont_name = 0
    cont = 0
    A_list = []
    for index, row in new_df.iterrows():
        if row["is_qualified"]:
            if cont_name == row["contig"]:
                countAs = get_As(cont, row["pos"], args)
            else:
                for seq_record in SeqIO.parse(args.reference, "fasta"):
                    if seq_record.name == row["contig"]:
                        countAs = get_As(seq_record, row["pos"], args)
                        cont_name = seq_record.name
                        cont = seq_record
                        break
            if row["rsum"] > 1:
                multiplier = args.multiplier
            else:
                multiplier = 0
            limit = dictionary["limit"].get(countAs)
            # this avoids division by 0 in the next lines:
            if row["coverage_before"] == 0:
                zero = 1
            else:
                zero = 0
            if ((row["rsum"] + row["fsum"]) * 100 > row["coverage_before"] 
                    - row["coverage_after"]) and ((row["rsum"] 
                    * multiplier >= row["fsum"]) or (((row["rsum"] 
                    + row["fsum"])/(row["coverage_before"] + zero)) > limit)):
                if args.feature[1] == "3":
                    feat = "tes"
                else:
                    feat = "tss"
            else:
                if args.feature[1] == "3":
                    feat = "Template-switching"
                else:
                    feat = "tss Template-switching"
        else:
            feat = None
            countAs = -2
        feat_list.append(feat)
        A_list.append(countAs)
    new_df["A_list"] = A_list
    new_df["feature"] = feat_list
    if args.feature[1] == "3":
        feat = "_tes"
    elif args.feature[1] == "5":
        feat = "_tss"
    else:
        feat = "tron"
    new_df.to_csv(args.feature_file.replace(".tsv", "{}.tsv".format(feat)),
                  index=False,
                  sep="\t")

def Stats(args):
    """
    Sets argument types and runs stat functions for features.
    """
    if not args.feature:
        args.feature = args.feature_file[-6:-4]
    print("Calculating {} feature statistics...".format(args.feature))
    if args.feature == "r5" or args.feature == "l3":
        args.strand = -1
    else:
        args.strand = 1
    if args.feature[1] == "3":
        args.distance = abs(args.distance) * args.strand * (-1)
        args.cov_sample = abs(args.cov_sample) * args.strand * (-1)
    else:
        args.distance = abs(args.distance) * args.strand
        args.cov_sample = abs(args.cov_sample) * args.strand
    if os.stat(args.feature_file).st_size == 0:
        print("Feature file {} is empty. There is nothing to do here.".format(
              args.feature_file))
    else:
        find_features(args)

###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    Stats(args)

def parsing():
    parser = ArgumentParser(description="This is the deal_with_ts module of \
                            LoRTIA, a Long-read RNA-Seq Transcript Isofom \
                            Annotator. This module is designed to sort out \
                            3' end template switching.")
    parser.add_argument("coverage_file",
                        help="The tsv file which contains the coverages.\
                        The tsv file should contain 3 columns: contig, \
                        position and coverage.",
                        metavar="coverage_file")
    parser.add_argument("feature_file",
                        help="A tab-separated values file containing feature\
                        statistics produced by the Samprocessor.",
                        metavar="feature_file")
    parser.add_argument("-m", "--minimum",
                        dest="minimum", 
                        help="The minimal number of reads for the feature to\
                        be accepted.",
                        type=int,
                        default=2,
                        metavar="[integer]")
    parser.add_argument("-f", "--feature",
                        dest="feature",
                        help="The feature that is examined. Options are \
                        'r5' for reverse strand 5' ends, 'l3' for \
                        reverse strand 3' ends, 'l5' for forward strand 5'\
                        ends, 'r3' for forward strand 3' ends and 'in' for \
                        introns. By default the tsv file's last two characters\
                        before the .tsv extension are considered.",
                        default=False,
                        metavar="[string]")
    parser.add_argument("-b", "--wobble",
                        dest="wobble",
                        help="The window, in which only one of each feature \
                        is expected, and locations with lesser support are \
                        considered to be derivatives of the major. The default\
                        value is 10, which means that only one feature of a \
                        kind can be described in a 21 nt bin (location +/-10 \
                        nt). This only applies to TSSs and TESs.",
                        type=int,
                        default=10,
                        metavar="[integer]")
    parser.add_argument("-r", "--reference",
                        dest="reference",
                        help="The reference fasta file. Template-switching \
                        in the case of putative introns is going to be checked\
                        according to this file.",
                        default="/mnt/c/Work/LT907985.2/Ref/LT907985.2.fasta",
                        metavar="[reference_fasta]")
    parser.add_argument("-t", "--ratio",
                        dest = "ratio",
                        help = "The minimal ratio of the coverage that a \
                        feature has to reach to be accepted. The default value\
                        is 0.001.",
                        type=float,
                        default=0.001,
                        metavar="[float]")
    parser.add_argument("-l", "--multiplier",
                        dest="multiplier",
                        help="This variable defines how many times more \
                        reads ought to end in the +/- [check_surroundings] nt vicinity of \
                        the given nucleotide that could not have been produced\
                        by template switching, compared to reads that could.\
                        The default value is 2.0, that means there has to be \
                        at least as many reads without template switching as \
                        with in a +/- [check_surroundings] nt vicinity in order for the\
                        [feature] to be called. Otherwise 'Template switching'\
                        is called.",
                        type=float,
                        default=1.0,
                        metavar="[float]")
    parser.add_argument("-d", "--distance",
                        dest="distance",
                        help="The distance from the feature position where coverage \
                        should be calculated. The default value is 15. A positive value\
                        should be given, the inward direction is calculated by the \
                        program automatically.",
                        type=int,
                        default=15,
                        metavar="[integer]")
    parser.add_argument("-s", "--cov_sample",
                        dest="cov_sample",
                        help="The number of nucleotides where the coverage \
                        should be averaged. This many consecutive nucleotides\
                        will be considered from the 'distance' towards the \
                        feature. Its absolute value has to be smaller than \
                        or equal to the value of 'distance'. The default value\
                        is 5.",
                        type=int,
                        default=5,
                        metavar="[integer]")
    parser.add_argument("-y", "--dictionary",
                        dest="dictionary",
                        help="The .tsv file that contains the percentage of \
                        overlapping reads that have to support a potential \
                        polyA site in order for it to be accepted as a TES.\
                        Leave it empty if you want to use the dict.tsv in your\
                        working directory. Provide a full path if you want to \
                        use a different file.",
                        type=str,
                        default=(str(os.path.abspath(__file__)).replace(
                                "deal_with_ts.py","")) + "dict.tsv",
                        metavar="[string]")
    parser.add_argument("-g", "--check_surroundings",
                        dest="check_surroundings",
                        help="The number of nucleotides surrounding the \
                        feature where the number of real and ts hits are \
                        to be compared. The default value is 10.",
                        type=int,
                        default=10,
                        metavar="[integer]")
    if not parser.parse_args().feature:
        parser.parse_args().feature = parser.parse_args().feature_file[-6:-4]
    return parser.parse_args()
    
if __name__== "__main__":
    main()