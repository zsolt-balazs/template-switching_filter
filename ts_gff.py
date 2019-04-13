#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd

def line_end(df, n_df, feature, sign):
    """
    Prepares gff line for transcript ends (TSS and TES).
    """
    for index, row in df.iterrows():
        i = len(n_df)
        if row["feature"] == feature:
            n_df.loc[i] = [row["contig"],
                             "LoRTIA", 
                             feature,
                             row["pos"],
                             row["pos"],
                             row["count"],
                             sign,
                             ".",
                             row["count"]]
            i += 1

def ts_gff(args):
    """
    Creates template switching gffs from stats files.
    """
    print("Creating {} template-switching gff files...".format(args.feature))
    cols = ["contig",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "info"]
    new_df = pd.DataFrame(columns=cols)
    if args.feature == "tss":
        filepos = "{}_ts_l5_{}.tsv".format(args.prefix, args.feature)
        fileneg = "{}_ts_r5_{}.tsv".format(args.prefix, args.feature)
    else:
        filepos = "{}_ts_r3_{}.tsv".format(args.prefix, args.feature)
        fileneg = "{}_ts_l3_{}.tsv".format(args.prefix, args.feature)
    dfpos = pd.read_csv(filepos, sep = "\t")
    dfneg = pd.read_csv(fileneg, sep = "\t")
    line_end(dfpos, new_df, args.feature, "+")
    line_end(dfneg, new_df, args.feature, "-")
    ts_df = pd.DataFrame(columns=cols)
    line_end(dfpos, ts_df, "Template-switching", "+")
    line_end(dfneg, ts_df, "Template-switching", "-")
    feat_gff = pd.read_csv("{}_{}.gff3".format(args.prefix, args.feature),
                           sep="\t",
                           header=None,
                           names=cols)
    alldfs = pd.concat([new_df, ts_df, feat_gff])
    summary = pd.DataFrame(columns=cols)
    for strand in set(alldfs["strand"]):
        stranddf = alldfs.loc[alldfs.strand == strand].copy()
        stranddf = stranddf.sort_values(by="start")
        counter = 0
        ID_list = []
        prev = -999
        for index, row in stranddf.iterrows():
            if prev not in range(row["start"] - args.wobble, row["start"] 
                             + args.wobble + 1):
                counter += 1
            ID_list.append(strand+str(counter))
            prev = row["start"]
        stranddf["ID"] = ID_list
        is_greatest_list = []
        for index, row in stranddf.iterrows():
            maximum = stranddf.loc[stranddf.ID == row["ID"]]["score"].max()
            is_greatest = row["score"] == maximum
            is_greatest_list.append(is_greatest)
        stranddf["is_greatest"] = is_greatest_list
        if (strand == "-" and args.feature == ("tes")) or (strand == "+" 
           and args.feature == ("tss")):
            leftmost = True
        else:
            leftmost = False
        is_picked_list = []
        if leftmost:
            for index, row in stranddf.iterrows():
                if row["is_greatest"]:
                    subdf = stranddf.loc[(stranddf.ID == row["ID"]) & (
                            stranddf.is_greatest == True)]
                    is_picked = row["start"] == subdf["start"].min()
                else:
                    is_picked = False
                is_picked_list.append(is_picked)
        else:
            for index, row in stranddf.iterrows():
                if row["is_greatest"]:
                    subdf = stranddf.loc[(stranddf.ID == row["ID"]) & (
                            stranddf.is_greatest == True)]
                    is_picked = row["start"] == subdf["start"].max()
                else:
                    is_picked = False
                is_picked_list.append(is_picked)
        stranddf["is_picked"] = is_picked_list
        chart = stranddf.loc[stranddf.is_picked == True].copy()
        summary = pd.concat([summary, chart], ignore_index=True, sort=False)
    summary = summary.sort_values(by=['start'])
    summary = summary.sort_values(by=['contig'])
    summary = summary.drop(["ID", "is_greatest", "is_picked"], 1)
    ts = summary.loc[summary.feature == "Template-switching"]
    ts.to_csv("{}_ts_{}.gff3".format(args.prefix, args.feature + "w" 
              + str(args.wobble)),
              index=False,
              header=False,
              sep="\t")
    feature = summary.loc[summary.feature == args.feature]
    feature.to_csv("{}_not_ts_{}.gff3".format(args.prefix, args.feature),
                   index=False,
                   header=False,
                   sep="\t")


###############################################################################
###                             Main function                               ###
###############################################################################

def main():
    args = parsing()
    ts_gff(args)

def parsing():
    parser = ArgumentParser(description="This is the third module of \
                            LoRTIA, a Long-read RNA-Seq Transcript Isofom \
                            Annotator. This module creates gff files from\
                            the feature statistics.")
    parser.add_argument("prefix",
                        help="The path and the prefix of the statistics, \
                        which are to be used for the gff.",
                        metavar="prefix")
    parser.add_argument("feature",
                        help="The type of feature for which the gff is \
                        generated. Options include 'tss' for transcriptional \
                        start sites, 'tes' for transcriptional end sites and \
                        'intron' for introns.",
                        metavar="feature")
    parser.add_argument("-w", "--wobble",
                        dest="wobble",
                        help="The wobble that is to be used to merge TSS and \
                        TES feature positions. The default is +/- 10nt.",
                        type=int,
                        default=10,
                        metavar="[integer]")
    return parser.parse_args()


if __name__== "__main__":
    main()
