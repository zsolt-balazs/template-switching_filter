#!/usr/bin/env python3
import os
from argparse import ArgumentParser
import ts_gff
import deal_with_ts

def main():
    args = parsing()


    args.feature_file = args.prefix + "_ts_l3.tsv"
    args.feature = args.feature_file[-6:-4]
    args.coverage_file = args.prefix + "_out_allcov.tsv"
    deal_with_ts.Stats(args)

    args.feature_file = args.prefix + "_ts_r3.tsv"
    args.feature = args.feature_file[-6:-4]
    args.coverage_file = args.prefix + "_out_allcov.tsv"
    deal_with_ts.Stats(args)
    
    args.feature = "tes"
    ts_gff.ts_gff(args)

def parsing():
    parser = ArgumentParser(description="This is the deal_with_ts module of \
                            LoRTIA, a Long-read RNA-Seq Transcript Isofom \
                            Annotator. This module is designed to sort out \
                            3' end template switching.")
    parser.add_argument("prefix",
                        help="The path and the prefix of the statistics, \
                        which are to be used for the gff.",
                        metavar="prefix")
    parser.add_argument("-r", "--reference",
                        dest="reference",
                        help="The reference fasta file. Template-switching \
                        in the case of putative introns is going to be checked\
                        according to this file.",
                        default="/mnt/c/Work/LT907985.2/Ref/LT907985.2.fasta",
                        metavar="[reference_fasta]")
    parser.add_argument("-m", "--minimum",
                        dest="minimum", 
                        help="The minimal number of reads for the feature to\
                        be accepted.",
                        type=int,
                        default=2,
                        metavar="[integer]")
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
    parser.add_argument("-t", "--ratio",
                        dest = "ratio",
                        help = "The minimal ratio of the coverage that a \
                        feature has to reach to be accepted. The default value\
                        is 0.001.",
                        type=float,
                        default=0.001,
                        metavar="[float]")
    parser.add_argument("-d", "--distance",
                        dest="distance",
                        help="The distance from the feature position where \
                        coverage should be calculated. The default value is \
                        15. A positive value should be given, the inward \
                        direction is calculated by the program automatically.",
                        type=int,
                        default=15,
                        metavar="[integer]")
    parser.add_argument("--cov_sample",
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
                        default=(str(os.path.abspath(__file__)).replace("ts_launcher.py",""))+"dict.tsv",
                        metavar="[file_with_the_limit_values]")
    parser.add_argument("-g", "--check_surroundings",
                        dest="check_surroundings",
                        help="The number of nucleotides surrounding the \
                        feature where the number of real and ts hits are \
                        to be compared. The default value is 10.",
                        type=int,
                        default=10,
                        metavar="[integer]")
    return parser.parse_args()    

if __name__== "__main__":
    main()