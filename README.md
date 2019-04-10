# Template switching filtering algorithm
These are the scripts which were used during the study highlighting the importance of template-switching artefacts in the analysis of polyadenylation (currently under review).
The filtering algorithm can only be used to filter artefactual polyA sites from long-read cDNA sequencing data.

## Contents

- [Importance](#importance)
- [The filtering algorithm](#algorithm)
- [Usage](#usage)
- [Advanced options](#options)

## <a name="importance"></a>Importance
Long-read cDNA sequencing is now accessible for many research labs. However, cDNA sequencing, especially long-read cDNA sequencing is susceptible to template switching artefacts. Template switching is known to produce spurious [splice junctions] and [chimeric reads], but its effect on analysis of alternative polyadenylation is largely neglected. Artefactual poly(A) sites are generally attributed to [internal priming] and a common approach is to exclude potential poly(A) sites which contain six or more consecutive As or more than 12 As in a 20-nt window. We have found that artefacts can also arise at regions which contain fewer As, at times as few as three. We argue that such artefacts are not likely to be produced by internal priming, rather by template switching. We have developed a filtering algorithm that is more specific and more sensitive than the conventional filtering method used to decrease the impact of internal priming. The filtering algorithm is based on the [LoRTIA] toolkit and can therefore handle any kind of cDNA sequencing platform that can be processed by it (PacBio Iso-Seq or MinION cDNA sequencing data).

## <a name="algorithm"></a>The filtering algorithm
The algorithm uses the [LoRTIA] output TSV files and first generates a list of potential polyA sites. The LoRTIA toolkit creates a list of TESs from the sites, that were not marked as template-switching artefacts (the length of the short-homologous sequence is controlled by the `shs_for_ts` option). However, the toolkit discards every site that was only supported by reads which could have arisen through template switching. These sites can be reevaluated and filtered by this algorithm. The user has several options to determine which sites should be considered as potential polyA sites. The options and their default settings are henceforward described as `option [default]`. The polyA site has to be supported by at least `--minimum [2]` number of reads ending at the same genomic position and by at least `--ratio [0.001]` proportion of the overlapping reads. In a 2*`--wobble[10]`+1 nucleotide long window the site with the most supporting reads is considered the potential polyA site. If two sites in the window have the same highest number of supporting reads, the most downstream site is selected. It is recommended to use the same settings that were used when running the LoRTIA toolkit.

From the potential polyA sites, 
- the ones that are supported at least as many correct reads as reads marked as "potential template switching"  
OR 
- the ones that have more supporting reads than a certain `limit` proportion of the overlapping reads
are considered as TESs and the rest are marked as template-switching artefacts. 
In actuality, the number of correct reads can be modified, by the `--multiplier [1.0]` option. E.g. etting the `--multiplier` 2.0 will mean that the number of correct will be counted as double. That means that if there are 4 artefactual and 2 correct reads at a potential polyA site, the site is discarded at default settings, but accepted if the `--multiplier` is set 2.0. The `limit` is dependent on the number of adenines counted in the region around the poly(A) site. The default values are calculated by the following formula: 0.8/(1+2^(-100*(1/(20-A_count)-0.08))), where A_count is the number of As upstream of a polyA site. The values can be modified by changing the `dict.tsv` file. A different dict file can also be specified using the `--dictionary` option.

The adenine content of a region is determined using the reference FASTA file. The 20 nucleotides immediately upstream of the polyA site are iterated one-by-one. Each A counts as +1, while other nucleotides count -1. A counter is run over the iteration until the counter reaches -1 or till the end of the 20 nucleotides. The highest value of the counter is the adenine content of a region.
Examples:
```txt
+1 +1 +1 +1 -1 +1 +1 -1 +1 +1 -1 -1 -1 +1 -1 -1 -1 +1 -1 -1
 A  A  A  A  C  A  A  G  A  A  C  G  T  A  C  T  G  A  G  T
 1  2  3  4  3  4  5  4  5  6  5  4  3  4  3  2  1  2  1  0
A count: 6

+1 +1 +1 -1 -1 +1 -1 -1 +1 +1 +1 -1 +1 +1 -1 +1 -1 +1 -1 -1
 A  A  A  G  T  A  C  T  A  A  A  G  A  A  T  A  C  A  T  G
 1  2  3  2  1  2  1  0  1  2  3  2  3  4  3  4  3  4  3  2
A count: 4
 
+1 +1 +1 -1 -1 
 A  A  C  C  T  G  A  A  A  A  T  C  T  A  C  G  C  A  C  A
 1  2  1  0 -1
A count: 2
```

## <a name="usage"></a>Usage
Apply the scripts on the outputs of the [LoRTIA] toolkit. 
Download all the files from this repository, add them to your path, and run:
```sh
ts_launcher.py /path/to/LoRTIA-output/prefix -r /path/to/reference.fasta
```
The `prefix` is the prefix of the output files of the LoRTIA toolkit, e.g. the output BAM file would be `prefix_out_sorted.bam`. The reference file has to be the same that was used for the mapping. The outputs of the algorithm are saved to the same folder where the LoRTIA TSV files are. The outputs are `prefix_ts_l3_tes.tsv` and `prefix_ts_r3_tes.tsv`, which are the summary tables of the polyA sites in A-rich regions and `prefix_ts_tes.gff3` and `prefix_not_ts_tes.gff3`, which are the GFF files of the template-switching artefacts and the genuine TESs, respectively.

## <a name="options"></a>Advanced options
Apart from the options detailed in the The filtering algorithm chapter, several other options can be specified the user:
- The argument `distance [15]` specifies the distance upstream of the polyA site, where the coverage value is to be calculated.
- The coverage is averaged over a `cov_sample [5]` number of nucleotides. The coverage value is used as the number of reads overlapping a certain polyA site. The default settings mean that the coverages of the nucleotides 19 to 15 nucleotides upstream of a TES are averaged to form the coverage value.

[LoRTIA]: https://github.com/zsolt-balazs/LoRTIA
[splice junctions]: https://www.sciencedirect.com/science/article/pii/S0888754305003770
[chimeric reads]: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012271
[internal priming]: https://www.pnas.org/content/99/9/6152
