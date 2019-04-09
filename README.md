# template-switching_filter
These are the scripts which were used during the study highlighting the importance of template-switching artefacts in the analysis of polyadenylation #currently under review
The filtering algorithm can only be used to filter long-read cDNA sequencing data.

### Importance
Long-read cDNA sequencing is now accessible for many research labs. However, cDNA sequencing, especially long-read cDNA sequencing is susceptible to template switching artefacts. Template switching is known to produce spurious [splice junctions] and [chimeric reads]. A common approach is to exclude 


### Usage
Apply the scripts on the outputs of the [LoRTIA] toolkit. 
```sh
ts_launcher.py /path/to/LoRTIA-output/prefix -r /path/to/reference.fasta
```
The `prefix` is the prefix of the output files of the LoRTIA toolkit, e.g. the output BAM file would be `prefix_out_sorted.bam`. The reference file has to be the same that was used for the mapping. 

[LoRTIA]: https://github.com/zsolt-balazs/LoRTIA
[splice junctions]: https://www.sciencedirect.com/science/article/pii/S0888754305003770
[chimeric reads]: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0012271
