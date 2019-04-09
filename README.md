# template-switching_filter
These are the scripts which were used during the template-switching study #currently under review

Apply the scripts on the outputs of the LoRTIA toolkit [LoRTIA]. 

###Usage
```sh
ts_launcher.py /path/to/LoRTIA-output/prefix -r /path/to/reference.fasta
```
The `prefix` is the prefix of the output files of the LoRTIA toolkit, so the output BAM file would be `prefix_out_sorted.bam`. The reference file has to be the same that was used for the mapping. 

[LoRTIA]: https://github.com/zsolt-balazs/LoRTIA
