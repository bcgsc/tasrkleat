# tasrkleat

### Version v0.2.2 (2018-XX-XX)

- added --dry-run option
- removed broken code on `R.pipeline_printout_graph`

### Version v0.2.1 (2018-07-11)

- added --batch 4 to gmap to speed up contigs2genome alignment
- added -t to bwa mem to speed up reads2contigs alignment

### Version v0.2 (2018-07-09)

- removed transfer step, leave it to Google Genomics API. This step is creating
  trouble somehow for large input files (e.g. with >2k genes), see the description of
  the problem on [GGP
  forum](https://groups.google.com/forum/#!topic/google-genomics-discuss/RQBscD6YSjk).
- added --batch 4 to gsnap to speed up reads2genome alignment

	
### Version v0.1.2 (2017-12-12)

- Updated Readme.md and ChangeLog.md, no code chagne.

### Version v0.1.1 (2017-12-11)

- Adjusted Gsnap arguments for better read2genome alignment (add --novelsplicing=1), for use of benchmarking DaPars
- Don't remove filtered reads by Biobloomfilter, for later use of benchmarking ContextMap2

### Version v0.1 (2017-01-06)

First release.
