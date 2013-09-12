sv-annotate
===========

A small pipeline to annotate and filter candidate structural variation calls, particularly
those made by BreakDancer.

Genomic structural variation callers like BreakDancer, DELLY, and HYDRA produce many 
false positive calls in addition to true positive calls. The goal of this pipeline
is to filter calls to produce a smaller set of higher-quality candidate calls for
wet-lab validation.

Much of the inspiration for this came from Aaron Quinlan's [description of a
typical workflow for use with his HYDRA SV caller](https://code.google.com/p/hydra-sv/wiki/TypicalWorkflow), and from his paper describing HYDRA.

Most of the processing takes place in the `annotate_bedpe.py` script. This script takes input in BEDPE format. If you are starting with BreaKDancer calls, you can convert them to BEDPE by using the `convertBreakDancerBedToBedpe.py` script.
The annotation script then categorizes the calls into their own separate files based on their implied SV type and overlap with annotation features. 
The script uses bedpe format so that it can use the "anchoring regions" of each call - the regions in which the supporting reads aligned - as well as
the actual called affected region (the "inner span" of the bedpe).

Calls are first divided based on orientation, length and inter/intra chromosome status to separate out candidate deletions,
insertions, tranlocations, and inversions. Then based on the SV event size, score, and overlaps with annotated features they are further
categorized. The remaining calls of each type in the "Stringent" category should hypothetically contain the most true positives.

The annotation files are:

`te_file`: This is a track of transposable elements. Currently we are using the repeatmasker track downloaded from UCSC for hg19. This file
is used in the following ways:

1. Apparent translocations in which one anchoring region overlaps a TE are more likely to be TE insertions in the sample as opposed to the reference, and the mapping to a distant TE is likely an alignment artifact.
2. Deletions where the deleted portion overlaps a TE are more likely TE insertions in the reference than actual deletions in the sample.

Optional annotation files:

* `common_deletions`: This is meant to be a set of variants commonly occuring in the population (and therefore less likely to be implicated in cancer) Currently taken from the 1000 Genomes project SV pilot data.
* `seg_dups`: segmental duplications track from UCSC
* `cent_tel`: centromeric and telomeric regions, currently just taken from UCSC by adding 100kb of flanking sequence to their centromere/telomere annotations

Other parameters:

* `insert_size`: expected insert size of the library
* `output_dir`: directory to be populated with annotated files
* `sample_name`: name of the sample to be added to track names in bed files
