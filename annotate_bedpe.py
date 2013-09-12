#!/usr/bin/env python

# This script takes input in bedpe format, ie from converting breakdancer output to bedpe format with convertBreakdancerBedToBedpe.py,
# and categorizes the calls into their own separate files based on their implied SV type and overlap with annotation features. The script
# uses bedpe format so that it can use the "anchoring regions" of each call - the regions in which the supporting reads aligned - as well as
# the actual called affected region (the "inner span" of the bedpe).

# Calls are first divided based on orientation, length and inter/intra chromosome status to separate out candidate deletions,
# insertions, tranlocations, and inversions. Then based on the SV event size, score, and overlaps with annotated features they are further
# categorized. The remaining calls of each type in the "Stringent" category should hypothetically contain the most true positives.

# The annotation files are:

# te_file: This is a track of transposable elements. Currently we are using the repeatmasker track downloaded from UCSC for hg19. This file
# is used in the following ways:
#
# 1) Apparent translocations in which one anchoring region overlaps a TE are more likely to be TE insertions in the sample as opposed to the
# reference, and the mapping to a distant TE is likely an alignment artifact.
# 2) Deletions where the deleted portion overlaps a TE are more likely TE insertions in the reference than actual deletions in the sample.
#
# Optional annotation files:
#
# common_deletions: This is meant to be a set of variants commonly occuring in the population (and therefore less likely to be implicated in
# cancer) Currently taken from the 1000 Genomes project SV pilot data.
#
# seg_dups: segmental duplications track from UCSC
#
# cent_tel: centromeric and telomeric regions, currently just taken from UCSC by adding 100kb of flanking sequence to their centromere/telomere annotations
#
# Other parameters:
#
# insert_size: expected insert size of the library
# output_dir: directory to be populated with annotated files
# sample_name: name of the sample to be added to track names in bed files

import sys
import os
import pybedtools
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("bedpe_file", help="BEDPE file of sv call anchoring regions")
parser.add_argument("insert_size", type=int, help="Insert size of the library")
parser.add_argument("output_dir", help="directory into which to put annotations")
parser.add_argument("sample_name", help="name of the sample being analyzed")
parser.add_argument("te_file", help="BED file of transposable element annotations")
parser.add_argument("--seg_dups_file", help="BED file of segmental duplications")
parser.add_argument("--cent_tel_file", help="BED file of centromeric/telomeric regions")
parser.add_argument("--common_deletions_file", help="BED file of common deletions")

args = parser.parse_args()

bedpe_file = args.bedpe_file
te_file = args.te_file
common_deletions_file = args.common_deletions_file
insert_size = args.insert_size
output_dir = args.output_dir
seg_dups_file = args.seg_dups_file
cent_tel_file = args.cent_tel_file
sample_name = args.sample_name

print 'analyzing bedpe file: ' + bedpe_file

log = open(output_dir + "/annotate.log", "w")
log.write("input file: {0}\n".format(bedpe_file))
log.write("sample name: ".format(sample_name))
log.write("te file: {0}\n".format(te_file))
log.write("common deletions file: {0}\n".format(common_deletions_file))
log.write("insert size: {0}\n".format(insert_size))
log.write("segmental duplications file: {0}\n".format(seg_dups_file))
log.write("centromeres and telomeres file: {0}\n".format(cent_tel_file))

common_deletion_overlap_pct = 0.5

# filter out calls that are longer than length
def bedpe_lt_length_filter(feature, length):
    if int(feature[4]) - int(feature[2]) < length:
        return True
    return False

# filter out calls that are shorter than length
def bedpe_gt_length_filter(feature, length):
    if int(feature[4]) - int(feature[2]) > length:
        return True
    return False

# filter out calls that aren't transchromosomal
def inter_chr_filter(feature):
    if feature[0] == feature[3]:
        return False
    return True

# filter out calls that are transchromosomal
def intra_chr_filter(feature):
    if feature[0] == feature[3]:
        return True
    return False

# filter out calls that have a score less than score
def score_gte_filter(feature, score):
    if int(feature[7]) >= score:
        return True
    return False

# filter out calls that have a score greater than score
def score_lt_filter(feature, score):
    if int(feature[7]) < score:
        return True
    return False

# filter out calls where the orientation of the reads in the clusters isn't the expected PE
def expected_orientation_filter(feature, matches):
    expected = False
    if (feature[8] == '+' and feature[9] == '-') or (feature[8] == '-' and feature[9] == '+'):
        expected = True
    if matches:
        return expected
    else:
        return not expected

# filter out calls by name
def name_not_in_set_filter(feature, name_set):
    if feature[6] in name_set:
        return False
    else:
        return True

# filter out calls where the ends match a te feature (added to the feature by pair_to_bed)
def bedpe_reciprocal_overlap_ends_filter(feature, overlap_pct):
    te_chr = feature[22]
    te_start = int(feature[23])
    te_end = int(feature[24])
    te_length = te_end - te_start
    return overlaps_by(te_chr, te_start, te_end, feature[0], int(feature[1]), int(feature[2]), overlap_pct) or overlaps_by(te_chr, te_start, te_end, feature[3], int(feature[4]), int(feature[5]), overlap_pct)

# filter out calls where the inner span overlaps a TE
def bedpe_reciprocal_overlap_ispan_filter(feature, overlap_pct):
    te_chr = feature[22]
    te_start = int(feature[23])
    te_end = int(feature[24])
    te_length = te_end - te_start
    return overlaps_by(te_chr, te_start, te_end, feature[0], int(feature[2]), int(feature[4]), overlap_pct) 

# Determine if the the two features overlap by overlap_pct
def overlaps_by(chr1, start1, end1, chr2, start2, end2,  overlap_pct):
    if chr1 != chr2:
        return False
    max_start = max(start1,start2)
    min_end = min(end1,end2)
    overlap_len = min_end - max_start
    len1 = end1 - start1
    len2 = end2 - start2
    return float(overlap_len) / len1 >= overlap_pct and float(overlap_len) / len2 >= overlap_pct

# write out a bed file
def write_bed(call, fh):
    for f in call.fields:
        fh.write(str(f) + "\t")
    fh.write("\n")

# find duplicate calls by looking for pairs of calls where the both ends match with a little slop
def merge_duplicate_breaks(calls, slop):
    dups = calls.pair_to_pair(calls, type="both", rdn=True, slop=slop).saveas()
    low_scoring_dups = set()
    for call in dups:
        n1 = call[6]        
        n2 = call[28]
        s1 = int(call[7])
        s2 = int(call[29])
        low_score_name = ""
        if s1 < s2:
            low_score_name = n1
        else:
            low_score_name = n2
        low_scoring_dups.add(low_score_name)
    return calls.filter(name_not_in_set_filter, low_scoring_dups).saveas()

# call script to convert bedpe to blocked bed12 format for a visualizable track
def convert_bedpe_to_bed12(bedpe_file, track_name):
    bed12_file = open(bedpe_file + ".bed", 'w')
    subprocess.call("bedpeToBed12.py -i {0} -d 1000000000 -n \"{1}\"".format(bedpe_file, track_name), shell=True, stdout=bed12_file)

# use sort -u to get rid of duplicate bed entries
def uniqify(bedtool, output_dir, file_name):
    ufile = open(output_dir + "/" + file_name, 'w') 
    subprocess.call("sort -u " + output_dir + "/tmp." + file_name, shell=True, stdout=ufile) 
    os.remove(output_dir + "/tmp." + file_name)
    return pybedtools.BedTool(ufile.name)

# break down the calls of type sv_type by region, overlaps
def save_output(master_out_bed, calls, output_dir, file_name, sample_name, sv_type, seg_dups, cent_tel):
    if len(calls) == 0:
        log.write("Zero calls of type " + sv_type + "\n")
        return
    track_name = sample_name + "_" + sv_type
    calls.saveas(output_dir + "/" + file_name + '.bedpe')
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '.bedpe', track_name)

    # get rid of duplicate calls
    calls = merge_duplicate_breaks(calls, 1000)
    print sv_type + "\tNON_DUPLICATE\t" + str(len(calls))
    log.write(sv_type + "\tNON_DUPLICATE\t" + str(len(calls)) + "\n")
    calls.saveas(output_dir + "/" + file_name + '_dedup.bedpe')
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup.bedpe', track_name)            

    # find calls overlapping segmental duplications
    if seg_dups is not None:
        seg_dup_overlap = calls.pair_to_bed(seg_dups, f=1, type="either").cut(xrange(0,22)).saveas(output_dir + "/tmp." + file_name + "_dedup_segdup.bedpe")
        seg_dup_overlap = uniqify(seg_dup_overlap, output_dir, file_name + "_dedup_segdup.bedpe")
        convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_segdup.bedpe', track_name + "_IN_SEG_DUPS")    
        subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_segdup.bedpe.bed', shell=True, stdout=master_out_bed)

    # find calls overlapping peri- centromeric and telomeric regions
    if cent_tel is not None:
        cent_tel_overlap = calls.pair_to_bed(cent_tel, f=1, type="either").cut(xrange(0,22)).saveas(output_dir + "/tmp." + file_name + "_dedup_cent_tel.bedpe")
        cent_tel_overlap = uniqify(cent_tel_overlap, output_dir, file_name + "_dedup_cent_tel.bedpe")
        convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_cent_tel.bedpe', track_name + "_IN_CENT_TEL")    
        subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_cent_tel.bedpe.bed', shell=True, stdout=master_out_bed)

    # subract SD and C/T calls from the stringent set
    if seg_dups is not None and len(seg_dup_overlap) > 0:
        stringent_minus_sd = calls.pair_to_pair(seg_dup_overlap, type="notboth").saveas()
    else:
        stringent_minus_sd = calls
    if cent_tel is not None and len(cent_tel_overlap) > 0:
        stringent_minus_ct = stringent_minus_sd.pair_to_pair(cent_tel_overlap, type="notboth").saveas()
    else:
        stringent_minus_ct = stringent_minus_sd
    
    print sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_ct))
    log.write(sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_ct)) + "\n")

    if len(stringent_minus_ct) == 0:
        return

    # filter out very short (< 1kb) and short (< 5kb) calls into their own category
    very_short_stringent = stringent_minus_ct.filter(bedpe_lt_length_filter, 1000).saveas()
    very_short_stringent.saveas(output_dir + "/" + file_name + "_dedup_stringent_very_short.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_very_short.bedpe', track_name + "_STRINGENT_LT_1KB")
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_very_short.bedpe.bed', shell=True, stdout=master_out_bed)
    if len(very_short_stringent) > 0:        
        stringent_minus_vs = stringent_minus_ct.pair_to_pair(very_short_stringent, type="notboth").saveas()
    else:
        stringent_minus_vs = stringent_minus_ct.saveas()
    print sv_type + "\tTOTAL_STRINGENT_MINUS_VERY_SHORT\t" + str(len(stringent_minus_vs))
    log.write(sv_type + "\tTOTAL_STRINGENT_MINUS_VERY_SHORT\t" + str(len(stringent_minus_vs)) + "\n")
    if len(stringent_minus_vs) == 0:
        return

    short_stringent = stringent_minus_vs.filter(bedpe_lt_length_filter, 5000).saveas()
    short_stringent.saveas(output_dir + "/" + file_name + "_dedup_stringent_short.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_short.bedpe', track_name + "_STRINGENT_LT_5KB")
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_short.bedpe.bed', shell=True, stdout=master_out_bed)
    if len(short_stringent) > 0:
        stringent_minus_vss = stringent_minus_vs.pair_to_pair(short_stringent, type="notboth").saveas()
    else:
        stringent_minus_vss = stringent_minus_vs.saveas()


    # from the remaining calls separate those with high score (breakdancer score = 99) from those with lower scores
    stringent_high_score = stringent_minus_vss.filter(score_gte_filter, 99).saveas(output_dir + "/" + file_name + "_dedup_stringent_high_score.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_high_score.bedpe', track_name + "_STRINGENT_HIGH_SCORE")    
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_high_score.bedpe.bed', shell=True, stdout=master_out_bed)
    
    stringent_low_score = stringent_minus_vss.filter(score_lt_filter, 99).saveas(output_dir + "/" + file_name + "_dedup_stringent_low_score.bedpe")
    convert_bedpe_to_bed12(output_dir + "/" + file_name + '_dedup_stringent_low_score.bedpe', track_name + "_STRINGENT_LOW_SCORE")    
    subprocess.call("cat " + output_dir + "/" + file_name + '_dedup_stringent_low_score.bedpe.bed', shell=True, stdout=master_out_bed)    

    if seg_dups is not None:
        print sv_type + "\tSEGMENTAL_DUPLICATION\t" + str(len(seg_dup_overlap))
    if cent_tel is not None:
        print sv_type + "\tIN_PERI_CENTROMERE_TELOMERE\t" + str(len(cent_tel_overlap))
    print sv_type + "\tSTRINGENT_VERY_SHORT\t" + str(len(very_short_stringent))
    print sv_type + "\tSTRINGENT_SHORT\t" + str(len(short_stringent))
    print sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_vss))
    print sv_type + "\tSTRINGENT_LOW_SCORE\t" + str(len(stringent_low_score))
    print sv_type + "\tSTRINGENT_HIGH_SCORE\t" + str(len(stringent_high_score))

    if seg_dups is not None:
        log.write( sv_type + "\tSEGMENTAL_DUPLICATION\t" + str(len(seg_dup_overlap)) + "\n")
    if cent_tel is not None:
        log.write( sv_type + "\tIN_PERI_CENTROMERE_TELOMERE\t" + str(len(cent_tel_overlap)) + "\n")
    log.write( sv_type + "\tSTRINGENT_VERY_SHORT\t" + str(len(very_short_stringent)) + "\n")
    log.write( sv_type + "\tSTRINGENT_SHORT\t" + str(len(short_stringent)) + "\n")
    log.write( sv_type + "\tTOTAL_STRINGENT\t" + str(len(stringent_minus_vss)) + "\n")
    log.write( sv_type + "\tSTRINGENT_LOW_SCORE\t" + str(len(stringent_low_score)) + "\n")
    log.write( sv_type + "\tSTRINGENT_HIGH_SCORE\t" + str(len(stringent_high_score)) + "\n")


bedpe_calls = pybedtools.BedTool(bedpe_file)

tes = pybedtools.BedTool(te_file)

if seg_dups_file is not None:
    seg_dups = pybedtools.BedTool(seg_dups_file)
else:
    seg_dups = None

if cent_tel_file is not None:
    cent_tel = pybedtools.BedTool(cent_tel_file)
else:
    cent_tel = None


if common_deletions_file is not None:
    common_deletions = pybedtools.BedTool(common_deletions_file)
else:
    common_deletions = None

master_out_bed = open(output_dir + "/" + sample_name + "_svs.bed", 'a')

num_calls = len(bedpe_calls)
print "TOTAL\tALL\t" + str(num_calls)
log.write("TOTAL\tALL\t" + str(num_calls) + "\n")

# divide calls by type:

# translocations: inter-chromosomal
inter_calls = bedpe_calls.filter(inter_chr_filter).saveas()
print "TRANSLOCATIONS\tALL\t" + str(len(inter_calls))
log.write("TRANSLOCATIONS\tALL\t" + str(len(inter_calls)) + "\n")
save_output(master_out_bed, inter_calls, output_dir, "translocations", sample_name, "TRANSLOCATIONS", seg_dups, cent_tel)

# translocations where one end maps to a TE are more like TE insertions in the sample
possible_te_insertions = inter_calls.pair_to_bed(tes, f=.75).saveas()
filtered_possible_te_insertions = possible_te_insertions.filter(bedpe_reciprocal_overlap_ends_filter, 0.75).saveas()
print "TRANSLOCATIONS-POSSIBLE_TE_INSERTIONS\tALL\t" + str(len(filtered_possible_te_insertions))
log.write("TRANSLOCATIONS-POSSIBLE_TE_INSERTIONS\tALL\t" + str(len(filtered_possible_te_insertions)) + "\n")
save_output(master_out_bed, filtered_possible_te_insertions, output_dir, "translocations_possible_te_insertions", sample_name, "TRANSLOCATIONS-POSSIBLE_TE_INSERTIONS", seg_dups, cent_tel)

# now working with intra chromosomal calls
intra_calls = bedpe_calls.filter(intra_chr_filter).saveas()

# expected paired end fragment orientation
expected_orientation = intra_calls.filter(expected_orientation_filter, matches=True).saveas()

# deletions <- right orientation, long fragment size
long_indel_intra_calls = expected_orientation.filter(bedpe_gt_length_filter, insert_size).saveas()
print "DELETIONS\tALL\t" + str(len(long_indel_intra_calls))
log.write("DELETIONS\tALL\t" + str(len(long_indel_intra_calls)) + "\n")
save_output(master_out_bed, long_indel_intra_calls, output_dir, "deletions", sample_name, "DELETIONS", seg_dups, cent_tel)

# special case of deletions - if the deleted area ("ispan") is a TE, more likely a TE insertion in the reference
possible_te_reference_insertions = long_indel_intra_calls.pair_to_bed(tes, type="ispan", f=.75).saveas()
filtered_possible_te_reference_insertions = possible_te_reference_insertions.filter(bedpe_reciprocal_overlap_ispan_filter, 0.75).saveas()
print "POSSIBLE_TE_INSERTIONS_IN_REFERENCE\tALL\t" + str(len(filtered_possible_te_reference_insertions))
log.write("POSSIBLE_TE_INSERTIONS_IN_REFERENCE\tALL\t" + str(len(filtered_possible_te_reference_insertions)) + "\n")
save_output(master_out_bed, filtered_possible_te_reference_insertions, output_dir, "possible_te_reference_insertions", sample_name, "POSSIBLE_TE_INSERTIONS_IN_REFERENCE", seg_dups, cent_tel)

# deletions that match population variants
if common_deletions is not None:
    common_deletions = long_indel_intra_calls.pair_to_bed(common_deletions, type="ispan", f=common_deletion_overlap_pct).saveas()
    filtered_possible_common_deletions = common_deletions.filter(bedpe_reciprocal_overlap_ispan_filter, common_deletion_overlap_pct).saveas()
    print "COMMON_DELETIONS\tALL\t" + str(len(filtered_possible_common_deletions))
    log.write("COMMON_DELETIONS\tALL\t" + str(len(filtered_possible_common_deletions)) + "\n")
    save_output(master_out_bed, filtered_possible_common_deletions, output_dir, "possible_common_deletions", sample_name, "COMMON_DELETIONS", seg_dups, cent_tel)

# insertions are shorter than the fragment size
short_indel_intra_calls = expected_orientation.filter(bedpe_lt_length_filter, insert_size).saveas()
print "INSERTIONS\tALL\t" + str(len(short_indel_intra_calls))
log.write("INSERTIONS\tALL\t" + str(len(short_indel_intra_calls)) + "\n")
save_output(master_out_bed, short_indel_intra_calls, output_dir, "insertions", sample_name, "INSERTIONS", seg_dups, cent_tel)

# inversions are what's left
unexpected_orientation = intra_calls.filter(expected_orientation_filter, matches=False).saveas()
print "INVERSION\tALL\t" + str(len(unexpected_orientation))
log.write("INVERSION\tALL\t" + str(len(unexpected_orientation)) + "\n")
save_output(master_out_bed, unexpected_orientation, output_dir, "inversions", sample_name,  "INVERSIONS", seg_dups, cent_tel)

pybedtools.cleanup()
log.close()
