#!/usr/bin/env python

import sys
import argparse

# converts the output of breakdancer (the breakdancer results .txt file and the *.bd.bed file generated with the -r option
# that has the actual supporting reads for each call) into a bedpe file.

parser = argparse.ArgumentParser()
parser.add_argument("bd_bed_file", help="Breakdancer bed file of supporting reads")
parser.add_argument("bd_out_file", help="Breakdancer output file")

args = parser.parse_args()
    
bd_bed_file = open(args.bd_bed_file, 'r')
bd_out_file = open(args.bd_out_file, 'r')

name = ''
c1 = ''
s1 = 100000000000000
e1 = 0
c2 = ''
s2 = 100000000000000
e2 = 0

num_svs = 0
num_pairs = 0

unused_fields = "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"

for line in bd_bed_file.readlines():
    if line.startswith("track"):        
        if name != '':
            if s2 < s1:
                ct = c1
                st = s1
                et = e1
                ot = o1
                c1 = c2
                s1 = s2
                e1 = e2
                o1 = o2
                c2 = ct
                s2 = st
                e2 = et
                o2 = ot
            print c1 + "\t" + str(s1) + "\t" + str(e1) + "\t" + c2 + "\t" + str(s2) + "\t" + str(e2) + "\t" + name + "\t" + str(score) + "\t" + o1 + "\t" + o2 + unused_fields
        name = ''
        c1 = ''
        s1 = 100000000000000
        e1 = 0
        c2 = ''
        o1 = ''
        s2 = 100000000000000
        e2 = 0
        o2 = ''
        track_fields = line.split()
        name = track_fields[1].split("=")[1]
        description = line.split("\"")[1]
        description_fields = description.split()
        desc_chrom = description_fields[1]
        desc_start = description_fields[2]
        desc_type = description_fields[3]
        bd_out_line = bd_out_file.readline()
        while bd_out_line.startswith('#'):
            bd_out_line = bd_out_file.readline()
        bd_out_fields = bd_out_line.split()
        bd_out_chrom = bd_out_fields[0]
        bd_out_start = bd_out_fields[1]
        bd_out_type = bd_out_fields[6]
#        print "comparing " + bd_out_chrom + "-" + bd_out_start + "-" + bd_out_type + " to " + desc_chrom + "-" + desc_start + "-" + desc_type
        assert(desc_chrom == bd_out_chrom and desc_start == bd_out_start and desc_type == bd_out_type)
        score = bd_out_fields[8]
        read1 = True
        num_svs = num_svs + 1
        num_pairs = 0
    else:
        fields = line.split()
        if read1:
            c1 = fields[0]
            if long(fields[1]) < s1:
                s1 = long(fields[1])
            if long(fields[2]) > e1:
                e1 = long(fields[2])
            o1 = fields[5]
            num_pairs = num_pairs + 1
        else:
            c2 = fields[0]
            if long(fields[1]) < s2:
                s2 = long(fields[1])
            if long(fields[2]) > e2:
                e2 = long(fields[2])
            o2 = fields[5]
        read1 = not read1

if s2 < s1:
    ct = c1
    st = s1
    et = e1
    ot = o1
    c1 = c2
    s1 = s2
    e1 = e2
    o1 = o2
    c2 = ct
    s2 = st
    e2 = et
    o2 = ot

print c1 + "\t" + str(s1) + "\t" + str(e1) + "\t" + c2 + "\t" + str(s2) + "\t" + str(e2) + "\t" + name + "\t" + str(num_pairs) + "\t" + o1 + "\t" + o2 + unused_fields
