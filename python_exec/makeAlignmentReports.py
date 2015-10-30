#!/usr/bin/env python

import sys
import os
import pysam

def Usage ():
    print("""USAGE:
    makeAlignmentReports.py -i <input_file> -o <output_dir> [-f -b]

Author: Adam Struck
""")

def makeAlignmentReports(input_file, file_format, file_basename, output_dir):
    if (file_format == "BAM"):
        samfile = pysam.Samfile(input_file, "rb")
    else:
        samfile = pysam.Samfile(input_file, "r")

    lastAlignedBase   = {}
    innerMateDistance = {}

    def recordToIMD(plus_rname, plus_end, minus_start, strand):
        if not (plus_rname, plus_end, minus_start, strand) in innerMateDistance:
            innerMateDistance[(plus_rname, plus_end, minus_start, strand)] = 1
        else:
            innerMateDistance[(plus_rname, plus_end, minus_start, strand)] += 1

    def recordToLAB(rname, pos, strand):
        if not (rname, pos, strand) in lastAlignedBase:
            lastAlignedBase[(rname, pos, strand)] = 1
        else:
            lastAlignedBase[(rname, pos, strand)] += 1

    for read in samfile.fetch( until_eof = True ):
        # look at the primary alignment of properly paired reads
        # properly paired definition depends on aligner used - this script infers this from the FLAG field. 
        if (read.is_proper_pair == True and read.is_read1 == True and read.is_secondary == False and read.tid == read.rnext):
            rname = samfile.getrname(read.tid)
            if (read.is_reverse == False and read.mate_is_reverse == True): 
                strand = "+"  # strand that read1 aligns to

                plus_start = int(read.pos) + 1
                plus_end = int(read.aend) + 1
                
                minus_start = int(read.pnext) + 1
                ## minus_end = ?

                if (plus_start <= minus_start):
                    recordToIMD(rname, plus_end, minus_start, strand)
        
            elif (read.is_reverse == True and read.mate_is_reverse == False):
                strand = "-"

                plus_start = int(read.pnext) + 1
                plus_end = "?"

                minus_start = int(read.pos) + 1
                minus_end = int(read.aend) + 1

                if (plus_start <= minus_start):
                    recordToIMD(rname, plus_end, minus_start, strand)                    

        # look at the primary alignment of unpaired reads
        elif (read.is_proper_pair == False and read.is_secondary == False and int(read.mapq) > 0):
            rname = samfile.getrname(read.tid)
            pos = int(read.aend) + 1
            
            if (read.is_reverse == False):
                strand = "+"
            else:
                strand = "-"

            recordToLAB(rname, pos, strand)

    samfile.close()
    
    print "Writing output to file..."

    if not os.path.isdir(output_dir):
        print "Making directory: %s" % (output_dir)
        os.makedirs(output_dir)
    
    output_file1 = os.path.join(output_dir, "%s.inner_mate_distances" % (file_basename))
    output_file2 = os.path.join(output_dir, "%s.last_aligned_base" % (file_basename))
    print "Writing results to %s & %s" % (output_file1, output_file2)
    
    fh_imd = open(output_file1, 'w')
    column_names = "%s%s%s%s%s%s%s%s%s" % ("chr", "\t", "plus", "\t", "minus", "\t", "strand", "\t", "n")
    fh_imd.write(column_names)
    fh_imd.write("\n")
    for key in sorted(innerMateDistance.keys()):
        output = "%s%s%s%s%s%s%s%s%s" % (key[0], "\t", key[1], "\t", key[2], "\t", key[3], "\t", innerMateDistance[key])
        fh_imd.write(output)
        fh_imd.write("\n")
    fh_imd.close()

    fh_lab = open(output_file2, 'w')
    column_names = "%s%s%s%s%s%s%s" % ("chr", "\t", "pos", "\t", "strand", "\t", "n")
    fh_lab.write(column_names)
    fh_lab.write("\n")
    for key in sorted(lastAlignedBase.keys()):
        output = "%s%s%s%s%s%s%s" % (key[0], "\t", key[1], "\t", key[2], "\t", lastAlignedBase[key])
        fh_lab.write(output)
        fh_lab.write("\n")
    fh_lab.close()

    print "Done"

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--input-file", dest="input_file", nargs=1, default=None, help="alignemnt file.")
    parser.add_option("-f", "--format", dest="file_format", nargs=1, default="BAM", help="is the alignment file in SAM or BAM format?")
    parser.add_option("-o", "--output-dir", dest="output_dir", nargs=1, default=".", help="output directory.")
    parser.add_option("-b", "--basename", dest="basename", nargs=1, default=None, help="basename to use for the output file.")
    (options, args) = parser.parse_args()

    if (options.input_file == None):
        Usage()
        print "Error: need to specify -i."
        return

    if (options.file_format != "BAM" and options.file_format != "SAM"):
        Usage()
        print "Error: only alignment files in SAM or BAM format are supported."
        return

    input_file  = os.path.abspath(os.path.expanduser(options.input_file))
    file_format = options.file_format
    output_dir  = os.path.abspath(os.path.expanduser(options.output_dir))
    
    if (options.basename == None):
        basename = os.path.basename(input_file)
        basename = basename.split(".")
        file_basename = basename[0]
    else:
        file_basename = options.basename
    
    makeAlignmentReports(input_file, file_format, file_basename, output_dir)

if __name__ == "__main__":
    main()
