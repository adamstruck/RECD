#!/usr/bin/env python

import sys
import os

def Usage ():
    print  """USAGE:                                                                                                                                                                                                                                                         
getIMDstatsFromMISO.py --gff-file <file.gff> --sample <file.bam> --output-dir <path>
gff_file        GFF file
sample          BAM file for sample
output_dir      Output directory path
"""

def getIMDstats(gff_file, gff_basename, sample, sample_basename, read_length, output_dir):
    if not os.path.isdir(output_dir):
        print "Making directory: %s" % (output_dir)
        os.makedirs(output_dir)
    if not os.path.isfile(os.path.join(output_dir, gff_basename + ".min_1000.const_exons.gff")):
        os.system("exon_utils.py --get-const-exons %s  --min-exon-size 1000 --output-dir %s" % (gff_file, output_dir))
        os.system("pe_utils.py --compute-insert-len %s %s.min_1000.const_exons.gff --sd-max 5 --output-dir %s" % (sample, os.path.join(output_dir, gff_basename), output_dir))
        os.system("python /gne/home/strucka/lib/R/RepeatAnalyzer/exec/getInnerMateDistances_RNAseq.py -i %s.insert_len -r %s -o %s" % (os.path.join(output_dir, sample), read_length, output_dir))
    else:
        os.system("pe_utils.py --compute-insert-len %s %s.min_1000.const_exons.gff --sd-max 5 --output-dir %s" % (sample, os.path.join(output_dir, gff_basename), output_dir))
        os.system("python /gne/home/strucka/lib/R/RepeatAnalyzer/exec/getInnerMateDistances_RNAseq.py -i %s.insert_len -r %s -o %s" % (os.path.join(output_dir, sample_basename), read_length, output_dir))

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--gff-file", dest="gff_file", nargs=1, default=None, help="GFF File")
    parser.add_option("--sample", dest="sample", nargs=1, default=None, help="Sample BAM File")
    parser.add_option("--output-dir", dest="output_dir", nargs=1, default=None, help="Output Directory")
    parser.add_option("--read-length", dest="read_length", nargs=1, default=75, type= "int" , help="Read Length.")
    (options, args) = parser.parse_args()

    if options.gff_file == None:
        Usage()
        print "Error: need to specify --gff-file."
        return

    if options.sample == None:
        Usage()
        print "Error: need to specify --sample."
        return

    if options.output_dir == None:
        Usage()
        print "Error: need to specify --output-dir."
        return

    read_length     = options.read_length
    gff_file        = os.path.abspath(os.path.expanduser(options.gff_file))
    gff_split       = os.path.splitext(os.path.basename(gff_file))
    gff_basename    = gff_split[0]
    sample          = os.path.abspath(os.path.expanduser(options.sample))
    sample_basename = os.path.basename(sample)
    output_dir      = os.path.abspath(os.path.expanduser(options.output_dir))

    getIMDstats(gff_file, gff_basename, sample, sample_basename, read_length, output_dir)

if __name__ == "__main__":
    main()
