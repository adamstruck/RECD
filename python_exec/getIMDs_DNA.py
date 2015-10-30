#!/usr/bin/env python

import sys
import os

def Usage ():
    print("""USAGE:
    getIMDs_DNA.py -s <sample_file> -o <output_dir> [-b]

Author: Adam Struck
""")

def getInnerMateDistances(sample_file, sample_basename, output_dir):
    fh = open(sample_file, "r")
    next(fh) # skip the header line
 
    imds = {}
    nline = 0
    for line in fh:
        field = line.strip("\n")
        field = line.split("\t")
        nline = nline+1
        if nline % 1000000 == 0:
            print nline
        plus  = int(field[1])
        minus = int(field[2])
        n     = int(field[4])
        dist   = minus - plus + 1
        if dist in imds:
            imds[dist] = imds[dist] + n
        else:
            imds[dist] = n

    fh.close()
    
    print "Writing output to file..."

    if not os.path.isdir(output_dir):
        print "Making directory: %s" % (output_dir)
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, "%s.imd_dist" % (sample_basename))
    print "Writing results to %s" % (output_file)

    fh_out = open(output_file, 'w')
    column_names = "%s%s%s" % ("inner.mate.distance", "\t", "n")
    fh_out.write(column_names)
    fh_out.write("\n")
    for dist in sorted(imds.keys()):
        output = "%s%s%s" % (dist, "\t", imds[dist])
        fh_out.write(output)
        fh_out.write("\n")

    fh_out.close()

    print "Done"

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-s", "--sample-file", dest="sample_file", nargs=1, default=None, help=".inner_mate_range Sample File From JUNKIE")
    parser.add_option("-o", "--output-dir", dest="output_dir", nargs=1, default=None, help="Output Directory.")
    parser.add_option("-b", "--basename", dest="basename", nargs=1, default=None, help="Basename for the Output File.")
    (options, args) = parser.parse_args()

    if options.sample_file == None:
        Usage()
        print "Error: need to specify -i."
        return

    if options.output_dir == None:
        Usage()
        print "Error: need to specify -o."
        return

    sample_file = os.path.abspath(os.path.expanduser(options.sample_file))
    output_dir  = os.path.abspath(os.path.expanduser(options.output_dir))
    
    if options.basename == None:
        basename = os.path.basename(sample_file)
        basename = basename.split(".")
        sample_basename = basename[0]
    else:
        sample_basename = options.basename
    
    getInnerMateDistances(sample_file, sample_basename, output_dir)

if __name__ == "__main__":
    main()
