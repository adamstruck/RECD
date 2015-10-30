#!/usr/bin/env python

import sys
import os

def Usage ():
    print  """USAGE:                                                                                                                                                                                                                                                         
getIMDs_RNA.py --gff-file <file.gff> --sample_file <file.inner_mate_range> --output-dir <path> [--basename]

Author: Adam Struck
"""

def getIMDstats(gff_file, sample_file, basename, output_dir):
    if not os.path.isdir(output_dir):
        print "Making directory: %s" % (output_dir)
        os.makedirs(output_dir)
    os.system("sed '1d' %s | perl -lane 'print unless($F[1] > $F[2])' > %s.tmp" % (sample_file, os.path.join(output_dir, basename)))
    os.system("bedtools intersect -f 1.0 -a %s -b %s > %s.imr_overlaps" % (os.path.join(output_dir, basename + ".tmp"), gff_file, os.path.join(output_dir, basename)))
    print "Found first set of intersections..."
    os.system("sed '1d' %s | perl -lane 'print unless($F[1] <= $F[2])' | awk ' { t = $2; $2 = $3; $3 = t; print; } ' | sed 's/\s/\t/g' > %s.tmp" % (sample_file, os.path.join(output_dir, basename)))
    os.system("bedtools intersect -f 1.0 -a %s -b %s | awk ' { t = $2; $2 = $3; $3 = t; print; } ' | sed 's/\s/\t/g' >> %s.imr_overlaps" % (os.path.join(output_dir, basename + ".tmp"), gff_file, os.path.join(output_dir, basename)))
    print "Found second set of intersections..."
    os.system("rm %s" % (os.path.join(output_dir, basename + ".tmp")))

    fh = open(os.path.join(output_dir, basename + ".imr_overlaps"), "r")

    print "Calulating inner mate distances..."

    imds = {}
    nline = 0
    for line in fh:
        field = line.strip("\n")
        field = line.split("\t")
        nline = nline+1
        if nline % 10000 == 0:
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
    output_file = os.path.join(output_dir, "%s.imd_dist" % (basename))
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

    os.system("rm %s" % (os.path.join(output_dir, basename + ".imr_overlaps")))


def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-g", "--gff-file", dest="gff_file", nargs=1, default=None, help="Constitutive Exons GFF File.")
    parser.add_option("-s", "--sample-file", dest="sample_file", nargs=1, default=None, help=".inner_mate_range Sample File From JUNKIE.")
    parser.add_option("-o", "--output-dir", dest="output_dir", nargs=1, default=None, help="Output Directory.")
    parser.add_option("-b", "--basename", dest="basename", nargs=1, default=None, help="Basename for the Output File.")
    (options, args) = parser.parse_args()

    if options.gff_file == None:
        Usage()
        print "Error: need to specify --gff-file."
        return

    if options.sample_file == None:
        Usage()
        print "Error: need to specify --sample-file."
        return

    if options.output_dir == None:
        Usage()
        print "Error: need to specify --output-dir."
        return

    gff_file    = os.path.abspath(os.path.expanduser(options.gff_file))
    sample_file = os.path.abspath(os.path.expanduser(options.sample_file))
    output_dir  = os.path.abspath(os.path.expanduser(options.output_dir))

    if options.basename == None:
        basename = os.path.basename(sample_file)
        basename = basename.split(".")
        basename = basename[0]
    else:
        basename = options.basename


    getIMDstats(gff_file, sample_file, basename, output_dir)

if __name__ == "__main__":
    main()
