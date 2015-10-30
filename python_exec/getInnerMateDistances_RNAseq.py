#!/usr/bin/env python

import sys
import os
import numpy

def Usage ():
    print  """USAGE:                                                                                                                                                                                                                                                         
getInnerMateDistances_RNAseq.py -i <input_file> -o <output_dir>                                                                                                                                                                                                                        
input_file       .insert_len file produced by pe_utils.py --compute-insert-len script included in the MISO package (Katz et al. Nat Methods. 2010)  
output_dir       Output directory.
"""

def _compare_keys(x, y):
    try:
        x = int(x)
    except ValueError:
        xint = False
    else:
        xint = True
    try:
        y = int(y)
    except ValueError:
        if xint:
            return -1
        return cmp(x.lower(), y.lower())
        # or cmp(x, y) if you want case sensitivity.
    else:
        if xint:
            return cmp(x, y)
        return 1

def ComputeInsertLengthDist(input_file, input_basename, read_length, output_dir):
    fh = open(input_file, 'r')
    ### Implement a dictionary that will keep track of inner mate range frequencies
    imrs = {}
    print "Building insert length distribution table..."
    for line in fh:
        if line[0] != "#":
            line = line.strip("\n")
            line = line.split("\t")
            values = line[1]
            values = values.split(",")
            for value in values:
                distance = (int(value) - read_length * 2)
                # distance = int(value)
                if distance in imrs:
                    imrs[distance] += 1
                else:
                    imrs[distance] = 1
        else:
            continue
    fh.close()

    if not os.path.isdir(output_dir):
        print "Making directory: %s" % (output_dir)
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, "%s_dist" % (input_basename))
    print "Writing results to %s" % (output_file)
    fh_out = open(output_file, 'w')
    column_names = "%s%s%s" % ("inner_mate_distance", "\t", "n")
    fh_out.write(column_names)
    fh_out.write("\n")
    sorted_insert_lengths = sorted(imrs.keys(), cmp=_compare_keys)
    for ins_len in sorted_insert_lengths:
        output = "%s%s%s" % (ins_len, "\t", imrs[ins_len])
        fh_out.write(output)
        fh_out.write("\n")
    
    fh_out.close()

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", dest="input_file", nargs=1, default=None, help="Input File.")
    parser.add_option("-o", dest="output_dir", nargs=1, default=None, help="Output Directory.")
    parser.add_option("-r", dest="read_length", nargs=1, default=75, type= "int" , help="Read Length.")
    (options, args) = parser.parse_args()

    if options.input_file == None:
        Usage()
        print "Error: need to specify -i."
        return

    if options.output_dir == None:
        Usage()
        print "Error: need to specify -o."
        return

    input_file     = os.path.abspath(os.path.expanduser(options.input_file))
    input_basename = os.path.basename(input_file)
    output_dir     = os.path.abspath(os.path.expanduser(options.output_dir))
    read_length    = options.read_length
    ComputeInsertLengthDist(input_file, input_basename, read_length, output_dir)

if __name__ == "__main__":
    main()
