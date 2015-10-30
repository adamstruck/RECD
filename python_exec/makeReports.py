#!/usr/bin/env python

import sys
import os
import pysam

def Usage ():
    print("""USAGE:
    makeAlignmentReports.py -i <input_file> -o <output_dir> [-f -b]

""")

def makeAlignmentReports(input_file, file_basename, output_dir):
    
    input_reads = pysam.Samfile(input_file, "rb")
    
    lastAlignedBase   = {}
    innerMateDistance = {}
    
    def recordToIMD(plus_rname, plus_end, minus_start, strand):
        try:
            innerMateDistance[(plus_rname, plus_end, minus_start, strand)] += 1
        except KeyError:
            innerMateDistance[(plus_rname, plus_end, minus_start, strand)] = 1

    def recordToLAB(rname, pos, strand):
        try:
            lastAlignedBase[(rname, pos, strand)] += 1
        except KeyError:
            lastAlignedBase[(rname, pos, strand)] += 1

    counter = 0
    for read in input_reads:

        counter += 1
        if (counter % 1000 == 0):
            print counter

        # check if read is paired and a primary alignment
        if (read.is_paired == True and read.is_secondary == False):
            # prevent doing the computation for both reads in the pair
            if (read.is_read1): 
                pos = input_reads.tell()
                try: 
                    mate = input_reads.mate(read)
                except ValueError:
                    continue
                finally:
                    input_reads.seek(pos)

                # are they concordant?
                if (read.tid == mate.tid):
                    rname = input_reads.getrname(read.tid)
                    # Check read orientation
                    if (read.is_reverse == False and read.mate_is_reverse == True): # FR reads
                        strand = "+"  # read1 strand

                        plus_start = int(read.pos) + 1
                        plus_end = int(read.aend) + 1
                
                        minus_start = int(mate.pos) + 1
                        minus_end = int(mate.aend) + 1

                        # make sure reads face inwards
                        if (plus_start <= minus_start):
                            recordToIMD(rname, plus_start, plus_end, minus_start, minus_end, strand)
                        else: 
                            continue

                    elif (read.is_reverse == True and read.mate_is_reverse == False): # RF reads
                        strand = "-" #read1 strand
                        
                        plus_start = int(mate.pos) + 1
                        plus_end = int(mate.aend) + 1
                        
                        minus_start = int(read.pos) + 1
                        minus_end = int(read.aend) + 1

                        # make sure reads face inwards
                        if (plus_start <= minus_start):
                            recordToIMD(rname, plus_start, plus_end, minus_start, minus_end, strand)
                        else:
                            continue

                else: # discordant read pairs
                    read_rname = input_reads.getrname(read.tid)
                    read_start = int(read.pos) + 1
                    read_end = int(read.aend) + 1

                    mate_rname = input_reads.getrname(mate.tid)
                    mate_start = int(mate.pos) + 1
                    mate_end = int(mate.aend) + 1
                    
                    if (read.is_reverse == False):
                        recordToLAB(read_rname, read_end, "+")
                    else:
                        recordToLAB(read_rname, read_start, "-")

                    if (mate.is_reverse == False):
                        recordToLAB(mate_rname, mate_end, "+")
                    else:
                        recordToLAB(mate_rname, mate_start, "-")

            else: # if read2, skip and continue loop
                continue

        elif (read.is_secondary == False): # primary alignments of unpaired reads
            rname = input_reads.getrname(read.tid)
            pos = int(read.aend) + 1
            
            if (read.is_reverse == False):
                strand = "+"
            else:
                strand = "-"

            recordToLAB(rname, pos, strand)

        else: # discard secondary alignments of unpaired reads
            continue

    input_reads.close()

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
    parser.add_option("-i", "--input-file", dest="input_file", nargs=1, default=None, help="BAM alignment file.")
    parser.add_option("-o", "--output-dir", dest="output_dir", nargs=1, default=".", help="output directory.")
    parser.add_option("-b", "--basename", dest="basename", nargs=1, default=None, help="basename to use for the output file.")
    (options, args) = parser.parse_args()

    if (options.input_file == None):
        Usage()
        print "Error: need to specify -i."
        return

    input_file  = os.path.abspath(os.path.expanduser(options.input_file))
    output_dir  = os.path.abspath(os.path.expanduser(options.output_dir))
    
    if (options.basename == None):
        basename = os.path.basename(input_file)
        basename = basename.split(".")
        file_basename = basename[0]
    else:
        file_basename = options.basename
    
    makeAlignmentReports(input_file, file_basename, output_dir)

if __name__ == "__main__":
    main()
