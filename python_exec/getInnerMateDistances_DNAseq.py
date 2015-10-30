#!/usr/bin/env python

from __future__ import print_function
import sys
import pydoc
import os
import re
import fileinput
import math
import argparse

def Usage ():
    print("""USAGE:                                                                                                                                                                                                                                                         
getinsertsize.py <SAM file> [options...] 
samtools view <BAM file> | getinsertsize.py - [options...] 

Author: Wei Li
Modified by: Adam Struck
""")

def getStats(dic,maxbound=-1):
  nsum = 0  
  n    = 0
  for (k,v) in dic.items():
    if maxbound!=-1 and k>maxbound:
      continue
    nsum = nsum+k*v
    n    = n+v
  
  meanv = (nsum*1.0)/n
  nsum  = 0
  n     = 0
  for (k,v) in dic.items():
    if maxbound!=-1 and k>maxbound:
      continue
    nsum = nsum+(k-meanv)*(k-meanv)*v
    n    = n+v
  
  varv = math.sqrt(nsum*1.0/(n-1))
  return (meanv,varv)

def getInnerMateDistances(SAMFILE, span_distribution_file, read_distribution_file):
    plrdlen  = {}
    plrdspan = {}

    objmrl = re.compile('([0-9]+)M$')

    nline = 0
    for lines in SAMFILE:
        field = lines.strip().split()
        nline = nline+1
        if nline%1000000 == 0:
            print(str(nline/1000000)+'M...',file=sys.stderr)
        if len(field)<12:
            continue
        try:
            mrl = objmrl.match(field[5])
            if mrl == None: # ignore non-perfect reads                                                                                                                                                              
                continue
            readlen = int(mrl.group(1))
            if readlen in plrdlen.keys():
                plrdlen[readlen] = plrdlen[readlen] + 1
            else:
                plrdlen[readlen] = 1
            if field[6] != '=':
                continue
            dist = int(field[8])
            if dist <= 0: # ignore neg dist                                                                                                                                                                         
                continue
            dist = (int(field[8]) - 2 * readlen)
            if dist < -readlen: # ignore outward facing reads
                continue
            if dist in plrdspan.keys():
                plrdspan[dist] = plrdspan[dist]+1
            else:
                plrdspan[dist] = 1
        except ValueError:
            continue

    readlenval = getStats(plrdlen)
    print('Read length: mean ' + str(readlenval[0]) + ', STD=' + str(readlenval[1]))
            
    if span_distribution_file != "NA":
        print('inner.mate.distance' + '\t' + 'n', file = span_distribution_file)
        for k in sorted(plrdspan.keys()):
            print(str(k) + '\t' + str(plrdspan[k]), file = span_distribution_file)

    if read_distribution_file != "NA":
        print('read.length' + '\t' + 'n', file = read_distribution_file)
        for k in sorted(plrdlen.keys()):
            print(str(k) + '\t' + str(plrdlen[k]), file = read_distribution_file)
    if len(plrdspan) == 0:
        print('No qualified paired-end reads found. Are they single-end reads?')
    else:
        maxv = max(plrdspan, key = plrdspan.get)
        spanval = getStats(plrdspan, maxbound = maxv*3)
        print('Read span: mean ' + str(spanval[0]) + ', STD=' + str(spanval[1]))

def main():
    parser = argparse.ArgumentParser(description='Automatically estimate the insert size of the paired-end reads for a given SAM/BAM file.')
    parser.add_argument('SAMFILE',type=argparse.FileType('r'),help='Input SAM file (use - from standard input)')
    parser.add_argument('--span-distribution-file','-s',type=argparse.FileType('w'),help='Write the distribution of the paired-end read span into a text file with name SPAN_DISTRIBUTION_FILE. This text file is tab-delimited, each line containing two numbers: the span and the number of such paired-end reads.')
    parser.add_argument('--read-distribution-file','-r',type=argparse.FileType('w'),help='Write the distribution of the paired-end read length into a text file with name READ_DISTRIBUTION_FILE. This text file is tab-delimited, each line containing two numbers: the read length and the number of such paired-end reads.')

    args = parser.parse_args()
    
    if args.SAMFILE == None:
        Usage()
        print("Error: Need to input a SAM file")
        return
    
    SAMFILE = args.SAMFILE

    if (args.read_distribution_file != None):
        read_distribution_file =  args.read_distribution_file
    else:
        read_distribution_file = "NA"

    if (args.span_distribution_file != None):
        span_distribution_file =  args.span_distribution_file
    else:
        span_distribution_file = "NA"

    getInnerMateDistances(SAMFILE, span_distribution_file, args.read_distribution_file)

if __name__ == "__main__":
    main()
