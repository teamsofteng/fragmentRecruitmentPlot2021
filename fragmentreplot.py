#!/usr/bin/env/python3
import argparse
import re
import os
import matplotlib.pyplot as fragment
# This script takes a sam file input and determines the genome position and percent identity of mapped reads to produce a scatter plot.
def main():
    # Input
    parser = argparse.ArgumentParser(description="Analysis")
    parser.add_argument('--source=','--source',type=str, required=True, dest="input")
    args = parser.parse_args()
    # Holding genome positions
    x_values = []
    # Holding percent identity scores
    y_values = []
    totalmatches = 0
    toTAL = 0
    name = ""
    # Parse out sample name for plot
    sample = args.input.split("_")[0]
    sp2 = ""
    if args.input.endswith(".sam"):
        for line in open(args.input,'r'):
            if not line.startswith("@"):
                # specifically looking at mapped reads through the MD:Z 
                if "MD:Z" in line:
                    breakdown = line.split("\t")
                    # Specifically grabbing the START position for the genome position
                    pos = breakdown[3]
                    x_values.append(int(pos))
                    name = breakdown[2]
                    cigar = breakdown[5]
                    # Searching for and collecting all MD:Z lines
                    MD = re.findall(r"MD:Z:[0-9A-Z]*",line)
                    # Determining matches through cigar
                    matches = re.findall(r"[0-9]*M",cigar)
                    finalmatch = [i.split("M")[0] for i in matches]
                    for e in range(0,len(finalmatch)):
                        totalmatches += int(finalmatch[int(e)])
                        toTAL = totalmatches
                    totalmatches = 0
                    # Splitting MDZ line to determine the number of mismatches
                    for m in MD:
                        spt = m.split(":")
                        sp2 = spt[2]
                    mismatches = len(re.findall('[ATCG\^]',sp2))
                    # Calculating percent identity
                    identity = 100 * (toTAL/(toTAL + mismatches))
                    y_values.append(identity)
    # Creating the scatter plot using matplotlib
    fragment.scatter(x_values,y_values,label=sample,color='g',s= 1)
    fragment.xlabel('Genome Position')
    fragment.ylabel('Percent Identity')
    #fragment.xticks(np.arange(1,max(x_values)))
    fragment.title(name)
    fragment.legend()
    fragment.show()

if __name__ == '__main__':
    main()


