import argparse
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
sys.path.insert(1, WEEK2_DIR)

from find_neighbors import neighbors
from distanceBetweenPatternAndStrings import distanceBetweenPatternAndStrings


def medianString(dna, k) :
    medianStr = ("", float("inf"))
    curr_kmer = "".join(["A"] * k)

    # Get all neighbors of 'curr_kmer' into a dictionary
    curr_kmer_neighbors = neighbors(curr_kmer, k)
    
    # Find distance between each neighbor against the DNA set, and get the minimum
    for kmer_neighbor in curr_kmer_neighbors :
        dist = distanceBetweenPatternAndStrings(kmer_neighbor, dna)

        if dist < medianStr[1] :
            medianStr = (kmer_neighbor, dist)
        
    return medianStr[0]


def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file will return the median string of the collection of DNAs.

                    Input File format :
                    ---------------------------------------
                    3
                    AAATTGACGCAT
                    GACGACCACGTT
                    CGTCAGCGCCTG
                    GCTGAGCACCGG
                    AGTTCGGGACAG
                    
                    Expected output :
                    ---------------------------------------
                    GAC
                '''
            )
    
    parser.add_argument('-f', '--file', required=False, help="Input file path.")

    return parser.parse_args()

def main() :
    args = parseArgs()
    dataset_path = "{}/{}_dataset.txt".format(DATASET_DIR, os.path.splitext(sys.argv[0])[0])
    
    # Default to the dataset folder, if not provided
    if not args.file :
        args.file = dataset_path
    
    with open(args.file, 'r') as f :
        k = f.readline().split()[0]
        dna = [_.replace("\n", "") for _ in f.readlines()]

    answer = medianString(dna, int(k))

    print(answer)

if __name__ == "__main__":
    main()