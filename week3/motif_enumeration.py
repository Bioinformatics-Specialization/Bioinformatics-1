import argparse
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
sys.path.insert(1, WEEK2_DIR)

from find_neighbors import neighbors


def motifEnumeration(dna, k, d) :
    pattern_dict = {}
    final_pattern = []

    for index, each_dna in enumerate(dna) :
        patterns = []
        for i in range(len(each_dna)-k+1) :
            kmer = each_dna[i:i+k]
            kmer_neighbors = neighbors(kmer, d)
            patterns = patterns + kmer_neighbors
        
        patterns = list(set(patterns))
        pattern_dict[index] = patterns

    for curr_kmer in pattern_dict[0] :
        flag = True
        for i in range(len(dna)) :
            if curr_kmer not in pattern_dict[i] :
                flag = False
                break
        
        if flag :
            final_pattern.append(curr_kmer)
    
    return " ".join(final_pattern)

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file will return the list of motifs from a given DNA strings.
                    Motifs can be mutated, but they all have to be seen in each DNA strings.

                    Input File format :
                    ---------------------------------------
                    3 1
                    ATTTGGC
                    TGCCTTA
                    CGGTATC
                    GAAAATT
                    
                    Expected output :
                    ---------------------------------------
                    ATA ATT GTT TTT
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
        k, d = f.readline().split()
        dna = [_.replace("\n", "") for _ in f.readlines()]

    motifs = motifEnumeration(dna, int(k), int(d))

    print(motifs)

if __name__ == "__main__":
    main()