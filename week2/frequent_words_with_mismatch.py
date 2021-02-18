import argparse
import os
import sys
from pathlib import Path
WEEK1_DIR = str(Path(__file__).resolve().parents[1]) + "/week1"
sys.path.insert(1, WEEK1_DIR)

from find_neighbors import neighbors
from reverse_complement import reverseComplement
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def frequentWordsWithMismatchesRC(text, k, d) :
    patterns = []
    freqMap = {}

    for i in range(len(text)-k+1):
        kmer = reverseComplement(text[i:i+k])
        neighborhood = neighbors(kmer, d)
        
        # Append all the reverse complement of the neighbors
        rc_neighborhood = neighborhood + [reverseComplement(neighbor) for neighbor in neighborhood]
        
        for neighbor in rc_neighborhood :
            if neighbor in freqMap :
                freqMap[neighbor] = freqMap[neighbor] + 1
            else :
                freqMap[neighbor] = 1
    
    m = max(freqMap.values())

    for key, val in freqMap.items() :
        if val == m :
            patterns.append(key)
    
    return " ".join(patterns)

def frequentWordsWithMismatches(text, k, d) :
    patterns = []
    freqMap = {}

    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        neighborhood = neighbors(kmer, d)
        for neighbor in neighborhood :
            if neighbor in freqMap :
                freqMap[neighbor] = freqMap[neighbor] + 1
            else :
                freqMap[neighbor] = 1
    
    m = max(freqMap.values())

    for key, val in freqMap.items() :
        if val == m :
            patterns.append(key)
    
    return " ".join(patterns)

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return most frequent kmers with 'd' mismatches in text.

                    Input File format :
                    ---------------------------------------
                    ACGTTGCATGTCGCATGATGCATGAGAGCT
                    4 1
                    
                    Expected output :
                    ---------------------------------------
                    GATG ATGC ATGT
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
        text = f.readline().strip()
        k, d = f.readline().rstrip().split(" ")

    frequent_kmers = frequentWordsWithMismatches(text, int(k), int(d))

    print(frequent_kmers)

    ''' Uncomment this for getting patterns with mismatch AND reverse complement '''
    #frequent_rc_kmers = frequentWordsWithMismatchesRC(text, int(k), int(d))
    #print(frequent_rc_kmers)

if __name__ == "__main__":
    main()