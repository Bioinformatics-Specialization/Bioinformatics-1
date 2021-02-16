import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
import time


def frequencyTable(text, k) :
    kmer_dict = {}

    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        
        if kmer in kmer_dict :
            kmer_dict[kmer] = kmer_dict[kmer] + 1
        else :
            kmer_dict[kmer] = 1
        
    return kmer_dict

def findClumps(genome, k, L, t) :
    clumped_kmers = list()
    for i in range(len(genome)-L+1):
        freqMap = frequencyTable(genome[i:i+L], k)
        
        [clumped_kmers.append(key) for key, val in freqMap.items() if t <= val]
                
    clumped_kmers = list(set(clumped_kmers))    

    return " ".join(clumped_kmers)

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="find_clumps.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return all kmers that clumps in the genome.

                    Input File format :
                    ---------------------------------------
                    CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
                    5 50 4
                    
                    Expected output :
                    ---------------------------------------
                    CGACA GAAGA
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
        genome = f.readline().strip()
        k, L, t = f.readline().split(" ")

    start = time.time()
    indices = findClumps(genome, int(k), int(L), int(t))
    end = time.time()

    print("Time taken...{}s".format(end-start))

    print(indices)

if __name__ == "__main__":
    main()