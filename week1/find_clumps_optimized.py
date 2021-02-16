import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
import time



def findClumps(genome, k, L, t) :
    clumped_kmers = []
    kmers_occurence_dict = {}
    
    # Populate kmers-indices dictionary by iterating the genome once.
    for i in range(len(genome)-k+1):
        kmer = genome[i:i+k]
        if kmer in kmers_occurence_dict :
            kmers_occurence_dict[kmer].append(i)
        else :
            kmers_occurence_dict[kmer] = [i]

    
    for key, indices in kmers_occurence_dict.items():
        # Base Case :
        if len(indices) < t :
            continue
        
        # Figure out if indices fall within the window and has at least 't' occurences.
        for i in range(len(indices)) :
            curr_idx = indices[i]
            idx_max = curr_idx + L - k            
            kmer_counter = 0

            for j in range(i,len(indices)):
                idx = indices[j]
                
                if (curr_idx <= idx) and (idx <= idx_max) :
                    kmer_counter = kmer_counter + 1

            if t <= kmer_counter :
                clumped_kmers.append(key)
                break
    
    return clumped_kmers

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="find_clumps_optimized.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return number of kmers (that clumps) from the genome.

                    Input File format :
                    ---------------------------------------
                    CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA
                    5 50 4
                    
                    Expected output :
                    ---------------------------------------
                    2
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
    clumps = findClumps(genome, int(k), int(L), int(t))
    end = time.time()

    print("Time taken...{}s".format(end-start))

    print(len(clumps))

if __name__ == "__main__":
    main()