import argparse
import random
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
WEEK3_DIR = str(Path(__file__).resolve().parents[1]) + "/week3"
sys.path.insert(1, WEEK2_DIR)
sys.path.insert(1, WEEK3_DIR)



def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    Uses Gibbs Sampler, which is similar to randomizedMotifSearch, but uses
                    slightly different RNG in order to retrieve the best motifs.

                    Input File format :
                    ---------------------------------------
                    8 5 100
                    CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA
                    GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
                    TAGTACCGAGACCGAAAGAAGTATACAGGCGT
                    TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
                    AATCCACCAGCTCCACGTGCAATGTTGGCCTA

                    Expected output :
                    ---------------------------------------
                    TCTCGGGG
                    CCAAGGTG
                    TACAGGCG
                    TTCAGGTG
                    TCCACGTG
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
        k, t, N = f.readline().split(" ")
        
        dna_strings = []
        [dna_strings.append(row.strip()) for i, row in enumerate(f.readlines())]
    
    

if __name__ == "__main__":
    main()