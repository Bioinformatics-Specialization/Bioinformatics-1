import argparse
import os
import sys
from hamming_distance import hammingDistance
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def neighbors(pattern, d) :
    if d == 0 : 
        return pattern

    if len(pattern) == 1 :
        return ['A', 'C', 'G', 'T']

    neighborhood = []

    suffixNeighbors = neighbors(pattern[1:], d)

    for sNeighbor in suffixNeighbors :
        if hammingDistance(pattern[1:], sNeighbor) < d :
            for nucleotide in ['A', 'C', 'G', 'T'] :
                neighborhood.append(nucleotide + sNeighbor)
        else :
            neighborhood.append(pattern[0] + sNeighbor)
            
    return neighborhood

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return all 'd' mismatches.

                    Input File format :
                    ---------------------------------------
                    ACG
                    1
                    
                    Expected output :
                    ---------------------------------------
                    CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG
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
        pattern = f.readline().strip()
        d = f.readline().strip()

    neighborhood = neighbors(pattern, int(d))

    print(" ".join(neighborhood))


if __name__ == "__main__":
    main()