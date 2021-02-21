import argparse
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
sys.path.insert(1, WEEK2_DIR)

from find_neighbors import neighbors


def medianString(dna, k) :

    return ""


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
        k = f.readline().split()
        dna = [_.replace("\n", "") for _ in f.readlines()]

    medianString(dna, k)

    print(dna)

if __name__ == "__main__":
    main()