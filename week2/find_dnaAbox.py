import argparse
import os
import sys
from pathlib import Path
WEEK1_DIR = str(Path(__file__).resolve().parents[1]) + "/week1"
sys.path.insert(1, WEEK1_DIR)

from minimum_skew import minimumSkew
from frequent_words_with_mismatch import frequentWordsWithMismatchesRC

from find_neighbors import neighbors
from reverse_complement import reverseComplement
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')



def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file will return the predicted dnaA box of a given genome.

                    Input File format :
                    ---------------------------------------
                    ACGTTGCATGTCGCATGATGCATGAGAGCT
                    
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
        genome_header = f.readline().strip()
        genome = f.readlines()

    # Remove newline and concatenate all lines into one single str
    genome = ''.join([line.replace('\n','') for line in genome])
    
    # Find where the minimum skew is : This will give where the oriC is.
    minimum_skew = minimumSkew(genome)
    
    WINDOW = 500
    k = 9
    d = 2
    location_ind1, location_ind2 = minimum_skew.split(" ")

    # Find frequent neighbors and it's reverse complement neighbors with mismatch
    dnaA_boxes = frequentWordsWithMismatchesRC(genome[int(location_ind1):int(location_ind1) + WINDOW], k, d)
    
    print(dnaA_boxes)


if __name__ == "__main__":
    main()