import argparse
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
sys.path.insert(1, WEEK2_DIR)

from hamming_distance import hammingDistance
from find_neighbors import neighbors


def distanceBetweenPatternAndStrings(pattern, dna) :
    k = len(pattern)
    distance = 0

    for each_dna in dna :
        hamming_distance = float('inf')

        for i in range(len(each_dna)-k+1):
            kmer = each_dna[i:i+k]
            h_dist = hammingDistance(pattern, kmer)

            if  h_dist < hamming_distance :
                hamming_distance = h_dist

        distance = distance + hamming_distance
    
    return distance


def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file will return the Score(Motifs) or sum of distances between Pattern and all strings in Dna.

                    Input File format :
                    ---------------------------------------
                    AAA
                    TTACCTTAAC GATATCTGTC ACGGCGTTCG CCCTAAAGAG CGTCAGAGGT
                    
                    Expected output :
                    ---------------------------------------
                    5
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
        dna = f.readline().split()

    distance = distanceBetweenPatternAndStrings(pattern, dna)

    print(distance)

if __name__ == "__main__":
    main()