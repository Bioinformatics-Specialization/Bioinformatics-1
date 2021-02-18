import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def hammingDistance(text_a, text_b) :
    dist = 0

    for i in range(len(text_a)):
        if text_a[i] != text_b[i]:
            dist = dist + 1

    return dist

def mismatchPatternCount(dna, pattern, d) :
    counter = 0

    for i in range(len(dna)-len(pattern)+1) :
        kmer = dna[i:i+len(pattern)]

        if hammingDistance(pattern, kmer) <= d :
            counter = counter + 1
    
    return counter


def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file will count the number of occurences of a given pattern with at most
                    'd' mismatches within a DNA string.

                    Input File format :
                    ---------------------------------------
                    GAGG
                    TTTAGAGCCTTCAGAGG
                    2

                    Expected output :
                    ---------------------------------------
                    4
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
        dna = f.readline().strip()
        d = f.readline().strip()

    num_matches = mismatchPatternCount(dna, pattern, int(d))

    print(num_matches)

if __name__ == '__main__' :
    main()