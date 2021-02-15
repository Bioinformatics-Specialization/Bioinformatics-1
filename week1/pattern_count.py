import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def patternCount(text, pattern):
    count = 0

    for i in range(len(text)-len(pattern)+1) :
        if text[i : i+len(pattern)] == pattern :
            count = count + 1
    
    return count

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="pattern_count.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file will count the number of occurences of a given pattern within a DNA string.

                    Input File format :
                    ---------------------------------------
                    GCGCG
                    GCG

                    Expected output :
                    ---------------------------------------
                    2
                '''
            )
    
    parser.add_argument('-f', '--file', required=False, help="Input file path.")

    return parser.parse_args()

def main() :
    args = parseArgs()
    text = ""
    pattern = ""

    dataset_path = "{}/{}_dataset.txt".format(DATASET_DIR, os.path.splitext(sys.argv[0])[0])
    
    # Default to the dataset folder, if not provided
    if not args.file :
        args.file = dataset_path
    
    with open(args.file, 'r') as f :
        text = f.readline().strip()
        pattern = f.readline().strip()

    num_matches = patternCount(text, pattern)

    print(num_matches)

if __name__ == '__main__' :
    main()