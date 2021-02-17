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

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="hamming_distance.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return number of mismatches between two given strings.

                    Input File format :
                    ---------------------------------------
                    GGGCCGTTGGT
                    GGACCGTTGAC
                    
                    Expected output :
                    ---------------------------------------
                    3
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
        text_a = f.readline().strip()
        text_b = f.readline().strip()

    distance = hammingDistance(text_a, text_b)

    print(distance)


if __name__ == "__main__":
    main()