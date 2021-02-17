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

def approximatePatternMatching(pattern, text, d):
    indices = []

    for i in range(len(text)-len(pattern)+1):
        kmer = text[i:i+len(pattern)]

        if hammingDistance(pattern, kmer) <= d :
            indices.append(str(i))

    return " ".join(indices)

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return number of mismatches between two given strings.

                    Input File format :
                    ---------------------------------------
                    ATTCTGGA
                    CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT
                    3
                    
                    Expected output :
                    ---------------------------------------
                    6 7 26 27
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
        text = f.readline().strip()
        d = f.readline().strip()

    indices = approximatePatternMatching(pattern, text, int(d))

    print(indices)


if __name__ == "__main__":
    main()