import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def reverseComplement(text):
    rc_text = []
    
    for letter in text[::-1]:
        if letter.upper() == "A" :
            rc_text.append("T")
        elif letter.upper() == "C" :
            rc_text.append("G")
        elif letter.upper() == "G" :
            rc_text.append("C")
        else :
            rc_text.append("A")
    
    return "".join(rc_text)

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="pattern_count.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return the most frequently count k-mer from a DNA string.

                    Input File format :
                    ---------------------------------------
                    ACGTTGCATGTCGCATGATGCATGAGAGCT
                    4

                    Expected output :
                    ---------------------------------------
                    CATG GCAT
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

    rc_pattern = reverseComplement(text)

    print(rc_pattern)

if __name__ == "__main__":
    main()