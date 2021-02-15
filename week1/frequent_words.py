import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def frequentWords(text, k) :
    kmer_dict = {}
    words = []

    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        
        if kmer in kmer_dict :
            kmer_dict[kmer] = kmer_dict[kmer] + 1
        else :
            kmer_dict[kmer] = 1

    frequent_count = max(kmer_dict.values())
    [words.append(k) for k, v in kmer_dict.items() if v == frequent_count]
            
    return ' '.join(words)

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
        k = f.readline().strip()

    frequent_kmers = frequentWords(text, int(k))

    print(frequent_kmers)

if __name__ == "__main__":
    main()