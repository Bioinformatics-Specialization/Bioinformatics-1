import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def patternMatching(pattern, genome):

    occurences = list()
    [occurences.append(str(i)) for i in range(len(genome)-len(pattern)+1) if pattern == genome[i:i+len(pattern)]]
            
    return ' '.join(occurences)


def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="reverse_complement.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return all the indices that the pattern is found in the genome.

                    Input File format :
                    ---------------------------------------
                    ATAT
                    GATATATGCATATACTT
                    
                    Expected output :
                    ---------------------------------------
                    1 3 9
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
        pattern = f.readline().strip()
        genome = f.readline().strip()


    indices = patternMatching(pattern, genome)

    print(indices)

if __name__ == "__main__":
    main()