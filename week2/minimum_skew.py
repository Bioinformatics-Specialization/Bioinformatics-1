import argparse
import os
import sys
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


def minimumSkew(genome):
    counter = 0
    minimum_val = 0
    skew_list = [0]

    for nucleotide in genome :
        if nucleotide == "A" :
            counter = counter + 0
        elif nucleotide == "T" :
            counter = counter + 0
        elif nucleotide == "G" :
            counter = counter + 1
        else :
            counter = counter - 1
        
        skew_list.append(counter)
    
    minimum_val = min(skew_list)

    indices = [str(i) for i, v in enumerate(skew_list) if v == minimum_val]

    return " ".join(indices)

def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="find_clumps_optimized.py",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    This file return GC skew (#G - #C) of the genome, and then
                    locate all the indices with the lowest value.

                    Input File format :
                    ---------------------------------------
                    CATGGGCATCGGCCATACGCC
                    
                    Expected output :
                    ---------------------------------------
                    0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2
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
        genome = f.readline().strip()

    minimum_indices = minimumSkew(genome)

    print(minimum_indices)


if __name__ == "__main__":
    main()