import argparse
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
sys.path.insert(1, WEEK2_DIR)




def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    Returns a collection of strings BestMotifs resulting from applying 
                    GreedyMotifSearch(Dna, k, t).

                    Input File format :
                    ---------------------------------------
                    3 5
                    GGCGTTCAGGCA
                    AAGAATCAGTCA
                    CAAGGAGTTCGC
                    CACGTCAATCAC
                    CAATAATATTCG

                    Expected output :
                    ---------------------------------------
                    CAG
                    CAG
                    CAA
                    CAA
                    CAA
                '''
            )
    
    parser.add_argument('-f', '--file', required=False, help="Input file path.")

    return parser.parse_args()

def main() :
    args = parseArgs()
    dataset_path = "{}/{}_dataset.txt".format(DATASET_DIR, os.path.splitext(sys.argv[0])[0])
    nucleotides = ["A", "C", "G", "T"]

    # Default to the dataset folder, if not provided
    if not args.file :
        args.file = dataset_path
    
    with open(args.file, 'r') as f :
        text = f.readline().strip()
        k = f.readline().strip()

        profile_matrix = {}
        for i, row in enumerate(f.readlines()) : 
            profile_matrix[nucleotides[i]] = [float(_) for _ in row.strip().split(" ")]

    kmer = profileMostProbableKmer(text, int(k), profile_matrix)

    print(kmer)

if __name__ == "__main__":
    main()