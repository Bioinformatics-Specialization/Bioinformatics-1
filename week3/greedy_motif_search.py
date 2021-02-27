import argparse
import os
import sys
from profile_most_probable import profileMostProbableKmer
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
sys.path.insert(1, WEEK2_DIR)
from hamming_distance import hammingDistance


def find_consensus_string(kmers):

    consensus = []

    for i in range(len(kmers[0])):
        nucleotides_dict = {"A": 0, "C": 0, "G": 0, "T": 0}

        for j in range(len(kmers)):
            curr_kmer = kmers[j]
            nucleotides_dict[curr_kmer[i]] = nucleotides_dict[curr_kmer[i]] + 1

        consensus.append( max(nucleotides_dict, key=nucleotides_dict.get) )

    return "".join(consensus)

def build_profile_matrix(kmers):
    k = len(kmers[0])
    # Initialize profile matrix
    profileMatrix = { "A": [0] * k, "C": [0] * k, "G": [0] * k, "T": [0] * k }    
    
    for kmer in kmers :
        for i in range(len(kmer)) :
            profileMatrix[kmer[i]][i] = profileMatrix[kmer[i]][i] + 1

    for key, val in profileMatrix.items() :
        profileMatrix[key] = [_/len(kmers) for _ in val]

    return profileMatrix

def greedyMotifSearch(dnas, k, t) :
    # Set best motif as the first kmers for each dna strings
    best_motifs = []
    [best_motifs.append(dna[:k]) for dna in dnas]
    best_motif_score = float('inf')

    count = 0
    for i in range(len(dnas[0])-k+1):
        motifs = []
        consensus_motif = ""
        curr_kmer = dnas[0][i:i+k]
        motifs.append(curr_kmer)
        
        for j in range(1, t) :
            # Get profile matrix
            profile = build_profile_matrix(motifs)

            # Find most probable kmer given the current profile
            motif = profileMostProbableKmer(dnas[j], k, profile)
            motifs.append(motif)
            
        consensus_motif = find_consensus_string(motifs)

        motif_score = 0
        for motif in motifs : 
            motif_score = motif_score + hammingDistance(consensus_motif, motif)
        
        if motif_score <= best_motif_score :
            best_motif_score = motif_score
            best_motifs = motifs

    return "\n".join(best_motifs)

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
        k, t = f.readline().split(" ")
        
        dna_strings = []
        [dna_strings.append(row.strip()) for i, row in enumerate(f.readlines())]
            
    motif = greedyMotifSearch(dna_strings, int(k), int(t))

    print(motif)

if __name__ == "__main__":
    main()