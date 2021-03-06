import argparse
import random
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')
WEEK2_DIR = str(Path(__file__).resolve().parents[1]) + "/week2"
WEEK3_DIR = str(Path(__file__).resolve().parents[1]) + "/week3"
sys.path.insert(1, WEEK2_DIR)
sys.path.insert(1, WEEK3_DIR)
from hamming_distance import hammingDistance
from profile_most_probable import profileMostProbableKmer
from distanceBetweenPatternAndStrings import distanceBetweenPatternAndStrings
from greedy_motif_search import build_profile_matrix


def find_motifs_from_profile(profile, dnas) :
    k = len(profile["A"])
    motifs = [profileMostProbableKmer(dna, k, profile) for dna in dnas]
    
    return motifs

def find_consensus_string(kmers):
    consensus = []

    for i in range(len(kmers[0])):
        nucleotides_dict = {"A": 0, "C": 0, "G": 0, "T": 0}

        for j in range(len(kmers)):
            curr_kmer = kmers[j]
            nucleotides_dict[curr_kmer[i]] = nucleotides_dict[curr_kmer[i]] + 1

        consensus.append( max(nucleotides_dict, key=nucleotides_dict.get) )

    return "".join(consensus)


def calculate_probability(kmer, profile_matrix) :
    probability = 1

    for i in range(len(kmer)) :
        nucleotide = kmer[i]
        probability = probability * profile_matrix[nucleotide][i]
    
    return probability

def randomizedMotifSearch(dnas, k, t) :
    best_motifs = []
    best_motifs_score = None

    # Randomly select one motif per each dna as a starting point
    for dna in dnas :
        rand_idx = random.randint(0, len(dna)-k)
        random_kmer = dna[rand_idx:rand_idx+k]
        best_motifs.append(random_kmer)
        
    motifs = best_motifs.copy()

    while True :
        profile = build_profile_matrix(motifs, pseudocounts=True)
        motifs = find_motifs_from_profile(profile, dnas)
        
        # Score all kmers in each dna against the profile
        next_motifs = []
        for dna in dnas :
            dna_kmers_scores_dict = {}

            for i in range(len(dna)-k+1) :
                kmer = dna[i:i+k]
                kmer_probability_score = calculate_probability(kmer, profile)
                dna_kmers_scores_dict[kmer] = kmer_probability_score
                
            # After scoring all kmers in one dna strand, get the max
            next_motif = max(dna_kmers_scores_dict, key=dna_kmers_scores_dict.get)


        consensus_motif = find_consensus_string(motifs)
        
        motifs_score = 0
        for motif in motifs : 
            motifs_score = motifs_score + hammingDistance(consensus_motif, motif)
        
        # Calculate score for the first best_motifs because best_score is 0 due to initialization
        if not best_motifs_score :
            consensus_best_motif = find_consensus_string(best_motifs)
            
            best_motifs_score = 0
            for best_motif in best_motifs : 
                best_motifs_score = best_motifs_score + hammingDistance(consensus_best_motif, best_motif)
        
        if motifs_score < best_motifs_score :
            best_motifs = motifs[:]
            best_motifs_score = motifs_score
        else :
            return best_motifs, best_motifs_score


def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    Uses randomness (Monte Carlo Algorithm) to get the best motifs quickly.

                    Input File format :
                    ---------------------------------------
                    8 5
                    CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA
                    GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG
                    TAGTACCGAGACCGAAAGAAGTATACAGGCGT
                    TAGATCAAGTTTCAGGTGCACGTCGGTGAACC
                    AATCCACCAGCTCCACGTGCAATGTTGGCCTA

                    Expected output :
                    ---------------------------------------
                    TCTCGGGG
                    CCAAGGTG
                    TACAGGCG
                    TTCAGGTG
                    TCCACGTG
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
    
    # Read from dataset file
    with open(args.file, 'r') as f :
        k, t = f.readline().split(" ")
        dna_strings = [row.strip() for i, row in enumerate(f.readlines())]
    
    motifs_score = float("inf")
    best_motifs = None
    for i in range(1000) :
        motifs, sc = randomizedMotifSearch(dna_strings, int(k), int(t))

        if sc < motifs_score :
            motifs_score = sc
            best_motifs = motifs
    
    for best_motif in best_motifs :
        print(best_motif)
    

if __name__ == "__main__":
    main()