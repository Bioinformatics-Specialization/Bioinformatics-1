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

def build_profile_matrix(kmers, pseudocounts=False):
    k = len(kmers[0])
    # Initialize profile matrix
    if pseudocounts : profileMatrix = { "A": [1] * k, "C": [1] * k, "G": [1] * k, "T": [1] * k }    
    else : profileMatrix = { "A": [0] * k, "C": [0] * k, "G": [0] * k, "T": [0] * k }    
    
    for kmer in kmers :
        for i in range(len(kmer)) :
            profileMatrix[kmer[i]][i] = profileMatrix[kmer[i]][i] + 1

    for key, val in profileMatrix.items() :
        profileMatrix[key] = [_/len(kmers) for _ in val]

    return profileMatrix


def find_motifs_from_profile2(profile, dnas) :
    motifs = []
    k = len(profile["A"])

    counter = 1
    for dna in dnas :
        # print("Looking into {}th DNA : {}".format(counter, dna))
        
        probable_kmer = profileMostProbableKmer(dna, k, profile)

        motifs.append(probable_kmer)
        counter = counter + 1
    
    return motifs



def find_motifs_from_profile(profile, dnas) :
    motifs = []
    k = len(profile["A"])

    counter = 1
    for dna in dnas :
        # print("Looking into {}th DNA : {}".format(counter, dna))
        curr_motif = ""
        curr_motif_score = float("-inf")

        for i in range(len(dna)-k+1) :
            curr_kmer = dna[i:i+k]
            curr_kmer_score = 1

            for j in range(len(curr_kmer)) :
                nucleotide = curr_kmer[j]
                curr_kmer_score = curr_kmer_score * profile[nucleotide][j]
            # print("{} has score of {}".format(curr_kmer, curr_kmer_score))

            if curr_motif_score < curr_kmer_score :
                curr_motif = curr_kmer
                curr_motif_score = curr_kmer_score

        motifs.append(curr_motif)
        counter = counter + 1
    
    return motifs


def find_consensus_string(kmers):

    consensus = []

    for i in range(len(kmers[0])):
        '''
        TODO : Need to change this dictionary to maybe list or some other ordered data store.
        It is biased because if occurences of "C" and "T" are equal, it will take "C", even if
        the kmers tell that "T" should be first.
        '''
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
    # random.seed(1)

    # Randomly select one motif per each dna as a starting point
    for dna in dnas :
        rand_idx = random.randint(0, len(dna)-k)
        random_kmer = dna[rand_idx:rand_idx+k]
        best_motifs.append(random_kmer)
        # print("{} -> {}".format(dna, random_kmer))
    motifs = best_motifs.copy()

    counter = 0
    while True :
        profile = build_profile_matrix(motifs, pseudocounts=True)
        # print("\nCalculated Profile : ")
        # for key, val in profile.items() : print("{} : {}".format(key, val))

        motifs = find_motifs_from_profile2(profile, dnas)
        # print("Motifs selected through MostProbableKmer()...")
        # print(motifs)

        next_motifs = []
        # Score all kmers in each dna against the profile
        for dna in dnas :
            # print(dna)
            dna_kmers_scores_dict = {}
            for i in range(len(dna)-k+1) :
                kmer = dna[i:i+k]
                kmer_probability_score = calculate_probability(kmer, profile)
                dna_kmers_scores_dict[kmer] = kmer_probability_score
                # print("{} : {}".format(kmer, kmer_probability_score))
            # After scoring all kmers in one dna strand, get the max
            next_motif = max(dna_kmers_scores_dict, key=dna_kmers_scores_dict.get)
            # print("Current dna's new motif is {}".format(next_motif))


        # print("**********")
        # print("**********")
        # print("**********")
        # print("**********")

        # for motif in motifs :
        #     print(motif)
        
        # print("\nBased off of profile, motifs are now ...\n{}".format(motifs))
        consensus_motif = find_consensus_string(motifs)
        # print("------")
        # print("consensus : {}".format(consensus_motif))

        motifs_score = 0
        for motif in motifs : 
            motifs_score = motifs_score + hammingDistance(consensus_motif, motif)
        
        # print("Motif score : {}".format(motifs_score))
        # Calculate score for the first best_motifs because best_score is 0 due to initialization
        if not best_motifs_score :
            consensus_best_motif = find_consensus_string(best_motifs)
            
            best_motifs_score = 0
            for best_motif in best_motifs : 
                best_motifs_score = best_motifs_score + hammingDistance(consensus_best_motif, best_motif)
        
        # print("\nBest Motifs is...\n{}".format(best_motifs))
        # print("Best motif score : {}".format(best_motifs_score))
        # print(".........")
        # print("Motifs (score: {}) vs. Best Motifs (score: {})".format(motifs_score, best_motifs_score))
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
    nucleotides = ["A", "C", "G", "T"]
    pseudo_counts = False

    # Default to the dataset folder, if not provided
    if not args.file :
        args.file = dataset_path
    
    with open(args.file, 'r') as f :
        k, t = f.readline().split(" ")
        
        dna_strings = []
        [dna_strings.append(row.strip()) for i, row in enumerate(f.readlines())]
    
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