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

def find_consensus_string(kmers):
    consensus = []

    for i in range(len(kmers[0])):
        nucleotides_dict = {"A": 0, "C": 0, "G": 0, "T": 0}

        for j in range(len(kmers)):
            curr_kmer = kmers[j]
            nucleotides_dict[curr_kmer[i]] = nucleotides_dict[curr_kmer[i]] + 1

        consensus.append( max(nucleotides_dict, key=nucleotides_dict.get) )

    return "".join(consensus)

def most_probable_kmer(dna, k, profile):
    kmer_prob = []

    for i in range(len(dna)-k+1):
        kmer = dna[i:i+k]
        probability = 1
        for j in range(len(kmer)) :
            probability = probability * profile[kmer[j]][j]
        
        kmer_prob.append((kmer, probability))
    
    return max(kmer_prob, key=lambda x: x[1])[0]

def build_profile_matrix(kmers, pseudocounts=False):
    k = len(kmers[0])
    # Initialize profile matrix
    if pseudocounts : profileMatrix = { "A": [1] * k, "C": [1] * k, "G": [1] * k, "T": [1] * k }    
    else : profileMatrix = { "A": [0] * k, "C": [0] * k, "G": [0] * k, "T": [0] * k }    
    
    for kmer in kmers :
        for i in range(len(kmer)) :
            profileMatrix[kmer[i]][i] = profileMatrix[kmer[i]][i] + 1

    for key, val in profileMatrix.items() :
        profileMatrix[key] = [_/(len(kmers)*2) for _ in val]

    return profileMatrix


def gibbsSampler(dnas, k, t, sample_times) :
    # Define all local vars
    motifs = []
    best_motifs = []
    best_motifs_score = float("inf")
    choices = list(range(0,len(dnas)-1))

    # Randomly pick kmers from each of the dna strings
    for dna in dnas :
        rind = random.randint(0, len(dna)-k)
        kmer = dna[rind:rind+k]
        motifs.append(kmer)
    
    # Perform sampling random dna to ignore / profiling motifs for N times
    for i in range(0, sample_times) :
        rind = random.choice(choices)
        motifs.pop(rind)

        profile = build_profile_matrix(motifs, pseudocounts=True)

        # Get probability sum of all the kmers from the excluded DNA        
        excluded_dna = dnas[rind]
        probable_kmer = most_probable_kmer(excluded_dna, k, profile)

        motifs.insert(rind, probable_kmer)

        consensus_motif = find_consensus_string(motifs)

        motifs_score = 0
        for motif in motifs : 
            motifs_score = motifs_score + hammingDistance(consensus_motif, motif)

        # Calculate score for the first best_motifs because best_score is 0 due to initialization
        if best_motifs :
            consensus_best_motif = find_consensus_string(best_motifs)
            
            best_motifs_score = 0
            for best_motif in best_motifs : 
                best_motifs_score = best_motifs_score + hammingDistance(consensus_best_motif, best_motif)

        if motifs_score < best_motifs_score :
            best_motifs = motifs[:]
            best_motifs_score = motifs_score

        # Redfining choices. Need to exclude the previously excluded index
        choices = list(range(0,len(dnas)-1))
        choices.pop(rind)
        
    return best_motifs, best_motifs_score


def parseArgs() :
    parser = argparse.ArgumentParser(
                prog="{}".format(__file__),
                formatter_class=argparse.RawDescriptionHelpFormatter,
                description='''\
                    Uses Gibbs Sampler, which is similar to randomizedMotifSearch, but uses
                    slightly different RNG in order to retrieve the best motifs.

                    Input File format :
                    ---------------------------------------
                    8 5 100
                    CGCCCCTCTCGGGGGTGTTCAGTAACCGGCCA
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
    if not args.file : args.file = dataset_path

    # Read from dataset file
    with open(args.file, 'r') as f :
        k, t, N = f.readline().split(" ")
        dna_strings = [row.strip() for i, row in enumerate(f.readlines())]
    
    motifs_score = float("inf")
    best_motifs = None
    for i in range(150) :
        motifs, sc = gibbsSampler(dna_strings, int(k), int(t), int(N))

        if sc < motifs_score :
            motifs_score = sc
            best_motifs = motifs
    
    for best_motif in best_motifs :
        print(best_motif)    
    print("Score : {}".format(motifs_score))
if __name__ == "__main__":
    main()