import argparse
import random
import os
import sys
from pathlib import Path
DATASET_DIR = os.path.join(os.getcwd(), 'datasets')


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

    # Randomly pick kmers from each of the dna strings
    for dna in dnas :
        print(dna)
        rind = random.randint(0, len(dna)-k)
        kmer = dna[rind:rind+k]
        motifs.append(kmer)
    
    print("\nSelected motifs are...")
    for m in motifs : print(m)

    # Perform sampling random dna to ignore / profiling motifs for N times
    for i in range(0, sample_times) :
        rind = random.randint(0, len(dnas)-1)
        motifs.pop(rind)
        print("\nExcluding {}th DNA string...".format(rind))
        for m in motifs : print(m)
        profile = build_profile_matrix(motifs, pseudocounts=True)
        for key, val in profile.items() : print(key, val)
        break



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
    
    for i in range(1000) :
        gibbsSampler(dna_strings, int(k), int(t), int(N))
        break
    

if __name__ == "__main__":
    main()