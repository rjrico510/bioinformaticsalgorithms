# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 02 Problem 2F
   Randomized motif search

"""
import argparse
import collections
import random

# random.seed(0) # make deterministic


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(description="Greedy motif search")
    parser.add_argument("data_file", help="input - 1st line - k t; rest - t strings")
    parser.add_argument(
        "-n",
        "--iterations",
        help="# iterations",
        type=int,
        default=1000,
        required=False,
    )
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line - k t; rest - t strings

    Returns:
        tuple: (list of DNA strings, k. t)
    """

    with open(filename) as f:
        lines = f.readlines()
        (k, t) = (int(i) for i in lines[0].strip().split())
        dna = [line.strip().upper() for line in lines[1:]]
    return (dna, k, t)


def compute_probability(kmer: str, profile: dict) -> float:
    """Compute k-mer probability from profile

    Args:
        kmer (str): k-mer
        profile (dict): profile matrix

    Returns:
        float: probability of the given k-mer
    """

    result = 1.0
    for i in range(0, len(kmer)):
        result *= profile[kmer[i]][i]

    return result


def profile_most_probable_first(txt: str, k: int, profile: dict) -> str:
    """Compute median string

    If more than one of the same probability is encountered, return only the 1st

    Args:
        txt (str): text to search (DNA string)
        k (int): k-mer size
        profile (dict): profile matrix: key={ACGT}, value = list

    Returns:
        set: most probable k-mer
    """
    result = None
    probability = -1

    for i in range(0, len(txt) - k + 1):
        kmer = txt[i : i + k]
        this_probability = compute_probability(kmer, profile)
        if this_probability > probability:
            probability = this_probability
            result = kmer

    return result


def generate_profile_with_pseudocounts(motifs: list) -> dict:
    """Generate profile matrix
    Initialize to 1 for every value to ensure no 0 entries

    Args:
        motifs (list): list of motifs (uppercase)

    Returns:
        dict: profile; keys = [ACGT]
    """
    BASES = ["A", "C", "G", "T"]
    n_motifs = len(motifs)

    count = dict(zip(BASES, [[], [], [], []]))
    for i in range(0, len(motifs[0])):
        this_count = dict(zip(BASES, [1] * 4))
        for j in range(0, n_motifs):
            key = motifs[j][i]
            this_count[key] += 1
        for base in BASES:
            count[base].append(this_count[base])

    profile = {}
    for key in count:
        profile[key] = [v / (2 * n_motifs) for v in count[key]]

    return profile


def compute_score(motifs: list) -> int:
    """compute score for a list of motifs of identical length

    For each position, find the most common base.
    The score for that position = (# motifs) - (# motifs w/ most common base)
    The overall score = sum of scores for each position

    Args:
        motifs (list): DNA strings of identical length, uppercase

    Returns:
        int: score
    """
    BASES = ["A", "C", "G", "T"]

    score = 0

    n_motifs = len(motifs)

    for i in range(0, len(motifs[0])):
        count = dict(zip(BASES, [0] * 4))
        for j in range(0, n_motifs):
            key = motifs[j][i]
            count[key] += 1
        max_count = max(count.values())
        score += n_motifs - max_count

    return score


def random_motif_search(dna: list, k: int, t: int) -> list:
    """Random motif search

    Randomly select k-mers from dna -> initial motifs
    - Compute score
    - Create a profile
    Iteratively produce a new set of motifs from previous profile
    - update motifs, score if new score < old score
    - exit otherwise

    Args:
        dna (list): list of equal length DNA sequences (must be > 1)
        k (int): kmer size
        t (int): # DNA sequences
    Returns:
        list: motif matrix of lowest score from kmers
    """
    len_motif = len(dna[0])
    best_motifs = []
    for dnai in dna:
        start = random.randrange(0, len_motif - k + 1)
        best_motifs.append(dnai[start : start + k])
    best_score = compute_score(best_motifs)
    profile = generate_profile_with_pseudocounts(best_motifs)

    while True:
        motifs = []
        for dnai in dna:
            motifs.append(profile_most_probable_first(dnai, k, profile))
        score = compute_score(motifs)

        if score < best_score:
            best_motifs = motifs
            best_score = score
            profile = generate_profile_with_pseudocounts(best_motifs)
        else:
            break

    return best_motifs


def main():
    """main"""
    args = parse_arguments()
    (dna, k, t) = parse_file(args.data_file)

    all_motifs = collections.defaultdict(int)
    for i in range(args.iterations):
        motifs = tuple(random_motif_search(dna, k, t))
        all_motifs[motifs] += 1

    result = max(all_motifs, key=all_motifs.get)
    print(all_motifs[result])
    print("\n".join(result))


if __name__ == "__main__":
    main()
