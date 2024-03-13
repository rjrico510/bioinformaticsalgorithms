#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 02 Problem 2C
   Profile-most probable k-mer

   Given a profile matrix, find the probability of each k-mer in
   a string & report the most probable

"""
import argparse

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description=" Profile-most probable k-mer problem"
    )
    parser.add_argument("data_file", help="input - 1st line - string; 2nd - k; next 4 - profile [ACGT]")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line - k; remaining lines - DNA

    Returns:
        tuple: (k, list of DNA strings)
    """
    BASES = ["A", "C", "G", "T"]

    with open(filename) as f:
        lines = f.readlines()
        txt = lines[0].strip().upper()
        k = int(lines[1].strip())
        probabilities = []
        for line in lines[2:6]:
            probabilities.append([float(n) for n in line.strip().split()])
        profile = dict(zip(BASES, probabilities))
    return (txt, k, profile)


def compute_probability(kmer: str, profile: dict) -> float:
    """Compute k-mer probability from profile

    Args:
        kmer (str): k-mer
        profile (dict): profile matrix

    Returns:
        float: probability of the given k-mer
    """

    result = 1.0
    for i in range(0,len(kmer)):
        result *= profile[kmer[i]][i]

    return result


def profile_most_probable(txt: str, k: int, profile: dict) -> set:
    """Compute median string

    Args:
        txt (str): text to search (DNA string)
        k (int): k-mer size
        profile (dict): profile matrix: key={ACGT}, value = list

    Returns:
        set: most probable k-mer(s)
    """
    result = set()
    probability = 0
    
    for i in range(0, len(txt) - k + 1):
        kmer = txt[i:i+k]
        this_probability = compute_probability(kmer, profile)
        if this_probability > probability:
            probability = this_probability
            result = set([kmer])
        elif this_probability == probability:
            result.add(kmer)
            
    return result


def main():
    """main
    """
    args = parse_arguments()
    (txt, k, profile) = parse_file(args.data_file)

    result = profile_most_probable(txt, k, profile)
    print(" ".join(result))

if __name__ == "__main__":
    main()