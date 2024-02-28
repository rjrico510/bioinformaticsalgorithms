#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1I
   Find the most frequent k-mers with mismatches <= d in a string

"""
import argparse
import collections
import copy

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Find the most frequent k-mers with mismatches <= d in a string"
    )
    parser.add_argument("data_file", help="input - 1st line - string; 2nd line k d")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - string, k, d

    Returns:
        tuple: (txt, k, d)
    """
    with open(filename) as f:
        lines = f.readlines()
        txt = lines[0].strip()
        tokens = lines[1].strip().split()
        (k, d) = (int(t) for t in tokens)
    return (txt, k, d)

def hamming(str1: str, str2: str) -> int:
    """Find the Hamming distance between 2 strings of equal length

    Args:
        str1 (str): string 1
        str2 (str): string 2

    Returns:
        int: Hamming distance
    """
    result = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            result += 1

    return result

def neighbors(pattern: str, d: int) -> set:
    """Given a pattern, find all strings matching with Hamming distance <= d
    Recursive function

    Args:
        pattern (str): base pattern
        d (int): Max Hamming distance

    Returns:
        set: all strings of Hamming distance <=d from pattern
    """
    NUCLEOTIDES = {"A", "C", "G", "T"}

    if d == 0: # only an exact match - just return the pattern
        return pattern
    if len(pattern) == 1: # the pattern is only 1 base long and d > 0
        return copy.deepcopy(NUCLEOTIDES)
    
    # recursively find suffix strings with Hamming distance <= d
    # For each suffix string 
    # - prefix w/ all nucleotides if Hamming distance < d
    # - prefix w/ only the 1st symbol of the original pattern id Hamming distance = d
    neighborhood = set()
    first_symbol = pattern[0]
    suffix_pattern = pattern[1:]
    suffix_neighbors = neighbors(suffix_pattern, d)
    for suffix_neighbor in suffix_neighbors:
        if hamming(suffix_pattern, suffix_neighbor) < d:
            neighborhood.update(("".join([n, suffix_neighbor]) for n in NUCLEOTIDES))
        else:
            neighborhood.add("".join([first_symbol, suffix_neighbor]))
    return neighborhood

def find_most_frequent_words_mismatches(txt: str, k: int, d: int) -> set:
    """Find the most frequent k-mers in text with <= d mismatches

    Args:
        txt (str): text to search
        k (int): k-mer length
        d (int): maximum allowed Hamming distance

    Returns:
        list: most frequent strings (w/ mismatches)
    """
    len_txt = len(txt)
    kmers = collections.defaultdict(int)

    for i in range(0, len_txt - k + 1):
        current_neighbors = neighbors(txt[i:i+k], d)
        for c in current_neighbors:
            kmers[c] += 1

    max_count = max(kmers.values())
    result = {key for key, value in kmers.items() if value == max_count}

    return result


def main():
    """main
    """
    args = parse_arguments()
    (txt, k, d) = parse_file(args.data_file)

    result = find_most_frequent_words_mismatches(txt, k, d)
    print(" ".join(result))

if __name__ == "__main__":
    main()