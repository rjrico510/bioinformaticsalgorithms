#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1N
   d-neighborhood of a string
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
    parser.add_argument("data_file", help="input - 1st line - string; 2nd line d")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - string, d

    Returns:
        tuple: (txt, d)
    """
    with open(filename) as f:
        lines = f.readlines()
        txt = lines[0].strip()
        d = int(lines[1].strip())
    return (txt, d)

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


def main():
    """main
    """
    args = parse_arguments()
    (txt, d) = parse_file(args.data_file)

    result = neighbors(txt, d)
    print("\n".join(result))

if __name__ == "__main__":
    main()