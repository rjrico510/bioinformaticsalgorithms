#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 02 Problem 2A
   Motif enumeration

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
        description="find all (k,d)-motifs from a list of DNA strings"
    )
    parser.add_argument("data_file", help="input - 1st line - k d; remaining lines - DNA")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line - k d; remaining lines - DNA

    Returns:
        tuple: (list of DNA strings, k, d)
    """
    with open(filename) as f:
        lines = f.readlines()
        tokens = lines[0].strip().split()
        (k, d) = (int(t) for t in tokens)
        dna = [line.strip() for line in lines [1:]]
    return (dna, k, d)

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

def motif_enumeration(dna: list, k: int, d: int) -> set:
    """Find all (k,d) motifs that appear in every DNA string
    (k,d)-motif = k-mer w/ maximum d allowed mismatches

    Args:
        dna (list): list of DNA sequences to search
        k (int): k-mer length
        d (int): maximum allowed Hamming distance

    Returns:
        list: all (k,d)-motifs in dna
    """
    result = set()

    for i in range(0, len(dna[0]) - k + 1):
        pattern = dna[0][i:i+k]
        current_neighbors = neighbors(pattern, d)
        for current_neighbor in current_neighbors: # iterate through all the neighbors
            match_all = True # assume the current neighbor is a match in all subsequent strings

            for current_dna in dna[1:]:
                match_current = False # assume there is no match for the current neighbor
                
                for i in range(0, len(current_dna) - k + 1):
                    hamming_distance = hamming(current_neighbor, current_dna[i:i+k])
                    if hamming_distance <= d:
                        match_current = True # found a match - can stop searching this string
                        break
                if not match_current:  # current neighbor not in current DNA
                    match_all = False
                    break # no need to continue searching the remaining DNA for this neighbor

            if match_all:
                result.add(current_neighbor)
                
    return result


def main():
    """main
    """
    args = parse_arguments()
    (dna, k, d) = parse_file(args.data_file)

    result = motif_enumeration(dna, k, d)
    print(" ".join(result))

if __name__ == "__main__":
    main()