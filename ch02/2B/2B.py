#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 02 Problem 2B
   Median string problem

   minimize d(Pattern, DNA) = sum(d(pattern, DNA(i)))
   DNA = list of DNA strings

"""
import argparse
import sys

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="find k-mer pattern minimizing d(pattern, dna)"
    )
    parser.add_argument("data_file", help="input - 1st line - k; remaining lines - DNA")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line - k; remaining lines - DNA

    Returns:
        tuple: (k, list of DNA strings)
    """
    with open(filename) as f:
        lines = f.readlines()
        k = int(lines[0].strip())
        dna = [line.strip() for line in lines [1:]]
    return (k, dna)

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

def number_to_symbol(i:int) -> str:
    """convert int [0-3] to [ACGT]
    
    Args:
        i (int): [0-3]

    Returns:
        str: single char in [ACGT]
    """
    SYMBOLS = ["A", "C", "G", "T"]
    return SYMBOLS[i]

def number_to_pattern(index: int, k:int) -> str:
    """Convert integer index to pattern

    Args:
        index (int): index corresponding to pattern
        k (int): k-mer length

    Returns:
        str: pattern
    """
    if k == 1:
        return number_to_symbol(index)
    
    prefix_index = index // 4
    remainder = index % 4
    symbol = number_to_symbol(remainder)
    prefix_pattern = number_to_pattern(prefix_index, k - 1)
    return  "".join([prefix_pattern, symbol])


def distance_between_strings(pattern, dna):
    """Total minimum hamming distance between a pattern and a list of DNA strings

    Args:
        pattern (str): pattern to match
        dna (list): list of DNA strings

    Returns:
        int: minimum hamming distance
    """
    result = 0
    k = len(pattern)
    for this_dna in dna:
        this_hd = sys.maxsize
        for i in range(0, len(this_dna) - k):
            this_d = hamming(pattern, this_dna[i:i+k])
            # note - seems if this_d is ever 0 you could skip the rest of the kmers in this_dna
            if this_d < this_hd:
                this_hd = this_d
        result += this_hd

    return result


def median_string(k: int, dna: list) -> set:
    """Compute median string

        minimize d(pattern, dna) = sum(d(pattern, DNA(i)))
        dna = list of DNA strings
        pattern = string of length k
        d = hamming distance between patter & DNA(i)

        Note - if more than 1 hit - return only the 1st

    Args:
        k (int): k-mer length
        dna (list): list of DNA strings

    Returns:
        set: pattern(s) minimizing hamming distance
    """
    result = None
    distance = sys.maxsize

    for i in range(0, (4**k)-1):
        pattern = number_to_pattern(i, k)
        this_distance  = distance_between_strings(pattern, dna)
        if distance > this_distance:
            result = set([pattern])
            distance = this_distance
        elif distance == this_distance:
            result.add(pattern)

    return result


def main():
    """main
    """
    args = parse_arguments()
    (k, dna) = parse_file(args.data_file)

    result = median_string(k, dna)
    print(" ".join(result))

if __name__ == "__main__":
    main()