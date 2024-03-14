#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 02 Problem 2H
   Distance between pattern & string

   Note - this is essentially lifted from 2B
"""
import argparse
import sys

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="find sum of hamming distances between a pattern & list of strings"
    )
    parser.add_argument("data_file", help="input - 1st line - pattern; 2nd line - strings")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line - pattern; 2nd line - strings

    Returns:
        tuple: (pattern, list of DNA strings)
    """
    with open(filename) as f:
        lines = f.readlines()
        pattern = lines[0].strip()
        dna = lines[1].strip().split()
    return (pattern, dna)

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

def main():
    """main
    """
    args = parse_arguments()
    (pattern, dna) = parse_file(args.data_file)

    result = distance_between_strings(pattern, dna)
    print(result)

if __name__ == "__main__":
    main()