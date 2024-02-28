#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1C
   DNA reverse complement 

"""
import argparse
import collections

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Given a DNA sequence, generate the reverse complement"
    )
    parser.add_argument("data_file", help="1 line DNA sequence")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> str:
    """Parse file

    Args:
        filename (str): file - DNA sequence

    Returns:
        str: dna_sequence
    """
    with open(filename) as f:
        lines = f.readlines()
        dna_sequence = lines[0].strip()
    return dna_sequence

def reverse_complement_dna(dna_sequence: str) -> str:
    """compute the reverse complement of a DNA sequence
        Assumes only ATCG in the string

    Args:
        dna_sequence (str): DNA sequence

    Returns:
        str: reverse complement of dna_sequence (upper case)
    """

    complement_map = {
        "A": "T",
        "C": "G",
        "G": "C", 
        "T": "A"
    }

    result = ""
    for base in reversed(dna_sequence.upper()):
        result = "".join([result, complement_map[base]]) 
    return result

def main():
    """main
    """
    args = parse_arguments()
    dna_sequence = parse_file(args.data_file)

    result = reverse_complement_dna(dna_sequence)
    print(result)

if __name__ == "__main__":
    main()