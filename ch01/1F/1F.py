#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1F
   Given a genome, find all positions that minimize the skew

"""
import argparse

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Given a DNA sequence, compute indices of min skew"
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

def min_skew(sequence: str) -> list:
    """compute the positions of minimum skew  (#G - #C)

    Args:
        sequence (str): DNA sequence

    Returns:
        list: indices of all positions minimizing skew
    """

    result = []
    skew = [0]
    end = 0
    for i in range(len(sequence)):
        base = (sequence[i].upper())
        if base == "C":
            delta = -1
        elif base == "G":
            delta = +1
        else:
            delta = 0
        skew.append(skew[end] + delta)
        end += 1

    min_value = min(skew)
    result = [i for i in range(len(skew)) if skew[i] == min_value]

    return result

def main():
    """main
    """
    args = parse_arguments()
    sequence = parse_file(args.data_file)

    result = min_skew(sequence)
    print(" ".join(str(e) for e in result))

if __name__ == "__main__":
    main()