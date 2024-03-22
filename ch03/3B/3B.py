# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3B
   String spelled by genome path
"""
import argparse


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(description="String spelled by genome path")
    parser.add_argument("data_file", help="input - k-mers in order")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 1st line - k; 2nd line - string

    Returns:
        tuple: (k, txt)
    """
    with open(filename) as f:
        lines = f.readlines()
        kmers = [line.strip() for line in lines]
    return kmers


def spell_string(kmers: list) -> str:
    """given sequence of n k-mers where
    last k-1 symbols in kmer(i) = 1st k-1 symbols in kmer(i+1)
    Output text of length k+n-1 such that ith kmer = ith kmer in text

    Args:
        kmer (list): list of k-mers in order

    Returns:
        str: text from overlapping k-mers
    """
    result = kmers[0]
    for i in range(1, len(kmers)):
        result += kmers[i][-1]
    return result


def main():
    """main"""
    args = parse_arguments()
    (kmers) = parse_file(args.data_file)

    result = spell_string(kmers)
    print(result)


if __name__ == "__main__":
    main()
