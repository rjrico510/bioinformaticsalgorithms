# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1E
   Given a genome, and the following parameters:
   k = k-mer size
   L = length of an interval in the genome
   t = minimum frequency of any k-mer in L

   Return all k-mers that appear in any L at least t times

"""
import argparse
import collections


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Find all distinct k-mers forming (L,t) clumps"
    )
    parser.add_argument("data_file", help="input - 1st line - genome; 2nd line k L t")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - genome, k, L, t

    Returns:
        tuple: (genome, k, L, t)
    """
    with open(filename) as f:
        lines = f.readlines()
        genome = lines[0].strip()
        tokens = lines[1].strip().split()
        (k, interval_length, min_frequency) = (int(t) for t in tokens)
    return (genome, k, interval_length, min_frequency)


def find_min_frequency_kmer(txt: str, k: int, min_freq: int) -> set:
    """find strings of length k (k-mers) appearing at least a min # times

    Args:
        text (str): text to search
        k (int): k-mer length
        min_freq (int): minimum frequency

    Returns:
        set: most frequent k-mers
    """
    kmers = collections.defaultdict(int)
    len_txt = len(txt)
    for i in range(0, len_txt - k + 1):
        kmers[txt[i : i + k]] += 1

    result = {key for key, value in kmers.items() if value >= min_freq}
    return result


def find_clumps(genome: str, k: int, interval_length: int, min_frequency: int) -> set:
    """Find all distinct k-mers forming (L,t) clumps

    Args:
        genome (str): genome
        k (int): k-mer length (k)
        interval_length(int): interval length (L)
        min_freq (int): minimum frequency (t)

    Returns:
        set: all k-mers in any interval that appear in frequency >= t
    """

    result = set()
    for i in range(0, len(genome) - interval_length + 1):
        this_result = find_min_frequency_kmer(
            genome[i : i + interval_length], k, min_frequency
        )
        result.update(this_result)

    return result


def main():
    """main"""
    args = parse_arguments()
    (genome, k, interval_length, min_frequency) = parse_file(args.data_file)

    result = find_clumps(genome, k, interval_length, min_frequency)
    result = sorted(result)
    print(" ".join(result))


if __name__ == "__main__":
    main()
