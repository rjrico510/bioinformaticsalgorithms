# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3C
   Overlap graph (as adjacency list)
"""
import argparse
import collections


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(description=" Overlap graph (as adjacency list)")
    parser.add_argument("data_file", help="input - list of k-mers")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - list of k-mers

    Returns:
        list: str k-mers
    """
    with open(filename) as f:
        lines = f.readlines()
        kmers = tuple([line.strip() for line in lines])
    return kmers


def construct_overlap_graph(kmers: list) -> dict:
    """given sequence of n k-mers
    Construct an adjacency list where for k-mer(i),
    another k-mer(j) is adjacent if suffix(k-mer(i)) = prefix(k-mer(j))

    Note that a kmer may be adjacent to itself

    Args:
        kmer (list): list of k-mers in order

    Returns:
        dict: key=kmer, value=list of ajacent kmers
    """
    result = collections.defaultdict(list)
    for this_kmer in kmers:
        for kmer in kmers:
            if this_kmer[1:] == kmer[0:-1]:
                result[this_kmer].append(kmer)
    return result


def print_adjacency_list(adjacency_list: dict) -> None:
    """Print adjacency list
    In lexographic order
    k-mer -> k-mer1

    Args:
        adjacency_list (dict): adjacency list as a dictionary
    """
    keys = list(adjacency_list.keys())
    keys.sort()
    for key in keys:
        values = adjacency_list[key]
        values.sort()
        for value in values:
            print(f"{key} -> {value}")


def main():
    """main"""
    args = parse_arguments()
    kmers = parse_file(args.data_file)

    result = construct_overlap_graph(kmers)
    print_adjacency_list(result)


if __name__ == "__main__":
    main()
