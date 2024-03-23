# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3E
   DeBruijn graph (as adjacency list) from k-mers
"""
import argparse
import collections


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description=" DeBruijn graph (as adjacency list) from k-mers"
    )
    parser.add_argument("data_file", help="input - list of kmers")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> list:
    """Parse file

    Args:
        filename (str): file - list of k-mers

    Returns:
       list: list of k-mer strings
    """
    with open(filename) as f:
        lines = f.readlines()
        kmers = [line.strip() for line in lines]
    return kmers


def construct_debruijn_graph_kmers(kmers: list) -> dict:
    """given a list of k-mers
    Construct an adjacency list where for k-mer(i),
    another k-mer(j) is adjacent if suffix(k-mer(i)) = prefix(k-mer(j))

    Note that a kmer may be adjacent to itself

    Args:
        kmer (list): list of k-mers in order

    Returns:
        dict: key=kmer, value=list of adjacent kmers
    """
    result = collections.defaultdict(list)
    for kmer in kmers:
        result[kmer[0:-1]].append(kmer[1:])
    return result


def print_adjacency_list(adjacency_list: dict) -> None:
    """Print adjacency list
    In lexographic order
    k-mer -> k-mer1, k-mer2, ...

    Args:
        adjacency_list (dict): adjacency list as a dictionary
    """
    keys = list(adjacency_list.keys())
    keys.sort()
    for key in keys:
        values = adjacency_list[key]
        values.sort()
        print(f"{key} -> {','.join(values)}")


def main():
    """main"""
    args = parse_arguments()
    kmers = parse_file(args.data_file)

    result = construct_debruijn_graph_kmers(kmers)
    print_adjacency_list(result)


if __name__ == "__main__":
    main()
