# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3D
   DeBruijn graph (as adjacency list)
"""
import argparse
import collections


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(description=" DeBruijn graph (as adjacency list)")
    parser.add_argument("data_file", help="input - k, txt")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - list of k-mers

    Returns:
       tuple: kmer size, txt
    """
    with open(filename) as f:
        lines = f.readlines()
        k = int(lines[0].strip())
        txt = lines[1].strip()
    return (k, txt)


def construct_debruijn_graph(k: int, txt: str) -> dict:
    """given k-mer size & a text string
    Construct an adjacency list where for k-mer(i),
    another k-mer(j) is adjacent if suffix(k-mer(i)) = prefix(k-mer(j))

    Note that a kmer may be adjacent to itself

    Args:
        kmer (list): list of k-mers in order

    Returns:
        dict: key=kmer, value=list of ajacent kmers
    """
    result = collections.defaultdict(list)
    for i in range(0, len(txt) - k + 1):
        result[txt[i : i + k - 1]].append(txt[i + 1 : i + k])
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
    (k, txt) = parse_file(args.data_file)

    result = construct_debruijn_graph(k, txt)
    print_adjacency_list(result)


if __name__ == "__main__":
    main()
