# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3I
   k-universal circular binary string
   (derived from 3H)
"""
import argparse
import collections
import copy


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(description="k-universal circular binary string")
    parser.add_argument("data_file", help="input - line 1 - k")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> int:
    """Parse file

    Args:
        filename (str): file - single line - k

    Returns:
       int: k
    """
    with open(filename) as f:
        lines = f.readlines()
        k = int(lines[0].strip())
    return k


def construct_kmers(k: int) -> list:
    """construct all binary kmers of length k
    e.g. k=1: [0,1]; k=2: [00,01,10,11]

    Args:
        k (int): k-mer length

    Returns:
        list: list of binary k-mers as strings
    """
    if k == 1:
        return ["0", "1"]
    else:
        return ["0" + n for n in construct_kmers(k - 1)] + [
            "1" + n for n in construct_kmers(k - 1)
        ]

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


def find_eulerian_cycle(graph: dict) -> list:
    """Find Eulerian cycle

    Args:
        graph (dict): Eulerian directed graph as an adjacency list

    Returns:
        list: sequence of nodes which traverse the cycle
    """
    nodes = list(graph.keys())
    num_edges = sum([len(graph[node]) for node in nodes])
    nodes.sort()
    start = nodes[0]
    cycle = [start]

    this_graph = copy.deepcopy(graph)
    unexplored_in_cycle = []  # nodes in cycle w/ unexplored edges

    while True:
        current_node = cycle[-1]
        while this_graph.get(current_node) is not None:
            cycle.append(this_graph[current_node].pop(0))
            unexplored_in_cycle = update_unexplored(
                unexplored_in_cycle, this_graph, current_node
            )
            current_node = cycle[-1]

        if len(cycle) == num_edges + 1:
            break
        else:  # save cycle but start at the new position
            start = unexplored_in_cycle[0]
            start_index = cycle.index(start)
            cycle = cycle[start_index:-1] + cycle[0 : start_index + 1]

    return cycle


def update_unexplored(
    unexplored_in_cycle: list, this_graph: dict, current_node: str
) -> list:
    """Update list of unexplored nodes

    Args:
        unexplored_in_cycle (list): current list of nodes w/ 1+ unexplored edge
        this_graph (dict): graph
        current_node (str): node to update

    Returns:
        list: updated list of unexplored nodes
    """
    if len(this_graph[current_node]) > 0:
        unexplored_in_cycle.append(current_node)
    else:
        del this_graph[current_node]
        unexplored_in_cycle = [u for u in unexplored_in_cycle if u != current_node]
    return unexplored_in_cycle


def cycle_to_string(cycle: list, k: int) -> None:
    """Convert cycle of (k-1)mers of length k-1 to a circular string contaning all binary k-mers

    e.g. 00, 01, 11, 10, 00 -> 0011
    Drop the last kmer as it is a repeat (part of the cycle)

    Args:
        cycle (list): Eulerian cycle - each entry is a node
        k (int): kmer length
    """
    result = ""
    if len(cycle) > 0:
        result = cycle[0]
    if len(cycle) > 1:
        for i in range(1, len(cycle)):
            result += cycle[i][-1]
    result = result[0 : -k + 1]  # drop the last k-mer
    return result


def construct_k_universal_circular_binary_string(k: int) -> str:
    """Construct k-universal circular binary string

    Args:
        k (int): length of kmer

    Returns:
        str: string formed from k-mers
    """
    kmers = construct_kmers(k)
    graph = construct_debruijn_graph_kmers(kmers)
    cycle = find_eulerian_cycle(graph)
    result = cycle_to_string(cycle, k)
    return result


def main():
    """main"""
    args = parse_arguments()
    k = parse_file(args.data_file)

    result = construct_k_universal_circular_binary_string(k)
    print(result)


if __name__ == "__main__":
    main()
