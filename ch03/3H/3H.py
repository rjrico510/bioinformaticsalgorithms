# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3H
   Reconstruct a string from its k-mer composition
   (assembled from 3G, 3B, 3E)
"""
import argparse
import collections
import copy


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Reconstruct a string from its k-mer composition"
    )
    parser.add_argument(
        "data_file", help="input - line 1 - k-mer length; remaining - list of k-mers"
    )
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - line 1 - k-mer length; remaining - list of k-mers

    Returns:
       tuple: (k, list of kmers)
    """
    with open(filename) as f:
        lines = f.readlines()
        k = int(lines[0].strip())
        kmers = [line.strip() for line in lines[1:]]
    return (k, kmers)


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


def find_unbalanced_nodes(graph: dict) -> tuple:
    """find unbalanced nodes (indegree != outdegree)
    in eulerian graph

    There should be exactly 1 node where indegree < outdegree (by 1)
    - this corresponds to the start node
    There should be exactly 1 node where outdegree < outdegree (by 1)
    - this corresponds to the end node

    Args:
        graph (dict): Eulerian directed graph as an adjacency list

    Returns:
        tuple: (node w in < out; node w/ in > out)
    """
    indegree = compute_indegree(graph)
    outdegree = compute_outdegree(graph)

    # unpack keys of each dictionary & take the union to get all keys
    all_nodes = set().union(*[indegree, outdegree])

    start_node = None
    end_node = None

    for node in all_nodes:
        node_indegree = indegree.get(node, 0)
        node_outdegree = outdegree.get(node, 0)
        if node_outdegree == node_indegree + 1:
            start_node = node
        elif node_outdegree + 1 == node_indegree:
            end_node = node

    return (start_node, end_node)


def compute_indegree(graph: dict) -> dict:
    """compute indegree for nodes in a graph

    Args:
        graph (dict): Eulerian directed graph as an adjacency list

    Returns:
        dict: key=node, value=# edges leading into the node
    """
    result = collections.defaultdict(int)
    for node in graph:
        for out_node in graph[node]:
            result[out_node] += 1
    return result


def compute_outdegree(graph: dict) -> dict:
    """compute outdegree for nodes in a graph

    Args:
        graph (dict): Eulerian directed graph as an adjacency list

    Returns:
        dict: key=node, value=# edges exiting the node
    """
    result = {k: len(v) for (k, v) in graph.items()}
    return result


def find_eulerian_path(graph: dict) -> list:
    """Find Eulerian path

    Args:
        graph (dict): Eulerian directed graph as an adjacency list

    Returns:
        list: sequence of nodes which traverse the path
    """
    (start_node, end_node) = find_unbalanced_nodes(graph)

    # balance the graph
    graph[end_node].append(start_node)  # assumes defaultdict(list)

    # now find the cycle
    nodes = list(graph.keys())
    num_edges = sum([len(graph[node]) for node in nodes])
    nodes.sort()
    start = start_node
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

    # reset the cycle back to a graph by
    # finding the edge end_node -> start_node
    # and reset
    start_test = 0
    while True:
        end_index = cycle.index(end_node, start_test)
        if cycle[end_index + 1] == start_node:
            start_index = end_index + 1
            cycle = cycle[start_index:-1] + cycle[0:start_index]
            break
        else:
            start_test = end_index + 1

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


def path_to_string(path: list, k: int) -> None:
    """Convert path of kmers of length k-1 to a string

    e.g. TTC, TCA, CAG -> TTCAG

    Args:
        path (list): Eulerian cycle - each entry is a node
        k (int): kmer length
    """
    result = ""
    if len(path) > 0:
        result = path[0]
        for i in range(1, len(path)):
            result += path[i][-1]
    return result


def reconstruct_string_from_kmers(kmers: list, k: int) -> str:
    """Construct string from k-mers

    Args:
        kmers (list): list of k-mers of equivalent length
        k (int): length of kmer

    Returns:
        str: string formed from k-mers
    """
    graph = construct_debruijn_graph_kmers(kmers)
    path = find_eulerian_path(graph)
    result = path_to_string(path, k)
    return result


def main():
    """main"""
    args = parse_arguments()
    (k, kmers) = parse_file(args.data_file)

    result = reconstruct_string_from_kmers(kmers, k)
    print(result)


if __name__ == "__main__":
    main()
