# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 3G
   Find an Eulerian path in an Eulerian directed graph
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
        description="Find an Eulerian path in an Eulerian directed graph"
    )
    parser.add_argument("data_file", help="input - graph as adjacency list")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> dict:
    """Parse file

    Args:
        filename (str): file - graph as adjacency list - lines of form i -> j,k,...

    Returns:
       dict: adjacency list as a dictionary
    """
    with open(filename) as f:
        lines = f.readlines()
        graph_entries = [line.strip() for line in lines]
        graph = generate_adjacency_list(graph_entries)
    return graph


def generate_adjacency_list(graph_entries: list) -> dict:
    """Generate an adjacency list from file input

    Args:
        graph_entries (list): entries of form key -> value1,value2,...

    Returns:
        dict: adjacency_list (each value is a list)
    """
    result = collections.defaultdict(list)
    for graph_entry in graph_entries:
        tokens = graph_entry.split()
        key = tokens[0]
        value = tokens[2].split(",")
        result[key] = value
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


def print_cycle(cycle: list) -> None:
    """Print cycle as
    node(1)->node(2)->...->node(n-1)->node(n)

    Args:
        cycle (list): Eulerian cycle - each entry is a node
    """
    print(f"{'->'.join(cycle)}")


def main():
    """main"""
    args = parse_arguments()
    graph = parse_file(args.data_file)

    result = find_eulerian_path(graph)
    print_cycle(result)


if __name__ == "__main__":
    main()
