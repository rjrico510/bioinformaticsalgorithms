# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1M
   Implement number_to_pattern

"""
import argparse


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(description="pattern_to_number")
    parser.add_argument("data_file", help="2-line file - index & k")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 2-line file - index & k

    Returns:
        tuple: (index, k)
    """
    with open(filename) as f:
        lines = f.readlines()
        index = int(lines[0].strip())
        k = int(lines[1].strip())
    return (index, k)


def number_to_symbol(i: int) -> str:
    """convert int [0-3] to [ACGT]

    Args:
        i (int): [0-3]

    Returns:
        str: single char in [ACGT]
    """
    SYMBOLS = ["A", "C", "G", "T"]
    return SYMBOLS[i]


def number_to_pattern(index: int, k: int) -> str:
    """Convert integer index to pattern

    Args:
        index: index corresponding to pattern
        k: k-mer length

    Returns:
        str: pattern
    """
    if k == 1:
        return number_to_symbol(index)

    prefix_index = index // 4
    remainder = index % 4
    symbol = number_to_symbol(remainder)
    prefix_pattern = number_to_pattern(prefix_index, k - 1)
    return "".join([prefix_pattern, symbol])


def main():
    """main"""
    args = parse_arguments()
    (index, k) = parse_file(args.data_file)

    result = number_to_pattern(index, k)
    print(f"{result}")


if __name__ == "__main__":
    main()
