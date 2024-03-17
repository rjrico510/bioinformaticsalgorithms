# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 03 Problem 2A
   String composition problem (all k-mers in a string)
"""
import argparse
import sys


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="String composition problem (all k-mers in a string)"
    )
    parser.add_argument("data_file", help="input - 1st line - k; 2nd line - string")
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
        k = int(lines[0].strip())
        txt = lines[1].strip()
    return (k, txt)


def composition_k(k: int, txt: str) -> list:
    """find all k-mers in a string

    Args:
        k (int): k-mer size
        txt (str):text string

    Returns:
        list: k-mers present in lexographic order
    """
    result = []
    for i in range(0, len(txt) - k + 1):
        result.append(txt[i : i + k])
    result.sort()
    return result


def main():
    """main"""
    args = parse_arguments()
    (k, txt) = parse_file(args.data_file)

    result = composition_k(k, txt)
    print("\n".join(result))


if __name__ == "__main__":
    main()
