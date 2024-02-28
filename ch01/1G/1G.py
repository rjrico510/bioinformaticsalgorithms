#/usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1G
   Find the Hamming distance between 2 strings

"""
import argparse

def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Find the Hamming distance between 2 strings - assume same length"
    )
    parser.add_argument("data_file", help="2-line file - 1 per string")
    args = parser.parse_args()
    return args

def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 2 lines

    Returns:
        tuple: (str1, str2)
    """
    with open(filename) as f:
        lines = f.readlines()
        str1 = lines[0].strip()
        str2 = lines[1].strip()
    return (str1, str2)

def hamming(str1: str, str2: str) -> int:
    """Find the Hamming distance between 2 strings of equal length

    Args:
        str1 (str): string 1
        str2 (str): string 2

    Returns:
        int: Hamming distance
    """
    result = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            result += 1

    return result

def main():
    """main
    """
    args = parse_arguments()
    (str1, str2) = parse_file(args.data_file)

    result = hamming(str1, str2)
    print(f"{result}")

if __name__ == "__main__":
    main()