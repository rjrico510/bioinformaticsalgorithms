# /usr/bin/env python3
"""Bioinformatics Algorithms Ch 01 Problem 1H
   Find all approximate occurrences of a pattern in a string

"""
import argparse


def parse_arguments() -> argparse.Namespace:
    """parse arguments

    Returns:
        argparse.Namespace: argument object
    """
    parser = argparse.ArgumentParser(
        description="Find all approximate occurrences of a pattern in a string"
    )
    parser.add_argument("data_file", help="3 line file - pattern, text, int")
    args = parser.parse_args()
    return args


def parse_file(filename: str) -> tuple:
    """Parse file

    Args:
        filename (str): file - 3 lines

    Returns:
        tuple: (txt, pattern, distance)
    """
    with open(filename) as f:
        lines = f.readlines()
        pattern = lines[0].strip()
        txt = lines[1].strip()
        distance = int(lines[2].strip())
    return (txt, pattern, distance)


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


def approx_match(txt: str, pattern: str, distance: int) -> list:
    """Find all approximate occurrences of a pattern in a string

    Args:
        txt (str): text to search
        pattern (str): pattern to match
        distance (int): Hamming distance (d)

    Returns:
        list: start positions of all instances of pattern in txt with < d mismatches
    """
    result = []

    len_txt = len(txt)
    len_pattern = len(pattern)

    for i in range(0, len_txt - len_pattern + 1):
        if hamming(pattern, txt[i : i + len_pattern]) <= distance:
            result.append(i)

    return result


def main():
    """main"""
    args = parse_arguments()
    (txt, pattern, distance) = parse_file(args.data_file)

    result = approx_match(txt, pattern, distance)
    print(f"{' '.join(str(e) for e in result)}")


if __name__ == "__main__":
    main()
